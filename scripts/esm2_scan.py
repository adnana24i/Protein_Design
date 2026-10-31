"""ESM2 zero-shot single-mutant scan of the trastuzumab scFv (README step 4, source 1).

WHAT THIS DOES
    For every position of the 242 aa parent scFv and every one of the 19 alternative amino
    acids, ask ESM2 how much more (or less) plausible the mutant residue is than the parent
    residue in this sequence context. No training and no structure are used.

        score(pos, wt -> mut) = log p(mut | context) - log p(wt | context)      (log-likelihood ratio, LLR)

    LLR > 0 : ESM2 finds the mutant at least as plausible as the parent residue.
    LLR << 0: ESM2 considers the change implausible (typically buried/structural or conserved residues).

SCORING MODES (both computed, Meier et al. 2021, "Language models enable zero-shot prediction of
the effects of mutations on protein function")
    masked : mask position i, run the model, read the distribution at i.  242 forward passes.
             Slower, usually the more accurate mode. This is the primary score (column llr_masked).
    wt     : one forward pass on the unmasked parent, read every position's distribution at once.
             Fast, slightly less accurate (column llr_wt). Kept as a cross-check.

CONTEXT
    scfv  : the full 242 aa single chain (what the lab actually expresses). Default.
    split : VH (120 aa) and VL (107 aa) scored separately without the linker. Useful to check that
            conclusions do not depend on the linker / on the model seeing the whole scFv.

WHAT THE SCORE DOES *NOT* TELL YOU
    ESM2 learned "looks like a natural protein", not "binds HER2". High LLR means likely stable and
    foldable, not higher affinity. It is biased to common residues. Use the output as ONE vote in the
    multi-source consensus (README step 5), and never mutate HER2-contact residues on LLR alone.

BUILT-IN SANITY CHECKS (printed at the end of every run)
    * The disulfide cysteines (VH C22/C96, VL C23/C88) should have strongly negative LLRs for all
      substitutions. If they do not, something is wrong (indexing, tokenisation, wrong sequence).
    * llr_masked and llr_wt should correlate (Spearman typically > 0.6).
    * The parent sequence is checked against the structure-derived sequence in data/interface_contacts.csv.

INPUTS
    data/interface_contacts.csv   (from scripts/interface_contacts.py) for region/contact annotations

OUTPUTS (default paths)
    data/scores/<model>_scan.csv        all 242 x 19 = 4,598 single mutants, ranked
    data/scores/<model>_positions.csv   per-position summary (wt probability, entropy, mean LLR)
    data/scores/<model>_shortlist.csv   mutations passing the safe-to-propose filters (see below)

SUSPECT WT
    The model sometimes assigns a tiny probability to the PARENT residue (e.g. VL F83 = scFv 218, p about 1e-4,
    where ESM2-650M predicts Glu with 74 %). Then EVERY substitution gets a huge positive LLR and floods the top
    of the ranking, although the residue is a conserved framework Phe. This is treated as a model artefact: such
    positions are flagged (`suspect_wt`, threshold --suspect-p) and excluded from the shortlist. Do not mutate
    them on ESM2 evidence alone; require agreement from an independent model (AbLang2, ProFam-1) and the structure.

SHORTLIST FILTER (a starting point for README step 5, not a final decision)
    keep : LLR(masked) >= --min-llr (default 0.0), the mutation is not in the linker
    drop : suspect-wt positions (see above), the linker, the four disulfide cysteines, any mutation TO cysteine (unpaired-Cys liability),
           and HER2-contact residues / positions closer than --min-dist (default 4.5 A) to HER2
           (override with --allow-contacts).
    Proline substitutions are not filtered here; leave them to the structure-based and QC steps.

USAGE
    python3 scripts/esm2_scan.py --model facebook/esm2_t12_35M_UR50D      # quick smoke test
    python3 scripts/esm2_scan.py --model facebook/esm2_t33_650M_UR50D     # main run (about 2.6 GB download)
    python3 scripts/esm2_scan.py --context split

Mutation naming: <wt><scfv_index><mut>, 1-based over the 242 aa scFv (e.g. Y33A), plus a
domain-local label (VH:Y33A or VL:N30A, PDB numbering) in column `label`.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# Parent from the hackathon deck (slide 15)
VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAED"
      "TAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
LINKER = "GGGGSGGGGSGGGGS"
VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTT"
      "PPTFGQGTKVEIK")
PARENT = VH + LINKER + VL
AAS = "ACDEFGHIKLMNPQRSTVWY"
# scFv indices (1-based) of the cysteines forming the variable-domain disulfides: VH 22-96, VL 23-88 (+120+15)
DISULFIDE_CYS = {22, 96, 158, 223}


def load_annotations(path):
    """Per-residue region / HER2 distance from the structure; also verifies the parent sequence."""
    df = pd.read_csv(path)
    ann = df.set_index("scfv_index").sort_index()
    seq = "".join(ann.loc[i, "aa"] for i in range(1, 121)) + LINKER + "".join(ann.loc[i, "aa"] for i in range(136, 243))
    if seq != PARENT:
        sys.exit("ERROR: parent sequence in the script differs from the structure-derived sequence")
    return ann


def domain_of(i):
    return "VH" if i <= 120 else ("linker" if i <= 135 else "VL")


def label_of(i, wt, mut):
    d = domain_of(i)
    local = i if d == "VH" else i - 135
    return f"{d}:{wt}{local}{mut}" if d != "linker" else f"linker:{wt}{i}{mut}"


@torch.no_grad()
def log_probs(model, tok, seqs, device, batch_size):
    """log-softmax over the vocabulary for each token of each sequence. Returns list of [L+2, V] cpu tensors."""
    out = []
    for s in range(0, len(seqs), batch_size):
        enc = tok(seqs[s:s + batch_size], return_tensors="pt", padding=True).to(device)
        logits = model(**enc).logits
        out.extend(torch.log_softmax(logits.float(), dim=-1).cpu())
    return out


def score_segment(model, tok, seq, device, batch_size):
    """Return (wt_marginal[L,V-restricted], masked_marginal[L,20]) log-prob tables over the 20 standard AAs."""
    aa_ids = torch.tensor([tok.convert_tokens_to_ids(a) for a in AAS])
    # (a) wild-type marginals: one pass on the unmasked sequence; token i+1 is residue i (<cls> first)
    lp_wt = log_probs(model, tok, [seq], device, 1)[0][1:len(seq) + 1][:, aa_ids]
    # (b) masked marginals: one pass per position with that residue replaced by <mask>
    masked = []
    for i in range(len(seq)):
        masked.append(seq[:i] + tok.mask_token + seq[i + 1:])
    lp_m = log_probs(model, tok, masked, device, batch_size)
    lp_masked = torch.stack([lp_m[i][i + 1][aa_ids] for i in range(len(seq))])
    return lp_wt.numpy(), lp_masked.numpy()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--context", choices=["scfv", "split"], default="scfv")
    ap.add_argument("--device", default=None, help="cpu | mps | cuda (default: auto)")
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--annotations", default=os.path.join(ROOT, "data", "interface_contacts.csv"))
    ap.add_argument("--outdir", default=os.path.join(ROOT, "data", "scores"))
    ap.add_argument("--min-llr", type=float, default=0.0, help="shortlist: minimum masked LLR")
    ap.add_argument("--min-dist", type=float, default=4.5, help="shortlist: minimum distance (A) to HER2")
    ap.add_argument("--suspect-p", type=float, default=0.01,
                    help="flag positions where the model gives the PARENT residue p below this (see SUSPECT WT)")
    ap.add_argument("--allow-contacts", action="store_true", help="shortlist: keep HER2-contact positions too")
    args = ap.parse_args()

    device = args.device or ("cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu")
    tag = args.model.split("/")[-1].replace("_UR50D", "") + ("" if args.context == "scfv" else "_split")
    os.makedirs(args.outdir, exist_ok=True)
    ann = load_annotations(args.annotations)
    print(f"parent OK ({len(PARENT)} aa) | model {args.model} | context {args.context} | device {device}")

    tok = AutoTokenizer.from_pretrained(args.model)
    model = EsmForMaskedLM.from_pretrained(args.model).to(device).eval()

    # segments: (start index in scFv, sequence)
    if args.context == "scfv":
        segments = [(0, PARENT)]
    else:
        segments = [(0, VH), (135, VL)]

    wt_lp = np.full((len(PARENT), 20), np.nan)
    ms_lp = np.full((len(PARENT), 20), np.nan)
    for start, seq in segments:
        print(f"scoring segment of {len(seq)} aa starting at scFv index {start + 1} ...", flush=True)
        a, b = score_segment(model, tok, seq, device, args.batch_size)
        wt_lp[start:start + len(seq)] = a
        ms_lp[start:start + len(seq)] = b

    rows = []
    for i, wt in enumerate(PARENT):
        idx = i + 1
        if np.isnan(ms_lp[i, 0]):  # linker in split mode
            continue
        w = AAS.index(wt)
        pw = np.exp(ms_lp[i]); pw = pw / pw.sum()  # restricted to the 20 standard AAs
        for j, mut in enumerate(AAS):
            if mut == wt:
                continue
            a = ann.loc[idx] if idx in ann.index else None
            rows.append({
                "scfv_index": idx, "domain": domain_of(idx), "wt": wt, "mut": mut,
                "mutation": f"{wt}{idx}{mut}", "label": label_of(idx, wt, mut),
                "region": "linker" if a is None else a["region"],
                "min_dist_to_her2": np.nan if a is None else a["min_dist_to_her2"],
                "her2_contact": False if a is None else bool(a["contact_4.5A"]),
                "p_wt_masked": pw[w], "suspect_wt": bool(pw[w] < args.suspect_p),
                "llr_masked": ms_lp[i, j] - ms_lp[i, w],
                "llr_wt": wt_lp[i, j] - wt_lp[i, w],
            })
    df = pd.DataFrame(rows)
    df["rank_masked"] = df["llr_masked"].rank(ascending=False, method="min").astype(int)
    df = df.sort_values("rank_masked")
    scan_path = os.path.join(args.outdir, f"{tag}_scan.csv")
    df.to_csv(scan_path, index=False, float_format="%.4f")

    # ---- per-position summary: how "locked" is each residue according to the model ----
    pos = []
    for i, wt in enumerate(PARENT):
        if np.isnan(ms_lp[i, 0]):
            continue
        p = np.exp(ms_lp[i]); p = p / p.sum()  # restrict to the 20 standard AAs
        sub = df[df.scfv_index == i + 1]
        pos.append({"scfv_index": i + 1, "domain": domain_of(i + 1), "wt": wt,
                    "region": sub.region.iloc[0], "her2_contact": bool(sub.her2_contact.iloc[0]),
                    "min_dist_to_her2": sub.min_dist_to_her2.iloc[0],
                    "p_wt_masked": p[AAS.index(wt)], "entropy": float(-(p * np.log(p + 1e-12)).sum()),
                    "mean_llr": sub.llr_masked.mean(), "n_llr_ge0": int((sub.llr_masked >= 0).sum()),
                    "best_mut": sub.sort_values("llr_masked").iloc[-1]["mut"], "best_llr": sub.llr_masked.max()})
    pdf = pd.DataFrame(pos)
    pos_path = os.path.join(args.outdir, f"{tag}_positions.csv")
    pdf.to_csv(pos_path, index=False, float_format="%.4f")

    # ---- shortlist (first-pass filter only, see module docstring) ----
    keep = (df.llr_masked >= args.min_llr) & (df.domain != "linker") & (df.mut != "C") & (~df.scfv_index.isin(DISULFIDE_CYS)) & (~df.suspect_wt)
    if not args.allow_contacts:
        keep &= (~df.her2_contact) & (df.min_dist_to_her2 > args.min_dist)
    sl = df[keep]
    sl_path = os.path.join(args.outdir, f"{tag}_shortlist.csv")
    sl.to_csv(sl_path, index=False, float_format="%.4f")

    # ---- sanity checks ----
    print("\n=== sanity checks ===")
    cys = df[df.scfv_index.isin(DISULFIDE_CYS)]
    print(f"disulfide Cys: mean masked LLR over all substitutions = {cys.llr_masked.mean():.2f} "
          f"(max {cys.llr_masked.max():.2f}); expected strongly negative")
    rho = df[["llr_masked", "llr_wt"]].corr(method="spearman").iloc[0, 1]
    print(f"Spearman(masked, wt-marginal) = {rho:.3f}; typically > 0.6")
    print("fraction of mutations with LLR >= 0, by region:")
    print(df.assign(ok=df.llr_masked >= 0).groupby("region").ok.mean().round(3).to_string())
    sus = df[df.suspect_wt].drop_duplicates("scfv_index")
    print(f"\nSUSPECT WT: {len(sus)} position(s) where ESM2 gives the parent residue p < {args.suspect_p}: "
          + ", ".join(f"{r.wt}{r.scfv_index}({r.label.split(':')[0]}, p={r.p_wt_masked:.4f})" for r in sus.itertuples()))
    print(f"\ncontact positions: mean LLR {df[df.her2_contact].llr_masked.mean():.2f} vs "
          f"non-contact {df[~df.her2_contact & (df.domain != 'linker')].llr_masked.mean():.2f}")
    print(f"\nTop 15 single mutants by masked LLR (shortlist filter applied):")
    print(sl.head(15)[["label", "mutation", "region", "min_dist_to_her2", "llr_masked", "llr_wt"]].to_string(index=False))
    print(f"\nwrote:\n  {scan_path}  ({len(df)} rows)\n  {pos_path}\n  {sl_path}  ({len(sl)} rows)")


if __name__ == "__main__":
    main()
