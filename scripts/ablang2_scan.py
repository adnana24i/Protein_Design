"""AbLang2 zero-shot single-mutant scan of the trastuzumab scFv (README step 4, source 2).

AbLang2 is an antibody-specific paired heavy/light language model, so unlike ESM2 it sees VH and VL
together as an antibody rather than as a generic protein. Scoring is the same masked-marginal LLR:

    score(pos, wt -> mut) = log p(mut | masked context) - log p(wt | masked context)

INPUT FORMAT
    AbLang2 takes the pair as "VH|VL" (separator token '|'), NOT the scFv with its (GGGGS)3 linker.
    The linker is therefore not scored by this model (no rows for scFv 121-135), which is correct:
    the linker is fixed by the hackathon rules anyway.

POSITION MAPPING (verified by the disulfide check below)
    token 0..119    -> VH 1..120    -> scFv index 1..120
    token 120       -> '|' separator, skipped
    token 121..227  -> VL 1..107    -> scFv index 136..242

KNOWN BIAS
    AbLang2 is trained on natural antibody repertoires, so it pulls sequences towards germline /
    common residues. That is useful for expression and stability, and is NOT evidence of better
    affinity. Trastuzumab is a humanised mouse antibody, so some of its binding-critical CDR residues
    are deliberately non-germline and AbLang2 is expected to dislike them. Treat a high score as
    "looks like a normal antibody", never as "binds HER2 better".

OUTPUT  data/scores/ablang2_scan.csv  (standard schema read by scripts/evaluate.py)

USAGE   python3 -I scripts/ablang2_scan.py
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch

import ablang2

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from validate import DISULFIDE_CYS, LINKER, PARENT, VH, VL  # noqa: E402

AAS = "ACDEFGHIKLMNPQRSTVWY"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--device", default="cpu", help="cpu or cuda (AbLang2 is small; cpu is fine)")
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--annotations", default=os.path.join(ROOT, "data", "interface_contacts.csv"))
    ap.add_argument("--outdir", default=os.path.join(ROOT, "data", "scores"))
    ap.add_argument("--suspect-p", type=float, default=0.01,
                    help="flag positions where the model gives the PARENT residue less than this probability")
    args = ap.parse_args()

    ann = pd.read_csv(args.annotations).set_index("scfv_index").sort_index()
    model = ablang2.pretrained("ablang2-paired", device=args.device)
    model.freeze()
    tok = model.tokenizer
    mask_id = tok.aa_to_token["*"]
    aa_ids = torch.tensor([tok.aa_to_token[a] for a in AAS])

    paired = VH + "|" + VL
    enc = tok([paired], pad=True, w_extra_tkns=False, device=args.device)[0]
    assert len(enc) == len(VH) + 1 + len(VL), "unexpected tokenisation length"
    sep = len(VH)
    assert int(enc[sep]) == tok.aa_to_token["|"], "separator not where expected"

    # token index -> scFv index (1-based); the separator has no scFv position
    tok2scfv = {i: i + 1 for i in range(len(VH))}
    tok2scfv.update({i: len(VH) + len(LINKER) + (i - sep) for i in range(sep + 1, len(enc))})

    # masked marginals: one forward pass per scored position
    positions = sorted(tok2scfv)
    lp = {}
    for s in range(0, len(positions), args.batch_size):
        chunk = positions[s:s + args.batch_size]
        batch = enc.repeat(len(chunk), 1)
        for row, t in enumerate(chunk):
            batch[row, t] = mask_id
        with torch.no_grad():
            logits = model.AbLang(batch)
        for row, t in enumerate(chunk):
            lp[t] = torch.log_softmax(logits[row, t].float(), -1)[aa_ids].numpy()
        print(f"  scored {min(s + args.batch_size, len(positions))}/{len(positions)} positions", end="\r", flush=True)
    print()

    rows = []
    for t, idx in tok2scfv.items():
        wt = PARENT[idx - 1]
        assert wt == tok.token_to_aa[int(enc[t])], f"sequence mismatch at scFv {idx}"
        w = AAS.index(wt)
        p = np.exp(lp[t]); p = p / p.sum()
        a = ann.loc[idx]
        for j, mut in enumerate(AAS):
            if mut == wt:
                continue
            rows.append({"scfv_index": idx, "domain": a.domain, "wt": wt, "mut": mut,
                         "mutation": f"{wt}{idx}{mut}", "region": a.region,
                         "min_dist_to_her2": a.min_dist_to_her2, "her2_contact": bool(a["contact_4.5A"]),
                         "p_wt_masked": p[w], "suspect_wt": bool(p[w] < args.suspect_p),
                         "score": lp[t][j] - lp[t][w]})
    df = pd.DataFrame(rows)
    df["rank"] = df.score.rank(ascending=False, method="min").astype(int)
    df = df.sort_values("rank")
    os.makedirs(args.outdir, exist_ok=True)
    path = os.path.join(args.outdir, "ablang2_scan.csv")
    df.to_csv(path, index=False, float_format="%.4f")

    cys = df[df.scfv_index.isin(DISULFIDE_CYS)]
    print("\n=== sanity checks ===")
    print(f"disulfide Cys mean LLR {cys.score.mean():.2f} (max {cys.score.max():.2f}); expected strongly negative")
    sus = df[df.suspect_wt].drop_duplicates("scfv_index")
    print(f"suspect-wt positions (p(parent) < {args.suspect_p}): "
          + (", ".join(f"{r.wt}{r.scfv_index}" for r in sus.itertuples()) if len(sus) else "none"))
    print("fraction with LLR >= 0 by region:")
    print(df.assign(ok=df.score >= 0).groupby("region").ok.mean().round(3).to_string())
    print(f"\ncontact mean LLR {df[df.her2_contact].score.mean():.2f} vs non-contact {df[~df.her2_contact].score.mean():.2f}")
    print("\nTop 10 (no filtering):")
    print(df.head(10)[["mutation", "region", "min_dist_to_her2", "score"]].to_string(index=False))
    print(f"\nwrote {path} ({len(df)} rows; linker positions 121-135 are not scored by this model)")


if __name__ == "__main__":
    main()
