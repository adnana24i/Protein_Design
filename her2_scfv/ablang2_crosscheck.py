#!/usr/bin/env python3
"""AbLang2 antibody-specific cross-check of the ESM-picked CDR substitutions.

Run with the ablang-env python:
    ~/Desktop/protein-design/ablang-env/bin/python ablang2_crosscheck.py

Loads the paired AbLang2 model, splits the parent scFv at the (GGGGS)3 linker,
and scores every CDR-position substitution with AbLang2 likelihood mode
(stepwise/masked-marginal). Cross-references the 10 model variants' (V06-V15)
3 model-picked substitutions each against the earlier ESM2 masked-marginal scan,
classifies each pick CONSENSUS / NEUTRAL / CONFLICT, and proposes an
antibody-preferred alternative for every CONFLICT. Writes ablang2_crosscheck.md.
Reports only — does NOT rebuild variants.
"""
import csv, glob, re, math
from pathlib import Path
import numpy as np
import ablang2

HERE = Path(__file__).resolve().parent
PARENT = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRF"
          "TISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSGGGGSGGGGSGGGGS"
          "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRS"
          "GTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
assert len(PARENT) == 242
VH = PARENT[:120]            # construct 1-120
VL = PARENT[135:]            # construct 136-242
assert PARENT[120:135] == "GGGGSGGGGSGGGGS"

CDRS = {"H1": (26, 35), "H2": (50, 66), "H3": (99, 109),
        "L1": (159, 169), "L2": (185, 191), "L3": (224, 232)}
PROTECTED = {50, 59, 99, 105, 226}
CYS = {22, 96, 158, 223}
AA20 = "ACDEFGHIKLMNPQRSTVWY"
HARD = ["NG", "NS", "NT", "DG", "DS", "DT"]

def cdr_of(pos):
    for n, (a, b) in CDRS.items():
        if a <= pos <= b:
            return n
    return None

def construct_to_idx(pos):
    """construct position (1-indexed) -> index in the AbLang2 processed string
    '<VH>|<VL>' (verified empirically)."""
    return pos if pos <= 120 else pos - 12

def hard_sites(seq):
    return {m.start() for mo in HARD for m in re.finditer("(?=" + mo + ")", seq)}

PARENT_HARD = hard_sites(PARENT)

def new_hard(pos, mut):
    s = PARENT[:pos-1] + mut + PARENT[pos:]
    return sorted(hard_sites(s) - PARENT_HARD)

# ---------------- AbLang2 scoring ----------------
def ablang_scan():
    model = ablang2.pretrained(model_to_use="ablang2-paired", device="cpu")
    tok = model.tokenizer
    aa2tok = {a: tok.aa_to_token[a] for a in AA20}
    # masked-marginal (stepwise_masking) logits, no ANARCI alignment (needs pandas)
    out = model([[VH, VL]], mode="likelihood", align=False, stepwise_masking=True)
    logits = np.asarray(out[0])            # [L_proc, vocab]
    # sanity: processed string length
    proc = f"<{VH}>|<{VL}>"
    assert logits.shape[0] == len(proc), (logits.shape, len(proc))
    tok_ids = [aa2tok[a] for a in AA20]
    scan = {}                              # pos -> {aa: (logit, prob20)}
    for name, (a, b) in CDRS.items():
        for pos in range(a, b + 1):
            idx = construct_to_idx(pos)
            assert proc[idx] == PARENT[pos-1], (pos, proc[idx], PARENT[pos-1])
            row = logits[idx]
            sub = np.array([row[t] for t in tok_ids], dtype=float)
            p20 = np.exp(sub - sub.max()); p20 = p20 / p20.sum()
            scan[pos] = {aa: (float(row[aa2tok[aa]]), float(p20[i]))
                         for i, aa in enumerate(AA20)}
    return scan

# ---------------- ESM scan ----------------
def esm_scan():
    prob = {}          # (pos, aa) -> p
    parent_p = {}      # pos -> p(parent aa)
    for fn in glob.glob(str(HERE / "esm_scan/out/chunk*/*.tsv")):
        with open(fn) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                lab = row["label"]; aa = row["aa"]
                if aa not in AA20:
                    continue
                pos = int(lab.split("_")[0][1:]); par = lab.split("_")[2]
                p = float(row["prob"])
                prob[(pos, aa)] = p
                if aa == par:
                    parent_p[pos] = p
    return prob, parent_p

# ---------------- parse model picks from variants.md ----------------
def parse_picks():
    md = (HERE / "variants.md").read_text()
    picks = {}
    for vid in [f"V{n:02d}" for n in range(6, 16)]:
        m = re.search(rf"### {vid}\b.*?\*\*Mutations \(\d+\):\*\*\s*([^\n]+)", md, re.S)
        muts = [x.strip() for x in m.group(1).split(",")]
        picks[vid] = [mu for mu in muts if mu not in ("D31E", "D102W")]
    return picks

def parse(mut):
    return mut[0], int(mut[1:-1]), mut[-1]

# ---------------- classification ----------------
# ESM favours every pick by construction (p_mut > p_wt). AbLang2 verdict drives
# the class, via the masked-marginal log-likelihood ratio llr = logit_mut-logit_wt.
CONS_T = 0.25     # llr >= +0.25 : AbLang2 also favours  -> CONSENSUS
CONF_T = -0.75    # llr <= -0.75 : AbLang2 clearly disfavours -> CONFLICT

def classify(llr):
    if llr >= CONS_T:
        return "CONSENSUS"
    if llr <= CONF_T:
        return "CONFLICT"
    return "NEUTRAL"

def best_alt(ab, esm_prob, esm_par, exclude_pos, prefer_cdr=None):
    """Antibody-preferred alternative for a CONFLICT pick: the substitution
    AbLang2 most favours (highest llr>=CONS_T) at an allowed CDR position that is
    not protected/cys, not already mutated in that variant, and introduces no new
    N/D-[GST] motif. To avoid trading a conflict for the reverse, candidates ESM
    does NOT disfavour (Δ >= 0) are preferred; only if none exist do we fall back
    to AbLang2-only picks. Prefers the conflicting pick's own CDR, then any CDR."""
    search = ([prefer_cdr] + [c for c in CDRS if c != prefer_cdr]) if prefer_cdr else list(CDRS)
    def collect(cdr_list):
        out = []
        for cdr in cdr_list:
            a, b = CDRS[cdr]
            for pos in range(a, b + 1):
                if pos in PROTECTED or pos in CYS or pos in exclude_pos:
                    continue
                wt = PARENT[pos-1]
                wl = ab[pos][wt][0]
                for aa in AA20:
                    if aa == wt or new_hard(pos, aa):
                        continue
                    llr = ab[pos][aa][0] - wl
                    if llr < CONS_T:
                        continue
                    e = esm_prob.get((pos, aa), 0.0) - esm_par.get(pos, 0.0)
                    out.append((cdr, pos, wt, aa, llr, ab[pos][aa][1], e))
        return out
    # try same CDR first, then all CDRs
    for cdr_list in ([[prefer_cdr]] if prefer_cdr else []) + [search]:
        cands = collect(cdr_list)
        if not cands:
            continue
        consensus = [c for c in cands if c[6] >= 0.0]   # ESM-neutral-or-better
        pool = consensus if consensus else cands
        pool.sort(key=lambda x: -x[4])
        return pool[0]
    return None

def main():
    ab = ablang_scan()
    esm_prob, esm_par = esm_scan()
    picks = parse_picks()

    # build rows
    rows = []
    for vid, ms in picks.items():
        # positions fixed in this variant: the two anchors + its 3 picks
        used_pos = {31, 102} | {parse(mu)[1] for mu in ms}
        for mu in ms:
            wt, pos, mut = parse(mu)
            cdr = cdr_of(pos)
            e_mut = esm_prob.get((pos, mut), 0.0)
            e_wt = esm_par.get(pos, 0.0)
            e_delta = e_mut - e_wt
            a_mut_logit, a_mut_p = ab[pos][mut]
            a_wt_logit, a_wt_p = ab[pos][wt]
            llr = a_mut_logit - a_wt_logit
            cls = classify(llr)
            alt = None
            if cls == "CONFLICT":
                alt = best_alt(ab, esm_prob, esm_par, exclude_pos=used_pos, prefer_cdr=cdr)
            rows.append(dict(vid=vid, mut=mu, pos=pos, cdr=cdr, wt=wt, mt=mut,
                             e_mut=e_mut, e_delta=e_delta,
                             a_p=a_mut_p, a_wt_p=a_wt_p, llr=llr, cls=cls, alt=alt))

    # ---- write markdown ----
    def alt_str(alt):
        if not alt:
            return "—"
        cdr, pos, wt, aa, llr, ap, e = alt
        return (f"{wt}{pos}{aa} (CDR-{cdr}, AbLang2 p={ap:.2f}, llr={llr:+.2f}; "
                f"ESM Δ={e:+.2f})")

    L = []
    w = L.append
    w("# AbLang2 antibody-specific cross-check of the model-picked CDR substitutions\n")
    w("**Scope.** The 10 model variants **V06-V15** each add 3 ESM-picked CDR "
      "substitutions on top of the fixed D102W + D31E anchors. Here each of those "
      "30 picks is re-scored with the **paired AbLang2** antibody language model "
      "(likelihood mode, masked-marginal / `stepwise_masking=True`) and compared "
      "with the original ESM2-650M masked-marginal scan. This is a **report only** "
      "— no variants were rebuilt.\n")
    w("## Method\n")
    w("- Parent scFv split at the (GGGGS)3 linker into **VH (construct 1-120)** and "
      "**VL (construct 136-242)**; fed to AbLang2 as a paired `[VH, VL]` input. "
      "VH construct position `p` maps to processed-string index `p`; VL position "
      "`p` to index `p-12` (verified against the parent residues).\n")
    w("- **AbLang2 score** = masked-marginal probability `p` of the residue at that "
      "CDR position (softmax over the 20 AAs), plus the log-likelihood ratio "
      "`llr = logit(mut) - logit(parent)`. `llr > 0` means AbLang2 prefers the "
      "mutant over the parent residue; `llr < 0` means it prefers the parent.\n")
    w("- **ESM score** = the mutant's probability `p_mut` from the earlier ESM2 "
      "masked scan, with `Δ = p_mut - p_parent` (positive for every pick, since the "
      "picks were chosen because ESM favoured them).\n")
    w(f"- **Classification** (ESM favours all picks by construction; AbLang2 is the "
      f"independent antibody vote): **CONSENSUS** if AbLang2 `llr >= {CONS_T:+.2f}` "
      f"(both favour), **CONFLICT** if `llr <= {CONF_T:+.2f}` (AbLang2 clearly "
      f"prefers the parent), **NEUTRAL** otherwise (AbLang2 roughly indifferent).\n")
    w("- **Conflict swaps** are the antibody-preferred substitution AbLang2 most "
      "favours (highest `llr`) at an allowed CDR position — never a protected "
      "residue (W99, R50, R59, Y105, H226) or conserved cysteine (C22, C96, C158, "
      "C223), and introducing no new N-[GST]/D-[GST] motif. Preference is given to "
      "the same CDR as the conflicting pick. **Suggestions only — not applied.**\n")

    # summary counts
    from collections import Counter
    c = Counter(r["cls"] for r in rows)
    w(f"## Summary\n")
    w(f"- **{c['CONSENSUS']}** CONSENSUS, **{c['NEUTRAL']}** NEUTRAL, "
      f"**{c['CONFLICT']}** CONFLICT (of {len(rows)} picks).\n")

    w("## All picks (V06-V15)\n")
    w("| Variant | Pick | CDR | ESM p_mut | ESM Δ | AbLang2 p | AbLang2 llr | Class | Suggested swap (conflicts) |")
    w("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        w(f"| {r['vid']} | {r['mut']} | {r['cdr']} | {r['e_mut']:.2f} | "
          f"{r['e_delta']:+.2f} | {r['a_p']:.2f} | {r['llr']:+.2f} | "
          f"**{r['cls']}** | {alt_str(r['alt']) if r['cls']=='CONFLICT' else ''} |")

    # per-variant conflict callout
    conflicts = [r for r in rows if r["cls"] == "CONFLICT"]
    w("\n## Conflicts and proposed antibody-preferred swaps\n")
    if not conflicts:
        w("None — AbLang2 did not clearly disfavour any model pick.\n")
    else:
        for r in conflicts:
            alt = r["alt"]
            w(f"- **{r['vid']} {r['mut']}** (CDR-{r['cdr']}): AbLang2 llr "
              f"{r['llr']:+.2f} (prefers parent {r['wt']}, p_parent={r['a_wt_p']:.2f} "
              f"vs p_mut={r['a_p']:.2f}). ESM still favours it (Δ={r['e_delta']:+.2f}). "
              f"→ antibody-preferred alternative: **{alt_str(alt)}**.")
    w("")
    w("## How to read this\n")
    w("CONSENSUS picks are the safest keeps (both a general PLM and an antibody PLM "
      "like the mutant). NEUTRAL picks are acceptable — AbLang2 neither endorses nor "
      "rejects them; keep them for diversity unless you want to be conservative. "
      "CONFLICT picks are where the antibody model pushes back on the ESM choice; the "
      "suggested swap is the antibody-native move if you choose to revise that variant. "
      "The user decides which (if any) swaps to apply.\n")

    (HERE / "ablang2_crosscheck.md").write_text("\n".join(L))

    # console echo
    print(f"{'Var':<4}{'pick':<7}{'cdr':<4}{'ESMp':>6}{'ESMd':>7}{'ABp':>6}{'llr':>7}  class")
    for r in rows:
        print(f"{r['vid']:<4}{r['mut']:<7}{r['cdr']:<4}{r['e_mut']:>6.2f}"
              f"{r['e_delta']:>+7.2f}{r['a_p']:>6.2f}{r['llr']:>+7.2f}  {r['cls']}"
              + (f"  -> {alt_str(r['alt'])}" if r['cls']=='CONFLICT' else ""))
    print(f"\n{c['CONSENSUS']} CONSENSUS / {c['NEUTRAL']} NEUTRAL / {c['CONFLICT']} CONFLICT")
    print("wrote ablang2_crosscheck.md")

if __name__ == "__main__":
    main()
