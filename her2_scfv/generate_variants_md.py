#!/usr/bin/env python3
"""Generate variants.md from the verified build (imports build_variants)."""
import re
from pathlib import Path
import build_variants as B

HERE = Path(__file__).resolve().parent
P = B.PARENT

def region(pos):
    c = B.cdr_of(pos)
    if c:
        return f"CDR-{c}"
    if pos in (85, 87):
        return "VH framework"
    return "framework"

HARD = ["NG", "NS", "NT", "DG", "DS", "DT"]

def hard_sites(seq):
    return {p for p, _ in B.motif_positions(seq, HARD)}

PARENT_HARD = hard_sites(P)

def mut_line(m):
    wt, pos, mu = B.parse(m)
    reg = region(pos)
    if m in B.LIT_INFO:
        src, reason = B.LIT_INFO[m]
    elif (pos, mu) in B.ESM_INFO:
        pm, pw = B.ESM_INFO[(pos, mu)]
        src = "ESM/MPNN scan"
        reason = (f"ESM2-650M masked-marginal scan prefers {mu} over parent {wt} "
                  f"(p={pm:.2f} vs {pw:.2f})")
        # annotate liability removal
        if wt == "N" and P[pos:pos+1] in ("G", "S", "T"):
            reason += f"; removes the parent N{pos}-{P[pos]}{pos+1} deamidation motif"
    else:
        src, reason = "ESM/MPNN scan", "model-favoured CDR substitution"
    return f"| {wt}{pos}{mu} | {wt} -> {mu} | {reg} | {src} | {reason} |"

def main():
    rows = []
    order = list(B.VARIANTS.keys())
    for vid in order:
        seq = B.apply_muts(B.VARIANTS[vid])
        _, _, errs, (level, total, liab) = B.verify(vid, B.VARIANTS[vid])
        assert not errs, (vid, errs)
        rows.append((vid, seq, level, total, liab))

    pl, pt, pliab = B.qc_scan(P)

    out = []
    w = out.append
    w("# Trastuzumab scFv — 15 affinity-matured variants\n")
    w("**Adaptyv × HerCode breast-cancer hackathon — HER2 scFv affinity maturation.**\n")
    w("## Strategy\n")
    w("This is **affinity maturation of an existing binder** (the trastuzumab scFv), "
      "not de-novo design. All 15 variants are built by rewriting a handful of CDR "
      "residues on the 242-aa parent and are anchored on the CDR-H3 **D102W** contact "
      "(Trp -> HER2 Ile613; Moon 2016) plus the CDR-H1 **D31E** edit (Zhou 2025).\n")
    w("- **V01-V05 (literature-grounded):** fixed sets drawn only from published, "
      "experimentally-validated anti-HER2 edits — the AH06 CDR-H3 quartet "
      "(D102W + M107F + D108A + Y109L, validated as a unit; Moon 2016), D31E "
      "(Zhou 2025), and the two permitted VH-framework stability edits S85N / R87T "
      "(Ng 2025).\n")
    w("- **V06-V15 (model-proposed):** each keeps the **D102W + D31E** anchors and "
      "adds **3 CDR substitutions** chosen by an **ESM2-650M masked-marginal scan** "
      "of every non-protected CDR position (`modal_esm2_predict_masked.py`). Picks are "
      "the substitutions the model scores above the parent residue, spread across "
      "different CDRs and kept diverse between variants.\n")
    w("## Guardrails (verified in code by `build_variants.py`, not by eye)\n")
    w("- >= 5 mutations from the parent; single chain, order **VH - (GGGGS)3 - VL**, "
      "length <= 260 aa (all variants are 242 aa).\n")
    w("- The **(GGGGS)3 linker** (positions 121-135) is never touched.\n")
    w("- Conserved disulfide cysteines **C22, C96 (VH), C158, C223 (VL)** are never "
      "mutated; total cysteine count stays even (4).\n")
    w("- Protected contacts never changed: **W99, R50, R59, Y105 (VH); H226 (VL)**.\n")
    w("- Mutations are confined to the six CDRs "
      "(H1 26-35, H2 50-66, H3 99-109, L1 159-169, L2 185-191, L3 224-232); the only "
      "framework exceptions are the documented S85N / R87T.\n")
    w("- **No new N-[GST] (deamidation) or D-[GST] (isomerization)** motif is created "
      "at any site that was clean in the parent. The anchors D31E and D102W actually "
      "**remove** the parent's two isomerization motifs (D31-T32, D102-G103), and the "
      "ESM picks N55S (V09) and N165S (V08) remove parent deamidation motifs "
      "(N55-G56, N165-T166).\n")
    w("## HER2 target — trimmed to subdomain IV\n")
    w("Trastuzumab binds HER2 **subdomain IV**. The co-fold target is therefore trimmed "
      "to **UniProt P04626 residues 510-652** (143 aa, `her2_domainIV.fasta`), which "
      "contains the trastuzumab epitope (~557-603). Using the full 630-aa ECD would be "
      "slower and noisier for the later cofolding step.\n")
    w("## Note on QC level\n")
    w(f"The parent scFv already scores **{pl} (severity {pt})** in the protein-qc "
      "sequence-liability scan — the trastuzumab framework carries several intrinsic "
      "deamidation/isomerization motifs. Every variant therefore inherits this baseline, "
      "but **none exceeds the parent severity** and several reduce it. The relevant "
      "hackathon constraint — *no newly introduced* N/D-[GST] liabilities — is met by "
      "all 15.\n")
    w("---\n")

    for vid, seq, level, total, liab in rows:
        muts = B.VARIANTS[vid]
        ndiff = sum(1 for i in range(len(seq)) if P[i] != seq[i])
        kind = "literature" if vid in ("V01","V02","V03","V04","V05") else "model-proposed (ESM)"
        w(f"### {vid} ({kind})\n")
        w(f"**Mutations ({len(muts)}):** {', '.join(muts)}  \n")
        w(f"**Mutation count vs parent:** {ndiff} | **Length:** {len(seq)} aa | "
          f"**protein-qc:** {level} (severity {total}; parent = {pt})\n")
        w("| Mutation | Parent -> Mutant | Region | Source | Reason |")
        w("|---|---|---|---|---|")
        for m in muts:
            w(mut_line(m))
        # per-variant liability delta
        removed = PARENT_HARD - hard_sites(seq)
        rem = sorted(hard_sites(seq))
        liab_strs = ", ".join(f"{t}({c})" for c, s, t in
                              sorted(liab, key=lambda x: -x[1]))
        w("")
        w(f"**N/D-[GST] hard sites:** {len(rem)} (parent has {len(PARENT_HARD)}); "
          f"removed vs parent: {len(removed)}. "
          f"**No new hard liability introduced.**")
        w(f"<sub>Full protein-qc flags: {liab_strs}</sub>\n")
        w("")

    (HERE / "variants.md").write_text("\n".join(out))
    print("wrote variants.md")

if __name__ == "__main__":
    main()
