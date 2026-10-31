#!/usr/bin/env python3
"""Build + verify the 15 affinity-matured trastuzumab scFv variants.

Affinity maturation of an EXISTING binder (NOT de-novo). Mutates only within the
six CDRs (construct numbering). Produces:
  - variants.fasta                 (parent + 15 variants)
  - her2_domainIV.fasta            (HER2 subdomain IV, UniProt P04626 510-652)
  - complexes/<ID>_complex.fasta   (2-chain Boltz inputs: >scfv + >her2, 16 files)
  - variants.md                    (full per-mutation documentation)
and prints a summary table. Every hard rule is checked in code; the script
aborts (AssertionError) if any variant fails.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

PARENT = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRF"
          "TISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSGGGGSGGGGSGGGGS"
          "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRS"
          "GTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
assert len(PARENT) == 242, len(PARENT)

# ---- construct layout (1-indexed) -------------------------------------------
VH = (1, 120)
LINKER = (121, 135)            # (GGGGS)3 — never mutate
VL = (136, 242)
assert PARENT[120:135] == "GGGGSGGGGSGGGGS"

CDRS = {"H1": (26, 35), "H2": (50, 66), "H3": (99, 109),
        "L1": (159, 169), "L2": (185, 191), "L3": (224, 232)}

def cdr_of(pos):
    for name, (a, b) in CDRS.items():
        if a <= pos <= b:
            return name
    return None

# protected residues (never change) + conserved disulfide cysteines
PROTECTED = {50: "R", 59: "R", 99: "W", 105: "Y", 226: "H"}
CYS = {22: "C", 96: "C", 158: "C", 223: "C"}
for p, aa in {**PROTECTED, **CYS}.items():
    assert PARENT[p - 1] == aa, (p, PARENT[p - 1], aa)

# ---- mutation sets -----------------------------------------------------------
# A mutation is written "<wt><pos><mut>" with 1-indexed construct position.
# source / reason metadata keyed by (pos, mut) for documentation.
LIT = "literature"
SCAN = "ESM/MPNN scan"

# literature citations/reasons for the fixed edits
LIT_INFO = {
    "D31E":  ("Zhou 2025 (J Immunol)", "CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif"),
    "D102W": ("Moon 2016 (Mol Cells)", "CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif"),
    "M107F": ("Moon 2016 (Mol Cells, AH06)", "CDR-H3; part of the AH06 combination, validated only together with D102W"),
    "D108A": ("Moon 2016 (Mol Cells, AH06)", "CDR-H3; part of the AH06 combination, validated only together with D102W"),
    "Y109L": ("Moon 2016 (Mol Cells, AH06)", "CDR-H3; part of the AH06 combination, validated only together with D102W"),
    "S85N":  ("Ng 2025 (mAbs)", "VH framework; binding-retained + stability gain (one of the two permitted framework edits)"),
    "R87T":  ("Ng 2025 (mAbs)", "VH framework; binding-retained + stability gain (one of the two permitted framework edits)"),
}

# ESM reasons: filled from the scan; delta/p recorded for traceability.
# (pos, mut) -> (p_mut, p_wt) from the masked-marginal ESM2-650M scan
ESM_INFO = {
    (32, "Y"): (0.798, 0.003), (168, "L"): (0.809, 0.035), (106, "Y"): (0.433, 0.026),
    (27, "Y"): (0.855, 0.098), (163, "S"): (0.537, 0.056), (188, "S"): (0.352, 0.003),
    (29, "F"): (0.665, 0.100), (165, "S"): (0.485, 0.035), (229, "Y"): (0.237, 0.052),
    (55, "S"): (0.350, 0.020), (167, "Y"): (0.365, 0.048), (190, "Q"): (0.367, 0.028),
    (54, "G"): (0.319, 0.033), (225, "G"): (0.240, 0.041), (104, "S"): (0.212, 0.017),
    (30, "S"): (0.238, 0.013), (169, "S"): (0.368, 0.139), (231, "L"): (0.212, 0.025),
    (52, "N"): (0.366, 0.135), (186, "S"): (0.301, 0.090), (166, "S"): (0.272, 0.103),
    (28, "T"): (0.407, 0.186), (162, "S"): (0.269, 0.116), (228, "Y"): (0.203, 0.035),
    (107, "Y"): (0.145, 0.017), (164, "L"): (0.282, 0.136), (232, "Y"): (0.190, 0.105),
    (100, "S"): (0.171, 0.086), (190, "E"): (0.236, 0.028), (231, "F"): (0.144, 0.025),
    (230, "A"): (0.119, 0.027),
}

ANCHOR = ["D31E", "D102W"]  # every model variant keeps these

VARIANTS = {
    # ---- 1-5 literature-grounded (exactly these sets) ----
    "V01": ["D31E", "D102W", "M107F", "D108A", "Y109L"],
    "V02": ["D102W", "M107F", "D108A", "Y109L", "S85N"],
    "V03": ["D102W", "M107F", "D108A", "Y109L", "R87T"],
    "V04": ["D102W", "M107F", "D108A", "Y109L", "S85N", "R87T"],
    "V05": ["D31E", "D102W", "M107F", "D108A", "Y109L", "S85N", "R87T"],
    # ---- 6-15 model-proposed: anchor D102W + D31E + 3 ESM-favoured CDR picks ----
    "V06": ANCHOR + ["T32Y", "V168L", "A106Y"],
    "V07": ANCHOR + ["F27Y", "D163S", "F188S"],
    "V08": ANCHOR + ["I29F", "N165S", "T229Y"],
    "V09": ANCHOR + ["N55S", "A167Y", "Y190Q"],
    "V10": ANCHOR + ["T54G", "Q225G", "F104S"],
    "V11": ANCHOR + ["K30S", "A169S", "P231L"],
    "V12": ANCHOR + ["Y52N", "A186S", "P230A"],
    "V13": ANCHOR + ["N28T", "Q162S", "T228Y"],
    "V14": ANCHOR + ["M107Y", "V164L", "T232Y"],
    "V15": ANCHOR + ["G100S", "Y190E", "P231F"],
}

# ---- mutation application / parsing -----------------------------------------
def parse(mut):
    wt, pos, mu = mut[0], int(mut[1:-1]), mut[-1]
    return wt, pos, mu

def apply_muts(muts):
    seq = list(PARENT)
    for m in muts:
        wt, pos, mu = parse(m)
        assert seq[pos - 1] == wt, f"{m}: parent has {seq[pos-1]} at {pos}, not {wt}"
        seq[pos - 1] = mu
    return "".join(seq)

# ---- liability scans (protein-qc expression-qc.md definitions) --------------
DEAMID = ["NG", "NS", "NT", "NN", "NH", "NA", "NQ"]       # N-deamidation
ISOMER = ["DG", "DS", "DT", "DD", "DH", "DN"]             # D-isomerization

def motif_positions(seq, motifs):
    out = set()
    for mo in motifs:
        for mm in re.finditer("(?=" + mo + ")", seq):
            out.add((mm.start(), mo))
    return out

def qc_scan(seq):
    """Full protein-qc sequence-liability scan -> (level, total, details)."""
    liab = []
    cys = seq.count("C")
    if cys % 2 != 0:
        liab.append(("cysteine", 20, f"odd cysteine count ({cys})"))
    for mm in re.finditer(r"[KR]{3,}", seq):
        liab.append(("polybasic", 8 + len(mm.group()), f"{mm.group()} @ {mm.start()}"))
    sev_de = {"NG": 8, "NS": 7, "NT": 6, "NN": 6, "NH": 5, "NA": 4, "NQ": 4}
    for mo, s in sev_de.items():
        for mm in re.finditer("(?=" + mo + ")", seq):
            liab.append(("deamidation", s, f"{mo} @ {mm.start()}"))
    sev_iso = {"DG": 6, "DS": 6, "DT": 5, "DD": 5, "DH": 5, "DN": 5}
    for mo, s in sev_iso.items():
        for mm in re.finditer("(?=" + mo + ")", seq):
            liab.append(("isomerization", s, f"{mo} @ {mm.start()}"))
    for mm in re.finditer(r"N[^P][ST]", seq):
        liab.append(("N-glyc", 3, f"{mm.group()} @ {mm.start()}"))
    # hydrophobic runs >= 6
    run = re.compile(r"[VILMFYW]{6,}")
    for mm in run.finditer(seq):
        liab.append(("hydrophobic", 8 + len(mm.group()), f"{mm.group()} @ {mm.start()}"))
    total = min(sum(s for _, s, _ in liab), 100)
    level = ("LOW" if total <= 15 else "MODERATE" if total <= 35
             else "HIGH" if total <= 60 else "CRITICAL")
    return level, total, liab

# ---- verification of a single variant ---------------------------------------
def verify(vid, muts):
    seq = apply_muts(muts)
    errs = []
    # >=5 mutations vs parent
    diff = [(i + 1, PARENT[i], seq[i]) for i in range(len(seq)) if PARENT[i] != seq[i]]
    if len(diff) < 5:
        errs.append(f"only {len(diff)} mutations (<5)")
    # declared mutations == actual diff
    declared = {(parse(m)[1], parse(m)[2]) for m in muts}
    actual = {(p, mu) for p, _, mu in diff}
    if declared != actual:
        errs.append(f"declared {sorted(declared)} != actual {sorted(actual)}")
    # length
    if len(seq) > 260:
        errs.append(f"length {len(seq)} > 260")
    # format VH-(GGGGS)3-VL, linker intact
    if seq[120:135] != "GGGGSGGGGSGGGGS":
        errs.append("linker altered")
    # every mutated position inside a CDR OR a permitted framework edit
    PERMITTED_FW = {(85, "N"), (87, "T")}
    for p, _, mu in diff:
        if cdr_of(p) is None and (p, mu) not in PERMITTED_FW:
            errs.append(f"mutation at {p} outside CDR/permitted-FW")
    # protected positions + cysteines unchanged
    for p, aa in {**PROTECTED, **CYS}.items():
        if seq[p - 1] != aa:
            errs.append(f"protected/cys {p} changed to {seq[p-1]}")
    # even cysteine count
    if seq.count("C") % 2 != 0:
        errs.append("odd cysteine count")
    # HARD RULE: no new N-[GST] or D-[GST] motif at a site that was clean in the
    # parent. Compared by START POSITION (an already-liable Asn/Asp is not "new").
    HARD = ["NG", "NS", "NT", "DG", "DS", "DT"]
    base = {p for p, _ in motif_positions(PARENT, HARD)}
    new = {(p, mo) for p, mo in motif_positions(seq, HARD) if p not in base}
    if new:
        errs.append(f"NEW N/D-[GST] motif(s): {sorted(new)}")
    level, total, liab = qc_scan(seq)
    return seq, diff, errs, (level, total, liab)

def main():
    results = {}
    print(f"{'ID':<5}{'#mut':>5}{'len':>5}  {'QC':<9}{'sev':>4}  mutations")
    all_ok = True
    for vid, muts in VARIANTS.items():
        seq, diff, errs, (level, total, liab) = verify(vid, muts)
        results[vid] = dict(seq=seq, muts=muts, diff=diff, errs=errs,
                            level=level, total=total, liab=liab)
        status = "PASS" if not errs else "FAIL:" + ";".join(errs)
        all_ok &= not errs
        print(f"{vid:<5}{len(diff):>5}{len(seq):>5}  {level:<9}{total:>4}  "
              f"{','.join(muts)}  [{status}]")
    # parent QC for reference
    pl, pt, _ = qc_scan(PARENT)
    print(f"\nparent QC: {pl} (sev {pt}); length {len(PARENT)}; cys {PARENT.count('C')}")
    if not all_ok:
        print("\n*** ONE OR MORE VARIANTS FAILED ***")
        sys.exit(1)
    print("\nAll 15 variants PASS every hard rule.")

    # ---------------- write FASTAs ----------------
    # variants.fasta (parent + 15)
    with open(HERE / "variants.fasta", "w") as f:
        f.write(f">parent\n{PARENT}\n")
        for vid in VARIANTS:
            f.write(f">{vid}\n{results[vid]['seq']}\n")

    # HER2 domain IV (UniProt P04626 residues 510-652)
    full = ""
    ff = HERE / "her2_full.fasta"
    for line in open(ff):
        if not line.startswith(">"):
            full += line.strip()
    assert len(full) > 1200, "her2_full.fasta looks wrong"
    her2_div = full[509:652]          # UniProt 510..652 inclusive, 1-indexed
    assert full[556:603]  # epitope ~557-603 present
    with open(HERE / "her2_domainIV.fasta", "w") as f:
        f.write(f">her2_domainIV_P04626_510-652\n{her2_div}\n")
    print(f"HER2 domain IV: {len(her2_div)} aa (P04626 510-652); "
          f"epitope 557-603 -> {full[556:603][:10]}...")

    # complexes/<ID>_complex.fasta for parent + 15
    cdir = HERE / "complexes"
    cdir.mkdir(exist_ok=True)
    def write_complex(name, scfv):
        with open(cdir / f"{name}_complex.fasta", "w") as f:
            f.write(f">scfv\n{scfv}\n>her2\n{her2_div}\n")
    write_complex("parent", PARENT)
    for vid in VARIANTS:
        write_complex(vid, results[vid]["seq"])
    print(f"wrote complexes/ (16 files) and variants.fasta (17 records)")

    # expose for the markdown writer
    import json
    meta = {vid: {
        "seq": results[vid]["seq"],
        "muts": results[vid]["muts"],
        "diff": results[vid]["diff"],
        "level": results[vid]["level"],
        "total": results[vid]["total"],
        "liab": results[vid]["liab"],
    } for vid in VARIANTS}
    meta["_parent"] = {"seq": PARENT, "level": pl, "total": pt}
    (HERE / "_variants_meta.json").write_text(json.dumps(meta, indent=0))

if __name__ == "__main__":
    main()
