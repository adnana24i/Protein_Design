"""Hackathon rule validator and mutation parser for scFv variants (README step 1).

Rules (deck slides 16 and 19):
    * single chain in the order VH, (GGGGS)3, VL
    * at least 5 mutations from the parent
    * up to 260 aa
Extra checks we add: only the 20 standard amino acids, linker intact, no duplicate sequences, and the
variable-domain disulfide cysteines (scFv 22, 96, 158, 223) untouched.

Mutations are counted as substitutions over the 242 aa parent. A variant of a different length (indel) is
reported as "indel" and cannot be scored additively by the evaluation pipeline; how the organisers count
mutations for indels is an open question (README section 8), so it is flagged as a warning.

Library use:   from validate import PARENT, mutations_of, validate_sequence
CLI:           python3 -I scripts/validate.py pool.csv        # CSV with columns id,sequence
"""
import sys

import pandas as pd

VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAED"
      "TAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
LINKER = "GGGGSGGGGSGGGGS"
VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTT"
      "PPTFGQGTKVEIK")
PARENT = VH + LINKER + VL
AAS = set("ACDEFGHIKLMNPQRSTVWY")
DISULFIDE_CYS = {22, 96, 158, 223}
MIN_MUTATIONS, MAX_LEN = 5, 260


def mutations_of(seq, parent=PARENT):
    """List of 'F218E'-style substitutions (1-based scFv index), or None if the length differs."""
    if len(seq) != len(parent):
        return None
    return [f"{w}{i}{m}" for i, (w, m) in enumerate(zip(parent, seq), 1) if w != m]


def validate_sequence(seq):
    """Return (errors, warnings, mutations). A sequence is valid for submission iff errors is empty."""
    errors, warns = [], []
    muts = mutations_of(seq)
    if not set(seq) <= AAS:
        errors.append(f"non-standard residues: {sorted(set(seq) - AAS)}")
    if len(seq) > MAX_LEN:
        errors.append(f"length {len(seq)} > {MAX_LEN}")
    if muts is None:
        warns.append(f"length {len(seq)} != 242 (indel): mutation count and scoring are ambiguous")
        if seq[:len(VH) // 2] != VH[:len(VH) // 2] or LINKER not in seq:
            errors.append("order/linker not VH-(GGGGS)3-VL")
    else:
        if seq[120:135] != LINKER:
            errors.append("linker (GGGGS)3 at 121-135 is not intact")
        if len(muts) < MIN_MUTATIONS:
            errors.append(f"only {len(muts)} mutations (< {MIN_MUTATIONS})")
        hit = [m for m in muts if int(m[1:-1]) in DISULFIDE_CYS]
        if hit:
            errors.append(f"disulfide cysteine mutated: {hit}")
    return errors, warns, muts or []


def validate_pool(df):
    """df has columns id, sequence. Adds columns valid, n_mut, mutations, errors, warnings."""
    rows = []
    dup = df.sequence.duplicated(keep=False)
    for r, d in zip(df.itertuples(), dup):
        e, w, m = validate_sequence(r.sequence)
        if d:
            e.append("duplicate sequence in pool")
        rows.append({"id": r.id, "sequence": r.sequence, "valid": not e, "n_mut": len(m),
                     "mutations": ";".join(m), "errors": " | ".join(e), "warnings": " | ".join(w)})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    out = validate_pool(pd.read_csv(sys.argv[1]))
    print(out.drop(columns="sequence").to_string(index=False))
    print(f"\n{int(out.valid.sum())}/{len(out)} valid")
    sys.exit(0 if out.valid.all() else 1)
