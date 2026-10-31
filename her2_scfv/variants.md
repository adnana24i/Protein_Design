# Trastuzumab scFv — 15 affinity-matured variants

**Adaptyv × HerCode breast-cancer hackathon — HER2 scFv affinity maturation.**

## Strategy

This is **affinity maturation of an existing binder** (the trastuzumab scFv), not de-novo design. All 15 variants are built by rewriting a handful of CDR residues on the 242-aa parent and are anchored on the CDR-H3 **D102W** contact (Trp -> HER2 Ile613; Moon 2016) plus the CDR-H1 **D31E** edit (Zhou 2025).

- **V01-V05 (literature-grounded):** fixed sets drawn only from published, experimentally-validated anti-HER2 edits — the AH06 CDR-H3 quartet (D102W + M107F + D108A + Y109L, validated as a unit; Moon 2016), D31E (Zhou 2025), and the two permitted VH-framework stability edits S85N / R87T (Ng 2025).

- **V06-V15 (model-proposed):** each keeps the **D102W + D31E** anchors and adds **3 CDR substitutions** chosen by an **ESM2-650M masked-marginal scan** of every non-protected CDR position (`modal_esm2_predict_masked.py`). Picks are the substitutions the model scores above the parent residue, spread across different CDRs and kept diverse between variants.

## Guardrails (verified in code by `build_variants.py`, not by eye)

- >= 5 mutations from the parent; single chain, order **VH - (GGGGS)3 - VL**, length <= 260 aa (all variants are 242 aa).

- The **(GGGGS)3 linker** (positions 121-135) is never touched.

- Conserved disulfide cysteines **C22, C96 (VH), C158, C223 (VL)** are never mutated; total cysteine count stays even (4).

- Protected contacts never changed: **W99, R50, R59, Y105 (VH); H226 (VL)**.

- Mutations are confined to the six CDRs (H1 26-35, H2 50-66, H3 99-109, L1 159-169, L2 185-191, L3 224-232); the only framework exceptions are the documented S85N / R87T.

- **No new N-[GST] (deamidation) or D-[GST] (isomerization)** motif is created at any site that was clean in the parent. The anchors D31E and D102W actually **remove** the parent's two isomerization motifs (D31-T32, D102-G103), and the ESM picks N55S (V09) and N165S (V08) remove parent deamidation motifs (N55-G56, N165-T166).

## HER2 target — trimmed to subdomain IV

Trastuzumab binds HER2 **subdomain IV**. The co-fold target is therefore trimmed to **UniProt P04626 residues 510-652** (143 aa, `her2_domainIV.fasta`), which contains the trastuzumab epitope (~557-603). Using the full 630-aa ECD would be slower and noisier for the later cofolding step.

## Note on QC level

The parent scFv already scores **HIGH (severity 54)** in the protein-qc sequence-liability scan — the trastuzumab framework carries several intrinsic deamidation/isomerization motifs. Every variant therefore inherits this baseline, but **none exceeds the parent severity** and several reduce it. The relevant hackathon constraint — *no newly introduced* N/D-[GST] liabilities — is met by all 15.

---

### V01 (literature)

**Mutations (5):** D31E, D102W, M107F, D108A, Y109L  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 43; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| M107F | M -> F | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| D108A | D -> A | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| Y109L | Y -> L | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V02 (literature)

**Mutations (5):** D102W, M107F, D108A, Y109L, S85N  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 47; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| M107F | M -> F | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| D108A | D -> A | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| Y109L | Y -> L | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| S85N | S -> N | VH framework | Ng 2025 (mAbs) | VH framework; binding-retained + stability gain (one of the two permitted framework edits) |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), NN @ 83(deamidation), DS @ 61(isomerization), DT @ 30(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V03 (literature)

**Mutations (5):** D102W, M107F, D108A, Y109L, R87T  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 48; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| M107F | M -> F | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| D108A | D -> A | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| Y109L | Y -> L | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| R87T | R -> T | VH framework | Ng 2025 (mAbs) | VH framework; binding-retained + stability gain (one of the two permitted framework edits) |

**N/D-[GST] hard sites:** 8 (parent has 9); removed vs parent: 1. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 30(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V04 (literature)

**Mutations (6):** D102W, M107F, D108A, Y109L, S85N, R87T  

**Mutation count vs parent:** 6 | **Length:** 242 aa | **protein-qc:** HIGH (severity 50; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| M107F | M -> F | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| D108A | D -> A | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| Y109L | Y -> L | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| S85N | S -> N | VH framework | Ng 2025 (mAbs) | VH framework; binding-retained + stability gain (one of the two permitted framework edits) |
| R87T | R -> T | VH framework | Ng 2025 (mAbs) | VH framework; binding-retained + stability gain (one of the two permitted framework edits) |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), NN @ 83(deamidation), DS @ 61(isomerization), DT @ 30(isomerization), DT @ 72(isomerization), DT @ 89(isomerization), NLT @ 84(N-glyc)</sub>


### V05 (literature)

**Mutations (7):** D31E, D102W, M107F, D108A, Y109L, S85N, R87T  

**Mutation count vs parent:** 7 | **Length:** 242 aa | **protein-qc:** HIGH (severity 45; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| M107F | M -> F | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| D108A | D -> A | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| Y109L | Y -> L | CDR-H3 | Moon 2016 (Mol Cells, AH06) | CDR-H3; part of the AH06 combination, validated only together with D102W |
| S85N | S -> N | VH framework | Ng 2025 (mAbs) | VH framework; binding-retained + stability gain (one of the two permitted framework edits) |
| R87T | R -> T | VH framework | Ng 2025 (mAbs) | VH framework; binding-retained + stability gain (one of the two permitted framework edits) |

**N/D-[GST] hard sites:** 6 (parent has 9); removed vs parent: 3. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), NN @ 83(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization), NLT @ 84(N-glyc)</sub>


### V06 (model-proposed (ESM))

**Mutations (5):** D31E, D102W, T32Y, V168L, A106Y  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 43; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| T32Y | T -> Y | CDR-H1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers Y over parent T (p=0.80 vs 0.00) |
| V168L | V -> L | CDR-L1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers L over parent V (p=0.81 vs 0.04) |
| A106Y | A -> Y | CDR-H3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers Y over parent A (p=0.43 vs 0.03) |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V07 (model-proposed (ESM))

**Mutations (5):** D31E, D102W, F27Y, D163S, F188S  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 43; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| F27Y | F -> Y | CDR-H1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers Y over parent F (p=0.85 vs 0.10) |
| D163S | D -> S | CDR-L1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers S over parent D (p=0.54 vs 0.06) |
| F188S | F -> S | CDR-L2 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers S over parent F (p=0.35 vs 0.00) |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V08 (model-proposed (ESM))

**Mutations (5):** D31E, D102W, I29F, N165S, T229Y  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 37; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| I29F | I -> F | CDR-H1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers F over parent I (p=0.67 vs 0.10) |
| N165S | N -> S | CDR-L1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers S over parent N (p=0.48 vs 0.04); removes the parent N165-T166 deamidation motif |
| T229Y | T -> Y | CDR-L3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers Y over parent T (p=0.24 vs 0.05) |

**N/D-[GST] hard sites:** 6 (parent has 9); removed vs parent: 3. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V09 (model-proposed (ESM))

**Mutations (5):** D31E, D102W, N55S, A167Y, Y190Q  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** MODERATE (severity 35; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| N55S | N -> S | CDR-H2 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers S over parent N (p=0.35 vs 0.02); removes the parent N55-G56 deamidation motif |
| A167Y | A -> Y | CDR-L1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers Y over parent A (p=0.36 vs 0.05) |
| Y190Q | Y -> Q | CDR-L2 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers Q over parent Y (p=0.37 vs 0.03) |

**N/D-[GST] hard sites:** 6 (parent has 9); removed vs parent: 3. **No new hard liability introduced.**
<sub>Full protein-qc flags: NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V10 (model-proposed (ESM))

**Mutations (5):** D31E, D102W, T54G, Q225G, F104S  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 43; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| T54G | T -> G | CDR-H2 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers G over parent T (p=0.32 vs 0.03) |
| Q225G | Q -> G | CDR-L3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers G over parent Q (p=0.24 vs 0.04) |
| F104S | F -> S | CDR-H3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers S over parent F (p=0.21 vs 0.02) |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V11 (model-proposed (ESM))

**Mutations (5):** D31E, D102W, K30S, A169S, P231L  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 46; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| K30S | K -> S | CDR-H1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers S over parent K (p=0.24 vs 0.01) |
| A169S | A -> S | CDR-L1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers S over parent A (p=0.37 vs 0.14) |
| P231L | P -> L | CDR-L3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers L over parent P (p=0.21 vs 0.03) |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization), NIS @ 27(N-glyc)</sub>


### V12 (model-proposed (ESM))

**Mutations (5):** D31E, D102W, Y52N, A186S, P230A  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 43; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| Y52N | Y -> N | CDR-H2 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers N over parent Y (p=0.37 vs 0.14) |
| A186S | A -> S | CDR-L2 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers S over parent A (p=0.30 vs 0.09) |
| P230A | P -> A | CDR-L3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers A over parent P (p=0.12 vs 0.03) |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V13 (model-proposed (ESM))

**Mutations (5):** D31E, D102W, N28T, Q162S, T228Y  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 43; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| N28T | N -> T | CDR-H1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers T over parent N (p=0.41 vs 0.19) |
| Q162S | Q -> S | CDR-L1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers S over parent Q (p=0.27 vs 0.12) |
| T228Y | T -> Y | CDR-L3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers Y over parent T (p=0.20 vs 0.04) |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V14 (model-proposed (ESM))

**Mutations (5):** D31E, D102W, M107Y, V164L, T232Y  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 43; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| M107Y | M -> Y | CDR-H3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers Y over parent M (p=0.14 vs 0.02) |
| V164L | V -> L | CDR-L1 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers L over parent V (p=0.28 vs 0.14) |
| T232Y | T -> Y | CDR-L3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers Y over parent T (p=0.19 vs 0.10) |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>


### V15 (model-proposed (ESM))

**Mutations (5):** D31E, D102W, G100S, Y190E, P231F  

**Mutation count vs parent:** 5 | **Length:** 242 aa | **protein-qc:** HIGH (severity 43; parent = 54)

| Mutation | Parent -> Mutant | Region | Source | Reason |
|---|---|---|---|---|
| D31E | D -> E | CDR-H1 | Zhou 2025 (J Immunol) | CDR-H1 Asp->Glu; lowered K_D ~1.5->0.2 nM in a germline-like anti-HER2 Fab; also removes the parent D31-T32 isomerization motif |
| D102W | D -> W | CDR-H3 | Moon 2016 (Mol Cells) | CDR-H3 anchor; Trp makes a new hydrophobic contact with HER2 Ile613; removes the parent D102-G103 isomerization motif |
| G100S | G -> S | CDR-H3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers S over parent G (p=0.17 vs 0.09) |
| Y190E | Y -> E | CDR-L2 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers E over parent Y (p=0.24 vs 0.03) |
| P231F | P -> F | CDR-L3 | ESM/MPNN scan | ESM2-650M masked-marginal scan prefers F over parent P (p=0.14 vs 0.03) |

**N/D-[GST] hard sites:** 7 (parent has 9); removed vs parent: 2. **No new hard liability introduced.**
<sub>Full protein-qc flags: NG @ 54(deamidation), NS @ 83(deamidation), NT @ 76(deamidation), NT @ 164(deamidation), DS @ 61(isomerization), DT @ 72(isomerization), DT @ 89(isomerization)</sub>

## In-silico ranking (ipSAE)

Ranked by **ipSAE_min** (min of the two asymmetric scFv<->HER2 directions; IPSAE 10/10) over the 16 Boltz co-folds; ipTM/pTM/pLDDT are Boltz complex confidence.

| Rank | Variant | ipSAE_min | ipTM | pTM | pLDDT | #mut |
|---|---|---|---|---|---|---|
| 1 | V12 | 0.821 (+0.059 vs parent) | 0.923 | 0.950 | 0.876 | 5 |
| 2 | V11 | 0.786 (+0.025 vs parent) | 0.913 | 0.945 | 0.907 | 5 |
| 3 | **parent** _(reference baseline)_ | 0.762 | 0.909 | 0.932 | 0.910 | 0 |
| 4 | V02 | 0.749 (-0.013 vs parent) | 0.891 | 0.935 | 0.884 | 5 |
| 5 | V04 | 0.742 (-0.019 vs parent) | 0.893 | 0.931 | 0.896 | 6 |
| 6 | V09 | 0.742 (-0.020 vs parent) | 0.901 | 0.928 | 0.904 | 5 |
| 7 | V08 | 0.731 (-0.031 vs parent) | 0.902 | 0.930 | 0.897 | 5 |
| 8 | V13 | 0.725 (-0.037 vs parent) | 0.905 | 0.931 | 0.908 | 5 |
| 9 | V07 | 0.720 (-0.041 vs parent) | 0.886 | 0.931 | 0.902 | 5 |
| 10 | V01 | 0.719 (-0.043 vs parent) | 0.880 | 0.925 | 0.897 | 5 |
| 11 | V10 | 0.708 (-0.054 vs parent) | 0.874 | 0.930 | 0.904 | 5 |
| 12 | V05 | 0.679 (-0.082 vs parent) | 0.861 | 0.915 | 0.888 | 7 |
| 13 | V15 | 0.661 (-0.101 vs parent) | 0.871 | 0.926 | 0.901 | 5 |
| 14 | V03 | 0.628 (-0.134 vs parent) | 0.849 | 0.916 | 0.887 | 5 |
| 15 | V14 | 0.600 (-0.161 vs parent) | 0.837 | 0.910 | 0.895 | 5 |
| 16 | V06 | 0.522 (-0.240 vs parent) | 0.846 | 0.900 | 0.903 | 5 |

> ipSAE_min is an interface-confidence rank, not a Kd — used to filter/rank, not to claim nM affinities. Measured Kd comes from the Adaptyv wet lab.
