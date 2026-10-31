# AbLang2 antibody-specific cross-check of the model-picked CDR substitutions

**Scope.** The 10 model variants **V06-V15** each add 3 ESM-picked CDR substitutions on top of the fixed D102W + D31E anchors. Here each of those 30 picks is re-scored with the **paired AbLang2** antibody language model (likelihood mode, masked-marginal / `stepwise_masking=True`) and compared with the original ESM2-650M masked-marginal scan. This is a **report only** — no variants were rebuilt.

## Method

- Parent scFv split at the (GGGGS)3 linker into **VH (construct 1-120)** and **VL (construct 136-242)**; fed to AbLang2 as a paired `[VH, VL]` input. VH construct position `p` maps to processed-string index `p`; VL position `p` to index `p-12` (verified against the parent residues).

- **AbLang2 score** = masked-marginal probability `p` of the residue at that CDR position (softmax over the 20 AAs), plus the log-likelihood ratio `llr = logit(mut) - logit(parent)`. `llr > 0` means AbLang2 prefers the mutant over the parent residue; `llr < 0` means it prefers the parent.

- **ESM score** = the mutant's probability `p_mut` from the earlier ESM2 masked scan, with `Δ = p_mut - p_parent` (positive for every pick, since the picks were chosen because ESM favoured them).

- **Classification** (ESM favours all picks by construction; AbLang2 is the independent antibody vote): **CONSENSUS** if AbLang2 `llr >= +0.25` (both favour), **CONFLICT** if `llr <= -0.75` (AbLang2 clearly prefers the parent), **NEUTRAL** otherwise (AbLang2 roughly indifferent).

- **Conflict swaps** are the antibody-preferred substitution AbLang2 most favours (highest `llr`) at an allowed CDR position — never a protected residue (W99, R50, R59, Y105, H226) or conserved cysteine (C22, C96, C158, C223), and introducing no new N-[GST]/D-[GST] motif. Preference is given to the same CDR as the conflicting pick. **Suggestions only — not applied.**

## Summary

- **14** CONSENSUS, **6** NEUTRAL, **10** CONFLICT (of 30 picks).

## All picks (V06-V15)

| Variant | Pick | CDR | ESM p_mut | ESM Δ | AbLang2 p | AbLang2 llr | Class | Suggested swap (conflicts) |
|---|---|---|---|---|---|---|---|---|
| V06 | T32Y | H1 | 0.80 | +0.79 | 0.49 | +3.81 | **CONSENSUS** |  |
| V06 | V168L | L1 | 0.81 | +0.77 | 0.56 | +2.20 | **CONSENSUS** |  |
| V06 | A106Y | H3 | 0.43 | +0.41 | 0.22 | +1.34 | **CONSENSUS** |  |
| V07 | F27Y | H1 | 0.86 | +0.76 | 0.01 | -4.43 | **CONFLICT** | K30S (CDR-H1, AbLang2 p=0.51, llr=+4.68; ESM Δ=+0.23) |
| V07 | D163S | L1 | 0.54 | +0.48 | 0.14 | -0.19 | **NEUTRAL** |  |
| V07 | F188S | L2 | 0.35 | +0.35 | 0.28 | +4.37 | **CONSENSUS** |  |
| V08 | I29F | H1 | 0.67 | +0.57 | 0.16 | -0.18 | **NEUTRAL** |  |
| V08 | N165S | L1 | 0.49 | +0.45 | 0.29 | +1.33 | **CONSENSUS** |  |
| V08 | T229Y | L3 | 0.24 | +0.19 | 0.27 | +0.21 | **NEUTRAL** |  |
| V09 | N55S | H2 | 0.35 | +0.33 | 0.07 | -0.37 | **NEUTRAL** |  |
| V09 | A167Y | L1 | 0.37 | +0.32 | 0.41 | +2.76 | **CONSENSUS** |  |
| V09 | Y190Q | L2 | 0.37 | +0.34 | 0.50 | +5.28 | **CONSENSUS** |  |
| V10 | T54G | H2 | 0.32 | +0.29 | 0.12 | +1.78 | **CONSENSUS** |  |
| V10 | Q225G | L3 | 0.24 | +0.20 | 0.00 | -5.23 | **CONFLICT** | T228S (CDR-L3, AbLang2 p=0.33, llr=+0.72; ESM Δ=+0.13) |
| V10 | F104S | H3 | 0.21 | +0.20 | 0.15 | +2.53 | **CONSENSUS** |  |
| V11 | K30S | H1 | 0.24 | +0.23 | 0.51 | +4.68 | **CONSENSUS** |  |
| V11 | A169S | L1 | 0.37 | +0.23 | 0.08 | +0.27 | **CONSENSUS** |  |
| V11 | P231L | L3 | 0.21 | +0.19 | 0.07 | -0.06 | **NEUTRAL** |  |
| V12 | Y52N | H2 | 0.37 | +0.23 | 0.09 | -1.42 | **CONFLICT** | T54D (CDR-H2, AbLang2 p=0.35, llr=+2.84; ESM Δ=+0.09) |
| V12 | A186S | L2 | 0.30 | +0.21 | 0.04 | -2.44 | **CONFLICT** | Y190Q (CDR-L2, AbLang2 p=0.50, llr=+5.28; ESM Δ=+0.34) |
| V12 | P230A | L3 | 0.12 | +0.09 | 0.01 | -3.76 | **CONFLICT** | T228S (CDR-L3, AbLang2 p=0.33, llr=+0.72; ESM Δ=+0.13) |
| V13 | N28T | H1 | 0.41 | +0.22 | 0.44 | +1.73 | **CONSENSUS** |  |
| V13 | Q162S | L1 | 0.27 | +0.15 | 0.00 | -7.30 | **CONFLICT** | A167S (CDR-L1, AbLang2 p=0.42, llr=+2.77; ESM Δ=+0.00) |
| V13 | T228Y | L3 | 0.20 | +0.17 | 0.01 | -2.76 | **CONFLICT** | P231W (CDR-L3, AbLang2 p=0.30, llr=+1.35; ESM Δ=-0.01) |
| V14 | M107Y | H3 | 0.15 | +0.13 | 0.05 | +2.63 | **CONSENSUS** |  |
| V14 | V164L | L1 | 0.28 | +0.15 | 0.04 | -0.52 | **NEUTRAL** |  |
| V14 | T232Y | L3 | 0.19 | +0.08 | 0.00 | -6.24 | **CONFLICT** | T228S (CDR-L3, AbLang2 p=0.33, llr=+0.72; ESM Δ=+0.13) |
| V15 | G100S | H3 | 0.17 | +0.08 | 0.07 | -0.81 | **CONFLICT** | M107F (CDR-H3, AbLang2 p=0.32, llr=+4.49; ESM Δ=+0.05) |
| V15 | Y190E | L2 | 0.24 | +0.21 | 0.36 | +4.94 | **CONSENSUS** |  |
| V15 | P231F | L3 | 0.14 | +0.12 | 0.00 | -2.95 | **CONFLICT** | T228S (CDR-L3, AbLang2 p=0.33, llr=+0.72; ESM Δ=+0.13) |

## Conflicts and proposed antibody-preferred swaps

- **V07 F27Y** (CDR-H1): AbLang2 llr -4.43 (prefers parent F, p_parent=0.56 vs p_mut=0.01). ESM still favours it (Δ=+0.76). → antibody-preferred alternative: **K30S (CDR-H1, AbLang2 p=0.51, llr=+4.68; ESM Δ=+0.23)**.
- **V10 Q225G** (CDR-L3): AbLang2 llr -5.23 (prefers parent Q, p_parent=0.50 vs p_mut=0.00). ESM still favours it (Δ=+0.20). → antibody-preferred alternative: **T228S (CDR-L3, AbLang2 p=0.33, llr=+0.72; ESM Δ=+0.13)**.
- **V12 Y52N** (CDR-H2): AbLang2 llr -1.42 (prefers parent Y, p_parent=0.38 vs p_mut=0.09). ESM still favours it (Δ=+0.23). → antibody-preferred alternative: **T54D (CDR-H2, AbLang2 p=0.35, llr=+2.84; ESM Δ=+0.09)**.
- **V12 A186S** (CDR-L2): AbLang2 llr -2.44 (prefers parent A, p_parent=0.50 vs p_mut=0.04). ESM still favours it (Δ=+0.21). → antibody-preferred alternative: **Y190Q (CDR-L2, AbLang2 p=0.50, llr=+5.28; ESM Δ=+0.34)**.
- **V12 P230A** (CDR-L3): AbLang2 llr -3.76 (prefers parent P, p_parent=0.64 vs p_mut=0.01). ESM still favours it (Δ=+0.09). → antibody-preferred alternative: **T228S (CDR-L3, AbLang2 p=0.33, llr=+0.72; ESM Δ=+0.13)**.
- **V13 Q162S** (CDR-L1): AbLang2 llr -7.30 (prefers parent Q, p_parent=0.64 vs p_mut=0.00). ESM still favours it (Δ=+0.15). → antibody-preferred alternative: **A167S (CDR-L1, AbLang2 p=0.42, llr=+2.77; ESM Δ=+0.00)**.
- **V13 T228Y** (CDR-L3): AbLang2 llr -2.76 (prefers parent T, p_parent=0.16 vs p_mut=0.01). ESM still favours it (Δ=+0.17). → antibody-preferred alternative: **P231W (CDR-L3, AbLang2 p=0.30, llr=+1.35; ESM Δ=-0.01)**.
- **V14 T232Y** (CDR-L3): AbLang2 llr -6.24 (prefers parent T, p_parent=0.52 vs p_mut=0.00). ESM still favours it (Δ=+0.08). → antibody-preferred alternative: **T228S (CDR-L3, AbLang2 p=0.33, llr=+0.72; ESM Δ=+0.13)**.
- **V15 G100S** (CDR-H3): AbLang2 llr -0.81 (prefers parent G, p_parent=0.17 vs p_mut=0.07). ESM still favours it (Δ=+0.08). → antibody-preferred alternative: **M107F (CDR-H3, AbLang2 p=0.32, llr=+4.49; ESM Δ=+0.05)**.
- **V15 P231F** (CDR-L3): AbLang2 llr -2.95 (prefers parent P, p_parent=0.08 vs p_mut=0.00). ESM still favours it (Δ=+0.12). → antibody-preferred alternative: **T228S (CDR-L3, AbLang2 p=0.33, llr=+0.72; ESM Δ=+0.13)**.

## How to read this

CONSENSUS picks are the safest keeps (both a general PLM and an antibody PLM like the mutant). NEUTRAL picks are acceptable — AbLang2 neither endorses nor rejects them; keep them for diversity unless you want to be conservative. CONFLICT picks are where the antibody model pushes back on the ESM choice; the suggested swap is the antibody-native move if you choose to revise that variant. The user decides which (if any) swaps to apply.
