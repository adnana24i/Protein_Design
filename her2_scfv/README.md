# HER2 scFv affinity maturation

Adaptyv × HerCode Breast-Cancer Hackathon (Zürich, 31 Oct 2026).
**Affinity maturation** of the trastuzumab scFv against HER2 — we improve an
existing binder by rewriting a few CDR residues. This is **not** de-novo design.

## What's here (deliverables)
- `variants.md` — full per-mutation documentation of all 15 variants (strategy,
  guardrails, per-mutation source + reason, QC result).
- `variants.fasta` — the parent + 15 variant scFvs (16 records).
- `her2_domainIV.fasta` — the co-fold target: HER2 subdomain IV, UniProt P04626
  residues **510-652** (143 aa), containing the trastuzumab epitope (~557-603).
- `complexes/<ID>_complex.fasta` — 2-chain Boltz inputs (`>scfv` + `>her2`
  domain IV) for the parent and all 15 variants (16 files), ready to co-fold.
- `build_variants.py` — builds + **verifies every hard rule in code** (aborts on
  any failure) and writes the FASTAs/complexes.
- `generate_variants_md.py` — renders `variants.md` from the verified build.
- `esm_scan/` — the ESM2-650M masked-marginal scan used to pick the model
  variants' substitutions (inputs + per-position probability TSVs).

## Design summary (15 variants, all 5-7 mutations, all 242 aa)
All anchored on **D102W** (CDR-H3; Trp→HER2 Ile613, Moon 2016) **+ D31E**
(CDR-H1, Zhou 2025).

- **V01-V05 — literature-grounded** fixed sets: AH06 quartet
  (D102W+M107F+D108A+Y109L, validated as a unit; Moon 2016), D31E (Zhou 2025),
  and the two permitted VH-framework edits S85N / R87T (Ng 2025).
- **V06-V15 — model-proposed**: D102W+D31E anchors + 3 ESM-favoured CDR
  substitutions each, spread across CDRs and kept diverse between variants.

## Guardrails (verified, not eyeballed)
≥5 mutations; single chain VH–(GGGGS)3–VL ≤260 aa; linker (121-135) untouched;
cysteines C22/C96/C158/C223 fixed and even (4); protected W99/R50/R59/Y105/H226
unchanged; mutations only in the six CDRs (plus documented S85N/R87T); **no new
N-[GST] / D-[GST]** liability at any site clean in the parent.

## Reproduce
```
python3 build_variants.py          # verify + write variants.fasta, her2_domainIV.fasta, complexes/
python3 generate_variants_md.py    # render variants.md

# ESM scan (already run; image patched to py3.11):
for c in 0 1 2; do
  GPU=A10G modal run ../modal_esm2_predict_masked.py \
    --input-faa esm_scan/masked_chunk${c}.faa --out-dir esm_scan/out --run-name chunk${c}
done
```

## LATER (user step — do not run as part of variant generation)
Co-fold each shortlist complex against HER2 domain IV, then rank with ipSAE.
```
GPU=L40S modal run ../modal_boltz.py --input-faa complexes/V01_complex.fasta \
  --params-str "--use_msa_server --seed 42" --out-dir out/V01
```
ipSAE is a **filter** to discard broken designs, not a fine-ranker — a small
ipSAE gap between two good variants is not a real affinity difference. After
cofolding, pick a diverse final set.

## Legacy seed files
`parent_complex.fasta`, `ah06_d31e_complex.fasta`, `her2_full.fasta`,
`build_inputs.sh`, and `out/parent/` are from the initial seed run and use the
**full 630-aa HER2 ECD**. The current campaign uses the trimmed domain-IV target
in `complexes/` + `her2_domainIV.fasta`.

## Sources
Moon et al. 2016 Mol Cells; Zhou et al. 2025 J Immunol; Ng et al. 2025 mAbs;
Mason et al. 2021 Nat Biomed Eng; Spanov et al. 2023 (deamidation liability).

## Setup note: biomodals Modal image build fix
Modal's default Python (3.15) breaks old pins. `modal_boltz.py` (micromamba
→ 3.11) and `modal_esm2_predict_masked.py` (micromamba → 3.11, torch 2.4.1 cu121)
were patched accordingly — do not revert.
