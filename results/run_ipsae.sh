#!/usr/bin/env bash
# Co-fold the merged pool against HER2 domain IV and rank by ipSAE.
# REPLICATES ARE THE POINT: the earlier her2_scfv ranking used a single seed, and gaps of
# 0.02-0.06 around the parent cannot be separated from seed noise. Run >=3 seeds and compare
# the mean +/- spread, not one number.
set -euo pipefail
cd "$(dirname "$0")/.."
for f in results/merged_complexes/*_complex.fasta; do
  id=$(basename "$f" _complex.fasta)
  for seed in 42 43 44; do
    GPU=L40S modal run her2_scfv/patched_biomodals/modal_boltz.py \
      --input-faa "$f" --params-str "--use_msa_server --seed $seed" \
      --out-dir "results/ipsae_out/${id}_s${seed}"
  done
done
