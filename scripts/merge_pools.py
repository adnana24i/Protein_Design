"""Merge the two independent 15-variant pools into one submission, using every evidence source.

WHY A MERGE
    Two pools were built without sharing assumptions:
      her2_scfv/   literature anchors (Moon 2016 AH06, Zhou 2025) + an ESM2-650M scan, with Boltz/ipSAE
      results/final_15.csv   measured yeast-display binding (Mason CDR-H3 library, ~36k variants)
    They agree on 9 mutations, reached independently. Only 15 slots are submitted in total, so this
    script composes one pool from the union of the evidence rather than picking a side.

EVIDENCE USED, STRONGEST FIRST
    1. MEASURED combination binder rate -- the library variants that actually carry a whole H3 core,
       compared against a mutation-count-matched baseline. Real data, but CDR-H3 only and yeast display.
    2. LITERATURE -- experimentally validated anti-HER2 edits (AH06 quartet; D31E).
    3. MODEL consensus votes across AbLang2 / ESM-1v / ESM2-650M / ESM2-150M.
    4. ipSAE -- carried forward for the two variants that beat the parent, but flagged PROVISIONAL:
       those runs were single-seed (--seed 42, one diffusion sample) and the her2_scfv README itself
       states a small ipSAE gap is not a real affinity difference. Must be re-run with replicates.

TIERS
    M  measured core        : highest measured enrichment, library-supported combinations
    X  cross-validated      : M107F (literature AH06 + best measured + 4/4 votes) combined with the
                              measured-best D102/D108 substitutions. NEW -- in neither original pool.
    L  literature arm       : AH06-based, keeps the D102W anchor. D102W has ZERO library observations,
                              so it is unverified rather than refuted; this arm preserves that bet.
    I  ipSAE survivors      : V11 / V12, the only two that beat the parent. Provisional.
    D  control              : model-only, no measured mutation. Tests whether the data-driven arms won.

OUTPUT  results/final_15_merged.csv, plus Boltz-ready complex FASTAs and a run script for ipSAE.
"""
import os, sys, itertools
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from validate import PARENT, validate_pool, DISULFIDE_CYS
from mason_affinity import load_mason, START as MSTART

ann = pd.read_csv(os.path.join(ROOT, "data", "interface_contacts.csv")).set_index("scfv_index")
cons = pd.read_csv(os.path.join(ROOT, "results", "eval", "consensus.csv"))
votes = {(r.scfv_index, r.mut): r.votes for r in cons.itertuples()}
lib, _ = load_mason(os.path.join(ROOT, "data", "mason"))
_base = lib.groupby("n_mut").label.mean(); _nm = lib.n_mut.values; _y = lib.label.values

def support(core):
    core = [(i, m) for i, m in core if MSTART <= i <= MSTART + 9]
    if not core: return 0, np.nan, np.nan
    m = np.ones(len(lib), bool)
    for i, mu in core:
        m &= np.array([s[i - MSTART] == mu for s in lib.AASeq])
    n = int(m.sum())
    if not n: return 0, np.nan, np.nan
    exp = float(np.mean([_base[k] for k in _nm[m]]))
    return n, float(_y[m].mean()), float(_y[m].mean()) / exp

def ok(i, m):
    a = ann.loc[i]
    return (not bool(a["contact_4.5A"])) and a.min_dist_to_her2 > 4.5 and m != "C" and i not in DISULFIDE_CYS

LIT = {(31,"E"),(102,"W"),(107,"F"),(108,"A"),(109,"L")}          # Moon 2016 AH06 + Zhou 2025 D31E
NONH3 = [(28,"T"),(79,"V"),(54,"G"),(29,"F"),(49,"G"),(168,"L"),(55,"S")]

# hand-specified cores per tier, each justified in the docstring
TIERS = {
 "M": [[(100,"A"),(101,"A"),(108,"A")], [(100,"A"),(101,"S"),(107,"F")],
       [(101,"A"),(102,"F"),(106,"V")], [(101,"A"),(106,"V"),(107,"Y")],
       [(100,"A"),(101,"A"),(107,"F")]],
 "X": [[(101,"A"),(102,"F"),(107,"F")], [(102,"F"),(107,"F"),(108,"A")],
       [(102,"A"),(107,"F"),(108,"A")], [(107,"F"),(108,"A")]],
 "L": [[(102,"W"),(107,"F"),(108,"A"),(109,"L")],
       [(31,"E"),(102,"W"),(107,"F"),(108,"A"),(109,"L")],
       [(102,"W"),(107,"F"),(108,"A"),(109,"L"),(85,"N")]],
 "I": [[(31,"E"),(102,"W"),(52,"N"),(186,"S"),(230,"A")],      # V12, best ipSAE
       [(31,"E"),(102,"W"),(30,"S"),(169,"S"),(231,"L")]],     # V11
 "D": [[(28,"T"),(29,"F"),(54,"G"),(55,"S"),(79,"V")]],
}

rows, used = [], []
for tier, cores in TIERS.items():
    for core in cores:
        muts = list(core)
        # top up to >=5 mutations with model-consensus changes outside the measured window
        for cand in NONH3:
            if len(muts) >= 5: break
            if cand[0] in {i for i, _ in muts}: continue
            if tier == "D" or ok(*cand): muts.append(cand)
        muts = tuple(sorted(muts))
        if muts in used: continue
        used.append(muts)
        n, rate, enr = support(muts)
        seq = list(PARENT)
        for i, m in muts: seq[i - 1] = m
        rows.append({"tier": tier, "muts": muts, "sequence": "".join(seq),
                     "n_lib": n, "binder_rate": rate, "enrichment": enr,
                     "n_literature": len(set(muts) & LIT),
                     "votes": sum(votes.get(k, 0) for k in muts)})

pool = pd.DataFrame(rows).head(15).reset_index(drop=True)
pool.insert(0, "id", [f"M{i:02d}_{t}" for i, t in zip(range(1, len(pool) + 1), pool.tier)])
pool["mutations"] = [";".join(f"{PARENT[i-1]}{i}{m}" for i, m in ms) for ms in pool.muts]
pool["parent_to_mutant"] = [", ".join(f"{PARENT[i-1]}{i}->{m}" for i, m in ms) for ms in pool.muts]
pool["h3_variant"] = [s[98:108] for s in pool.sequence]
pool["n_mut"] = [len(m) for m in pool.muts]
v = validate_pool(pool[["id", "sequence"]])
pool["valid"] = v.valid.values; pool["errors"] = v.errors.values

out = os.path.join(ROOT, "results", "final_15_merged.csv")
pool.drop(columns="muts").to_csv(out, index=False, float_format="%.3f")
print(f"{int(pool.valid.sum())}/{len(pool)} valid -> {out}\n")
print(pool[["id","tier","parent_to_mutant","h3_variant","n_lib","binder_rate",
            "enrichment","n_literature","votes","n_mut","valid"]].to_string(index=False))

# ---- Boltz-ready inputs so ipSAE is one command away ----
her2 = os.path.join(ROOT, "her2_scfv", "her2_domainIV.fasta")
cdir = os.path.join(ROOT, "results", "merged_complexes"); os.makedirs(cdir, exist_ok=True)
if os.path.exists(her2):
    tgt = "".join(l.strip() for l in open(her2) if not l.startswith(">"))
    for r in pool.itertuples():
        with open(os.path.join(cdir, f"{r.id}_complex.fasta"), "w") as f:
            f.write(f">scfv\n{r.sequence}\n>her2\n{tgt}\n")
    with open(os.path.join(cdir, "parent_complex.fasta"), "w") as f:
        f.write(f">scfv\n{PARENT}\n>her2\n{tgt}\n")
    with open(os.path.join(ROOT, "results", "run_ipsae.sh"), "w") as f:
        f.write("""#!/usr/bin/env bash
# Co-fold the merged pool against HER2 domain IV and rank by ipSAE.
# REPLICATES ARE THE POINT: the earlier her2_scfv ranking used a single seed, and gaps of
# 0.02-0.06 around the parent cannot be separated from seed noise. Run >=3 seeds and compare
# the mean +/- spread, not one number.
set -euo pipefail
cd "$(dirname "$0")/.."
for f in results/merged_complexes/*_complex.fasta; do
  id=$(basename "$f" _complex.fasta)
  for seed in 42 43 44; do
    GPU=L40S modal run her2_scfv/patched_biomodals/modal_boltz.py \\
      --input-faa "$f" --params-str "--use_msa_server --seed $seed" \\
      --out-dir "results/ipsae_out/${id}_s${seed}"
  done
done
""")
    os.chmod(os.path.join(ROOT, "results", "run_ipsae.sh"), 0o755)
    print(f"\nwrote {len(pool)+1} Boltz complexes -> {cdir}")
    print("wrote results/run_ipsae.sh (3 seeds per variant)")
