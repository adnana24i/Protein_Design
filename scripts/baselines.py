"""Trivial baseline "models" in the same format as the real ones, so the evaluation pipeline can tell
whether a deep-learning model adds anything beyond a substitution table or noise.

    blosum62 : score = BLOSUM62[wt][mut]. Knows nothing about this protein, only generic substitution
               similarity. Conservative substitutions score higher.
    random   : score ~ N(0, 1), fixed seed. Any metric where a real model is not clearly better than this
               is noise.

Writes data/scores/blosum62_scan.csv and data/scores/random_scan.csv with the standard schema
(scfv_index, wt, mut, mutation, domain, region, min_dist_to_her2, her2_contact, score).

Usage: python3 -I scripts/baselines.py
"""
import os

import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AAS = "ACDEFGHIKLMNPQRSTVWY"

ann = pd.read_csv(os.path.join(ROOT, "data", "interface_contacts.csv")).set_index("scfv_index")
B62 = substitution_matrices.load("BLOSUM62")
rng = np.random.default_rng(0)

rows = []
for idx, r in ann.sort_index().iterrows():
    for mut in AAS:
        if mut == r.aa:
            continue
        rows.append({"scfv_index": idx, "wt": r.aa, "mut": mut, "mutation": f"{r.aa}{idx}{mut}",
                     "domain": r.domain, "region": r.region, "min_dist_to_her2": r.min_dist_to_her2,
                     "her2_contact": bool(r["contact_4.5A"]), "blosum": B62[r.aa][mut]})
df = pd.DataFrame(rows)
out = os.path.join(ROOT, "data", "scores")
os.makedirs(out, exist_ok=True)
df.assign(score=df.blosum).drop(columns="blosum").to_csv(os.path.join(out, "blosum62_scan.csv"), index=False)
df.assign(score=rng.standard_normal(len(df))).drop(columns="blosum").to_csv(os.path.join(out, "random_scan.csv"), index=False)
print(f"wrote blosum62_scan.csv and random_scan.csv ({len(df)} rows each; VH+VL only, linker not mutated)")
