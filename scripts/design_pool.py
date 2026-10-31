"""Compose the 15-variant submission pool and rank it with our own evaluation pipeline.

SELECTION LOGIC (in priority order)
    1. MEASURED evidence first, and measured for the COMBINATION, not just the single mutations.
       An earlier version of this script ranked by the sum of single-mutation effects. The library then
       showed that no variant carried the resulting 5-6 mutation H3 combinations at all (0-1 observations),
       so those sums were pure extrapolation across untested epistasis. This version instead scores an H3
       core by the MEASURED BINDER RATE of the library variants that actually carry it, and requires at
       least --min-support of them. Library-wide binder rate is 0.29, so a core at 0.6+ is a real signal.
    2. MODEL consensus second, for positions the measured data cannot reach (H1, H2, framework, VL).
       These are plausibility votes, not affinity evidence -- see results/eval/report.md.
    3. HARD FILTERS, from the deck rules and our QC: at least 5 mutations, VH-(GGGGS)3-VL, <= 260 aa,
       never a HER2 contact residue (<4.5 A), never a disulfide cysteine, never a mutation TO cysteine.
    4. DIVERSITY: no two submitted variants may share more than `--max-overlap` of their mutations.

TIERS  (the pool is deliberately not 15 copies of the single best-scoring design)
    A  best-supported core : the H3 cores with the highest MEASURED combination binder rate, minimally
                             extended to reach 5 mutations. Highest evidence.
    B  supported + model   : a supported H3 core plus more model-backed changes outside the window.
    C  developability      : as B, but carries the N55S/N55G deamidation fix (N55-G56 motif).
    D  model-only control  : NO measured mutations at all. Adaptyv measures all 15, so these slots test
                             whether measured-data-guided design actually beat model-only design.

OUTPUT results/final_15.csv  (id, tier, sequence, mutations, scores) -- feed to scripts/evaluate.py
"""
import argparse, itertools, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from validate import PARENT, validate_pool, DISULFIDE_CYS

ap = argparse.ArgumentParser()
ap.add_argument("--max-overlap", type=float, default=0.60)
ap.add_argument("--min-support", type=int, default=30,
                help="minimum library variants that must carry an H3 core for it to be usable")
ap.add_argument("--out", default=os.path.join(ROOT, "results", "final_15.csv"))
args = ap.parse_args()

ann = pd.read_csv(os.path.join(ROOT, "data", "interface_contacts.csv")).set_index("scfv_index")
cons = pd.read_csv(os.path.join(ROOT, "results", "eval", "consensus.csv"))
emp = pd.read_csv(os.path.join(ROOT, "data", "scores", "mason_empirical_scan.csv"))
emp = emp[emp.n_obs >= 50].copy()
emp["pct"] = emp.score.rank(pct=True)
meas = {(r.scfv_index, r.mut): (r.score, r.pct) for r in emp.itertuples()}

def ok(i, m):
    a = ann.loc[i]
    return (not bool(a["contact_4.5A"])) and a.min_dist_to_her2 > 4.5 and m != "C" and i not in DISULFIDE_CYS

# ---- measured pool: best mutation(s) per H3 position ----
MEASURED = {}
for (i, m), (eff, pct) in meas.items():
    if pct >= 0.70 and ok(i, m):
        MEASURED.setdefault(i, []).append((pct, eff, m))
for i in MEASURED:
    MEASURED[i].sort(reverse=True)

# ---- model pool: consensus mutations OUTSIDE the measured window ----
model = cons[(cons.votes >= 3) & cons.eligible & ~cons.scfv_index.between(99, 108)]
MODEL = {}
for r in model.itertuples():
    if ok(r.scfv_index, r.mut):
        MODEL.setdefault(r.scfv_index, []).append((r.votes, -r.mean_rank_pct, r.mut))
for i in MODEL:
    MODEL[i].sort(reverse=True)

print(f"measured positions: {sorted(MEASURED)}")
print(f"model positions outside window: {sorted(MODEL)}")

def seq_of(muts):
    s = list(PARENT)
    for i, m in muts: s[i - 1] = m
    return "".join(s)

def score(muts):
    me = sum(meas.get(k, (0.0, 0.0))[0] for k in muts)
    mp = [meas[k][1] for k in muts if k in meas]
    vt = 0
    for i, m in muts:
        row = cons[(cons.scfv_index == i) & (cons.mut == m)]
        vt += int(row.votes.iloc[0]) if len(row) else 0
    return me, (float(np.mean(mp)) if mp else np.nan), len(mp), vt

# ---- H3 cores with real combination support in the library ----
sys.path.insert(0, HERE)
from mason_affinity import load_mason, START as MSTART
lib, _ = load_mason(os.path.join(ROOT, "data", "mason"))
lib_rate = lib.label.mean()
libseq = lib.AASeq.values
libY = lib.label.values

_base_by_nmut = lib.groupby("n_mut").label.mean()
_libnm = lib.n_mut.values

def support(core):
    """Measured binder rate of the library variants that actually carry this whole core.

    Returns (n carriers, binders, observed rate, mutation-count-matched expected rate).
    The carriers also hold other random mutations in the window (median 7 in total), so the rate is
    conditional on the core, not the rate of a clean 3-mutation variant. The matched baseline is the
    mean library binder rate at each carrier's own mutation count, which is the fair comparison --
    the raw library average of 0.29 would flatter these cores.
    """
    mask = np.ones(len(lib), bool)
    for i, m in core:
        mask &= np.array([s[i - MSTART] == m for s in libseq])
    n = int(mask.sum())
    if not n:
        return 0, 0, np.nan, np.nan
    exp = float(np.mean([_base_by_nmut[k] for k in _libnm[mask]]))
    return n, int(libY[mask].sum()), float(libY[mask].mean()), exp

cores = []
mpos = sorted(MEASURED)
for k in (2, 3, 4):
    for combo in itertools.combinations(mpos, k):
        for picks in itertools.product(*[MEASURED[i][:3] for i in combo]):
            core = tuple(sorted((i, p[2]) for i, p in zip(combo, picks)))
            n, nb, rate, exp = support(core)
            if n >= args.min_support:
                cores.append({"core": core, "n_support": n, "n_binders": nb,
                              "binder_rate": rate, "matched_base": exp})
cores = sorted(cores, key=lambda c: -c["binder_rate"])
print(f"library-wide binder rate {lib_rate:.3f}; {len(cores)} H3 cores with >= {args.min_support} observations")
print("top cores:")
for c in cores[:8]:
    print("   ", "+".join(f"{PARENT[i-1]}{i}{m}" for i, m in c["core"]),
          f"n={c['n_support']:>4} binders={c['n_binders']:>4} rate={c['binder_rate']:.3f} "
          f"matched_base={c['matched_base']:.3f} enrichment={c['binder_rate']/c['matched_base']:.2f}x")

cands = []
mod_sorted = sorted(MODEL, key=lambda i: -MODEL[i][0][0])
# A: best-supported cores, minimally extended (most measured-anchored). Generated FIRST so that an
# identical mutation tuple is labelled A rather than being deduped into B/C later.
for c in cores[:8]:
    core = c["core"]
    for extra in itertools.combinations(mod_sorted[:5], max(5 - len(core), 1)):
        muts = tuple(sorted(list(core) + [(j, MODEL[j][0][2]) for j in extra]))
        if any(j == 55 for j, _ in muts): continue          # that is tier C by definition
        cands.append(("A", muts, c))
# B / C: supported core plus a wider set of model-backed changes outside the window
for c in cores[:40]:
    core = c["core"]
    need = max(0, 5 - len(core))
    for extra in itertools.combinations(mod_sorted[:7], max(need, 2)):
        muts = tuple(sorted(list(core) + [(j, MODEL[j][0][2]) for j in extra]))
        if len(muts) > 7: continue
        tier = "C" if any(j == 55 for j, _ in muts) else "B"
        cands.append((tier, muts, c))
for combo in itertools.combinations(mod_sorted, 5):
    muts = tuple(sorted((j, MODEL[j][0][2]) for j in combo))
    cands.append(("D", muts, None))

rows = []
seen = set()
for tier, muts, c in cands:
    if muts in seen: continue
    seen.add(muts)
    me, mp, nm, vt = score(muts)
    rows.append({"tier": tier, "muts": muts, "n_mut": len(muts), "sum_measured": me,
                 "mean_measured_pctile": mp, "n_measured": nm, "votes": vt,
                 "n_support": (c or {}).get("n_support", 0),
                 "n_binders": (c or {}).get("n_binders", 0),
                 "binder_rate": (c or {}).get("binder_rate", np.nan),
                 "matched_base": (c or {}).get("matched_base", np.nan)})
df = pd.DataFrame(rows)
print(f"{len(df)} candidate variants generated")

QUOTA = {"A": 5, "B": 4, "C": 3, "D": 3}
chosen = []

def try_fill(tier, n, overlap):
    sub = df[df.tier == tier].sort_values(["binder_rate", "votes"], ascending=False)
    for r in sub.itertuples():
        if len([c for c in chosen if c["tier"] == tier]) >= n:
            break
        s_ = set(r.muts)
        if any(len(s_ & set(c["muts"])) / len(s_ | set(c["muts"])) > overlap for c in chosen):
            continue
        chosen.append({"tier": tier, "muts": r.muts, "sum_measured": r.sum_measured,
                       "mean_measured_pctile": r.mean_measured_pctile,
                       "n_measured": r.n_measured, "votes": r.votes,
                       "n_support": r.n_support, "n_binders": r.n_binders,
                       "binder_rate": r.binder_rate, "matched_base": r.matched_base})

for overlap in (args.max_overlap, 0.70, 0.80, 0.90):
    for tier, n in QUOTA.items():
        try_fill(tier, n, overlap)
    if len(chosen) >= 15:
        break
if len(chosen) < 15:
    print(f"WARNING: only {len(chosen)} variants met the diversity constraint")
print(f"selected {len(chosen)} variants (final overlap cap {overlap:.2f})")

out = []
for n, c in enumerate(sorted(chosen, key=lambda x: (x["tier"], -(x["binder_rate"] if x["binder_rate"] == x["binder_rate"] else -1))), 1):
    v = seq_of(c["muts"])
    rate, base = c["binder_rate"], c["matched_base"]
    out.append({"id": f"TZ{n:02d}_{c['tier']}", "tier": c["tier"], "sequence": v,
                "parent_to_mutant": ", ".join(f"{PARENT[i-1]}{i}->{m}" for i, m in c["muts"]),
                "h3_parent": PARENT[98:108], "h3_variant": v[98:108],
                "mutations": ";".join(f"{PARENT[i-1]}{i}{m}" for i, m in c["muts"]),
                "n_mut": len(c["muts"]), "n_measured": c["n_measured"],
                "sum_measured_effect": round(c["sum_measured"], 3),
                "mean_measured_pctile": round(c["mean_measured_pctile"], 3) if c["n_measured"] else np.nan,
                "model_votes": c["votes"],
                "h3_core_n_library": c["n_support"],
                "h3_core_binders": c["n_binders"],
                "h3_core_binder_rate": round(rate, 3) if rate == rate else np.nan,
                "h3_core_matched_base": round(base, 3) if base == base else np.nan,
                "h3_core_enrichment": round(rate / base, 2) if (rate == rate and base == base) else np.nan})
pool = pd.DataFrame(out)
val = validate_pool(pool[["id", "sequence"]])
pool["valid"] = val.valid.values; pool["errors"] = val.errors.values
pool["length"] = pool.sequence.str.len()
os.makedirs(os.path.dirname(args.out), exist_ok=True)
pool.to_csv(args.out, index=False)
print(f"\n{int(pool.valid.sum())}/{len(pool)} valid; wrote {args.out}")
print(pool[["id","tier","parent_to_mutant","h3_parent","h3_variant","h3_core_binders",
            "h3_core_n_library","h3_core_binder_rate","h3_core_matched_base","h3_core_enrichment",
            "model_votes","valid"]].to_string(index=False))
