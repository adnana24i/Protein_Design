"""Evaluation pipeline: run every model's scan through the same checks and compare them (README step 5).

WHY THIS EXISTS
    We have no measured affinity data for variants until the lab results arrive (about 3 weeks after the
    hackathon), so a model's real predictive performance cannot be measured today. This pipeline therefore
    has THREE layers, from "available now" to "available after the lab":

    Layer 1  INTRINSIC SANITY CHECKS (no labels needed). Each model is tested against things we already know:
             * disulfide cysteines (C22/C96/C158/C223) must rank among the most intolerant quarter of positions
               (a model that lets them mutate freely is broken or mis-indexed). Lenient on purpose.
             * how often the model keeps the parent residue as its top choice (wild-type recovery)
             * does the model see HER2-contact CDR residues as less tolerant than non-contact CDR residues?
               (AUC < 0.5 means it rates contacts as MORE tolerant, i.e. it is blind to binding.)
             * suspect-wt positions (model gives the parent residue p < 0.01, see esm2_scan.py)
    Layer 2  CROSS-MODEL AGREEMENT. Spearman between models (mutation level and position level), overlap
             of their top-100 mutations, and a vote table: which mutations several independent models
             rank in their top 20%. Trivial baselines (BLOSUM62, random) are included as references:
             a real model that does not beat BLOSUM62 on a check adds nothing for that check.
    Layer 3  EXTRINSIC PERFORMANCE (needs variants; labels optional). Given a CSV of variants it
             (a) validates the hackathon rules, (b) scores each variant under every model by summing the
             per-mutation scores (additive, ignores epistasis), (c) when measured KD is supplied, computes
             Spearman(score, -log10 KD), AUROC for "better than parent" and compares with the trivial
             predictor "fewer mutations is better". The same code works for any labelled set, e.g. a public
             trastuzumab mutational scan (use --no-rules) or our own 15 variants after the lab returns.

ADDING A NEW MODEL (AbLang2, DSM, ProFam-1, ...)
    Write data/scores/<name>_scan.csv with columns scfv_index, mut (1-letter) and either `score`
    or `llr_masked` (higher = more favourable). Optional: suspect_wt (bool). Nothing else changes: the
    next run of this script picks it up automatically. Names ending in `_split` and the baselines
    are excluded from the consensus vote.

USAGE
    python3 -I scripts/evaluate.py                                     # layers 1 and 2
    python3 -I scripts/evaluate.py --variants results/final_15.csv     # + layer 3 (rules and pool scores)
    python3 -I scripts/evaluate.py --variants labelled.csv --parent-kd 0.82   # + performance vs measured KD
    variants CSV: id, and either `sequence` or `mutations` (e.g. "Y33A;N55S", scFv indexing), optional `kd_nM`.

OUTPUTS  results/eval/   metrics_intrinsic.csv, spearman_mutations.csv, spearman_positions.csv,
                         overlap_top100.csv, consensus.csv, [pool_scores.csv, metrics_extrinsic.csv],
                         eval_summary.png, report.md
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from validate import DISULFIDE_CYS, PARENT, validate_pool  # noqa: E402

BASELINES = {"blosum62", "random"}
# Measured-affinity reference, NOT a model: it is ground truth to validate against, and it covers only the
# 190 mutations of the Mason CDR-H3 window, so including it in the inner join would shrink the comparison
# set from ~4,300 mutants to 190. Loaded separately by load_reference().
REFERENCE = "mason_empirical"


# ----------------------------------------------------------------------------- loading
def load_models(scores_dir, ann, only=None):
    """Wide table (index scfv_index, mut) with one score column per model, restricted to VH/VL positions."""
    cols, meta = {}, {}
    for path in sorted(glob.glob(os.path.join(scores_dir, "*_scan.csv"))):
        name = os.path.basename(path)[:-len("_scan.csv")]
        if name == REFERENCE or (only and name not in only):
            continue
        df = pd.read_csv(path)
        col = "score" if "score" in df else "llr_masked"
        df = df[df.scfv_index.isin(ann.index)]
        cols[name] = df.set_index(["scfv_index", "mut"])[col].rename(name)
        sus = df["suspect_wt"].astype(bool) if "suspect_wt" in df else pd.Series(False, index=df.index)
        meta[name] = {"is_llr": col == "llr_masked", "baseline": name in BASELINES,
                      "suspect_pos": set(df.loc[sus, "scfv_index"])}
    if not cols:
        sys.exit(f"no *_scan.csv found in {scores_dir}")
    wide = pd.concat(cols.values(), axis=1, join="inner")
    return wide, meta


def load_reference(scores_dir):
    """Measured per-mutation affinity effects (Mason CDR-H3 window), or None if not generated yet."""
    path = os.path.join(scores_dir, f"{REFERENCE}_scan.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    return df.set_index(["scfv_index", "mut"])[["score", "n_obs"]]


def measured_validation(wide, ref, min_obs=50):
    """Spearman of each model against MEASURED binding effects, on the mutations they share.

    Restricted to mutations observed in at least `min_obs` library variants, because an effect size
    estimated from a handful of observations is mostly noise.
    """
    keep = ref[ref.n_obs >= min_obs]
    common = wide.index.intersection(keep.index)
    out = {}
    for m in wide:
        v = wide.loc[common, m].values
        rho, p = spearmanr(v, keep.loc[common, "score"].values)
        out[m] = {"measured_rho": rho, "measured_p": p, "measured_n": len(common)}
    return pd.DataFrame(out).T


def annotate(wide, ann):
    pos = wide.index.get_level_values(0)
    a = ann.loc[pos]
    return pd.DataFrame({"wt": a.aa.values, "region": a.region.values, "contact": a["contact_4.5A"].values,
                         "dist": a.min_dist_to_her2.values}, index=wide.index)


# ----------------------------------------------------------------------------- layer 1
def intrinsic(wide, meta, ann, info):
    rows = []
    cdr_pos = ann.index[ann.region != "FR"]
    for m in wide:
        s = wide[m]
        pos = s.index.get_level_values(0)
        cys = s[pos.isin(DISULFIDE_CYS)]
        pm = s.groupby(level=0).mean()
        # percentile of the 4 Cys positions among all positions by mean score (0 = most intolerant position)
        cys_pct = float(pm.rank(pct=True).loc[sorted(DISULFIDE_CYS)].mean())
        cdr = pm.loc[pm.index.isin(cdr_pos)]
        non_contact = ~ann.loc[cdr.index, "contact_4.5A"].values
        auc = roc_auc_score(non_contact, cdr.values) if 0 < non_contact.sum() < len(non_contact) else np.nan
        fav = s > 0
        rows.append({
            "model": m,
            "n_mutants": len(s),
            "cys_mean": cys.mean(),
            "cys_pos_pctile": cys_pct,                               # expect < 0.25
            "wt_top1": (s.groupby(level=0).max() < 0).mean() if meta[m]["is_llr"] else np.nan,
            "frac_fav_FR": fav[info.region.values == "FR"].mean(),   # score > 0
            "frac_fav_CDR": fav[info.region.values != "FR"].mean(),
            "contact_sens_auc": auc,                                 # > 0.5: contacts rated less tolerant
            "n_suspect_pos": len(meta[m]["suspect_pos"]),
            "check_cys": "PASS" if cys_pct < 0.25 else "FAIL",
        })
    return pd.DataFrame(rows).set_index("model")


# ----------------------------------------------------------------------------- layer 2
def agreement(wide, k=100):
    rho_m = wide.corr(method="spearman")
    rho_p = wide.groupby(level=0).mean().corr(method="spearman")
    tops = {m: set(wide[m].nlargest(k).index) for m in wide}
    jac = pd.DataFrame({a: {b: len(tops[a] & tops[b]) / len(tops[a] | tops[b]) for b in wide} for a in wide})
    return rho_m, rho_p, jac


def consensus(wide, meta, info, top_frac=0.2):
    models = [m for m in wide if not meta[m]["baseline"] and not m.endswith("_split")]
    if not models:
        return pd.DataFrame(), models
    pct = wide[models].rank(pct=True, ascending=False)
    vote = pct <= top_frac
    vote &= wide[models] > 0
    out = wide[models].add_prefix("score_").copy()
    out["votes"] = vote.sum(axis=1)
    out["n_models"] = len(models)
    out["mean_rank_pct"] = pct.mean(axis=1)
    out = out.join(info)
    pos = out.index.get_level_values(0)
    out["suspect"] = [any(p in meta[m]["suspect_pos"] for m in models) for p in pos]
    out["eligible"] = (~out.contact & (out.dist > 4.5) & (out.index.get_level_values(1) != "C")
                       & ~pos.isin(DISULFIDE_CYS) & ~out.suspect)
    out = out.reset_index().sort_values(["eligible", "votes", "mean_rank_pct"], ascending=[False, False, True])
    out.insert(2, "mutation", out.wt + out.scfv_index.astype(str) + out.mut)
    return out, models


# ----------------------------------------------------------------------------- layer 3
def apply_mutations(mutstr):
    seq = list(PARENT)
    for m in [x for x in str(mutstr).split(";") if x]:
        wt, i, mut = m[0], int(m[1:-1]), m[-1]
        if PARENT[i - 1] != wt:
            sys.exit(f"mutation {m}: parent has {PARENT[i - 1]} at {i}")
        seq[i - 1] = mut
    return "".join(seq)


def load_variants(path):
    df = pd.read_csv(path)
    if "sequence" not in df:
        df["sequence"] = df["mutations"].map(apply_mutations)
    return df


def score_variants(pool, wide):
    out = {}
    for m in wide:
        sc = []
        for muts in pool.mutations:
            keys = [(int(x[1:-1]), x[-1]) for x in muts.split(";") if x]
            if not keys:
                sc.append(np.nan)
                continue
            vals = [wide[m].get(k, np.nan) for k in keys]
            sc.append(np.nan if np.isnan(vals).any() else float(np.sum(vals)))
        out[m] = sc
    return pd.DataFrame(out, index=pool.index)


def extrinsic(pool, scores, parent_kd):
    y = -np.log10(pool.kd_nM.values)
    improved = (pool.kd_nM.values < parent_kd) if parent_kd else None
    rows = []
    preds = scores.copy()
    preds["n_mutations_(fewer=better)"] = -pool.n_mut.values
    for m in preds:
        ok = ~preds[m].isna().values
        if ok.sum() < 5:
            continue
        rho, p = spearmanr(preds[m].values[ok], y[ok])
        auc = np.nan
        if improved is not None and 0 < improved[ok].sum() < ok.sum():
            auc = roc_auc_score(improved[ok], preds[m].values[ok])
        rows.append({"predictor": m, "n": int(ok.sum()), "spearman_vs_-log10KD": rho, "p_value": p,
                     "auroc_better_than_parent": auc})
    return pd.DataFrame(rows).set_index("predictor")


# ----------------------------------------------------------------------------- reporting
def figure(rho_m, intr, path):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    im = ax[0].imshow(rho_m.values, vmin=-1, vmax=1, cmap="RdBu_r")
    ax[0].set_xticks(range(len(rho_m)), rho_m.columns, rotation=45, ha="right")
    ax[0].set_yticks(range(len(rho_m)), rho_m.index)
    for i in range(len(rho_m)):
        for j in range(len(rho_m)):
            ax[0].text(j, i, f"{rho_m.values[i, j]:.2f}", ha="center", va="center", fontsize=8)
    ax[0].set_title("Spearman between models (mutation level)")
    fig.colorbar(im, ax=ax[0], fraction=0.046)
    ax[1].barh(intr.index, intr.cys_pos_pctile, color=["tab:green" if v < 0.25 else "tab:red" for v in intr.cys_pos_pctile])
    ax[1].axvline(0.25, color="k", ls="--", lw=1)
    ax[1].set_xlabel("Cys-position tolerance percentile (want < 0.25)")
    ax[1].set_title("Sanity check: disulfide cysteines")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores-dir", default=os.path.join(ROOT, "data", "scores"))
    ap.add_argument("--annotations", default=os.path.join(ROOT, "data", "interface_contacts.csv"))
    ap.add_argument("--models", nargs="*", help="restrict to these model names")
    ap.add_argument("--variants", help="CSV: id, sequence|mutations, [kd_nM]")
    ap.add_argument("--parent-kd", type=float, default=0.82, help="parent KD in nM (deck: 0.82)")
    ap.add_argument("--no-rules", action="store_true", help="skip hackathon rule validation (e.g. single-mutant sets)")
    ap.add_argument("--top-frac", type=float, default=0.2, help="consensus: a model votes for its top fraction")
    ap.add_argument("--outdir", default=os.path.join(ROOT, "results", "eval"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ann = pd.read_csv(args.annotations).set_index("scfv_index").sort_index()
    wide, meta = load_models(args.scores_dir, ann, args.models)
    info = annotate(wide, ann)
    print(f"models: {list(wide)} | mutants in common: {len(wide)}")

    ref = load_reference(args.scores_dir)
    intr = intrinsic(wide, meta, ann, info)
    if ref is not None:
        intr = intr.join(measured_validation(wide, ref))
    rho_m, rho_p, jac = agreement(wide)
    cons, cmodels = consensus(wide, meta, info, args.top_frac)
    intr.to_csv(os.path.join(args.outdir, "metrics_intrinsic.csv"), float_format="%.4f")
    rho_m.to_csv(os.path.join(args.outdir, "spearman_mutations.csv"), float_format="%.4f")
    rho_p.to_csv(os.path.join(args.outdir, "spearman_positions.csv"), float_format="%.4f")
    jac.to_csv(os.path.join(args.outdir, "overlap_top100.csv"), float_format="%.4f")
    if len(cons):
        cons.to_csv(os.path.join(args.outdir, "consensus.csv"), index=False, float_format="%.4f")
    figure(rho_m, intr, os.path.join(args.outdir, "eval_summary.png"))

    rep = ["# Evaluation report", "",
           f"Models: {', '.join(wide)}; {len(wide)} single mutants (VH and VL positions) in common.", "",
           "## Layer 1. Intrinsic sanity checks", "", "```", intr.round(3).to_string(), "```", "",
           "- `cys_pos_pctile` < 0.25 (`check_cys` PASS): the four disulfide cysteines rank in the most intolerant quarter of positions (0 = most intolerant, random is about 0.5). The criterion is deliberately lenient: it is a bug detector, not a quality score.",
           "- `wt_top1`: fraction of positions where the parent residue is the model's top choice (log-ratio models only).",
           "- `contact_sens_auc` > 0.5: among CDR positions the model rates HER2-contact residues as less tolerant. "
           "About 0.5 or below: the model is blind to binding. A crystal contact is not proof of mutational "
           "sensitivity, so read this as descriptive.", "",
           *( [] if ref is None else [
               "### Layer 1b. Validation against MEASURED binding",
               "",
               "`measured_rho` is the Spearman correlation between the model's score and the measured "
               "per-mutation binding effect from the Mason CDR-H3 library (scFv 99-108), over mutations "
               "observed at least 50 times. This replaces the crystal-contact proxy with real data, but "
               "covers only the CDR-H3 window — it says nothing about H1, H2, framework or VL. "
               "See `results/mason/report.md` for the variant-level AUROC, which is the stronger test.",
               "", "```",
               intr[["measured_rho", "measured_p", "measured_n"]].round(3).to_string(), "```", ""]),
           "## Layer 2. Agreement between models", "", "Spearman, mutation level:", "```", rho_m.round(2).to_string(), "```", "",
           "Spearman, position level (mean score per position):", "```", rho_p.round(2).to_string(), "```", "",
           "Jaccard overlap of top-100 mutations:", "```", jac.round(2).to_string(), "```", ""]
    if len(cons):
        el = cons[cons.eligible]
        rep += [f"Consensus vote ({', '.join(cmodels)}; top {args.top_frac:.0%} and score > 0). "
                f"Eligible = not HER2-contact, > 4.5 A from HER2, not Cys, not suspect-wt: {len(el)} mutations; "
                f"votes distribution: {el.votes.value_counts().sort_index().to_dict()}.", "", "```",
                el.head(15)[["mutation", "region", "dist", "votes", "mean_rank_pct"]].round(3).to_string(index=False), "```", ""]
        if len(cmodels) < 2:
            rep += ["**Only one model is in the vote, so `votes` carries no consensus information yet.** Add more models.", ""]
    if len(cmodels) and set(wide.columns) - set(cmodels) - BASELINES:
        rep += ["Note: `_split` variants are robustness checks of the same model and are not counted as separate votes.", ""]

    if args.variants:
        pool = load_variants(args.variants)
        val = validate_pool(pool[["id", "sequence"]]) if not args.no_rules else None
        pool["mutations"] = [";".join(f"{w}{i}{m}" for i, (w, m) in enumerate(zip(PARENT, s), 1) if w != m)
                             if len(s) == len(PARENT) else "" for s in pool.sequence]
        pool["n_mut"] = [len(x.split(";")) if x else 0 for x in pool.mutations]
        scores = score_variants(pool, wide)
        table = pd.concat([pool[["id", "n_mut"]], scores.add_prefix("sum_")], axis=1)
        if val is not None:
            table.insert(2, "valid", val.valid.values)
            table.insert(3, "errors", val.errors.values)
        pos_sets = [set(int(x[1:-1]) for x in s.split(";") if x) for s in pool.mutations]
        mut_sets = [set(s.split(";")) for s in pool.mutations]
        jacc = [len(a & b) / max(1, len(a | b)) for i, a in enumerate(mut_sets) for b in mut_sets[i + 1:]]
        rep += ["## Layer 3. Variant pool", "",
                f"{len(pool)} variants; " + (f"{int(val.valid.sum())} pass the hackathon rules; " if val is not None else "rules not checked; ")
                + f"{len(set().union(*pos_sets))} distinct positions mutated; mean pairwise mutation overlap "
                  f"(Jaccard) {np.mean(jacc):.2f}." if jacc else "", "", "```", table.round(2).to_string(index=False), "```", ""]
        table.to_csv(os.path.join(args.outdir, "pool_scores.csv"), index=False, float_format="%.4f")
        if "kd_nM" in pool and pool.kd_nM.notna().sum() >= 5:
            pool_l = pool.dropna(subset=["kd_nM"])
            ext = extrinsic(pool_l, scores.loc[pool_l.index], args.parent_kd)
            ext.to_csv(os.path.join(args.outdir, "metrics_extrinsic.csv"), float_format="%.4f")
            rep += ["### Performance against measured KD", "", "```", ext.round(3).to_string(), "```", "",
                    f"Best measured KD {pool_l.kd_nM.min():.2f} nM; {int((pool_l.kd_nM < args.parent_kd).sum())}/{len(pool_l)} "
                    f"better than the parent ({args.parent_kd} nM). With n this small the Spearman has wide uncertainty, "
                    "see p_value; compare models only through the gap to the `n_mutations` and BLOSUM62 baselines.", ""]
        elif "kd_nM" in pool:
            rep += ["_Fewer than 5 labelled variants: no performance metrics computed._", ""]
        else:
            rep += ["_No `kd_nM` column: variants were validated and scored but no performance can be measured yet._", ""]

    with open(os.path.join(args.outdir, "report.md"), "w") as f:
        f.write("\n".join(rep))
    print("\n".join(rep))
    print(f"\nwrote {args.outdir}/ (report.md, eval_summary.png, csv tables)")


if __name__ == "__main__":
    main()
