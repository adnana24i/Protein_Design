"""Test our sequence models against MEASURED HER2 binding (Mason et al. trastuzumab CDR-H3 library).

WHY
    Every model in data/scores/ scores "does this residue look normal here?". The evaluation pipeline
    showed they are anti-correlated with the HER2 interface, but that used crystal contact distance as a
    PROXY for "matters for binding". This script replaces the proxy with real binding labels, and then
    builds a model that targets affinity directly instead of naturalness.

THE DATA  (data/mason/, from github.com/dahjan/DMS_opt, Mason et al.)
    A yeast-display library of trastuzumab CDR-H3 variants sorted into antigen-positive / antigen-negative.
    The 10 randomised residues are the scFv positions 99-108 (WT = WGGDGFYAMD), verified by locating the
    constant 'CSR' motif that precedes CDR-H3. That window contains 4 of our 20 HER2-contact residues
    (W99, G103, F104, Y105) and 2 of our consensus picks (D102, A106).

THREE PROBLEMS WITH THE DATA, AND WHAT WE DO ABOUT THEM
    1. ~2,342 sequences appear in BOTH classes (sorting noise; even WT does, 3 positive vs 171 negative
       reads). We resolve each sequence by total read count and record how many were ambiguous;
       --drop-ambiguous removes them instead, as a robustness check.
    2. The library is combinatorial (median 7 mutations), not single mutants. Our per-mutation scores are
       therefore summed over a variant's mutations, which ignores epistasis BY CONSTRUCTION.
    3. **The dominant confound**: binder rate falls monotonically with mutation count (0.83 at 3 mutations
       to 0.07 at 10). "Fewer mutations" alone is a powerful classifier, and summed LLRs correlate with
       mutation count almost by definition. Every headline number here is therefore reported BOTH overall
       and STRATIFIED BY MUTATION COUNT. Only the stratified number answers the real question:
       *among variants with the same number of mutations, can the model tell binders from non-binders?*

WHAT IT PRODUCES
    1. Zero-shot: AUROC per model, overall and per mutation-count stratum, against three baselines
       (mutation count, BLOSUM62, random).
    2. An EMPIRICAL per-mutation affinity score for all 10x19 mutations in the window, estimated by
       logistic regression on one-hot mutations WITH mutation count as a covariate, so the effect sizes
       are confound-adjusted. This is in the same (scfv_index, mut) format as our model scans, so it can
       be compared with them directly and consumed by evaluate.py.
    3. The direct test: Spearman between each model's LLR and that empirical affinity score.
    4. A supervised affinity model (one-hot and, with --embeddings, ESM2 embeddings) trained on a held-out
       split, to show what targeting affinity instead of naturalness actually buys. A one-hot learning
       curve at matched training sizes is included so the embedding model is compared fairly rather than
       against a one-hot model trained on 18x more data.
       NOTE on the embedding features: residue embeddings are MEAN-POOLED over the 10-residue window,
       which discards position. That is a deliberate simplicity choice and it costs accuracy here;
       concatenating the 10 per-residue vectors (10 x 1280) would preserve position and is the obvious
       next thing to try if the embedding route is pursued.
    5. A verdict on our own consensus picks that fall in this window (D102, A106).

OUTPUTS  data/scores/mason_empirical_scan.csv     empirical affinity score per mutation (model-comparable)
         results/mason/zeroshot_auroc.csv         AUROC per model, overall and stratified
         results/mason/model_vs_empirical.csv     Spearman of each model against measured effects
         results/mason/supervised.csv             held-out performance of the affinity-trained models
         results/mason/consensus_verdict.csv      what the data says about our own picks
         results/mason/report.md

USAGE    python3 -I scripts/mason_affinity.py
         python3 -I scripts/mason_affinity.py --embeddings --embed-n 6000    # + ESM2-embedding model
         python3 -I scripts/mason_affinity.py --drop-ambiguous               # robustness check
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from validate import PARENT  # noqa: E402

AAS = "ACDEFGHIKLMNPQRSTVWY"
START = 99                      # scFv index (1-based) of the first randomised CDR-H3 residue
WT = "WGGDGFYAMD"               # parent residues at scFv 99-108
MIN_STRATUM = 200               # minimum variants in a mutation-count stratum to report its AUROC
MIN_OBS = 50                    # minimum library observations for a per-mutation effect to be trusted.
                                # The logistic regression shrinks rarely-seen mutations toward zero, so a
                                # near-zero effect at low n_obs means "no data", NOT "no effect".


# --------------------------------------------------------------------------- data
def load_mason(dirpath, drop_ambiguous=False):
    pos = pd.read_csv(os.path.join(dirpath, "mHER_H3_AgPos.csv"))
    neg = pd.read_csv(os.path.join(dirpath, "mHER_H3_AgNeg.csv"))
    df = pd.concat([pos.assign(label=1), neg.assign(label=0)])[["AASeq", "Count", "label"]]
    assert (df.AASeq.str.len() == 10).all(), "expected 10-residue CDR-H3 variants"
    assert PARENT[START - 1:START - 1 + 10] == WT, "WT window does not match the parent scFv"

    # one row per sequence; a sequence seen in both classes is resolved by total read count
    agg = df.groupby(["AASeq", "label"], as_index=False).Count.sum()
    piv = agg.pivot(index="AASeq", columns="label", values="Count").fillna(0.0)
    piv.columns = ["neg_reads", "pos_reads"]
    n_ambig = int(((piv.neg_reads > 0) & (piv.pos_reads > 0)).sum())
    if drop_ambiguous:
        piv = piv[~((piv.neg_reads > 0) & (piv.pos_reads > 0))]
    out = piv.reset_index()
    out["label"] = (out.pos_reads > out.neg_reads).astype(int)
    out = out[out.pos_reads != out.neg_reads]          # unresolvable ties
    out["n_mut"] = [sum(a != b for a, b in zip(s, WT)) for s in out.AASeq]
    out["mutations"] = [[(START + i, m) for i, (w, m) in enumerate(zip(WT, s)) if w != m] for s in out.AASeq]
    return out.reset_index(drop=True), n_ambig


def load_model_scores(scores_dir):
    """{model: Series indexed by (scfv_index, mut)} restricted to the randomised window."""
    import glob
    out = {}
    for path in sorted(glob.glob(os.path.join(scores_dir, "*_scan.csv"))):
        name = os.path.basename(path)[:-len("_scan.csv")]
        if name.startswith("mason_"):
            continue
        df = pd.read_csv(path)
        col = "score" if "score" in df else "llr_masked"
        df = df[df.scfv_index.between(START, START + 9)]
        out[name] = df.set_index(["scfv_index", "mut"])[col]
    return out


def score_variants(data, series):
    """Additive model score per variant; NaN if any constituent mutation is unscored."""
    vals = []
    for muts in data.mutations:
        s = [series.get(k, np.nan) for k in muts]
        vals.append(np.nan if (not s or np.isnan(s).any()) else float(np.sum(s)))
    return np.array(vals)


# --------------------------------------------------------------------------- analyses
def stratified_auroc(y, score, n_mut):
    """AUROC within each mutation-count stratum, plus a sample-weighted mean across strata."""
    rows = []
    for k in sorted(set(n_mut)):
        m = (n_mut == k) & ~np.isnan(score)
        if m.sum() < MIN_STRATUM or len(set(y[m])) < 2:
            continue
        rows.append({"n_mut": k, "n": int(m.sum()), "auroc": roc_auc_score(y[m], score[m])})
    t = pd.DataFrame(rows)
    mean = np.average(t.auroc, weights=t.n) if len(t) else np.nan
    return t, mean


def empirical_effects(data):
    """Confound-adjusted per-mutation effect sizes.

    Logistic regression: label ~ one-hot(mutation at each position) + n_mut.
    Including n_mut absorbs the 'fewer mutations bind better' confound, so each mutation's coefficient is
    its effect holding the number of mutations fixed. Positive = enriched among binders.
    """
    cols = [(START + i, a) for i in range(10) for a in AAS if a != WT[i]]
    idx = {c: j for j, c in enumerate(cols)}
    X = np.zeros((len(data), len(cols) + 1), dtype=np.float32)
    for r, muts in enumerate(data.mutations):
        for k in muts:
            X[r, idx[k]] = 1.0
    X[:, -1] = data.n_mut.values
    lr = LogisticRegression(max_iter=3000, C=1.0)
    lr.fit(X, data.label.values)
    coef = lr.coef_[0]
    eff = pd.DataFrame({"scfv_index": [c[0] for c in cols], "mut": [c[1] for c in cols],
                        "wt": [WT[c[0] - START] for c in cols], "effect": coef[:len(cols)]})
    # how often each mutation was actually observed, as a confidence proxy
    obs = pd.Series(0, index=pd.MultiIndex.from_tuples(cols))
    for muts in data.mutations:
        for k in muts:
            obs[k] += 1
    eff["n_obs"] = [obs[(r.scfv_index, r.mut)] for r in eff.itertuples()]
    eff["mutation"] = eff.wt + eff.scfv_index.astype(str) + eff.mut
    return eff, float(coef[-1])


def supervised(data, use_embeddings, embed_n, device):
    """Affinity-targeted models, evaluated on a held-out split."""
    tr, te = train_test_split(data, test_size=0.25, random_state=0, stratify=data.label)
    rows = []

    def onehot(d):
        X = np.zeros((len(d), 10 * 20), dtype=np.float32)
        for r, s in enumerate(d.AASeq):
            for i, a in enumerate(s):
                X[r, i * 20 + AAS.index(a)] = 1.0
        return X

    m = LogisticRegression(max_iter=3000).fit(onehot(tr), tr.label)
    rows.append({"model": "supervised: one-hot sequence", "n_train": len(tr), "n_test": len(te),
                 "auroc_test": roc_auc_score(te.label, m.predict_proba(onehot(te))[:, 1])})

    # same, but with mutation count removed as a feature source: trained within a single stratum
    big = data.n_mut.value_counts().idxmax()
    sub = data[data.n_mut == big]
    if len(sub) > 1000 and sub.label.nunique() > 1:
        s_tr, s_te = train_test_split(sub, test_size=0.25, random_state=0, stratify=sub.label)
        m2 = LogisticRegression(max_iter=3000).fit(onehot(s_tr), s_tr.label)
        rows.append({"model": f"supervised: one-hot, within n_mut={big} only",
                     "n_train": len(s_tr), "n_test": len(s_te),
                     "auroc_test": roc_auc_score(s_te.label, m2.predict_proba(onehot(s_te))[:, 1])})

    # learning curve: one-hot at matched training sizes, so the embedding model below is compared fairly
    for n in (500, 2000, 8000):
        sub = data.sample(min(n, len(data)), random_state=0)
        s_tr, s_te = train_test_split(sub, test_size=0.25, random_state=0, stratify=sub.label)
        mm = LogisticRegression(max_iter=3000).fit(onehot(s_tr), s_tr.label)
        rows.append({"model": f"supervised: one-hot, n_total={n}", "n_train": len(s_tr), "n_test": len(s_te),
                     "auroc_test": roc_auc_score(s_te.label, mm.predict_proba(onehot(s_te))[:, 1])})

    if use_embeddings:
        import torch
        from transformers import AutoTokenizer, EsmForMaskedLM
        name = "facebook/esm2_t33_650M_UR50D"
        tok = AutoTokenizer.from_pretrained(name)
        esm = EsmForMaskedLM.from_pretrained(name, output_hidden_states=True).to(device).eval()
        sub = data.sample(min(embed_n, len(data)), random_state=0)
        seqs = [PARENT[:START - 1] + s + PARENT[START - 1 + 10:] for s in sub.AASeq]
        embs = []
        with torch.no_grad():
            for i in range(0, len(seqs), 4):
                enc = tok(seqs[i:i + 4], return_tensors="pt", padding=True).to(device)
                h = esm(**enc).hidden_states[-1]
                # mean-pool the randomised window only (token j+1 is residue j)
                embs.append(h[:, START:START + 10].mean(1).float().cpu().numpy())
                print(f"  embedding {min(i + 4, len(seqs))}/{len(seqs)}", end="\r", flush=True)
        print()
        X = np.concatenate(embs)
        y = sub.label.values
        Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.25, random_state=0, stratify=y)
        m3 = LogisticRegression(max_iter=3000).fit(Xtr, ytr)
        rows.append({"model": "supervised: ESM2 embeddings", "n_train": len(Xtr), "n_test": len(Xte),
                     "auroc_test": roc_auc_score(yte, m3.predict_proba(Xte)[:, 1])})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.path.join(ROOT, "data", "mason"))
    ap.add_argument("--scores-dir", default=os.path.join(ROOT, "data", "scores"))
    ap.add_argument("--outdir", default=os.path.join(ROOT, "results", "mason"))
    ap.add_argument("--drop-ambiguous", action="store_true", help="drop sequences seen in both classes")
    ap.add_argument("--embeddings", action="store_true", help="also train on ESM2 embeddings (slow)")
    ap.add_argument("--embed-n", type=int, default=6000)
    ap.add_argument("--device", default=None)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    data, n_ambig = load_mason(args.data, args.drop_ambiguous)
    models = load_model_scores(args.scores_dir)
    y = data.label.values
    nm = data.n_mut.values
    print(f"{len(data)} usable variants ({y.sum()} binders); {n_ambig} sequences were in both classes "
          f"({'dropped' if args.drop_ambiguous else 'resolved by read count'})")
    print(f"models: {list(models)}")

    # ---- 1. zero-shot ----
    rows, strat_tables = [], {}
    scored = {m: score_variants(data, s) for m, s in models.items()}
    scored["baseline: fewer mutations"] = -nm.astype(float)
    for m, sc in scored.items():
        ok = ~np.isnan(sc)
        t, mean = stratified_auroc(y, sc, nm)
        strat_tables[m] = t
        rows.append({"model": m, "n_scored": int(ok.sum()),
                     "auroc_overall": roc_auc_score(y[ok], sc[ok]),
                     "auroc_stratified": mean})
    zs = pd.DataFrame(rows).sort_values("auroc_stratified", ascending=False).set_index("model")
    zs.to_csv(os.path.join(args.outdir, "zeroshot_auroc.csv"), float_format="%.4f")

    # ---- 2. empirical per-mutation effects ----
    eff, nmut_coef = empirical_effects(data)
    emp = eff.rename(columns={"effect": "score"})[["scfv_index", "wt", "mut", "mutation", "score", "n_obs"]]
    emp.to_csv(os.path.join(args.scores_dir, "mason_empirical_scan.csv"), index=False, float_format="%.4f")

    # ---- 3. model LLR vs measured effect ----
    rows = []
    for m, s in models.items():
        v = np.array([s.get((r.scfv_index, r.mut), np.nan) for r in eff.itertuples()])
        ok = ~np.isnan(v)
        rho, p = spearmanr(v[ok], eff.effect.values[ok])
        rows.append({"model": m, "n_mutations": int(ok.sum()), "spearman_vs_measured": rho, "p_value": p})
    mve = pd.DataFrame(rows).sort_values("spearman_vs_measured", ascending=False).set_index("model")
    mve.to_csv(os.path.join(args.outdir, "model_vs_empirical.csv"), float_format="%.4f")

    # ---- 4. supervised ----
    device = args.device or ("cuda" if __import__("torch").cuda.is_available() else
                             "mps" if __import__("torch").backends.mps.is_available() else "cpu") \
        if args.embeddings else "cpu"
    sup = supervised(data, args.embeddings, args.embed_n, device)
    sup.to_csv(os.path.join(args.outdir, "supervised.csv"), index=False, float_format="%.4f")

    # ---- 5. our own consensus picks that fall in this window ----
    cons_path = os.path.join(ROOT, "results", "eval", "consensus.csv")
    verdict = pd.DataFrame()
    missed = pd.DataFrame()
    if os.path.exists(cons_path):
        cons = pd.read_csv(cons_path)
        cons = cons[cons.scfv_index.between(START, START + 9)]
        verdict = cons.merge(eff[["scfv_index", "mut", "effect", "n_obs"]], on=["scfv_index", "mut"], how="inner")
        # percentile is only meaningful among well-observed mutations
        rel = eff[eff.n_obs >= MIN_OBS].copy()
        rel["pct"] = rel.effect.rank(pct=True)
        pct = dict(zip(zip(rel.scfv_index, rel.mut), rel.pct))
        verdict["reliable"] = verdict.n_obs >= MIN_OBS
        verdict["measured_pctile"] = [pct.get((r.scfv_index, r.mut), np.nan) for r in verdict.itertuples()]
        verdict = verdict[["mutation", "region", "dist", "votes", "eligible", "effect",
                           "measured_pctile", "n_obs", "reliable"]]
        verdict.sort_values(["reliable", "measured_pctile"], ascending=False).to_csv(
            os.path.join(args.outdir, "consensus_verdict.csv"), index=False, float_format="%.4f")
        # mutations the measured data likes that our FILTERS threw away -- the filters' false negatives
        missed = verdict[verdict.reliable & ~verdict.eligible & (verdict.measured_pctile >= 0.9)]
        missed.to_csv(os.path.join(args.outdir, "filtered_out_but_good.csv"), index=False, float_format="%.4f")

    # ---- report ----
    R = [f"# Measured-affinity benchmark (Mason et al. CDR-H3, scFv {START}-{START + 9})", "",
         f"{len(data)} variants, {int(y.sum())} binders ({y.mean():.1%}). "
         f"{n_ambig} sequences appeared in both sorted classes and were "
         f"{'dropped' if args.drop_ambiguous else 'resolved by read count'}.",
         f"Median mutations per variant: {int(np.median(nm))}. "
         f"Binder rate falls from {data[data.n_mut == 3].label.mean():.2f} at 3 mutations to "
         f"{data[data.n_mut == 10].label.mean():.2f} at 10 — the confound the stratified column controls for.",
         "", "## 1. Zero-shot: do our models predict measured binding?", "", "```",
         zs.round(3).to_string(), "```", "",
         "`auroc_overall` is confounded by mutation count. **`auroc_stratified` is the honest number**: "
         "AUROC computed within each mutation-count group and averaged, i.e. among variants with the same "
         "number of mutations, can the model separate binders from non-binders? 0.5 = no signal.", "",
         "## 2. Model scores vs. measured per-mutation effects", "", "```",
         mve.round(3).to_string(), "```", "",
         f"Measured effects come from logistic regression with mutation count as a covariate "
         f"(its own coefficient is {nmut_coef:.3f} per mutation, confirming the confound). "
         f"Spearman is over the {len(eff)} mutations in the window.", "",
         "## 3. Models trained to predict affinity instead of naturalness", "", "```",
         sup.to_string(index=False), "```", ""]
    if len(verdict):
        ours = verdict[(verdict.votes >= 3) & verdict.reliable].sort_values("measured_pctile", ascending=False)
        best = eff[eff.n_obs >= MIN_OBS].nlargest(10, "effect").merge(
            verdict[["mutation", "votes", "eligible"]], on="mutation", how="left")
        R += [f"## 4. What the measured data says about our own picks", "",
              f"Restricted to mutations observed at least {MIN_OBS} times; rarer ones are shrunk toward "
              f"zero by the regression and carry no information. `measured_pctile` ranks the mutation "
              f"among the {int((eff.n_obs >= MIN_OBS).sum())} reliably-measured mutations in this window.",
              "", "### Our 3+/4 consensus picks, graded", "", "```",
              ours.round(3).to_string(index=False), "```", "",
              "### The 10 best measured mutations in the window, and whether we would have picked them",
              "", "```",
              best[["mutation", "effect", "n_obs", "votes", "eligible"]].round(3).to_string(index=False),
              "```", ""]
        if len(missed):
            R += ["### Mutations our FILTERS discarded that the data ranks in the top 10%", "", "```",
                  missed.round(3).to_string(index=False), "```", "",
                  "These are false negatives of our own pipeline, not of the models: the models voted for "
                  "them and a filter removed them. Worth re-examining before the pool is fixed.", ""]
    R += ["## Caveats", "",
          "- **CDR-H3 only** (scFv 99-108). Says nothing about our H1, H2, framework or VL candidates.",
          "- **Yeast display, not SPR/BLI.** Apparent affinity with avidity effects, confounded by display "
          "level; Adaptyv measures cell-free protein. Correlated, not the same quantity.",
          "- **Binary bins, not KD.** No continuous affinity anywhere in this dataset.",
          "- **Additive scoring ignores epistasis**, and a combinatorial library is where epistasis matters "
          "most. This is a floor on model performance, not a ceiling.",
          "- Mutation-count stratification removes the dominant confound but not display-level effects."]
    with open(os.path.join(args.outdir, "report.md"), "w") as f:
        f.write("\n".join(R))
    print("\n".join(R))
    for m, t in strat_tables.items():
        if len(t):
            print(f"\n{m} by stratum:\n{t.round(3).to_string(index=False)}")
    print(f"\nwrote {args.outdir}/ and data/scores/mason_empirical_scan.csv")


if __name__ == "__main__":
    main()
