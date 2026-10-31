# Measured-affinity benchmark (Mason et al. CDR-H3, scFv 99-108)

36368 variants, 10554 binders (29.0%). 2342 sequences appeared in both sorted classes and were resolved by read count.
Median mutations per variant: 8. Binder rate falls from 0.93 at 3 mutations to 0.07 at 10 — the confound the stratified column controls for.

## 1. Zero-shot: do our models predict measured binding?

```
                           n_scored  auroc_overall  auroc_stratified
model                                                               
esm1v_t33_650M_UR90S_1        36367          0.727             0.637
blosum62                      36367          0.567             0.562
esm2_t33_650M_split           36367          0.655             0.562
esm2_t30_150M                 36367          0.626             0.556
ablang2                       36367          0.678             0.541
esm2_t33_650M                 36367          0.637             0.537
random                        36367          0.530             0.517
baseline: fewer mutations     36368          0.774             0.500
```

`auroc_overall` is confounded by mutation count. **`auroc_stratified` is the honest number**: AUROC computed within each mutation-count group and averaged, i.e. among variants with the same number of mutations, can the model separate binders from non-binders? 0.5 = no signal.

## 2. Model scores vs. measured per-mutation effects

```
                        n_mutations  spearman_vs_measured  p_value
model                                                             
esm1v_t33_650M_UR90S_1          190                 0.168    0.020
esm2_t33_650M_split             190                 0.164    0.024
esm2_t33_650M                   190                 0.130    0.075
random                          190                 0.115    0.115
blosum62                        190                 0.055    0.453
ablang2                         190                 0.048    0.508
esm2_t30_150M                   190                -0.046    0.531
```

Measured effects come from logistic regression with mutation count as a covariate (its own coefficient is -0.930 per mutation, confirming the confound). Spearman is over the 190 mutations in the window.

## 3. Models trained to predict affinity instead of naturalness

```
                                   model  n_train  n_test  auroc_test
            supervised: one-hot sequence    27276    9092    0.891423
supervised: one-hot, within n_mut=8 only     7172    2391    0.838557
        supervised: one-hot, n_total=500      375     125    0.874855
       supervised: one-hot, n_total=2000     1500     500    0.871238
       supervised: one-hot, n_total=8000     6000    2000    0.899790
             supervised: ESM2 embeddings     1500     500    0.851505
```

## 4. What the measured data says about our own picks

Restricted to mutations observed at least 50 times; rarer ones are shrunk toward zero by the regression and carry no information. `measured_pctile` ranks the mutation among the 99 reliably-measured mutations in this window.

### Our 3+/4 consensus picks, graded

```
mutation region  dist  votes  eligible  effect  measured_pctile  n_obs  reliable
   M107F     H3  9.36      4     False   2.082            1.000  11565      True
   M107Y     H3  9.36      3     False   1.194            0.929   5658      True
   D102S     H3  5.08      4      True   0.995            0.889   4376      True
    W99Y     H3  4.05      3     False   0.955            0.859   8733      True
   D102T     H3  5.08      4      True   0.877            0.848   1674      True
   D102G     H3  5.08      4      True   0.856            0.838   2759      True
   M107L     H3  9.36      3     False   0.763            0.828   6091      True
   G101S     H3  5.89      3      True   0.404            0.747   1226      True
   F104Y     H3  4.27      4     False   0.262            0.697   2223      True
   M107N     H3  9.36      3     False  -0.197            0.465   4311      True
   A106T     H3  6.83      3      True  -0.249            0.394   3968      True
   F104S     H3  4.27      4     False  -0.312            0.333   2886      True
   A106D     H3  6.83      3      True  -0.315            0.323    787      True
    W99S     H3  4.05      3     False  -0.675            0.131    272      True
   M107I     H3  9.36      3     False  -0.755            0.101   4034      True
```

### The 10 best measured mutations in the window, and whether we would have picked them

```
mutation  effect  n_obs  votes  eligible
   M107F   2.082  11565      4     False
   D102P   1.684   2087      0      True
   D102Y   1.333   1837      1      True
   D102F   1.286   2325      0      True
   D102A   1.214   2099      2      True
   D102V   1.210   2598      0      True
   D102H   1.199   1568      0      True
   M107Y   1.194   5658      3     False
   D102I   1.165   2316      0      True
   D102C   1.065   4745      0     False
```

### Mutations our FILTERS discarded that the data ranks in the top 10%

```
mutation region  dist  votes  eligible  effect  measured_pctile  n_obs  reliable
   M107F     H3  9.36      4     False   2.082            1.000  11565      True
   M107Y     H3  9.36      3     False   1.194            0.929   5658      True
   D102C     H3  5.08      0     False   1.065            0.909   4745      True
```

These are false negatives of our own pipeline, not of the models: the models voted for them and a filter removed them. Worth re-examining before the pool is fixed.

## Caveats

- **CDR-H3 only** (scFv 99-108). Says nothing about our H1, H2, framework or VL candidates.
- **Yeast display, not SPR/BLI.** Apparent affinity with avidity effects, confounded by display level; Adaptyv measures cell-free protein. Correlated, not the same quantity.
- **Binary bins, not KD.** No continuous affinity anywhere in this dataset.
- **Additive scoring ignores epistasis**, and a combinatorial library is where epistasis matters most. This is a floor on model performance, not a ceiling.
- Mutation-count stratification removes the dominant confound but not display-level effects.