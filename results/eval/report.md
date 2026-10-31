# Evaluation report

Models: ablang2, blosum62, esm1v_t33_650M_UR90S_1, esm2_t30_150M, esm2_t33_650M, esm2_t33_650M_split, random; 4313 single mutants (VH and VL positions) in common.

## Layer 1. Intrinsic sanity checks

```
                        n_mutants  cys_mean  cys_pos_pctile  wt_top1  frac_fav_FR  frac_fav_CDR  contact_sens_auc  n_suspect_pos check_cys  measured_rho  measured_p  measured_n
model                                                                                                                                                                           
ablang2                      4313    -8.499           0.089      NaN        0.008         0.155             0.248              6      PASS         0.086       0.396        99.0
blosum62                     4313    -2.105           0.081      NaN        0.108         0.104             0.481              0      PASS         0.079       0.437        99.0
esm1v_t33_650M_UR90S_1       4313   -10.445           0.064    0.784        0.011         0.118             0.178              5      PASS         0.290       0.004        99.0
esm2_t30_150M                4313   -10.294           0.035    0.537        0.052         0.220             0.296              8      PASS        -0.052       0.610        99.0
esm2_t33_650M                4313    -9.901           0.062    0.652        0.026         0.179             0.261              4      PASS         0.138       0.173        99.0
esm2_t33_650M_split          4313    -9.368           0.091    0.670        0.027         0.175             0.267              0      PASS         0.186       0.066        99.0
random                       4313     0.113           0.666      NaN        0.492         0.471             0.521              0      FAIL         0.224       0.026        99.0
```

- `cys_pos_pctile` < 0.25 (`check_cys` PASS): the four disulfide cysteines rank in the most intolerant quarter of positions (0 = most intolerant, random is about 0.5). The criterion is deliberately lenient: it is a bug detector, not a quality score.
- `wt_top1`: fraction of positions where the parent residue is the model's top choice (log-ratio models only).
- `contact_sens_auc` > 0.5: among CDR positions the model rates HER2-contact residues as less tolerant. About 0.5 or below: the model is blind to binding. A crystal contact is not proof of mutational sensitivity, so read this as descriptive.

### Layer 1b. Validation against MEASURED binding

`measured_rho` is the Spearman correlation between the model's score and the measured per-mutation binding effect from the Mason CDR-H3 library (scFv 99-108), over mutations observed at least 50 times. This replaces the crystal-contact proxy with real data, but covers only the CDR-H3 window — it says nothing about H1, H2, framework or VL. See `results/mason/report.md` for the variant-level AUROC, which is the stronger test.

```
                        measured_rho  measured_p  measured_n
model                                                       
ablang2                        0.086       0.396        99.0
blosum62                       0.079       0.437        99.0
esm1v_t33_650M_UR90S_1         0.290       0.004        99.0
esm2_t30_150M                 -0.052       0.610        99.0
esm2_t33_650M                  0.138       0.173        99.0
esm2_t33_650M_split            0.186       0.066        99.0
random                         0.224       0.026        99.0
```

## Layer 2. Agreement between models

Spearman, mutation level:
```
                        ablang2  blosum62  esm1v_t33_650M_UR90S_1  esm2_t30_150M  esm2_t33_650M  esm2_t33_650M_split  random
ablang2                    1.00      0.31                    0.74           0.58           0.65                 0.66    0.01
blosum62                   0.31      1.00                    0.36           0.36           0.35                 0.34   -0.00
esm1v_t33_650M_UR90S_1     0.74      0.36                    1.00           0.76           0.84                 0.84   -0.00
esm2_t30_150M              0.58      0.36                    0.76           1.00           0.86                 0.84    0.01
esm2_t33_650M              0.65      0.35                    0.84           0.86           1.00                 0.99    0.01
esm2_t33_650M_split        0.66      0.34                    0.84           0.84           0.99                 1.00    0.01
random                     0.01     -0.00                   -0.00           0.01           0.01                 0.01    1.00
```

Spearman, position level (mean score per position):
```
                        ablang2  blosum62  esm1v_t33_650M_UR90S_1  esm2_t30_150M  esm2_t33_650M  esm2_t33_650M_split  random
ablang2                    1.00      0.12                    0.76           0.57           0.65                 0.66    0.04
blosum62                   0.12      1.00                    0.19           0.19           0.22                 0.21   -0.05
esm1v_t33_650M_UR90S_1     0.76      0.19                    1.00           0.75           0.85                 0.85   -0.07
esm2_t30_150M              0.57      0.19                    0.75           1.00           0.85                 0.83   -0.01
esm2_t33_650M              0.65      0.22                    0.85           0.85           1.00                 0.98   -0.01
esm2_t33_650M_split        0.66      0.21                    0.85           0.83           0.98                 1.00    0.00
random                     0.04     -0.05                   -0.07          -0.01          -0.01                 0.00    1.00
```

Jaccard overlap of top-100 mutations:
```
                        ablang2  blosum62  esm1v_t33_650M_UR90S_1  esm2_t30_150M  esm2_t33_650M  esm2_t33_650M_split  random
ablang2                    1.00      0.03                    0.27           0.16           0.23                 0.22    0.01
blosum62                   0.03      1.00                    0.04           0.02           0.03                 0.04    0.01
esm1v_t33_650M_UR90S_1     0.27      0.04                    1.00           0.17           0.30                 0.32    0.01
esm2_t30_150M              0.16      0.02                    0.17           1.00           0.28                 0.31    0.01
esm2_t33_650M              0.23      0.03                    0.30           0.28           1.00                 0.72    0.01
esm2_t33_650M_split        0.22      0.04                    0.32           0.31           0.72                 1.00    0.01
random                     0.01      0.01                    0.01           0.01           0.01                 0.01    1.00
```

Consensus vote (ablang2, esm1v_t33_650M_UR90S_1, esm2_t30_150M, esm2_t33_650M; top 20% and score > 0). Eligible = not HER2-contact, > 4.5 A from HER2, not Cys, not suspect-wt: 3402 mutations; votes distribution: {0: 3092, 1: 214, 2: 60, 3: 25, 4: 11}.

```
mutation region  dist  votes  mean_rank_pct
   D102S     H3  5.08      4          0.016
   A106Y     H3  6.83      4          0.018
    N28T     H1 13.25      4          0.023
   D102G     H3  5.08      4          0.023
    A79V     FR 13.84      4          0.028
    T54G     H2  9.02      4          0.032
    N28S     H1 13.25      4          0.035
    T54D     H2  9.02      4          0.035
    T54Y     H2  9.02      4          0.037
   D102T     H3  5.08      4          0.041
    T54N     H2  9.02      4          0.043
    N55S     H2  5.31      3          0.025
    I29F     H1 10.99      3          0.026
    T54S     H2  9.02      3          0.029
    N55G     H2  5.31      3          0.034
```

Note: `_split` variants are robustness checks of the same model and are not counted as separate votes.

## Layer 3. Variant pool

15 variants; 15 pass the hackathon rules; 15 distinct positions mutated; mean pairwise mutation overlap (Jaccard) 0.21.

```
    id  n_mut  valid errors  sum_ablang2  sum_blosum62  sum_esm1v_t33_650M_UR90S_1  sum_esm2_t30_150M  sum_esm2_t33_650M  sum_esm2_t33_650M_split  sum_random
TZ01_A      5   True               -0.61          -2.0                       -2.77               2.01              -4.31                    -2.95       -4.72
TZ02_A      5   True               -3.43          -4.0                       -2.19               1.97              -1.87                    -0.19       -4.76
TZ03_A      5   True                6.68           0.0                        2.39               5.18               2.78                     3.08       -3.02
TZ04_A      5   True                3.86          -2.0                        2.96               5.13               5.22                     5.85       -3.05
TZ05_A      5   True               -0.18          -3.0                        0.42               1.33              -0.15                    -0.18       -2.93
TZ06_B      5   True               -3.30          -1.0                       -4.75               3.24              -0.78                     0.31       -2.83
TZ07_B      5   True                3.99           1.0                        0.41               6.41               6.31                     6.34       -1.12
TZ08_B      5   True               -0.08          -4.0                       -0.20              -1.50               3.52                     3.23       -1.87
TZ09_B      5   True               -5.80          -3.0                       -0.37               5.34               2.15                     2.43       -2.13
TZ10_C      5   True                0.72           0.0                       -0.65               6.48               4.58                     5.35        1.78
TZ11_C      5   True                3.07           1.0                        1.00               5.04               3.99                     4.28       -1.28
TZ12_C      5   True                2.14           0.0                       -2.92               6.65               5.25                     5.48        1.59
TZ13_D      5   True                3.92          -1.0                        8.04               8.31               8.80                     8.85       -2.10
TZ14_D      5   True                6.31          -1.0                        5.58               6.99               9.30                     8.66       -1.29
TZ15_D      5   True                5.02           2.0                        5.39               4.27               7.18                     7.35       -3.45
```

_No `kd_nM` column: variants were validated and scored but no performance can be measured yet._
