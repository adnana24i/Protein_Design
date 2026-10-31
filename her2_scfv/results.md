# HER2 scFv co-folds — ipSAE ranking

Ranked by **ipSAE_min** (min of the two asymmetric scFv<->HER2 directions; IPSAE, PAE/dist cutoffs 10/10) across the 16 Boltz co-folds (parent + V01-V15). ipTM / pTM / pLDDT are the Boltz complex confidence metrics (0-1).

| Rank | Variant | ipSAE_min | ipTM | pTM | pLDDT | #mut |
|---|---|---|---|---|---|---|
| 1 | V12 | 0.821 (+0.059 vs parent) | 0.923 | 0.950 | 0.876 | 5 |
| 2 | V11 | 0.786 (+0.025 vs parent) | 0.913 | 0.945 | 0.907 | 5 |
| 3 | **parent** _(reference baseline)_ | 0.762 | 0.909 | 0.932 | 0.910 | 0 |
| 4 | V02 | 0.749 (-0.013 vs parent) | 0.891 | 0.935 | 0.884 | 5 |
| 5 | V04 | 0.742 (-0.019 vs parent) | 0.893 | 0.931 | 0.896 | 6 |
| 6 | V09 | 0.742 (-0.020 vs parent) | 0.901 | 0.928 | 0.904 | 5 |
| 7 | V08 | 0.731 (-0.031 vs parent) | 0.902 | 0.930 | 0.897 | 5 |
| 8 | V13 | 0.725 (-0.037 vs parent) | 0.905 | 0.931 | 0.908 | 5 |
| 9 | V07 | 0.720 (-0.041 vs parent) | 0.886 | 0.931 | 0.902 | 5 |
| 10 | V01 | 0.719 (-0.043 vs parent) | 0.880 | 0.925 | 0.897 | 5 |
| 11 | V10 | 0.708 (-0.054 vs parent) | 0.874 | 0.930 | 0.904 | 5 |
| 12 | V05 | 0.679 (-0.082 vs parent) | 0.861 | 0.915 | 0.888 | 7 |
| 13 | V15 | 0.661 (-0.101 vs parent) | 0.871 | 0.926 | 0.901 | 5 |
| 14 | V03 | 0.628 (-0.134 vs parent) | 0.849 | 0.916 | 0.887 | 5 |
| 15 | V14 | 0.600 (-0.161 vs parent) | 0.837 | 0.910 | 0.895 | 5 |
| 16 | V06 | 0.522 (-0.240 vs parent) | 0.846 | 0.900 | 0.903 | 5 |

> ipSAE_min is an interface-confidence rank, not a Kd — used to filter/rank, not to claim nM affinities. Measured Kd comes from the Adaptyv wet lab.
