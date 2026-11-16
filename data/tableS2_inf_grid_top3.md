# Top-3 INF configurations per decile (CV TPR)

Decile 1 corresponds to targets with 0–5 known AEs (sparsest); decile 10 corresponds to 163–256 AEs (densest). 'config_idx' is the position of the parameter triple in `param_combos.pkl`. Source: ``analyse_parameter_slices_quantile`` applied per fold to ``cv_gs_results_l3.pkl`` (fold-mean estimator, t-based 95% half-width).

| decile | rank | config | deflator exp | INF aggregator | pooling | CV TPR | ± 95% CI |
|---|---|---|---|---|---|---|---|
| 1 | 1 | 12 | 0 | product | max | 0.42 | ± 0.11 |
| 1 | 2 | 13 | 0 | product | mean | 0.42 | ± 0.17 |
| 1 | 3 | 9 | 0 | sum | max | 0.40 | ± 0.12 |
| 2 | 1 | 7 | 0 | — | — | 0.66 | ± 0.15 |
| 2 | 2 | 27 | 0.5 | product | mean | 0.64 | ± 0.09 |
| 2 | 3 | 15 | 0.3 | sum | min | 0.64 | ± 0.07 |
| 3 | 1 | 15 | 0.3 | sum | min | 0.91 | ± 0.06 |
| 3 | 2 | 20 | 0.3 | product | mean | 0.89 | ± 0.07 |
| 3 | 3 | 27 | 0.5 | product | mean | 0.89 | ± 0.08 |
| 4 | 1 | 17 | 0.3 | sum | mean | 0.98 | ± 0.03 |
| 4 | 2 | 31 | 0.75 | sum | mean | 0.98 | ± 0.03 |
| 4 | 3 | 28 | 0.75 | — | — | 0.98 | ± 0.03 |
| 5 | 1 | 17 | 0.3 | sum | mean | 0.99 | ± 0.01 |
| 5 | 2 | 32 | 0.75 | product | min | 0.99 | ± 0.01 |
| 5 | 3 | 3 | — | sum | mean | 0.99 | ± 0.01 |
| 6 | 1 | 31 | 0.75 | sum | mean | 1.00 | ± 0.00 |
| 6 | 2 | 30 | 0.75 | sum | max | 1.00 | ± 0.00 |
| 6 | 3 | 29 | 0.75 | sum | min | 1.00 | ± 0.00 |
| 7 | 1 | 17 | 0.3 | sum | mean | 0.98 | ± 0.01 |
| 7 | 2 | 27 | 0.5 | product | mean | 0.98 | ± 0.01 |
| 7 | 3 | 26 | 0.5 | product | max | 0.98 | ± 0.01 |
| 8 | 1 | 0 | — | — | — | 0.98 | ± 0.04 |
| 8 | 2 | 23 | 0.5 | sum | max | 0.98 | ± 0.04 |
| 8 | 3 | 22 | 0.5 | sum | min | 0.98 | ± 0.04 |
| 9 | 1 | 34 | 0.75 | product | mean | 0.99 | ± 0.01 |
| 9 | 2 | 30 | 0.75 | sum | max | 0.99 | ± 0.01 |
| 9 | 3 | 5 | — | product | max | 0.99 | ± 0.01 |
| 10 | 1 | 25 | 0.5 | product | min | 0.97 | ± 0.02 |
| 10 | 2 | 20 | 0.3 | product | mean | 0.97 | ± 0.02 |
| 10 | 3 | 17 | 0.3 | sum | mean | 0.97 | ± 0.03 |
