# Result tables for the manuscript (generated from outputs/06_lumawarp_audit)

## Edge audit given the ordering (10 replicates; MGM-PC-Stable from the draft for reference)

| dataset             | learner                        | ordering   |   recall |   recall_mcse |   precision |   noise_kept |   edges_called |
|:--------------------|:-------------------------------|:-----------|---------:|--------------:|------------:|-------------:|---------------:|
| full_set            | LumaWarp                       | given      |    0.852 |         0.009 |       0.932 |        0.001 |           57.7 |
| full_set            | spline additive                | given      |    0.846 |         0.008 |       0.928 |        0.002 |           57.5 |
| full_set            | logistic/linear                | given      |    0.803 |         0.007 |       0.992 |        0     |           51   |
| full_set            | XGBoost                        | given      |    0.733 |         0.01  |       0.908 |        0.002 |           51   |
| full_set            | MGM-PC-Stable (draft, 50 reps) | none       |    0.78  |       nan     |       1     |      nan     |          nan   |
| reference_subsample | LumaWarp                       | given      |    0.476 |         0.02  |       0.909 |        0.004 |           33.7 |
| reference_subsample | spline additive                | given      |    0.459 |         0.012 |       0.823 |        0.007 |           35.5 |
| reference_subsample | logistic/linear                | given      |    0.478 |         0.006 |       0.912 |        0     |           33.1 |
| reference_subsample | XGBoost                        | given      |    0.414 |         0.015 |       0.965 |        0.005 |           27.2 |
| reference_subsample | MGM-PC-Stable (draft, 50 reps) | none       |    0.56  |       nan     |       0.99  |      nan     |          nan   |
| astronaut_set       | LumaWarp                       | given      |    0.459 |         0.015 |       0.843 |        0.003 |           35   |
| astronaut_set       | spline additive                | given      |    0.454 |         0.016 |       0.749 |        0.007 |           38.2 |
| astronaut_set       | logistic/linear                | given      |    0.471 |         0.014 |       0.926 |        0.001 |           32.2 |
| astronaut_set       | XGBoost                        | given      |    0.408 |         0.009 |       0.912 |        0.002 |           28.2 |
| astronaut_set       | MGM-PC-Stable (draft, 50 reps) | none       |    0.55  |       nan     |       0.97  |      nan     |          nan   |

## DAG-parent adjusted estimation of the per-unit contrast (50 replicates; the ordering-only route is in estimation/summary.csv)

| treatment       | dataset             | estimator              |   reps |   truth |   mean_estimate |    bias |   empirical_se |   mean_se |   coverage |   sign_right |
|:----------------|:--------------------|:-----------------------|-------:|--------:|----------------:|--------:|---------------:|----------:|-----------:|-------------:|
| k_citrate       | full_set            | Naive difference       |     50 | -0.0202 |          0.0322 |  0.0524 |         0.0103 |  nan      |     nan    |         0    |
| k_citrate       | full_set            | LW-TMLE                |     50 | -0.0202 |         -0.0202 | -0      |         0.0109 |    0.0096 |       0.92 |         0.92 |
| k_citrate       | full_set            | SL-TMLE                |     50 | -0.0202 |         -0.0204 | -0.0002 |         0.0106 |    0.0095 |       0.92 |         0.94 |
| k_citrate       | full_set            | Logistic g-computation |     50 | -0.0202 |         -0.0189 |  0.0013 |         0.008  |  nan      |     nan    |         1    |
| k_citrate       | reference_subsample | Naive difference       |     50 | -0.0202 |          0.0357 |  0.0559 |         0.0325 |  nan      |     nan    |         0.1  |
| k_citrate       | reference_subsample | LW-TMLE                |     50 | -0.0202 |         -0.0181 |  0.0021 |         0.0328 |    0.0311 |       0.96 |         0.78 |
| k_citrate       | reference_subsample | SL-TMLE                |     50 | -0.0202 |         -0.0207 | -0.0005 |         0.031  |    0.0306 |       0.94 |         0.76 |
| k_citrate       | reference_subsample | Logistic g-computation |     50 | -0.0202 |         -0.0178 |  0.0024 |         0.0279 |  nan      |     nan    |         0.82 |
| k_citrate       | astronaut_set       | Naive difference       |     50 | -0.0152 |          0.0189 |  0.0341 |         0.0348 |  nan      |     nan    |         0.3  |
| k_citrate       | astronaut_set       | LW-TMLE                |     50 | -0.0152 |         -0.0222 | -0.007  |         0.0305 |    0.0268 |       0.82 |         0.74 |
| k_citrate       | astronaut_set       | SL-TMLE                |     50 | -0.0152 |         -0.0268 | -0.0116 |         0.0296 |    0.0254 |       0.74 |         0.78 |
| k_citrate       | astronaut_set       | Logistic g-computation |     50 | -0.0152 |         -0.0188 | -0.0036 |         0.0266 |  nan      |     nan    |         0.8  |
| bisphosphonates | full_set            | Naive difference       |     50 | -0.0083 |         -0.0081 |  0.0003 |         0.009  |  nan      |     nan    |         0.82 |
| bisphosphonates | full_set            | LW-TMLE                |     50 | -0.0083 |         -0.004  |  0.0044 |         0.0132 |    0.0125 |       0.96 |         0.64 |
| bisphosphonates | full_set            | SL-TMLE                |     50 | -0.0083 |         -0.0036 |  0.0047 |         0.0131 |    0.0123 |       0.96 |         0.64 |
| bisphosphonates | full_set            | Logistic g-computation |     50 | -0.0083 |         -0.0076 |  0.0007 |         0.0089 |  nan      |     nan    |         0.82 |
| bisphosphonates | reference_subsample | Naive difference       |     50 | -0.0083 |         -0.0046 |  0.0037 |         0.0223 |  nan      |     nan    |         0.58 |
| bisphosphonates | reference_subsample | LW-TMLE                |     50 | -0.0083 |         -0.0012 |  0.0071 |         0.0439 |    0.0383 |       0.84 |         0.48 |
| bisphosphonates | reference_subsample | SL-TMLE                |     50 | -0.0083 |         -0.0012 |  0.0072 |         0.0448 |    0.0391 |       0.88 |         0.54 |
| bisphosphonates | reference_subsample | Logistic g-computation |     50 | -0.0083 |         -0.0063 |  0.0021 |         0.0255 |  nan      |     nan    |         0.58 |
| bisphosphonates | astronaut_set       | Naive difference       |     50 | -0.0063 |         -0.0036 |  0.0027 |         0.0293 |  nan      |     nan    |         0.5  |
| bisphosphonates | astronaut_set       | LW-TMLE                |     50 | -0.0063 |          0.0073 |  0.0137 |         0.0453 |    0.0382 |       0.84 |         0.46 |
| bisphosphonates | astronaut_set       | SL-TMLE                |     50 | -0.0063 |          0.0102 |  0.0165 |         0.0462 |    0.04   |       0.82 |         0.44 |
| bisphosphonates | astronaut_set       | Logistic g-computation |     50 | -0.0063 |          0.0046 |  0.011  |         0.0357 |  nan      |     nan    |         0.44 |

## Propensity-bound sensitivity, K+ citrate, DAG parents, LW-TMLE (50 replicates)

|   bound |   coverage_astronaut_set |   coverage_full_set |   coverage_reference_subsample |   bias_astronaut_set |   bias_full_set |   bias_reference_subsample |   clipped_astronaut_set |   clipped_full_set |   clipped_reference_subsample |
|--------:|-------------------------:|--------------------:|-------------------------------:|---------------------:|----------------:|---------------------------:|------------------------:|-------------------:|------------------------------:|
|   0.001 |                     0.8  |                0.92 |                           0.96 |               -0.007 |           0.001 |                      0.004 |                   0     |              0     |                         0     |
|   0.01  |                     0.82 |                0.92 |                           0.96 |               -0.008 |           0     |                      0.004 |                   0.029 |              0.031 |                         0.026 |
|   0.025 |                     0.82 |                0.92 |                           0.96 |               -0.007 |          -0     |                      0.002 |                   0.155 |              0.149 |                         0.109 |
|   0.05  |                     0.82 |                0.9  |                           0.94 |               -0.006 |          -0     |                      0.002 |                   0.4   |              0.331 |                         0.288 |

## Attribution with LumaWarp (lw) versus the Super Learner (sl) as the explained model (10 replicates, 32 background records)

| dataset             | method              |   reps |   tau_lw |   tau_sl |   diff |   diff_mcse |   top5_lw |   top5_sl |   auc_lw |   auc_sl |
|:--------------------|:--------------------|-------:|---------:|---------:|-------:|------------:|----------:|----------:|---------:|---------:|
| full_set            | interventional_shap |     10 |    0.814 |    0.811 |  0.003 |       0.014 |      1    |      1    |    0.82  |    0.822 |
| full_set            | ordering_only_shap  |     10 |    0.407 |    0.359 |  0.048 |       0.03  |      0.7  |      0.62 |    0.82  |    0.822 |
| full_set            | standard_shap       |     10 |    0.404 |    0.352 |  0.052 |       0.029 |      0.68 |      0.6  |    0.82  |    0.822 |
| reference_subsample | interventional_shap |     10 |    0.63  |    0.618 |  0.012 |       0.032 |      0.84 |      0.78 |    0.752 |    0.796 |
| reference_subsample | ordering_only_shap  |     10 |    0.34  |    0.307 |  0.033 |       0.042 |      0.6  |      0.58 |    0.752 |    0.796 |
| reference_subsample | standard_shap       |     10 |    0.347 |    0.29  |  0.057 |       0.054 |      0.64 |      0.6  |    0.752 |    0.796 |
| astronaut_set       | interventional_shap |     10 |    0.579 |    0.551 |  0.027 |       0.046 |      0.78 |      0.8  |    0.771 |    0.79  |
| astronaut_set       | ordering_only_shap  |     10 |    0.354 |    0.361 | -0.007 |       0.03  |      0.56 |      0.6  |    0.771 |    0.79  |
| astronaut_set       | standard_shap       |     10 |    0.345 |    0.345 |  0     |       0.023 |      0.54 |      0.58 |    0.771 |    0.79  |

## Curved-mechanism variant: hybrid recovery of the three curved edges (full set, 20 replicates)

| arm                         | method             |   reps |   edges |   precision |   recall |       f1 |   u_edges_found |   u_edges_oriented |   saturating_found |   false_edges |   added_false |   orient_acc |
|:----------------------------|:-------------------|-------:|--------:|------------:|---------:|---------:|----------------:|-------------------:|-------------------:|--------------:|--------------:|-------------:|
| off-centre U-shapes (arm A) | Hybrid, chain veto |     20 |   57.9  |    0.930147 | 0.853968 | 0.890086 |             3   |               3    |                  1 |          4.1  |          2.75 |     0.775878 |
| off-centre U-shapes (arm A) | Hybrid, no veto    |     20 |   59.6  |    0.903151 | 0.853968 | 0.877647 |             3   |               3    |                  1 |          5.8  |          2.75 |     0.775878 |
| off-centre U-shapes (arm A) | PC-stable alone    |     20 |   55.15 |    0.975983 | 0.853968 | 0.910677 |             3   |               1.95 |                  1 |          1.35 |        nan    |     0.755935 |
| symmetric U-shapes (arm B)  | Hybrid, chain veto |     20 |   59.15 |    0.912422 | 0.856349 | 0.88335  |             3   |               3    |                  1 |          5.2  |          3.25 |     0.764762 |
| symmetric U-shapes (arm B)  | Hybrid, no veto    |     20 |   60    |    0.899565 | 0.856349 | 0.877253 |             3   |               3    |                  1 |          6.05 |          3.25 |     0.764762 |
| symmetric U-shapes (arm B)  | PC-stable alone    |     20 |   53    |    0.963697 | 0.810317 | 0.880172 |             0.1 |               0.1  |                  1 |          1.95 |        nan    |     0.751008 |

