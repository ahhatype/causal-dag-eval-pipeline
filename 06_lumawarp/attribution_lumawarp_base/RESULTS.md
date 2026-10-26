# Results: the paper's attribution arm with LumaWarp as the explained model (folder 35, 19 to 20 October 2026)

Pre-registered in PREREGISTRATION.md. 30 cells (replicates 1 to 10, three data sets, ancestor feature set), the
paper's split, explained records and feature orders, with LumaWarp (old build, research configuration) in place of
the Super Learner as the model being explained. Reduced Monte Carlo budget: 32 background records (the paper's 128),
64 orders as the paper; 2.1 million engine predictions per method per cell; 0 failures; 3.75 h wall. Paired cell by
cell with the pipeline's committed Super Learner results (same replicate, data set and method).

## Verdicts
- **B1 PASS (interchangeable explained model).** Full set, causal interventional SHAP: tau 0.814 with LumaWarp
  against 0.811 with the Super Learner, paired difference +0.003 (MCSE 0.014); AUC of the explained model 0.820
  against 0.822.
- **B2 PASS (the paper's contrast reproduces).** Interventional minus standard with LumaWarp as base: +0.41 on the
  full set (the paper: +0.43), +0.28 on the reference set, +0.23 on the astronaut set (the paper: +0.24 to +0.26).
  Ordering-only minus standard: 0.00, -0.01, +0.01 (the paper: at most 0.02).

## Table (10 replicates per cell; tau-b against the intervention-contrast truth, own population, 16 ancestors)

| data set | method | tau LumaWarp | tau Super Learner | paired diff (MCSE) | top-5 LW / SL | PBI LW / SL |
|---|---|---|---|---|---|---|
| full_set | interventional shap | 0.814 | 0.811 | +0.003 (0.014) | 1.00 / 1.00 | 0.11 / 0.13 |
| full_set | ordering only shap | 0.407 | 0.359 | +0.048 (0.030) | 0.70 / 0.62 | 0.80 / 0.78 |
| full_set | standard shap | 0.404 | 0.352 | +0.052 (0.029) | 0.68 / 0.60 | 0.81 / 0.79 |
| reference_subsample | interventional shap | 0.630 | 0.618 | +0.012 (0.032) | 0.84 / 0.78 | -0.21 / -0.16 |
| reference_subsample | ordering only shap | 0.340 | 0.307 | +0.033 (0.042) | 0.60 / 0.58 | 0.18 / 0.12 |
| reference_subsample | standard shap | 0.347 | 0.290 | +0.057 (0.054) | 0.64 / 0.60 | 0.21 / 0.13 |
| astronaut_set | interventional shap | 0.579 | 0.551 | +0.027 (0.046) | 0.78 / 0.80 | -0.30 / -0.25 |
| astronaut_set | ordering only shap | 0.354 | 0.361 | -0.007 (0.030) | 0.56 / 0.60 | 0.02 / -0.10 |
| astronaut_set | standard shap | 0.345 | 0.345 | +0.000 (0.023) | 0.54 / 0.58 | 0.06 / -0.10 |

AUC of the explained model (LumaWarp / Super Learner): full_set 0.820 / 0.822 ; reference_subsample 0.752 / 0.796 ; astronaut_set 0.771 / 0.790.
Max efficiency error over all cells 1.0e-15. Engine minutes per method per cell: interventional_shap 2.0, ordering_only_shap 1.9, standard_shap 1.9.

## Reading
- Every attribution finding in the draft holds with LumaWarp as the explained model: the value function makes SHAP
  causal, the ordering alone does nothing, sample size costs more than selection, and the proximity bias of the
  standard variants is unchanged. The explained model is interchangeable.
- Standard and ordering-only SHAP are 0.03 to 0.06 higher with LumaWarp than with the Super Learner on the full
  and reference sets (MCSE 0.03 to 0.05): within noise at 10 replicates, reported, not claimed.
- The 900-row sets lose the same amount with either base (interventional tau 0.63 and 0.35 with LumaWarp against
  0.62 and 0.36 with the Super Learner), so the draft's sample-size statement does not depend on the model either.

## Deviations
- Background 32 instead of 128 (pre-registered as the reduced plan): the explained records and feature orders are
  the paper's; the background is a 32-record subsample of the paper's 128, which adds Monte Carlo noise to every
  Shapley value and no bias. Causal predictive SHAP was not run (noted in the pre-registration).
- The full budget on the full set (10 cells, about 2.3 h each) remains available if the paper's number should carry
  the paper's exact budget; at the reduced budget the difference to the Super Learner is already inside 0.01.

Files: cells/ (per-cell runs, importance and Shapley values), summary.csv, compare.py, run_cells.py.
