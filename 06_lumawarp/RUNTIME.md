# Runtime record for the Lumawarp analyses (what produced the committed tables)

All runs on one machine: Ubuntu Linux (kernel 7.0.0), 16 CPU cores, NVIDIA GeForce RTX 5090 Laptop GPU (driver
595.99.02). Engine fits ran on the GPU except where noted in the folder RESULTS files; inference is CPU only.

## Engine
| use | executable | version string | binary date | SHA-256 (first 16) |
|---|---|---|---|---|
| edge audit, estimation nuisances, attribution base model, nonlinear hybrid | Lumawarp engine, training executable | version 0.0.1, September 2026 build | 2026-09-08 | 333178de463caed7 |
| same | Lumawarp engine, inference executable | version 0.0.1, September 2026 build | 2026-09-08 | 689b1335b12f57fb |

The October 2026 build (version 0.0.1, 2026-10-12, f3b38fb9d7ce5602 / 178b878984a626cd) was used only for checks
that are not in the manuscript; every reported row comes from the September build.
Configuration files as used (keys and values): `edge_audit/lumawarp_research.cfg`, `estimation/lumawarp_playbook.cfg`
and `attribution_lumawarp_base/cfg_attribution.cfg`. The engine's public interface is the two executables, the
configuration file and the CSV inputs and outputs; the model file and its reader are not distributed, so the
monotonicity reading is reported from the committed shape_flags.csv.

## R (estimation, bound sensitivity, structure validation of the nonlinear variant)
R 4.3.3; tmle 2.1.1; SuperLearner 2.0.42; glmnet 5.1; earth 5.3.6; simcausal 0.5.7; rCausalMGM 1.0.1; data.table
1.18.6.1; dagitty 0.3.4. `tmle()` calls: family binomial, fluctuation logistic, gbound 0.025/0.975 (and 0.001, 0.01,
0.05 in the bound sensitivity), prescreenW.g FALSE, target.gwt FALSE.

## Python (job construction, reference learners, scoring, attribution driver)
Python 3.12.3; numpy 2.5.3; pandas 3.0.6; scikit-learn 1.9.1; xgboost 3.4.1; shap 0.52.0; scipy 1.18.1; networkx
3.7; causal-learn (PC-stable in the nonlinear hybrid). The attribution driver imports the pipeline's own
`04_attribution_validation` package at the pinned revision for split, records, orders, SCM and scoring.

## Data identity
Main-generator replicates were regenerated on this machine from the pipeline at commit 79a1b56 (R 4.3.3, simcausal
0.5.7) and carry configuration hash d2671ab6dc028eccc7a32b0cc8312afc, the hash recorded in the repository's committed
`02_data/provenance.txt` (generated on the authors' machine with R 4.4.1, simcausal 0.5.6); the committed
`truth_values.csv` is byte-identical to the regenerated one. Nonlinear variant: configuration hash
376d9da24a8c2cddef59263ea2578760, generated 2026-10-13 from the working tree at 74a7e72 with the recalibrated
`calibration.yaml` that was committed afterwards (faf959e on the historical branch) and is the calibration in
`config/variants/nonlinear_symmetric_u/`; regenerating from the opt-in configuration reproduces the replicate files
byte for byte. The 50-replicate MGM-PC-Stable outputs on the variant (`outputs/06_lumawarp/nonlinear_sensitivity/armB/
structure_validation/`, provenance commit 74a7e72, data config hash 376d9da2, 2026-10-14) and on the off-centre arm A
(`.../armA/`, data config hash 93b850f1, 2026-10-13) are the pipeline's own structure-validation outputs.

## Timing of pre-registrations and thresholds (file-system timestamps; the run folders are not under version control)
Per-learner edge thresholds (`thresholds.json`, Lumawarp 5, logistic 4, spline 2, XGBoost 3) and the tuning table
(`threshold_tuning_seeds100-102.csv`, tuning on a separate simulator's seeds 100 to 102): written 2026-10-05 06:30,
before that folder's first test fit (06:30:40) and nine days before the paper-data audit. The paper-data audit's
PREREGISTRATION.md and score.py: 2026-10-14 09:52:58; its first engine fit: 2026-10-14 09:54:21 (run log). Estimation
PREREGISTRATION.md preceded its first fit on 2026-10-14 likewise; attribution and bound-sensitivity pre-registrations
on 2026-10-19 preceded their runs the same day.

## Deposit (not in the repository)
`lumawarp_nuisance.csv.gz` (35.3 MB; the cross-fitted Lumawarp outcome and propensity predictions for all 600
estimation units, keyed by data set, treatment, route, replicate and record id) and `lumawarp_predictions.csv.gz`
(2.9 MB; also committed under `outputs/06_lumawarp/edge_audit/`), for the data archive accompanying the paper.
