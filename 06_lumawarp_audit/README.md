# 06_lumawarp_audit: ordering-informed edge audit, DAG-adjusted estimation and sensitivity checks with a learned-kernel ensemble

Secondary analyses on the pipeline's own replicate data (`02_data`, commit 79a1b56 generator), run 14 to 20 October
2026 by Lauren Morgenthaler (Lucidity Sciences). Each directory holds the pre-registration written before its first
fit, the scripts, and a RESULTS.md; result tables are in `outputs/06_lumawarp_audit/`, with manuscript-ready tables in
`outputs/06_lumawarp_audit/tables/TABLES.md`.

| directory | question | engine fits | headline |
|---|---|---|---|
| `edge_audit/` | given the DAG's topological order, how many of the 63 analysis-graph edges does a per-node flexible audit recover, at what precision, on the three data sets | 1,050 | full-set recall 0.85 at precision 0.93 (spline 0.85 at 0.93, logistic 0.80 at 0.99); K+ citrate -> urine chemistry 10 of 10 replicates; 900-row sets 0.46 to 0.48 for every learner; 0 of 160 shape flags on the full set |
| `estimation/` | does adjustment by the DAG's parents (or by the ordering alone) recover the per-unit contrast for K+ citrate and bisphosphonates, with a calibrated interval, with LumaWarp, a Super Learner, or a logistic model as the nuisance learner | 600 units x 6 fits | full-set K+ citrate: naive +0.032 (sign wrong 50 of 50) becomes -0.020 against a true -0.020, coverage 0.92, identical across learners and routes; at 900 rows bias removed, SE 0.03 exceeds the effect |
| `estimation_bounds/` | is the astronaut-set under-coverage (0.82) the propensity bounds or the sample | 0 (re-targeting only) | no bound from 0.001 to 0.05 changes coverage or bias by more than 0.002: the sample |
| `attribution_lumawarp_base/` | do the attribution findings hold when LumaWarp replaces the Super Learner as the explained model, with the same split, records and feature orders | 30 | causal interventional tau 0.814 vs 0.811 on the full set; interventional minus standard +0.41 (draft +0.43); ordering-only minus standard 0.00 |
| `nonlinear_sensitivity/` | what curvature does to shape validation inside this generator (branch `nonmonotone-variant`: three edges made U-shaped, one saturating) | 0 here (PC-stable and the hybrid are CPU) | symmetric U-shapes: MGM-PC-Stable recovers none of the three curved edges in 150 data sets and the linear implication test certifies their absence; the flexible two-way audit and a basis-expanded conditional-independence test recover all three |
| `shape_flag_checks/` | the monotonicity flag's dependence on the test's functional form and on the engine build | 144 | the flag is a kernel-slice degree-centroid reading (format-1 model files); the newer build's importance file does not carry it |

## Method in one paragraph (the ordering-informed audit)
For every analysis-graph variable with at least one predecessor, fit it on all earlier variables in the topological
order plus 12 pure-noise columns (6 continuous, 6 binary; `make_jobs.py`). Sweep each candidate across its range with
the others at their median (`score.csv`, 21 points per input); edge score = sweep span / largest noise-column span.
Keep an edge if the score exceeds the learner's threshold (LumaWarp 5x, logistic 4x, spline 2x, XGBoost 3x; tuned on
separate simulated seeds and never changed, see `edge_audit/PREREGISTRATION.md`). Orientation is given by the order.
Measures: adjacency recall and precision against the 63 edges, noise columns kept, shape flags (continuous parents of
continuous targets whose relevance-weighted degree centroid exceeds every continuous noise column's).

## Method in one paragraph (estimation)
Per replicate, three stratified folds; the learned-kernel ensemble (playbook profile, `lumawarp_playbook.cfg`) and a
Super Learner (SL.glm, SL.glmnet, SL.earth) fit the outcome on treatment plus covariates and the treatment on the
covariates on the other folds, predicting the held-out fold at treatment 0 and 1; `tmle()` (R package tmle; binomial,
logistic fluctuation, bounds 0.025 to 0.975) targets each; logistic g-computation and the naive difference alongside.
Adjustment sets: the exposure's analysis-graph parents, or every variable earlier in the order. Truth: the pipeline's
per-unit contrast for the population each data set represents.

## Requirements
- The LumaWarp engine binaries (`lumawarp_training`, `lumawarp_inference`; Lucidity Sciences, licensed). Scripts look
  for them on PATH or in `LUMAWARP_TRAIN` / `LUMAWARP_INFER`. The shape flag reads format-1 model files, which the
  September build writes and the October build does not (`shape_flag_checks/RESULTS_importance_file_flag.md`).
- Python with numpy, pandas, scipy, scikit-learn, xgboost, shap, causal-learn, networkx, pyyaml, threadpoolctl
  (`PYTHON` selects the interpreter; default `python3`). `attribution_lumawarp_base/` imports the pipeline's own
  `04_attribution_validation` package for the split, records, orders, SCM and scoring.
- R with tmle, SuperLearner, glmnet, earth.
- `DATA_DIR` as in `.env` (default `02_data`); the analysis graph is read from `outputs/03_shape_validation/`.

## Running
    cd 06_lumawarp_audit/edge_audit && python make_jobs.py && ./run_jobs.sh && python references.py && python score.py
    cd 06_lumawarp_audit/estimation && PLAN=full ./run_all.sh          # engine driver plus three R workers; then python summarize.py
    cd 06_lumawarp_audit/estimation_bounds && ./run_all.sh && python summarize.py
    cd 06_lumawarp_audit/attribution_lumawarp_base && BACKGROUND=32 ./run_all.sh && python compare.py
    nonlinear_sensitivity: generate data on branch nonmonotone-variant (make recalibrate), then ./run_after_data.sh
Each run script is resumable. Engine time on one RTX 5090 laptop: edge audit about 12 h (two streams), estimation
about 6 h, attribution about 4 h at 32 background records (about 24 h at the paper's 128).

## Outputs committed (`outputs/06_lumawarp_audit/`)
`edge_audit/` summary.csv, scores_by_rep.csv, shape_flags.csv, edges_all.csv.gz (every candidate score and decision);
`estimation/` summary.csv, results_all_reps.csv; `estimation_bounds/` summary.csv, results_all.csv;
`attribution_lumawarp_base/` summary.csv, per_run.csv, importance.csv; `nonlinear_sensitivity/armA|armB/`
hybrid_summary.csv, audit_summary.csv, hybrid_by_rep.csv, gaps.csv; `shape_flag_checks/` the flag's operating-range
map and the importance-file comparison; `tables/` the manuscript tables. Per-record Shapley values and fitted models
are not committed (available on request).

## Deviations recorded
- Attribution ran with 32 background records (the paper: 128); explained records and feature orders are the paper's.
  Causal predictive SHAP was not ported.
- The edge audit's full-set stage ran as two engine streams with a claim file; it was paused and resumed once.
- Replicates: edge audit and attribution 1 to 10; estimation 1 to 50.
