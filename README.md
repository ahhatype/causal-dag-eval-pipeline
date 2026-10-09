# causal-dag-eval-pipeline

## Overview

<!-- TO BE WRITTEN BY THE AUTHOR: purpose of the project, study design, and how the pieces fit together. -->

_(Overview to be written.)_

---

## Setup

R is needed for data generation and shape validation. Package versions (`simcausal`, `rCausalMGM`, `data.table`, `yaml`, `dagitty`, `testthat`) are pinned in `renv.lock`. Attribution validation is Python, managed with [uv](https://docs.astral.sh/uv/) and pinned in `04_attribution_validation/uv.lock`.

```bash
Rscript -e 'renv::restore()'   # install the pinned R packages
(cd 04_attribution_validation && uv sync)   # install the pinned Python packages
cp .env.example .env           # then add your API key; .env is git-ignored
make data                      # generate the data sets
make shape                     # structural validation on the generated data
make attribution               # attribution validation on the generated data
make test                      # run all tests
```

`DATA_DIR` in `.env` (default `./02_data`) sets where generated files are written; the full data are git-ignored (see [Data and results availability](#data-and-results-availability)).

## Data generation

All code lives in [01_data_generation/](01_data_generation/). Three synthetic data sets are drawn from one prespecified DAG with [simcausal](https://cran.r-project.org/package=simcausal), so the true causal structure and total effects are known. The DAG is the NASA Human System Risk Board renal-stone DAG (SA-07566; 51 nodes, 75 edges) with three added nodes (**cumulative mission duration** as the exposure, **mission era** as a root, **pre-flight fitness** as a root), giving 54 nodes and 90 edges. The outcome is **nephrolithiasis** (binary). NASA's DAG specifies structure only: all equations, coefficients, era values and detection parameters here are study design parameters, not estimates of real-world effects.

The data are generated with four further edges that are not in SA-07566, so that two drugs are confounded and confounding attribution can be evaluated: **individual factors** and **urine concentration** → **K+ citrate** (prescription by indication), and **cumulative mission duration** and **pre-flight fitness** → **bisphosphonates** (bone-loss countermeasure). Both drugs are also taken in the 1960s (5%; 15% in the 2000s), which is not historical, so that every record has treated and untreated counterparts. These edges are part of the true model and of the analysis graph.

| | Full-set | Reference subsample | Astronaut-set |
|---|---|---|---|
| n | 10,000 | 900 | 900 |
| Source | Simulated from the DAG | Random subsample of the full-set | The same replicate's full-set, selected by sampling |
| Selection | None | None | Retained with probability expit(α − 1.0·individual factors + 1.0·pre-flight fitness), α set for ~30% retention |
| Era mix | 25% 1960s / 75% 2000s | as full-set | not era-dependent |
| Seed (replicate 1) | 20261004 | 20261007 | full-set's |
| Role | Ideal data, no selection | Separates sample-size from selection effects | Small, selected, astronaut-like cohort |

Each data set is generated in 50 independent replicates. Replicate r uses the seeds above plus 1000·(r − 1). Ground truth is a population quantity, computed once with seed 20261006. Output files are the same set of 46 observed columns in every data set (7 latent nodes withheld; latent values go to separate `*_latent.csv` debug files).

### How it works

1. **Config** ([01_data_generation/config/](01_data_generation/config/)): the single source of truth.
   - `nodes.csv`, `edges.csv`: node roles/scales and every edge with its type and coefficient. Built from the DAG appendix (Supplementary Tables 1–2) by `tools/build_config_tables.py`.
   - `simulation_edges.csv`: the four edges into K+ citrate and bisphosphonates that are not in SA-07566 (moderate tier, ±0.8 on the log-odds scale).
   - `params.yaml`: seeds, sample sizes, era values, uptake probabilities, detection sensitivity/specificity, selection parameters, base rates. Entries marked PROVISIONAL are pending domain-expert confirmation.
   - `dag.txt`: the adapted DAGitty code used to validate the config.
2. **Validation** (`R/validate_dag.R`): `edges.csv` must reproduce `dag.txt` exactly: 54 nodes, 90 edges, same edge set, acyclic, 16 outcome ancestors, 19 other observed features, 7 latent nodes. The simulation edges add no node and change no ancestor set.
3. **Calibration** (`R/calibrate.R`): a pilot run (n = 200,000) freezes, in topological order, each continuous node's standardization and residual SD, the logit intercept of each binary node (e.g. 10% nephrolithiasis prevalence), the era-specific uptake intercepts of the two confounded drugs, and the selection intercept α. Results are stored in `config/calibration.yaml` with a hash of the config files. Generation stops if the config has changed since calibration; re-run with `make recalibrate`.
4. **Model** (`R/model_spec.R`, `R/dag.R`): each node becomes a formula built from `edges.csv` and `params.yaml`.
   - Mission era is Bernoulli (0.75 for the 2000s). The 10 observed and 2 latent mission-level (context) nodes are fixed by era, so edges among them are kept in the graph but not separately parameterized. Cumulative mission duration is drawn from truncated era-specific distributions (1960s; 2000s Shuttle/ISS mixture).
   - Continuous children: linear combination of parents plus noise, in SD units of the child. Continuous parents enter standardized; binary parents (drugs, procedures) enter as 0/1, so their coefficient is the shift for users vs non-users. Binary children: logistic. Detection nodes use sensitivity/specificity; near-deterministic edges use a 5% flip.
   - Every node has its own uniform noise node, so interventions change a node's formula without shifting any random draws (**common random numbers**).
5. **Data sets** (`R/pipeline.R`): in each replicate:
   - The full-set is simulated.
   - The reference subsample is a random 900 of that replicate's full-set.
   - The astronaut-set applies the selection rule to that replicate's 10,000 full-set records (about 3,000 retained) and keeps the first 900 retained (`R/select.R`), so astronaut and full-set comparisons are paired within a replicate.
6. **Ground truth** (`R/ground_truth.R`): the importance of each of the 35 features for nephrolithiasis, computed from the true model on 500,000 records. The model has no closed form for these quantities (logistic outcome, calibrated constants, mixture and threshold nodes), so they are evaluated by intervention with common random numbers; the only error is Monte Carlo error, reported for every value.
   - The node formulas are evaluated directly (the strings simcausal uses), so an intervention re-evaluates only the target's descendants; a test checks this reproduces simcausal exactly.
   - **Per-unit** (`per_unit`): risk difference in the outcome's expected probability between the 75th and 25th percentile of the source population (1 vs 0 for binary features).
   - **Population-scaled** (`pop`): within each era e, μ(x, e) is the mean outcome probability under do(X = x); the truth is Σ_e P(e) · E_{x∼p(x|e)} |μ(x, e) − m_e|, with m_e = E_{x∼p(x|e)} μ(x, e). Binary features are weighted by within-era prevalence; continuous features use the midpoints of 20 within-era equal-probability bins. For mission era itself there is no within-era variation, and the spread is taken across eras.
   - **Record-level** (`rec`, sensitivity): the same centre computed record by record with each record's own noise, then mean |P_j(observed) − centre_j|, mirroring how mean |SHAP| takes absolute values per record.
   - Each is computed for the **source** population (all records) and the **selected** population (records retained at baseline, membership held fixed, then intervened on).
   - Standard errors: paired across records for per-unit and record-level; batch means over 20 batches for population-scaled.
   - For binary features, the within-era population-scaled truth is compared with 2p(1 − p) times the observational contrast; a gap beyond 3 MCSE flags confounding (`truth_binary_checks.csv`). K+ citrate and bisphosphonates are flagged; thiazides is not.
7. **Checks** (`R/checks.R`):
   - Every replicate has the right sizes.
   - Prevalence, era share and retention are within tolerance in every replicate, and prevalence is on target averaged over replicates.
   - Selection induces the individual factors–fitness association on average.
   - 35 features are present and no latent columns are written.
   - Non-ancestors have exactly zero truth of every type. A failed check stops generation. The exception is the ancestor ranking: adjacent per-unit truths within 2 standard errors of each other are reported as a warning (near-tie).

### Run

```bash
make data                                   # or: Rscript 01_data_generation/generate_all.R
```

Takes about 1.5 minutes (ground truth about 1 minute). Outputs in `02_data/` (about 360 MB):

| File | Contents |
|---|---|
| `replicates/rNNN/full_set.csv`, `reference_subsample.csv`, `astronaut_set.csv` | Observed columns (46 incl. outcome), one folder per replicate |
| `replicates/rNNN/*_latent.csv` | Latent nodes, for debugging only |
| `truth_values.csv` | Long table: feature, truth_type (per_unit, pop, rec), population (source, selected), value, mcse; low/high contrast values for per-unit |
| `truth_binary_checks.csv` | Per binary feature, population and era: prevalence, population-scaled truth, its observational counterpart and the confounding flag |
| `scm_check.csv` | 500 simcausal records with their noise draws, for checking other copies of the equations |
| `feature_sets.csv` | The 35-feature (all) and 16-feature (ancestor) sets |
| `data_summary.csv` | Per replicate and data set: n, seed, era share, realized prevalence, individual factors–fitness correlation, retention |
| `provenance.txt` | Git commit, config hash, replicate count, R and package versions |

### Tests

`make test` runs `01_data_generation/tests/`: DAG fidelity, topological order, seed reproducibility, common-random-number behavior, calibration targets and staleness, duration truncation, confounded drug uptake, astronaut selection from the full-set, direct formula evaluation against simcausal, the population-scaled spread (2p(1 − p)|Δ| for binary features), ground truth, the checks themselves, an end-to-end run with two replicates at reduced size, and (when `references/` is present locally) that the config tables rebuild exactly from the appendix.

## Shape validation

Code in [03_shape_validation/](03_shape_validation/) asks whether the structure of the DAG is consistent with the data. It tests the DAG's conditional independencies directly, and learns a structure with causal discovery and compares the two. Every analysis is repeated on each of the 50 replicates of each data set, and results are reported as rates or means with Monte Carlo standard errors.

### How it works

1. **Analysis graph** (`R/analysis_graph.R`): the DAG restricted to observed variables.
   - The ten observed mission-level nodes vary only by era, so they are collapsed into mission era, and their edges are redirected to it.
   - Edges from the two mission-level nodes that are equal in both eras (altered gravity, humidity) are dropped.
   - Latent nodes are removed.
   - Result: 36 variables and 63 edges (59 from SA-07566 plus the four simulation edges).
   - Its equivalence class (CPDAG: 51 compelled, 12 reversible edges) is the reference for scoring orientations.
2. **Testable implications** (`R/implications.R`): Shipley's d-separation basis set, one implication per non-adjacent pair, X ⊥ Y | pa(X) ∪ pa(Y). That gives 567 implications.
   - X is the later of the two in topological order.
   - Each implication is tested with a likelihood-ratio test of Y in the regression of X on Z (linear for continuous X, logistic for binary X; continuous variables standardized). The effect is Y's coefficient per SD.
   - P values are adjusted with Benjamini–Hochberg at 0.05.
   - Calibration: in each replicate, the p values of the 566 implications not expected to fail under selection are compared with a uniform distribution (Kolmogorov–Smirnov), and their raw rejection rate at 0.05 is recorded, over all 566 and over those with an additive response.
   - Implications are flagged when X's generating mechanism is not additive linear/logistic, and when they are expected to fail under the astronaut-set's selection (X and Y d-connected given Z and a selection node whose parents are individual factors and pre-flight fitness). Exactly one is: individual factors ⊥ pre-flight fitness.
3. **Power** (`R/power.R`): a positive control.
   - Each true edge is deleted in turn, and the single false independence this creates is tested (Bonferroni over the 63 edges).
   - The detection rate per edge shows which parts of the graph the implication tests can check at each sample size.
4. **Causal discovery** (`R/discovery.R`): MGM-PC from rCausalMGM, with no background knowledge.
   - The MGM skeleton's sparsity penalties λ are chosen per data set by StEPS.
   - PC-Stable is run on that skeleton with four collider-orientation rules: MGM-PC-Stable (sepsets; primary), MGM-CPC-Stable (conservative), MGM-MPC-Stable (majority) and MGM-PC-Max-Stable (max-p).
   - Binary variables with fewer than 5 records in a category, which rCausalMGM cannot fit, are left out of discovery for that data set and recorded.
   - Edge stability: the primary variant is rerun on 100 subsamples of floor(0.632·n) records drawn without replacement (first replicate only).
5. **Comparison** (`R/compare.R`, `R/dsep.R`): each learned CPDAG is extended to a DAG.
   - A consistent extension is found with the Dor–Tarsi algorithm. All consistent extensions imply the same independencies.
   - If none exists (conflicting orientations), the comparison is summarized as the median, min and max over acyclic orientations of the undirected edges (enumerated, or sampled beyond 64).
   - Each graph's implications are checked for d-separation in the other. An analysis-graph implication that fails in the learned graph points to an added edge; a learned implication that fails in the analysis graph points to an omitted edge.
   - Adjacency precision/recall and orientation agreement with the true CPDAG are reported alongside.
6. **Checks** (`R/pipeline.R`): the run stops if the data were generated from a different config than the current calibration (stale data), or if the analysis graph does not have 36 nodes, 63 edges, 567 implications and one selection-sensitive implication.

Settings (FDR level, power α, PC-Stable α, variants, StEPS subsamples, stability subsamples, seeds, flagged nodes) are in `config/params.yaml`.

### Run

```bash
make shape                                  # or: Rscript 03_shape_validation/run_all.R
Rscript 03_shape_validation/run_all.R --resume                       # continue an interrupted run
Rscript 03_shape_validation/run_all.R --fresh                        # discard saved units and start again
Rscript 03_shape_validation/run_all.R --replicates=2 --subsamples=0   # quick check
Rscript 03_shape_validation/run_calibration.R                         # implication tests and calibration only
```

The full run (50 replicates × 3 data sets × 4 variants, with StEPS) takes about 2–3 hours when run alone, plus about 30 minutes for the stability analysis on replicate 1; the implication tests alone take a few minutes. Each replicate × data set is saved in `units/` when complete. A later run refuses to touch them unless given `--resume` (continue, only if settings, config and code are unchanged) or `--fresh` (discard). **Reproducibility:** rCausalMGM's `pcStable` has no seed, so MGM-PC-Stable (sepset rule) keeps its skeleton on repeat runs but can orient colliders differently; the other three variants repeat exactly. It reads `DATA_DIR` and writes to `OUTPUT_DIR/03_shape_validation/` (default `./outputs`, git-ignored):

| File | Contents |
|---|---|
| `summary.csv` | Per data set and variant: rates (selection violation detected, any unexpected rejection, comparison available) and mean ± MC SE of every per-replicate metric |
| `per_replicate.csv` | One row per replicate, data set and variant: rejections by flag, detectable edges, λ, variables left out of discovery, structure metrics, implication comparison |
| `implication_rejection_rates.csv` | Each of the 567 implications with its flags and rejection rate per data set |
| `power_by_edge.csv` | Per edge and data set: detection rate and median p in the positive control |
| `edge_recovery.csv` | Per data set and variant: how often each true edge is found and oriented as in the true CPDAG, and how often each false edge appears |
| `implication_calibration.csv`, `implication_calibration_summary.csv` | Per replicate and data set: KS p value and raw rejection rate of implications expected to hold; per data set: share of replicates rejecting uniformity and mean raw rejection rate |
| `edge_stability_r001_<dataset>.csv` | Subsampling stability of the primary variant, replicate 1 |
| `analysis_graph_edges.csv`, `analysis_graph_cpdag.csv` | Analysis graph (with the source edges behind each edge) and its CPDAG |
| `detail/` | Implication tests for every replicate; learned CPDAGs per variant for replicate 1 |
| `provenance.txt` | Git commit, data config hash, replicates analysed, R and package versions (with the R version each package was built under); `implication_calibration_provenance.txt` likewise for `run_calibration.R` |

### Tests

`make test` also runs `03_shape_validation/tests/`:
- **Settings:** they parse to the expected types.
- **Analysis graph:** size and collapsing.
- **Basis set:** 567 implications, all implied by the graph.
- **d-separation:** agreement with dagitty.
- **Selection labelling.**
- **Likelihood-ratio test:** on planted dependence and independence.
- **Calibration:** uniform p values pass, skewed ones fail, the selection implication is excluded.
- **Comparison:** detects an added and an omitted edge.
- **Consistent extension:** found, or correctly absent for an unchorded 4-cycle; when absent, the comparison reports a range and is deterministic.
- **Orientation scoring:** against the true CPDAG.
- **Positive control.**
- **Stale data:** refused.
- **Discovery wrapper:** recovers a small mixed-data chain.
- **End to end:** a run on two generated replicates.

## Attribution validation

Code in [04_attribution_validation/](04_attribution_validation/) asks whether attribution methods recover the relative importance of the DAG's nodes. Each method's ranking of features is compared with the ground truth from data generation, on every replicate of every data set.

### How it works

1. **Inputs** (`data.py`): the replicate data sets and ground truth from `DATA_DIR`. The analysis graph is built from `01_data_generation/config` with the same rules as shape validation; a test checks the two are identical. The run stops if the data were generated from a different config.
2. **Feature sets:** all 35 observed features, and the 16 ancestors of nephrolithiasis. The all-features set adds 15 outcome descendants and 4 variables that are neither ancestors nor descendants, all with a true effect of zero.
3. **Predictive model** (`superlearner.py`): one super learner per data set and feature set.
   - Library: main-terms logistic regression, L2-penalized logistic regression, random forest, gradient boosting.
   - Weights: a convex combination chosen to minimize 10-fold cross-validated log loss, then the learners are refit on all training data.
   - Records are split 70/30, stratified by outcome. Held-out AUC is reported next to the AUC of the true outcome probability. That bounds what a model can reach on the ancestor set; with all features, outcome descendants let a model exceed it.
4. **Explained records:** within a replicate and data set, both feature sets use the same split, and every method explains the same 64 held-out records, drawn at random, against the same 128 training records as background. Records are not swapped to guarantee outcome events; runs with fewer than 4 events among them are flagged (`explained_events_low`). All attributions are on the probability scale. Standard, ordering-only and interventional attributions are checked to sum to the model's prediction for each record; Ng-style attributions are rescaled to do so.
5. **Methods:**
   - **Standard SHAP** (`standard.py`): shap's permutation explainer, 64 feature orders (32 antithetic pairs), background features filled in independently of the graph.
   - **Ordering-only asymmetric SHAP** (`ordering.py`), a marginal value function with causal order: the same graph-consistent orders as interventional SHAP, but features outside a coalition are taken from the background sample, as in standard SHAP. Comparing it with the other two separates the effect of respecting the graph's order from the effect of propagating interventions.
   - **Ng-style causal SHAP** (`ng.py`), after Ng et al. (arXiv:2509.00846), with the analysis graph supplied in place of PC + IDA.
     - Each edge's strength is the |coefficient| of the parent in one linear regression of the child on all its parents.
     - A feature's causal weight is the normalized sum, over its directed paths to the outcome, of the product of edge strengths. Non-ancestors get zero.
     - The causal value of a coalition averages the model over 64 draws of the other features in graph order: roots from their empirical distribution, continuous features from a linear and binary features from a logistic regression on their feature parents.
     - Shapley values of this value function are estimated from 64 feature orders (32 antithetic pairs), with the same draws shared by every coalition along an order. Each is multiplied by the feature's causal weight, then all are rescaled to sum to f(x) − E f(X).
     - The weighted sum the rescaling divides by can be near zero. Each record's pre-rescaling total, rescaling factor and attributions before and after rescaling are saved. Importance excludes records whose |pre-rescaling total| is below 0.001; importance over all records is reported as `ng_style_shap_all_records`.
   - **Interventional asymmetric SHAP** (`interventional.py`): a coalition is valued by the model's mean prediction when its features are set by intervention.
     - A structural causal model is fitted to the training data: each node regressed on its analysis-graph parents, linear for continuous and logistic for binary nodes. It includes non-feature nodes such as the outcome, so interventions propagate to its descendants.
     - Background records are abducted to exogenous noise; binary noise is drawn uniformly from the range consistent with the observed value.
     - 64 feature orders are sampled uniformly from the orders consistent with ancestry in the analysis graph, using a Markov chain of adjacent swaps.
6. **Oracle reference** (`truesim.py`): on the ancestor set, for replicates 1–10 of every data set, each method is rerun with the true model in place of the fitted one, on the same explained records.
   - The predictor is the true outcome probability, P(nephrolithiasis | individual factors, mineralized renal material); both parents are observed members of the ancestor set (checked before the run).
   - Standard and ordering-only SHAP use 128 background records drawn from 20,000 simulated from the true model. Ng-style SHAP draws out-of-coalition features from the true structural equations, with path weights from that simulated population. Interventional SHAP uses the true structural equations and the true noise of those background records.
   - `truesim.py` is a Python copy of the generating equations; a test checks it reproduces simcausal's records node by node.
   - This separates the estimand gap (oracle against truth) from estimation error (fitted against oracle).
7. **Evaluation** (`evaluation.py`): mean |SHAP| over the explained records gives each feature's importance. On the 16 ancestors in every run, against every truth type (per-unit, population-scaled, record-level) and both populations (source, selected), with mission era included (primary) and excluded:
   - Kendall's τ_b with |truth|. Primary measures: τ_b against the per-unit and population-scaled truths (co-primary), for each data set's own population (source for the full-set and reference subsample, selected for the astronaut-set).
   - Ties: truths closer than 2% of the largest |truth| are tied (chained between neighbours; primary). Sensitivity: 5%, and ties within 1.96 combined Monte Carlo standard errors. Tied features share one value in τ_b; tie groups under every rule are written to `tie_groups.csv`.
   - Secondary: weighted τ (Vigna 2015, hyperbolic weights, scipy's `weightedtau`), and the L1 distance between importance and |truth|, each scaled to sum to 1 (0 to 2).
   - Descriptive: top-5 recovery (every feature tied with the fifth counts as a true top-5 feature).
   - The proximity bias index: importance-weighted mean distance to the outcome under the truth minus under the method. Positive means credit is pooled near the outcome.
   - All-features runs also report the non-ancestor credit share.
   - Paired contrasts within replicate: the order effect τ_b(ordering-only) − τ_b(standard) and the value-function effect τ_b(interventional) − τ_b(ordering-only). A contrast is **robust** if it has the same sign and a 95% interval excluding zero under both co-primary truths, otherwise **truth-dependent**.
8. **Ng-style convergence** (`run_ng_convergence.py`): on replicate 1 of each data set and feature set, two independent halves at 150 and 300 orders are compared by the correlation of per-feature mean |φ|, its largest difference and τ_b (criterion: τ_b ≥ 0.9 and correlation ≥ 0.95). If 300 fails, the budget is doubled until it passes or a run exceeds 4 times the production runtime.
9. **ConfoundingSHAP** (Brockschmidt et al. 2026): not yet run. The confounded drugs in the data-generating process are in place; settings are reserved in `config/params.yaml` (`confounding_shap`, disabled). <!-- TO BE ADDED: ConfoundingSHAP method and results -->

Settings (library folds, split, explained records, background size, permutations, Monte Carlo budgets, order sampler, seeds) are in `config/params.yaml`.

### Run

```bash
make attribution                            # or: cd 04_attribution_validation && uv run python run_all.py
uv run python run_all.py --resume           # continue an interrupted run
uv run python run_all.py --replicates 1 --datasets reference_subsample --feature-sets ancestor --workers 1   # one cell
uv run python run_ng_convergence.py         # Ng-style convergence study
```

Each replicate × data set × feature set is one cell, seeded on its own and run single-threaded, so results do not depend on the number of workers (7 by default). Finished cells are kept in `cells/`. A later run refuses to touch them unless given `--resume` (continue, only if settings and code are unchanged) or `--fresh` (discard). A failed cell is logged and the others continue; outputs are combined once every cell is complete. On 7 cores the full design (300 cells) takes about 5–7 hours; the slowest cell (full-set, all features) about 15 minutes. Outputs go to `OUTPUT_DIR/04_attribution_validation/`:

| File | Contents |
|---|---|
| `summary.csv` | Per data set, feature set, method, truth type, population and scope: mean and MC SE of each metric, AUCs and runtime; `primary` marks each data set's own population with mission era included |
| `per_run.csv` | Long format, one row per replicate, data set, feature set, method, truth type, population and scope: metrics, seed, super learner weights and cross-validated log loss per learner, AUCs, efficiency error, explained events, runtime; for Ng-style, the share of records excluded and the rescaling factor's median, maximum and share above 50 |
| `contrasts.csv` | Paired order and value-function effects, mean, MC SE, 95% interval and robust/truth-dependent classification |
| `oracle_decomposition.csv` | Per method, replicate (1–10), data set, truth type and population: τ_b of oracle and fitted against truth, and fitted against oracle (τ_b and L1) |
| `tie_groups.csv` | Tie groups of the truth under each rule |
| `ng_convergence.csv` | Ng-style convergence study |
| `importance.csv` | Each feature's importance in every run |
| `values.csv.gz` | Signed attribution of every feature for every explained record and run (oracle runs as `<method>_oracle`), with the record's outcome; for Ng-style, each record's Shapley total, weighted (pre-rescaling) total, rescaling factor and pre-rescaling attributions (`pre_<feature>`) |
| `provenance.txt` | Git commit, code hash, data config hash, package versions and settings |
| `cells/` | Per-cell results and the run's settings manifest |

### Tests

`make test` also runs `04_attribution_validation/tests/`:
- **Analysis graph:** matches the design and shape validation.
- **Order sampler:** uniform over consistent orders.
- **Ancestry constraints:** pass through non-feature nodes.
- **Structural model:** abduction and simulation round-trip, coefficients recovered.
- **Ordering-only SHAP:** efficiency, and no credit to an ancestor the model does not use.
- **Closed form:** on a chain x → m with a model that uses only m, interventional SHAP credits x with 0.8 (x − E x) and ordering-only SHAP credits it with zero.
- **Tied effects:** treated as tied in τ_b and top-k recovery; margin ties chain between neighbours, MCSE ties use the combined error.
- **Weighted τ:** disagreements at the top cost more than at the bottom. **L1 distance:** compares shares.
- **True model copy:** reproduces simcausal's records node by node; the true sampler fixes a coalition and propagates only downstream.
- **Paired contrasts:** robust only when both co-primary truths agree.
- **Interventional SHAP:** credits an ancestor acting through a model feature, gives zero to an irrelevant feature, and is efficient.
- **Ng weights:** zero for non-ancestors, local accuracy; pre-rescaling attributions sum to the weighted total and rescale to the reported values.
- **Shapley estimator:** averaging over all orders reproduces exact Shapley values.
- **Ng sampler:** binary features stay binary, draws are shared across coalitions, fixing a parent shifts its child.
- **Permutation SHAP:** efficiency.
- **Super learner:** weights form a convex combination.
- **Metrics.**
- **Stale data:** refused.
- **End to end:** one run on generated data with small budgets, including the oracle.
- **Parallel runs:** identical to sequential.
- **Cell directory:** never reused with changed settings or code, never cleared without `--fresh`.

## Data and results availability

**Committed.** Code, configuration, seeds and the frozen calibration (`01_data_generation/config`), plus the small files that are the source data for the paper's tables and figures:
- `02_data/`: `data_summary.csv` (realized prevalence, era share and selection-induced correlation per replicate), `feature_sets.csv`, `truth_values.csv`, `truth_binary_checks.csv` and `provenance.txt`.
- `outputs/03_shape_validation/` and `outputs/04_attribution_validation/`: the summary, per-replicate and per-run tables, and per-edge and per-implication results. Each folder's `provenance.txt` records the code commit and data config hash.

**Archived, not committed.** The 50 replicate data sets per data set (about 360 MB as CSV), the per-cell results and the per-record attributions (`values.csv.gz`) are deposited in a DOI-minting repository with the tagged release: <!-- TO BE WRITTEN BY THE AUTHOR: DOI and repository at submission -->.

**Regenerating.** The data are fully determined by the config and seeds; `make data` rebuilds them and every analysis refuses to run on data from a different config hash. Check `02_data/provenance.txt` against the archived copy before comparing results across machines, because random number generation can differ across platforms.

To prepare the archive from a finished run:

```bash
tar -czf causal-dag-eval-data.tar.gz 02_data outputs
```

---

## License

MIT, see [LICENSE](LICENSE).
