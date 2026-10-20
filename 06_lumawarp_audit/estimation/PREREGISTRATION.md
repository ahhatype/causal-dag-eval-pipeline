# Pre-registration: estimating the two confounded drug contrasts on the paper's data sets (written 14 October 2026, before any fit)

## Question
The draft (13 October) calls the intervention contrast "the quantity a countermeasure decision needs" and reports a
crude K+ citrate contrast of +0.031 against a true -0.020, but no method in the draft estimates the contrast. Does
adjustment by the DAG recover it, with a calibrated interval, on the paper's own three data sets, and does the
answer depend on the nuisance learner or on how much of the graph is used?

## Data (fixed)
The 50 replicates of the three data sets regenerated from the published pipeline (`ahhatype/causal-dag-eval-pipeline`,
commit 79a1b56; generator code and configuration identical to the local checkout). Treatments: K+ citrate (truth
-0.0202 source, -0.0152 selected) and bisphosphonates (-0.0083, -0.0063). Truth: the pipeline's per-unit
(0 vs 1) intervention contrast, source population for the full-set and reference-set, selected population for the
astronaut-set (its "own" truth, as in the draft).

## Adjustment routes
- DAG parents: K+ citrate {individual factors, mission era, urine concentration}; bisphosphonates {cumulative
  mission duration, mission era, pre-flight fitness}. The draft's own confounder sets.
- Ordering only: every analysis-graph variable earlier than the drug in the DAG's topological order (no graph beyond
  the order).

## Estimators (per replicate, as folders 25 and 32)
Three stratified folds; LumaWarp (playbook profile) and Super Learner (glm, glmnet, earth) nuisance fits for g and Q;
one `tmle()` call per learner (binomial, logistic fluctuation, gbound 0.025 to 0.975); logistic g-computation and
the naive difference alongside.

## Plan
Stage 1 (for the 14 October team meeting): K+ citrate, reference-set and astronaut-set, both routes, replicates 1 to
10. Stage 2: all treatments, data sets and routes, replicates 1 to 50.

## Measures and bars (stage 2, K+ citrate, full-set, DAG parents, primary)
Bias, empirical and mean SE, coverage. P1: |bias| of LW-TMLE <= 0.004. P2: coverage >= 0.90. P3: parity with
SL-TMLE (bias gap <= 0.002, coverage gap <= 0.05). Reported: the 900-row sets (coverage expected to hold; sign
recovery rate reported), the ordering route, bisphosphonates (weak confounding by design), the naive contrast.

## What each outcome lets us say
- P1 to P3 pass: the paper's countermeasure framing gets its payoff: the DAG's adjustment set turns the sign-flipped
  crude contrast into an unbiased estimate with a valid interval, on the paper's own data, with the learner
  interchangeable. Report beside ConfoundingSHAP as the second half of the confounding analysis.
- Coverage fails at 900 rows: the sample-size headline extends to estimation.
