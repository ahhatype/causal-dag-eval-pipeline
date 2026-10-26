# Pre-registration: the paper's attribution arm with LumaWarp as the base model (written 19 October 2026, before any fit)

## Question
The draft's attribution arm (standard, ordering-only, causal interventional and causal predictive SHAP) explains a
Super Learner fitted to each data set. LumaWarp appears nowhere in that arm. If LumaWarp replaces the Super Learner as
the model being explained, with everything else held fixed (the same split, the same 64 explained records, the same
128 background records, the same sampled feature orders, the same truth and tie rules), do the rankings of the 16
ancestors agree with the truth as well as the paper's, and does the draft's main contrast (interventional minus
standard) reproduce?

## Data (fixed)
The regenerated pipeline data (commit 79a1b56), replicates 1 to 10 of full-set, reference-set and astronaut-set,
ancestor feature set (16 inputs). Seeds, split, explained records and background exactly as the pipeline's
`prepare_run` draws them (the pipeline's own code is imported, not reimplemented).

## Base model
LumaWarp old build, research configuration (8 members, 7,200 s budget, axis dumps and report off), binary task,
fitted on the pipeline's 70 percent training split of each replicate. Predictions on the probability scale (p_1).

## Methods
Standard SHAP (permutation, 64 antithetic orders over the 128 background records), ordering-only SHAP and causal
interventional SHAP (64 graph-consistent orders from the pipeline's OrderSampler with the pipeline's seed; the
interventional value function uses the pipeline's fitted linear-logistic SCM on the training split). Each method is
run exactly as in the pipeline except that every row the value function needs is generated first and predicted in
one engine call per 2 million rows, since the engine is a command-line program. Causal predictive SHAP is not run
(its rescaling and per-record samplers would need a further port; see Deviations if added later).

## Measures (the paper's)
Kendall tau-b against the intervention-contrast truth (each data set's own population, mission era included, primary
tie rule delta2), top-5 recovery, proximity bias index, plus the AUC of the explained model. Paired with the
pipeline's committed per-replicate Super Learner results for the identical cells.

## Bars
- B1 (interchangeable base model): on the full-set, causal interventional SHAP with LumaWarp has mean tau within
  0.05 of the paper's Super Learner value (paired over the 10 replicates), and the LumaWarp model's AUC is within
  0.01 of the Super Learner's.
- B2 (the paper's contrast reproduces): interventional minus standard tau with LumaWarp as base is at least 0.30 on
  the full-set (the paper reports 0.43 with the Super Learner).
- Reported: the two 900-row sets; ordering-only minus standard (the paper's inconclusive 0.01 contrast).

## What each outcome lets us say
- B1 and B2 pass: the attribution findings do not depend on the explained model; LumaWarp can take the Super
  Learner's place in Table 2 with the same conclusions. A LumaWarp row for the attribution step.
- B1 fails with LumaWarp lower: the ranking depends on the fitted model's shape; report and do not tout.
- B1 fails with LumaWarp higher: a lead, to be checked at 50 replicates before any claim.
