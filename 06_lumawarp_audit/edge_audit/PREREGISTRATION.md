# Pre-registration: ordering-informed edge audit on the paper's data sets (written 14 October 2026, before any fit)

## Question
The draft's placeholder asks for "exploratory shape analysis with the commercial learned-kernel classifier, reported
against the same three data sets and the measures in Methods". Given the established DAG's ordering (not its edges),
how many of the analysis graph's 63 edges does a flexible per-node audit recover, at what precision, on the three
data sets, and does the engine do anything a spline or logistic model does not?

## Data (fixed)
The regenerated replicates 1 to 10 of full-set, reference-set and astronaut-set (published pipeline, commit 79a1b56).
Analysis graph: 36 variables, 63 edges (the pipeline's `analysis_graph_edges.csv`), topological order from it.

## Method (folder 13's rule, unchanged)
For every analysis-graph variable with at least one predecessor, fit it on all earlier variables in the order plus
12 pure-noise columns (6 continuous, 6 yes/no). Sweep each candidate over its range with the others at their
median; edge score = sweep span / largest noise-column span. Keep an edge if the score exceeds the learner's
pre-registered threshold: LumaWarp 5x, logistic/linear 4x, spline additive 2x, XGBoost 3x (folder 13's
thresholds, tuned on other seeds and never changed). Orientation is given by the supplied order.
Shape flag (folder 14's rule) read from every saved model: continuous parents of continuous targets flagged if the
centroid beats every continuous noise column. The generator is linear, so the correct count is zero.

## Measures (per replicate and data set, then mean, with Monte Carlo SE; the draft's own measures)
Adjacency recall and precision against the 63 edges; edges called; noise columns kept (false-keep rate); the share
of the 63 edges recovered by edge type; shape flags raised.

## Bars
- E1: LumaWarp recall on the full-set >= 0.80 at precision >= 0.90 (the draft's MGM-PC: 0.78 at 1.00 with no ordering).
- E2: shape flags on the full-set <= 5% of continuous-parent cases.
- E3 (reported): LumaWarp within 0.10 recall of the best reference at each data set.
## What each outcome lets us say
- E1, E2 pass and E3 "matches": the placeholder reads "given the DAG's ordering, a per-node flexible audit recovers
  X% of edges at Y% precision, the same with a spline model; the learned-kernel model's monotonicity reading fired
  zero times, as it should on this generator". Honest, and it fills the slot.
