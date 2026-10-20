# Pre-registration: a linear audit passes a wrong graph, a nonlinear audit rejects it (written 13 October 2026)

## Question
The paper's audit step tests a graph's d-separation implications (Shipley basis set) with linear and logistic
likelihood-ratio tests, as in Andy's pipeline (`03_shape_validation/R/implications.R`). If a graph is missing a
non-monotone edge, does that audit notice? Does a nonlinear test?

## Data (fixed)
nl_v2 (`sims_v2.py`), n = 10,000, seeds 400 to 407. 51 nodes, 75 edges; the four non-monotone edges are
Nutrients -> Bone Remodeling (inverted U), Hydration -> Urine Concentration (U), Microbiome -> Urine Chemistry (U),
Urine Chemistry -> Mineralized Renal Material (wave).

## Graphs audited
- (a) the true graph;
- (b) the true graph with the four non-monotone edges deleted: a wrong graph that a linear-only pipeline would
  produce, with four new implications (the deleted pairs given the parents of both);
- (c) the PC-stable skeleton from folder 27 (Fisher-z, alpha 0.01), oriented by the true ordering, as the
  discovery output a pipeline would actually audit.

## Tests (one implication = X independent of Y given pa(X) union pa(Y), X later in topological order)
- Linear/logistic LRT, the pipeline's test: regression of X on Z with and without Y (linear for continuous X,
  logistic for binary X), chi-square on one degree of freedom.
- GCM (Shah and Peters 2020), the nonlinear test: regress X on Z and Y on Z with flexible models (cubic splines,
  6 knots, on continuous covariates; linear in binary covariates; least squares for continuous targets, logistic
  for binary ones), then test the normalised mean of the product of residuals against N(0, 1).
Both with Benjamini-Hochberg at 0.05 over the implications of the graph being audited.

## Measures (per graph, test and seed, then mean over seeds)
- rejections among the implications; graph rejected = any BH rejection.
- For graph (b): how many of the four deletion implications are rejected (0 to 4).
- For graph (a): share of implications rejected (false-rejection rate; BH keeps it near 0.05 or below).
- For graph (c): rejections, and whether they sit on the non-monotone pairs.

## Pass bars
- A1 (the linear audit is blind): the linear LRT rejects at most 1 of the 4 deletion implications on average.
- A2 (the nonlinear audit sees them): the GCM rejects at least 3 of 4 on average.
- A3 (the nonlinear test is not just trigger-happy): on the true graph the GCM's false-rejection rate is at most
  0.10, and within 0.05 of the linear test's.

## What each outcome lets us say
- All three pass: the figure for the paper's second pitch. A pipeline whose audit is linear certifies a graph with
  the non-monotone edges missing; swapping the test for GCM rejects it at the same false-alarm rate.
- A2 fails: these particular mechanisms are not detectable by GCM at this size given the conditioning sets; the
  claim is weakened to "a nonlinear test is necessary but may not suffice".
- A3 fails: the nonlinear test buys detection with false alarms; report the trade-off.
