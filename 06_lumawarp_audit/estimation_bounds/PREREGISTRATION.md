# Pre-registration: propensity-bound sensitivity of the folder 34 estimates (written 19 October 2026, before any run)

## Question
Folder 34 clipped 8 to 21 percent of propensity fits at the 0.025/0.975 bounds, most on the 900-row sets, and LW-TMLE
covered the truth in only 0.82 of astronaut-set replicates for K+ citrate. Is the under-coverage caused by the bounds
(a tuning choice) or by the sample (900 rows, 68 to 93 events)?

## Data (fixed)
Folder 34's 600 `tmle_inputs.csv` files: the cross-fitted LumaWarp and Super Learner outcome and propensity
predictions already computed. No engine fits and no Super Learner refits; only the targeting step is rerun.

## Settings
`tmle()` as in folder 34 (binomial, logistic fluctuation, prescreenW.g = FALSE, target.gwt = FALSE) with gbound in
{0.025/0.975 (folder 34), 0.05/0.95, 0.01/0.99, 0.001/0.999}, for the LumaWarp and the Super Learner nuisances.

## Measures
Bias, mean SE, empirical SE, coverage, per data set, treatment, route, learner and bound.

## Bars
- S1 (bounds): if some bound raises K+ citrate coverage on both 900-row sets to >= 0.90 for LW-TMLE while full-set
  coverage stays >= 0.90 and |bias| <= 0.004, the under-coverage is a tuning matter and that bound is reported as a
  sensitivity. Otherwise the under-coverage is the sample's and the draft's sample-size statement stands for estimation.
- Reported: the clipping share at each bound.
