# Pre-registration: a non-monotone variant of the paper's 54-node generator (written 13 October 2026, before any run)

## Question
The paper's central experiment (`causal-dag-eval-pipeline`, 54 nodes, 90 edges, 50 replicates of three data sets)
has only linear and logistic mechanisms, so its audit and discovery steps cannot show the "linearity hides" half of
the argument on the paper's own testbed. This variant adds curved mechanisms to the same generator so that the
hybrid (folders 27, 29) and the nonlinear audit (folder 30) can be run inside the paper's experiment.

## The variant (branch `nonmonotone-variant` of the pipeline repository; one commit, generator untouched otherwise)
`edges.csv` gains a `mechanism` column (default `linear`) and `model_spec.R`'s `lin_expr()` reads it. Four edges
change, chosen where the physiology argues for a non-monotone or saturating dose-response and where the edge is
continuous-to-continuous so the shape is identifiable from the data:

| Edge | Mechanism | Coefficient (was) | Rationale |
|---|---|---|---|
| Nutrients -> Bone Remodeling | inverted U, centre 0.2 | 0.35 (0.25) | deficiency and excess both raise turnover (the vitamin D evidence in the rework note) |
| Microbiome -> Urine Chemistry | U, centre 0.2 | 0.35 (0.20) | departure from a healthy range in either direction |
| Hydration -> Urine Concentration | U, centre 0.2 | 0.45 (-0.60) | low hydration concentrates urine; the U adds an over-hydration arm so that the edge is invisible to a linear test (a strong monotone edge would be found anyway) |
| Urine Chemistry -> Mineralized Renal Material | saturating, 1.2 tanh | 0.60 (0.60) | monotone nonlinear control: should still be found by linear tests |
All curved terms are centred to mean zero under standardized parents, so the children's calibration (`make
recalibrate`) keeps their scales and the outcome prevalence. `nonadditive_nodes` in the shape-validation config
lists the four children so their implications are flagged as in the pipeline's own convention.

## Plan
1. `make recalibrate` then `make data` on the branch: 50 replicates of the three data sets, new ground truth.
2. `make shape`: the pipeline's own implication tests, power control and MGM-PC discovery, unchanged, on the variant.
   Prediction: the three U-shaped edges are found in under 10% of replicates by MGM-PC; the saturating edge in over
   90%; the linear implication tests do not flag the missing U-edges (they are not in the audited graph, so no
   implication corresponds to them) and their calibration stays uniform.
3. Hybrid (folder 29's rule, XGBoost gap rule, then LumaWarp if time allows) on the full-set replicates: prediction,
   the three U-shaped edges recovered and oriented in over 80% of replicates with at most 0.2 false edges per
   replicate after the chain veto.
4. Nonlinear audit (folder 30's GCM) on the MGM-PC graph and on the true graph with the U-edges deleted.
5. Attribution (`make attribution`) if time allows, to see what the curved edges do to the SHAP rankings.

## Bars
Same as folders 27, 29 and 30, applied per replicate on the full-set. Reported for the 900-row sets.

## Blockers recorded before the run
`dagitty` (and therefore the generator's validation step and the comparison code) needs the R package `V8`, which
needs `libcurl4-openssl-dev` on this machine; `rCausalMGM` and `testthat` are being installed. Until then the branch
is a patch, not a result.
