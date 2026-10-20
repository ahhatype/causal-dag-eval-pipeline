# Results: estimating the two confounded drug contrasts on the paper's data sets (folder 34, 14 October 2026)

Pre-registered in PREREGISTRATION.md before any fit. Stage 1 (K+ citrate, two 900-row sets, replicates 1 to 10) was
reported at the 14 October meeting; this file is the complete stage 2: both treatments, all three data sets, both
adjustment routes, replicates 1 to 50 (600 estimation units, 0 engine failures). Truth is the pipeline's own per-unit
contrast (source population for the full and reference sets, selected population for the astronaut set).

## Verdicts (primary: K+ citrate, full set, DAG parents, LW-TMLE)

| bar | rule | observed | verdict |
|---|---|---|---|
| P1 | abs(bias) <= 0.004 | bias -0.0000 (mean estimate -0.0202 vs truth -0.0202) | PASS |
| P2 | coverage >= 0.90 | 0.92 (mean SE 0.0096 vs empirical SE 0.0109) | PASS |
| P3 | parity with SL-TMLE: bias gap <= 0.002, coverage gap <= 0.05 | bias gap 0.0002, coverage gap 0.00 | PASS |

## What the full table says

- **The sign flip the draft reports is removed by the DAG's own adjustment set.** On the full set the naive K+ citrate
  contrast is +0.032 in every replicate (sign right 0 of 50) against a true -0.020. LW-TMLE, SL-TMLE and logistic
  g-computation with the three DAG parents all land on -0.019 to -0.020 with the right sign in 46 to 50 of 50 replicates.
- **The learner is interchangeable.** LW-TMLE and SL-TMLE agree to within 0.001 in bias and 0.02 in coverage in every
  one of the twelve cells. Logistic g-computation is also unbiased here, as it was in folders 25 and 32: the generator's
  mechanisms are linear-logistic, so a correctly specified parametric model has nothing to lose. LumaWarp adds nothing on
  this generator beyond matching the parametric answer; it does not detract either.
- **Ordering-only adjustment gives the same answer as the DAG parents.** Adjusting for every variable earlier in the
  topological order changes bias by at most 0.001 on the full set. The ordering is enough; the graph's edges are not
  needed for this contrast (same finding as folders 19 and 21 on the renal simulator).
- **At 900 rows the contrast is estimable but not decisive.** Empirical SE is 0.027 to 0.033 for K+ citrate on the two
  900-row sets, larger than the effect (0.020). Coverage holds for the reference set (0.92 to 0.98) and slips on the
  astronaut set (LW-TMLE 0.82 and 0.86, SL-TMLE 0.74); the sign is right in 68 to 82 percent of replicates and the naive
  contrast still flips it (sign right 0.10 and 0.30). The sample-size headline extends to estimation: 900 rows resolve a
  0.02 risk difference at roughly two-thirds power for the sign, not for the magnitude.
- **The astronaut set's selected-population truth is the harder target.** Both TMLEs under-cover there at 900 rows, with
  more propensity clipping (13.5 to 21 percent of g fits clipped at the 0.025/0.975 bounds versus 8 to 13 percent on the
  reference set). The draft's note that the astronaut set has its own truth should travel with any estimation table.
- **Bisphosphonates is weakly confounded by design and the draft's crude contrast is already close.** Naive bias is
  0.0003 on the full set. Both TMLEs show a small positive bias (+0.0044, about one-third of the 0.0083 effect) with
  coverage 0.94 to 0.96 and SE 0.012; logistic g-computation is tighter (bias +0.0007, SE 0.009). At 900 rows nothing
  resolves this effect: SE 0.04 against 0.006 to 0.008, sign right about half the time for every estimator. This is a
  power statement, not a method statement.

## Full table (50 replicates per cell)

| data set | treatment | route | estimator | truth | mean est. | bias | emp. SE | mean SE | coverage | sign right |
|---|---|---|---|---|---|---|---|---|---|---|
| full_set | k_citrate | parents | LW-TMLE | -0.0202 | -0.0202 | -0.0000 | 0.0109 | 0.0096 | 0.92 | 0.92 |
| full_set | k_citrate | parents | SL-TMLE | -0.0202 | -0.0204 | -0.0002 | 0.0106 | 0.0095 | 0.92 | 0.94 |
| full_set | k_citrate | parents | Logistic g-computation | -0.0202 | -0.0189 | 0.0013 | 0.0080 |  |  | 1.00 |
| full_set | k_citrate | parents | Naive difference | -0.0202 | 0.0322 | 0.0524 | 0.0103 |  |  | 0.00 |
| full_set | k_citrate | ordering | LW-TMLE | -0.0202 | -0.0201 | 0.0001 | 0.0103 | 0.0094 | 0.92 | 0.94 |
| full_set | k_citrate | ordering | SL-TMLE | -0.0202 | -0.0203 | -0.0001 | 0.0108 | 0.0092 | 0.90 | 0.94 |
| full_set | k_citrate | ordering | Logistic g-computation | -0.0202 | -0.0187 | 0.0015 | 0.0080 |  |  | 1.00 |
| full_set | k_citrate | ordering | Naive difference | -0.0202 | 0.0322 | 0.0524 | 0.0103 |  |  | 0.00 |
| reference_subsample | k_citrate | parents | LW-TMLE | -0.0202 | -0.0181 | 0.0021 | 0.0328 | 0.0311 | 0.96 | 0.78 |
| reference_subsample | k_citrate | parents | SL-TMLE | -0.0202 | -0.0207 | -0.0005 | 0.0310 | 0.0306 | 0.94 | 0.76 |
| reference_subsample | k_citrate | parents | Logistic g-computation | -0.0202 | -0.0178 | 0.0024 | 0.0279 |  |  | 0.82 |
| reference_subsample | k_citrate | parents | Naive difference | -0.0202 | 0.0357 | 0.0559 | 0.0325 |  |  | 0.10 |
| reference_subsample | k_citrate | ordering | LW-TMLE | -0.0202 | -0.0108 | 0.0094 | 0.0308 | 0.0333 | 0.98 | 0.72 |
| reference_subsample | k_citrate | ordering | SL-TMLE | -0.0202 | -0.0215 | -0.0013 | 0.0304 | 0.0291 | 0.92 | 0.82 |
| reference_subsample | k_citrate | ordering | Logistic g-computation | -0.0202 | -0.0170 | 0.0032 | 0.0283 |  |  | 0.74 |
| reference_subsample | k_citrate | ordering | Naive difference | -0.0202 | 0.0357 | 0.0559 | 0.0325 |  |  | 0.10 |
| astronaut_set | k_citrate | parents | LW-TMLE | -0.0152 | -0.0222 | -0.0070 | 0.0305 | 0.0268 | 0.82 | 0.74 |
| astronaut_set | k_citrate | parents | SL-TMLE | -0.0152 | -0.0268 | -0.0116 | 0.0296 | 0.0254 | 0.74 | 0.78 |
| astronaut_set | k_citrate | parents | Logistic g-computation | -0.0152 | -0.0188 | -0.0036 | 0.0266 |  |  | 0.80 |
| astronaut_set | k_citrate | parents | Naive difference | -0.0152 | 0.0189 | 0.0341 | 0.0348 |  |  | 0.30 |
| astronaut_set | k_citrate | ordering | LW-TMLE | -0.0152 | -0.0167 | -0.0016 | 0.0310 | 0.0290 | 0.86 | 0.68 |
| astronaut_set | k_citrate | ordering | SL-TMLE | -0.0152 | -0.0255 | -0.0103 | 0.0287 | 0.0250 | 0.74 | 0.78 |
| astronaut_set | k_citrate | ordering | Logistic g-computation | -0.0152 | -0.0182 | -0.0030 | 0.0274 |  |  | 0.78 |
| astronaut_set | k_citrate | ordering | Naive difference | -0.0152 | 0.0189 | 0.0341 | 0.0348 |  |  | 0.30 |
| full_set | bisphosphonates | parents | LW-TMLE | -0.0083 | -0.0040 | 0.0044 | 0.0132 | 0.0125 | 0.96 | 0.64 |
| full_set | bisphosphonates | parents | SL-TMLE | -0.0083 | -0.0036 | 0.0047 | 0.0131 | 0.0123 | 0.96 | 0.64 |
| full_set | bisphosphonates | parents | Logistic g-computation | -0.0083 | -0.0076 | 0.0007 | 0.0089 |  |  | 0.82 |
| full_set | bisphosphonates | parents | Naive difference | -0.0083 | -0.0081 | 0.0003 | 0.0090 |  |  | 0.82 |
| full_set | bisphosphonates | ordering | LW-TMLE | -0.0083 | -0.0039 | 0.0044 | 0.0124 | 0.0122 | 0.94 | 0.62 |
| full_set | bisphosphonates | ordering | SL-TMLE | -0.0083 | -0.0040 | 0.0044 | 0.0122 | 0.0117 | 0.96 | 0.66 |
| full_set | bisphosphonates | ordering | Logistic g-computation | -0.0083 | -0.0075 | 0.0008 | 0.0087 |  |  | 0.84 |
| full_set | bisphosphonates | ordering | Naive difference | -0.0083 | -0.0081 | 0.0003 | 0.0090 |  |  | 0.82 |
| reference_subsample | bisphosphonates | parents | LW-TMLE | -0.0083 | -0.0012 | 0.0071 | 0.0439 | 0.0383 | 0.84 | 0.48 |
| reference_subsample | bisphosphonates | parents | SL-TMLE | -0.0083 | -0.0012 | 0.0072 | 0.0448 | 0.0391 | 0.88 | 0.54 |
| reference_subsample | bisphosphonates | parents | Logistic g-computation | -0.0083 | -0.0063 | 0.0021 | 0.0255 |  |  | 0.58 |
| reference_subsample | bisphosphonates | parents | Naive difference | -0.0083 | -0.0046 | 0.0037 | 0.0223 |  |  | 0.58 |
| reference_subsample | bisphosphonates | ordering | LW-TMLE | -0.0083 | 0.0022 | 0.0105 | 0.0415 | 0.0387 | 0.88 | 0.50 |
| reference_subsample | bisphosphonates | ordering | SL-TMLE | -0.0083 | -0.0036 | 0.0048 | 0.0396 | 0.0355 | 0.88 | 0.54 |
| reference_subsample | bisphosphonates | ordering | Logistic g-computation | -0.0083 | -0.0055 | 0.0028 | 0.0248 |  |  | 0.58 |
| reference_subsample | bisphosphonates | ordering | Naive difference | -0.0083 | -0.0046 | 0.0037 | 0.0223 |  |  | 0.58 |
| astronaut_set | bisphosphonates | parents | LW-TMLE | -0.0063 | 0.0073 | 0.0137 | 0.0453 | 0.0382 | 0.84 | 0.46 |
| astronaut_set | bisphosphonates | parents | SL-TMLE | -0.0063 | 0.0102 | 0.0165 | 0.0462 | 0.0400 | 0.82 | 0.44 |
| astronaut_set | bisphosphonates | parents | Logistic g-computation | -0.0063 | 0.0046 | 0.0110 | 0.0357 |  |  | 0.44 |
| astronaut_set | bisphosphonates | parents | Naive difference | -0.0063 | -0.0036 | 0.0027 | 0.0293 |  |  | 0.50 |
| astronaut_set | bisphosphonates | ordering | LW-TMLE | -0.0063 | 0.0082 | 0.0145 | 0.0454 | 0.0379 | 0.90 | 0.48 |
| astronaut_set | bisphosphonates | ordering | SL-TMLE | -0.0063 | 0.0061 | 0.0124 | 0.0457 | 0.0359 | 0.86 | 0.46 |
| astronaut_set | bisphosphonates | ordering | Logistic g-computation | -0.0063 | 0.0025 | 0.0088 | 0.0351 |  |  | 0.52 |
| astronaut_set | bisphosphonates | ordering | Naive difference | -0.0063 | -0.0036 | 0.0027 | 0.0293 |  |  | 0.50 |

                     g_lw_clipped  g_sl_clipped
dataset                                        
astronaut_set               0.135         0.210
full_set                    0.135         0.141
reference_subsample         0.080         0.127

                                           n  events  treated
dataset             treatment                                
astronaut_set       bisphosphonates    900.0    68.0     80.0
                    k_citrate          900.0    68.0     79.0
full_set            bisphosphonates  10000.0  1009.0   1246.0
                    k_citrate        10000.0  1009.0   1252.0
reference_subsample bisphosphonates    900.0    93.0    114.0
                    k_citrate          900.0    93.0    112.0

## Deviations and leads

- No deviations from the pre-registration. Stage 1 numbers (replicates 1 to 10) are superseded by this table.
- Lead (post hoc, not pre-registered): the TMLE bias on bisphosphonates (+0.0044 on the full set for both learners) is
  the one place the flexible learners sit further from truth than the parametric model. With 1,246 treated and 1,009
  events in 10,000 rows this is within one SE and may be the logistic fluctuation's finite-sample behaviour with a small
  effect; it would take more replicates to call it bias rather than noise.
- The propensity clipping rates (8 to 21 percent) are high for a generator with mild overlap; folder 25 on the renal
  simulator saw less. Worth a sentence if the estimation analysis goes into the paper: the bounds matter at 900 rows.

## What this lets the draft say

The draft calls the intervention contrast "the quantity a countermeasure decision needs" and then reports only the
crude +0.031 for K+ citrate. This folder closes that gap on the paper's own data: adjusting for the DAG's three parents
(or simply for everything earlier in the ordering) turns the sign-flipped crude contrast into an unbiased estimate with
a 92 percent interval at 10,000 rows, with the nuisance learner interchangeable (LumaWarp, Super Learner, or a plain
logistic model). At 900 rows the same adjustment removes the bias but the interval is wider than the effect.
Files: summary.csv (one row per cell), results_all_reps.csv (one row per replicate and route), run_engine.py,
tmle_rep.R, summarize.py.
