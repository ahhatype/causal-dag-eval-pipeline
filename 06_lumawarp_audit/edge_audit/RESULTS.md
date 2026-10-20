# Results: ordering-informed edge audit on the paper's data sets (folder 33, 14 to 15 October 2026)

Pre-registered in PREREGISTRATION.md before any fit. 1,050 LumaWarp fits (old build, format-1 models kept for the
shape flag), 0 failures; logistic/linear, spline additive and XGBoost on the identical files. Replicates 1 to 10 of
the three regenerated data sets (pipeline commit 79a1b56): full-set 10,000 rows, reference-set and astronaut-set
900 rows. 36 analysis-graph variables, 63 edges, orientation given by the DAG's topological order.

## Verdicts

| bar | rule | observed | verdict |
|---|---|---|---|
| E1 | full-set recall >= 0.80 at precision >= 0.90 | recall 0.852 (MCSE 0.009), precision 0.932 | PASS |
| E2 | full-set shape flags <= 5% of continuous-parent cases | 0 of 160 | PASS |
| E3 | within 0.10 recall of the best reference at each data set | +0.006 (full), -0.002 (reference), -0.013 (astronaut) | matches |

## Summary (mean of 10 replicates; MCSE of recall 0.01 to 0.02)

| data set | learner | edges called | recall | precision | noise kept |
|---|---|---|---|---|---|
| full-set (10,000) | LumaWarp | 57.7 | 0.852 | 0.932 | 0.001 |
| | logistic/linear | 51.0 | 0.803 | 0.992 | 0.000 |
| | spline additive | 57.5 | 0.846 | 0.928 | 0.002 |
| | XGBoost | 51.0 | 0.733 | 0.908 | 0.002 |
| | draft's MGM-PC, no ordering (50 reps) | | 0.78 | 1.00 | |
| reference-set (900) | LumaWarp | 33.7 | 0.476 | 0.909 | 0.004 |
| | logistic/linear | 33.1 | 0.478 | 0.912 | 0.000 |
| | spline additive | 35.5 | 0.459 | 0.823 | 0.007 |
| | XGBoost | 27.2 | 0.414 | 0.965 | 0.005 |
| | draft's MGM-PC, no ordering | | 0.56 | 0.99 | |
| astronaut-set (900) | LumaWarp | 35.0 | 0.459 | 0.843 | 0.003 |
| | logistic/linear | 32.2 | 0.471 | 0.926 | 0.001 |
| | spline additive | 38.2 | 0.454 | 0.749 | 0.007 |
| | XGBoost | 28.2 | 0.408 | 0.912 | 0.002 |
| | draft's MGM-PC, no ordering | | 0.55 | 0.97 | |

Per-replicate LumaWarp on the full-set: recall 0.81 to 0.89, precision 0.84 to 0.97 (replicate 1 is the only one
under 0.90). Shape flags: full-set 0 of 160; reference-set 6 of 160; astronaut-set 6 of 160 (3.75%, all on edges
whose centroid sits just above the noisiest noise column, 1.58 to 2.44 against 1.57 to 2.34).

## What the edge-level breakdown says

- **At 10,000 rows the audit recovers the countermeasure edges that MGM-PC never finds.** K+ citrate -> urine
  chemistry is kept in 10 of 10 replicates by every learner (the pipeline's MGM-PC: 0 of 50 in every data set).
  Thiazides -> urine chemistry is weaker: LumaWarp 3 of 10, logistic 6, spline 5, XGBoost 7. Bisphosphonates -> bone
  remodeling 10 of 10 for all but XGBoost (8).
- **The misses are the same edges for every learner.** Mission era -> bisphosphonates, mission era -> detect long-term
  outcomes, medications -> evacuation and medications -> loss of crew life are found 0 of 10 by all four learners
  at 10,000 rows; evacuation -> loss of crew life and mission era -> thiazides 0 to 1 of 10. These are the generator's
  smallest coefficients on rare binary outcomes, and they set the recall ceiling (about 0.90) for any per-node audit.
- **LumaWarp's extra recall over logistic comes with lower precision, and spline matches it exactly.** LumaWarp and
  spline call 57.5 to 57.7 edges and land within 0.006 in recall and 0.004 in precision. Logistic calls 51 with
  near-perfect precision. The engine does nothing on this linear generator that an additive spline does not.
- **At 900 rows nothing recovers the countermeasures.** K+ citrate -> urine chemistry 0 to 3 of 10, thiazides and
  bisphosphonates 0 of 10, for every learner on both 900-row sets. Recall halves (0.41 to 0.48) for all learners;
  LumaWarp equals logistic and does not beat it. The draft's MGM-PC without an ordering reports higher recall here
  (0.55 to 0.56) at a looser bar (any adjacency, no noise-column yardstick), so the two numbers are not comparable;
  what is comparable is that both fail on the countermeasure edges at 900 rows.
- **The recurring false edges are the same across learners and sets**: evacuation -> loss of mission (4 of 10 in all
  three sets; both are children of the same parents), cumulative mission duration -> ultrasound and -> hydration
  (the skewed exposure seen in folder 31), and resistive exercise -> tamsulosin on the full-set. Every one of them is
  a confounded sibling pair the ordering cannot separate, not an engine artefact.
- **The shape flag is clean at 10,000 rows and fires at 3.75% on 900 rows**, all marginal (centroid within 0.5 of the
  noise ceiling). The generator is linear, so these 12 are false positives and set the flag's small-sample floor.

## Deviations
- The full-set stage was run as two engine streams (run_jobs2.sh, forward and reverse order with a claim file) after
  the first single-stream fit took over two minutes; the engine's fit is unchanged, only the scheduling.
- Paused once (14 October, 80 of 350 full-set fits) at Lauren's request and resumed with the same scripts; the two
  interrupted fits were refit from scratch.

## Files
scores_by_rep.csv (one row per replicate, data set and learner), summary.csv, edges_all.csv (every candidate score
and keep decision), shape_flags.csv, score.py, references.py, make_jobs.py, run_jobs.sh, run_jobs2.sh, graph54.py.
