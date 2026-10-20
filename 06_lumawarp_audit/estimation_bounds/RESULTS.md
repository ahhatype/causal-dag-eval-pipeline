# Results: propensity-bound sensitivity (folder 36, 19 October 2026)

Pre-registered in PREREGISTRATION.md. All 600 of folder 34's units re-targeted at four propensity bounds with the
LumaWarp nuisances as stored and the Super Learner nuisances refitted with the same seeds (the stored LumaWarp
estimate at the original bound reproduces folder 34 to 1e-16; 0 failures).

## Verdict
S1: no bound restores K+ citrate coverage to 0.90 on both 900-row sets. LW-TMLE coverage on the astronaut set is
0.80 to 0.82 at every bound from 0.001 to 0.05, while the full set stays 0.90 to 0.92 and the reference set 0.94 to
0.96. The under-coverage is the sample's (900 rows, 68 events), not a tuning choice.

## What else the table says
- Bias is insensitive to the bounds: the K+ citrate estimates move by at most 0.002 between the tightest and the
  loosest bound on every data set and for both learners. The sign flip is removed at every setting.
- The bounds only change how many propensities are clipped (0 to 33 percent) and the estimated SE slightly; they
  do not change the answer. The 0.025/0.975 setting used in folders 25, 32 and 34 is not driving any result.
- The astronaut set's residual under-coverage comes with a small negative bias for both learners (LW -0.007, SL
  -0.012 at the original bound) relative to the selected-population truth, which is the harder target the draft
  itself notes.

## K+ citrate, DAG-parent route (50 replicates per cell)

                       bias                                    clipped                                   coverage                                    mean_se                             
dataset       astronaut_set full_set reference_subsample astronaut_set full_set reference_subsample astronaut_set full_set reference_subsample astronaut_set full_set reference_subsample
learner bound                                                                                                                                                                            
lw      0.001        -0.007    0.001               0.004         0.000    0.000               0.000          0.80     0.92                0.96         0.027    0.011               0.033
        0.010        -0.008    0.000               0.004         0.029    0.031               0.026          0.82     0.92                0.96         0.027    0.011               0.033
        0.025        -0.007   -0.000               0.002         0.155    0.149               0.109          0.82     0.92                0.96         0.027    0.010               0.031
        0.050        -0.006   -0.000               0.002         0.400    0.331               0.288          0.82     0.90                0.94         0.024    0.008               0.028
sl      0.001        -0.010    0.001               0.002         0.002    0.000               0.001          0.72     0.88                0.96         0.026    0.011               0.033
        0.010        -0.012    0.000               0.000         0.071    0.040               0.044          0.74     0.90                0.96         0.026    0.010               0.032
        0.025        -0.012   -0.000              -0.000         0.234    0.157               0.158          0.74     0.92                0.94         0.025    0.009               0.031
        0.050        -0.010   -0.000               0.000         0.449    0.332               0.335          0.76     0.92                0.96         0.023    0.008               0.027

## For the draft
The estimation paragraph can say that coverage at 900 rows was 0.82 on the selected set and that this did not
depend on the propensity truncation (0.001 to 0.05 tried). It is a sample-size statement, consistent with the
draft's headline. Files: results_all.csv, summary.csv, rerun_bounds.R, run_worker.sh.
