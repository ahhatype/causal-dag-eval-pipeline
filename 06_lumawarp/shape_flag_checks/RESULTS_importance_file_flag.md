# Results: the importance file does not carry the shape flag (folder 40, 19 October 2026)

Pre-registered in PREREGISTRATION.md. 144 new-build fits on folder 28's tables (0 failures, 3.8 h), four readings
from the feature-importance file scored under folder 14's noise rule against the 216 known cases.

## Verdict: F1 FAIL for every reading

| reading | non-monotone flagged (of 144) | monotone false (of 72) | AUC curved vs linear | kappa vs centroid flag |
|---|---|---|---|---|
| h0 variance | 144 | 72 | 0.30 | 0.00 |
| h1 variance | 144 | 69 | 0.53 | 0.03 |
| importance / theta (relative variance) | 143 | 72 | 0.33 | 0.01 |
| importance | 144 | 72 | 0.51 | 0.00 |
| centroid flag, old build (folder 28) | 105 | 0 | | |

Every reading sits above the noise columns for every real parent, curved or straight, at every sample size. They
separate "has an effect" from "has no effect" perfectly and "curved" from "straight" not at all (AUC 0.30 to 0.53,
chance 0.50). The operator variances measure how unevenly an input acts along its range once it matters, which a
linear parent does as much as a curved one; the degree centroid measures which basis degrees carry the effect,
which only the kernel slice exposes.

## What this settles
- The shape flag stays a kernel-slice reading from the old build's model file; the Methods wording as written is
  the correct one. The new build cannot supply the flag today, with or without the importance file.
- Relevance is the one quantity the two builds share (folder 26 agreement 0.96), and any Methods sentence that names
  the new build should be about relevance only.
- A lead for Lexi, not for the paper: exposing the per-axis kernel slice (or just its diagonal in the degree basis)
  in the new build's report or diagnostics would make the flag readable from the current build.
Files: flags_by_reading.csv, summary.csv, score.py, run_jobs.sh.
