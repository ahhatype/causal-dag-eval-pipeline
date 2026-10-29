# Pre-registration: can the new build's importance file replace the kernel-slice shape flag? (written 19 October 2026, before any fit)

## Question
The monotonicity flag used in folders 14, 28 and 33 reads the degree centroid of each input's kernel slice, which
only the old build's model file exposes. The new build (12 October) writes a feature-importance file with, per
input, the variance of the learned h0 and h1 operator coefficients along the axis, theta, and importance
(theta x mean relative variance). A curved effect usually needs a non-uniform operator, so one of these columns may
carry the same information. Does any of them reproduce the centroid flag's detections without raising its false
flags?

## Data (fixed)
Folder 28's 144 training tables, copied unchanged (nl_v2 with amplitude 0.5, 1.0, 1.5; innovation sd x1, x2;
2,000, 5,000, 10,000 rows; two seeds; four continuous targets each). Truth per case from folder 28's
`shape_flags.csv`: 144 non-monotone parent cases and 72 monotone (linear or threshold) controls, plus the old
build's centroid flag for every case.

## Fits (144)
New build, research configuration with `feature_importance = true` (folder 26's cfg_research.cfg).

## Readings (each applied with folder 14's rule: flag the parent if its value exceeds every continuous noise
column's value in the same model)
R1 h0 variance; R2 h1 variance; R3 relative variance = importance / theta (the operator's non-uniformity with the
relevance divided out); R4 importance itself (expected to fail: it is a relevance measure).

## Bars
- F1: a reading flags at least 80 percent as many non-monotone cases as the centroid flag does (the centroid
  flagged 106 of 144) with at most 4 of 72 monotone controls flagged (the centroid flagged 0).
- F2 (reported): case-level agreement with the centroid flag (Cohen's kappa) and the AUC of each reading for
  non-monotone versus monotone over all 216 cases.

## What each outcome lets us say
- F1 passes for some reading: the shape flag can be read from the new build's CSV; the Methods wording can describe
  it as a reading from the fitted model's reported coefficients, and the old-build dependency disappears.
- F1 fails: the h0/h1 variances measure non-uniformity, not curvature; the flag stays a kernel-slice reading from
  the old build, and the Methods wording stands as written.
