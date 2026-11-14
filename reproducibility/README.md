# Reproducibility checks

Two checks of what this repository reports. Both only read their inputs; neither reads `.env` or needs credentials or a network connection. Each prints a JSON result and writes a file only where `--json` points.

## Replay the stored results

`verify_saved_results.py` recomputes the reported summary tables from the committed per-run and per-replicate results and compares them with the stored tables. It needs Python with NumPy and pandas; the attribution module's environment has both.

```bash
python3 reproducibility/verify_saved_results.py
# or, in the attribution module's pinned environment:
cd 04_attribution_validation && uv run python ../reproducibility/verify_saved_results.py --repo ..
```

`--repo` defaults to the current directory. Stdout carries a compact JSON summary of the comparisons. `--json <file>` also writes the full report: the SHA-256 of each input file as checked out (so the hashes differ between LF and CRLF checkouts), the runtime versions, the comparison counts, any discrepancies and the headline values. The exit status is zero only if every comparison passes. Numeric comparisons use absolute and relative tolerances of `1e-10`; row keys, row counts, classifications and other text fields must match exactly.

The attribution tables are recomputed with the module's own aggregation functions, compiled from `04_attribution_validation/src/attribution_validation/pipeline.py` without importing the package, so the fitting libraries are not needed. The other tables are rebuilt from their formulas. The `06_lumawarp` checks read the committed combined CSVs, because the `results/` and `cells/` shards that the `06_lumawarp` summary scripts expect are not committed.

Coverage:

- Attribution summaries, paired method contrasts, and sample-size and selection contrasts.
- Structure-validation means and Monte Carlo standard errors.
- ConfoundingSHAP summaries from the per-run metrics.
- TMLE, g-computation and unadjusted summaries from the retained estimates, and the propensity-bound summaries.
- The paired Lumawarp and Super Learner explained-predictor comparison from the per-run files.
- Edge-audit retention decisions from the stored candidate scores, the per-replicate audit scores and the audit summary.
- Both nonlinear hybrid summary tables from the per-replicate metrics.

This replays the aggregation of saved output. It does not regenerate the data, refit a predictive or structural model, recompute SHAP from per-record credits, reproduce nuisance predictions or regenerate the Lumawarp monotonicity flags. The nonlinear implication-test summary is not rebuilt, because the per-test results behind it are not committed.

## Check regenerated data

`verify_generated_data.py` compares the SHA-256 of every replicate CSV in a data directory with `expected-data-sha256.csv`: 50 replicates, three data sets, observed and latent files, 300 files in all. It uses only the Python standard library.

```bash
make data
python3 reproducibility/verify_generated_data.py --data-dir 02_data
```

Generate with the frozen calibration (`make data`, not `make recalibrate`). To leave existing data in place, point `DATA_DIR` in `.env` at another relative directory first (for example `DATA_DIR=../regenerated-data`) and pass the same directory as `--data-dir`. The exit status is zero only if all 300 files are present and match; the result lists each missing or differing file, and `--json <file>` saves it.
