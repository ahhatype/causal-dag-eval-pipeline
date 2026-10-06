# CLAUDE.md

## Project
causal-dag-eval-pipeline: evaluates an established causal DAG (NASA renal-stone DAG SA-07566, adapted) in synthetic data — structural validation (testable implications, causal discovery) and attribution validation (standard and causal SHAP against simulated total effects).

## Layout
- `01_data_generation/` (R): simulates replicate data sets and ground truth from `config/`.
- `02_data/`: generated data (`DATA_DIR`); git-ignored.
- `03_shape_validation/` (R): implication tests, power control, MGM-PC discovery.
- `04_attribution_validation/` (Python, uv): attribution methods and evaluation.
- `outputs/` (`OUTPUT_DIR`): analysis outputs; git-ignored.
- `references/`: drafts and notes; git-ignored, never cite it from README or code.

## Conventions
- R: config-driven (`config/params.yaml`), testthat tests in each module's `tests/`, packages pinned in `renv.lock`.
- Python: uv-managed project with pinned lockfile, pytest tests.
- Seeds live in each module's `config/params.yaml`; never use unseeded randomness.
- Directories come from `.env` (`DATA_DIR`, `OUTPUT_DIR`); never hard-code keys.
- Comments are sparse; the README carries intent. Never describe the editing process in code, comments or README.

## Rules
- Never read, print, or commit `.env`.
- Ask before adding heavy dependencies.
- Add tests for any new metric.

## Commands
- Generate data: `make data` (`make recalibrate` after changing `01_data_generation/config`)
- Shape validation: `make shape`
- Tests: `make test`
