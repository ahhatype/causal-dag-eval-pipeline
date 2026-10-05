# CLAUDE.md

## Project
causal-dag-eval-pipeline: evaluates causal DAGs against reference graphs (structural metrics such as SHD, precision/recall on edges, orientation accuracy).

## Conventions
- Python 3.11+, type hints, `ruff` for lint/format, `pytest` for tests.
- Config comes from environment variables (see `.env.example`); never hard-code keys.
- Keep runs reproducible: seed from `RANDOM_SEED`, write results to `OUTPUT_DIR`.
- Raw data lives in `DATA_DIR` and is git-ignored; do not commit datasets or outputs.

## Rules
- Never read, print, or commit `.env`.
- Ask before adding heavy dependencies.
- Add tests for any new metric.

## Commands
- Install: `pip install -r requirements.txt`
- Test: `pytest`
- Lint: `ruff check . && ruff format .`
