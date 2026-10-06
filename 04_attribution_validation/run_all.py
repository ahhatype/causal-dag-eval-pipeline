"""Usage: uv run python run_all.py [--replicates N] [--methods m1,m2] [--datasets d1,d2] [--feature-sets f1,f2]"""

import argparse
from pathlib import Path

import yaml

from attribution_validation import data as D
from attribution_validation.pipeline import run

here = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("--replicates", type=int)
ap.add_argument("--methods")
ap.add_argument("--datasets")
ap.add_argument("--feature-sets")
a = ap.parse_args()
split = lambda s: s.split(",") if s else None  # noqa: E731
params = yaml.safe_load((here / "config" / "params.yaml").read_text())
out = D.env_dir("OUTPUT_DIR", "./outputs") / "04_attribution_validation"
summary = run(params, D.env_dir("DATA_DIR", "./02_data"), out, a.replicates, split(a.methods), split(a.datasets),
              split(a.feature_sets))
print(summary.to_string())
print(f"Wrote attribution outputs to {out}")
