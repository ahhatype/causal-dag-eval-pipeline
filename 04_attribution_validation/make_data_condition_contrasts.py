"""Usage: uv run python make_data_condition_contrasts.py

Sample-size and selection contrasts of the attribution metrics, from the stored per_run.csv (nothing is rerun):
data_condition_contrasts.csv.
"""

import pandas as pd

from attribution_validation import data as D
from attribution_validation.pipeline import data_condition_contrasts


def main() -> None:
    out = D.env_dir("OUTPUT_DIR", "./outputs") / "04_attribution_validation"
    contrasts = data_condition_contrasts(pd.read_csv(out / "per_run.csv"))
    contrasts.to_csv(out / "data_condition_contrasts.csv", index=False)
    print(f"wrote data_condition_contrasts.csv ({len(contrasts)} rows) to {out}")


if __name__ == "__main__":
    main()
