"""Usage: uv run python make_table1.py

Table 1 data and the ancestor-displacement check, from the stored results (nothing is rerun):
table1_ordering.csv, ancestor_displacement.csv, ancestor_displacement_summary.csv, then table1_rendered.md and
table1_wide.csv with the manuscript labels from config/labels.yaml.
"""

from pathlib import Path

import pandas as pd

from attribution_validation import data as D
from attribution_validation import table1 as T

HERE = Path(__file__).resolve().parent
DATASETS = ["full_set", "reference_subsample", "astronaut_set"]


def main() -> None:
    labels = T.load_labels(HERE / "config" / "labels.yaml")
    data_dir = D.env_dir("DATA_DIR", "./02_data")
    out = D.env_dir("OUTPUT_DIR", "./outputs") / "04_attribution_validation"
    importance = pd.read_csv(out / "importance.csv")
    truth = D.truth_values(data_dir)
    ties = pd.read_csv(out / "tie_groups.csv")
    fsets = pd.read_csv(data_dir / "feature_sets.csv")
    fsets = fsets[fsets["in_all_features_set"]].reset_index(drop=True)
    long = T.table1_long(importance, truth, ties, fsets, labels, DATASETS)
    disp, summ = T.displacement(importance, fsets, labels, DATASETS)
    long.to_csv(out / "table1_ordering.csv", index=False)
    disp.to_csv(out / "ancestor_displacement.csv", index=False)
    summ.to_csv(out / "ancestor_displacement_summary.csv", index=False)
    wide = T.wide_table(long, disp)
    pd.concat([t.assign(dataset=labels["datasets"][ds]) for ds, t in wide.items()], ignore_index=True).to_csv(
        out / "table1_wide.csv", index=False)
    (out / "table1_rendered.md").write_text(T.markdown(wide, labels, summ))
    print(f"wrote table1_ordering.csv ({len(long)} rows), ancestor_displacement.csv ({len(disp)}), "
          f"ancestor_displacement_summary.csv ({len(summ)}), table1_wide.csv, table1_rendered.md to {out}")


if __name__ == "__main__":
    main()
