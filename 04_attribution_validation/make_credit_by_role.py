"""Usage: uv run python make_credit_by_role.py

Share of attribution by feature role (ancestor, descendant, unrelated) in the all-features runs, from the stored
importance.csv (nothing is rerun): credit_by_role.csv. The shares are checked against the summary's non-ancestor share.
"""

import pandas as pd

from attribution_validation import data as D
from attribution_validation.evaluation import check_non_ancestor_share, credit_by_role, roles_from_feature_sets

METHODS = ["standard_shap", "ordering_only_shap", "ng_style_shap", "interventional_shap", "ng_style_shap_all_records"]


def main() -> None:
    data_dir = D.env_dir("DATA_DIR", "./02_data")
    out = D.env_dir("OUTPUT_DIR", "./outputs") / "04_attribution_validation"
    importance = pd.read_csv(out / "importance.csv")
    credit = credit_by_role(importance, roles_from_feature_sets(pd.read_csv(data_dir / "feature_sets.csv")),
                            [m for m in METHODS if m in set(importance["method"])])
    check = check_non_ancestor_share(credit, pd.read_csv(out / "summary.csv"), D.TRUTH_POPULATION)
    credit.to_csv(out / "credit_by_role.csv", index=False)
    print(f"wrote credit_by_role.csv ({len(credit)} rows) to {out}; largest difference from the summary's "
          f"non-ancestor share {check['difference'].abs().max():.2g}")


if __name__ == "__main__":
    main()
