import numpy as np
import pandas as pd
import pytest
import yaml
from pathlib import Path

from attribution_validation import table1 as T

LABELS = yaml.safe_load((Path(__file__).parents[1] / "config" / "labels.yaml").read_text())
PER_UNIT = LABELS["truth_type"]["per_unit"]


def test_labels_cover_every_key_used_by_the_tables():
    assert set(T.TRUTH_COLUMNS) <= set(LABELS["truth_type"])
    assert set(T.TABLE_METHODS) <= set(LABELS["methods"])
    assert {"full_set", "reference_subsample", "astronaut_set"} <= set(LABELS["datasets"])
    assert LABELS["truth_type"]["per_unit"] == "Intervention-contrast"
    assert "ng_style_shap_all_records" not in T.TABLE_METHODS


def test_competition_rank_shares_the_best_rank_within_a_group_and_skips():
    v = pd.Series({"a": 0.5, "b": 0.3, "c": 0.29, "d": 0.1})
    assert competition_rank_list(v, [["b", "c"]]) == [1, 2, 2, 4]
    assert competition_rank_list(v, None) == [1, 2, 3, 4]


def competition_rank_list(v, groups):
    return T.competition_rank(v, groups).tolist()


def toy_inputs():
    fsets = pd.DataFrame({"id": ["a", "b", "c", "n1", "n2"], "name": list("abc") + ["n1", "n2"],
                          "in_ancestor_set": [True, True, True, False, False], "in_all_features_set": True})
    rows = []
    for pop in ("source", "selected"):
        for tt, vals in (("per_unit", [0.30, -0.20, -0.19, 0.0, 0.0]), ("pop", [0.1, 0.2, 0.3, 0.0, 0.0])):
            for f, v in zip(fsets["id"], vals):
                rows.append({"feature": f, "name": f, "ancestor_set": f in "abc", "truth_type": tt, "population": pop,
                             "value": v, "mcse": 0.0, "low": np.nan, "high": np.nan})
    truth = pd.DataFrame(rows)
    ties = pd.DataFrame([{"truth_type": tt, "population": pop, "scope": "with_era", "rule": "delta2", "group": g, "size": len(m),
                          "features": "; ".join(m), "max_abs_truth": 0, "min_abs_truth": 0}
                         for tt, pop, groups in (("per_unit", "source", [["a"], ["b", "c"]]), ("per_unit", "selected", [["a"], ["b", "c"]]),
                                                 ("pop", "source", [["c"], ["b"], ["a"]]), ("pop", "selected", [["c"], ["b"], ["a"]]))
                         for g, m in enumerate(groups, 1)])
    return fsets, truth, ties


def toy_importance(reps=20, swap_in_all_run=True):
    """Ancestor run: a > b > c. All-features run: b and a swapped (a known displacement), c unchanged; n1 large."""
    rng = np.random.default_rng(0)
    rows = []
    for ds in ("full_set", "reference_subsample", "astronaut_set"):
        for m in T.TABLE_METHODS:
            for r in range(1, reps + 1):
                noise = rng.normal(0, 0.001, 5)
                anc = {"a": 0.30, "b": 0.20, "c": 0.10}
                allr = ({"a": 0.20, "b": 0.30, "c": 0.10, "n1": 0.25, "n2": 0.02} if swap_in_all_run
                        else {**anc, "n1": 0.25, "n2": 0.02})
                for f, v in anc.items():
                    rows.append({"method": m, "feature": f, "importance": v + noise[0], "replicate": r, "dataset": ds, "feature_set": "ancestor"})
                for f, v in allr.items():
                    rows.append({"method": m, "feature": f, "importance": v + noise[1], "replicate": r, "dataset": ds, "feature_set": "all"})
    return pd.DataFrame(rows)


def test_truth_ranks_share_tie_groups_and_non_ancestors_share_the_bottom_rank():
    fsets, truth, ties = toy_inputs()
    b = T.truth_block(truth, ties, fsets, "full_set", "per_unit").set_index("feature")
    assert b["rank"].to_dict() == {"a": 1, "b": 2, "c": 2, "n1": 4, "n2": 4}     # bottom rank follows the 3 ancestors
    assert b["share_pct"].sum() == pytest.approx(100) and b.loc["n1", "share_pct"] == 0
    p = T.truth_block(truth, ties, fsets, "astronaut_set", "pop").set_index("feature")
    assert p["rank"].to_dict() == {"c": 1, "b": 2, "a": 3, "n1": 4, "n2": 4}


def test_method_block_means_ranks_and_share_mcse():
    fsets, _, _ = toy_inputs()
    imp = toy_importance()
    b = T.method_block(imp, fsets, "full_set", "standard_shap").set_index("feature")
    assert b["n_replicates"].eq(20).all() and b["share_pct"].sum() == pytest.approx(100)
    assert b["rank"].to_dict() == {"a": 3, "b": 1, "c": 4, "n1": 2, "n2": 5}
    assert (b["mcse_share_pct"] > 0).all() and (b["mcse_share_pct"] < 0.5).all()


def test_known_swap_is_flagged_as_displaced_and_tau_is_negative_for_it():
    fsets, _, _ = toy_inputs()
    d, s = T.displacement(toy_importance(), fsets, LABELS, ["full_set"])
    one = d[(d["dataset"] == "full_set") & (d["method"] == "standard_shap")].set_index("feature")
    assert one["rank_ancestor_run"].to_dict() == {"a": 1, "b": 2, "c": 3}
    assert one["rank_all_run_among_ancestors"].to_dict() == {"a": 2, "b": 1, "c": 3}      # non-ancestors removed before ranking
    assert one["displaced"].to_dict() == {"a": True, "b": True, "c": False}
    row = s[(s["method"] == "standard_shap")].iloc[0]
    assert row["n_displaced"] == 2 and row["n_replicates"] == 20 and row["n_ancestors"] == 3
    assert row["tau_b_mean"] == pytest.approx(1 / 3, abs=1e-9)                            # one discordant pair of three
    assert row["tau_b_mcse"] == pytest.approx(0, abs=1e-9)


def test_no_displacement_when_the_orderings_agree():
    fsets, _, _ = toy_inputs()
    d, s = T.displacement(toy_importance(swap_in_all_run=False), fsets, LABELS, ["full_set"])
    assert not d["displaced"].any() and s["tau_b_mean"].to_numpy() == pytest.approx(1)


def test_wide_table_orders_by_first_truth_and_marks_ancestors_and_displaced_cells():
    fsets, truth, ties = toy_inputs()
    imp = toy_importance()
    long = T.table1_long(imp, truth, ties, fsets, LABELS, ["full_set"])
    disp, summ = T.displacement(imp, fsets, LABELS, ["full_set"])
    assert len(long) == 5 * 6 and set(long["column"]) == set(T.TRUTH_COLUMNS) | set(T.TABLE_METHODS)
    w = T.wide_table(long, disp)["full_set"]
    assert list(w.columns) == ["Feature", PER_UNIT, "Population-scaled", "Standard SHAP", "Ordering-only SHAP",
                               "Causal predictive SHAP", "Causal interventional SHAP"]
    assert w["Feature"].tolist()[:3] == ["a †", "b †", "c †"]                              # ordered by per-unit truth rank, ancestors marked
    assert w.loc[0, "Standard SHAP"].endswith("*") and w.loc[1, "Standard SHAP"].endswith("*") and not w.loc[2, "Standard SHAP"].endswith("*")
    assert not w.loc[0, PER_UNIT].endswith("*")                                  # truth columns are never starred
    assert not w.loc[3, "Standard SHAP"].endswith("*") and "†" not in w.loc[3, "Feature"]  # non-ancestor
    assert T.cell_text(3, 12.345, True) == "3 (12.3%)*" and T.cell_text(1, 0.0) == "1 (0.0%)"
    md = T.markdown({"full_set": w}, LABELS, summ)
    assert f"### {LABELS['datasets']['full_set']}" in md and f"| Feature | {PER_UNIT} |" in md and "†" in md
