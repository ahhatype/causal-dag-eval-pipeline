import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import supp_tables as S  # noqa: E402

N_DATASETS = len(S.style.DATASETS)


def attribution_rows(scope):
    """one row per truth x data set (own population, plus the astronaut-set on the source population) x method"""
    s = pd.read_csv(S.ATTR / "summary.csv")
    s = s[(s.scope == scope) & s.method.isin(S.METHOD_ORDER6)]
    own = s[s.population == s.dataset.map(S.OWN)]
    src = s[(s.dataset == "astronaut_set") & (s.population == "source")]
    return len(pd.concat([own, src]).drop_duplicates(["truth_type", "dataset", "population", "method"]))


def reported_contrasts():
    paired, cond, structure = S.contrast_files()
    p = paired[paired.truth_type.isin(S.CO_PRIMARY) & (paired.population == paired.dataset.map(S.OWN))
               & (paired.feature_set == "ancestor")]
    c = cond[cond.truth_type.isin(S.CO_PRIMARY) & cond.metric.isin(S.COND_METRICS)
             & cond.contrast.isin(["sample_size_effect", "selection_effect"])]
    return len(p) + len(c) + len(structure)


def test_supplementary_data_keep_every_row():
    n_edges, n_contrasts = S.write_supplementary_data()
    assert n_edges == len(pd.read_csv(S.STRUCTURE / "edge_recovery.csv"))
    assert n_contrasts == sum(len(f) for f in S.contrast_files())


def expected_counts():
    nodes = len(pd.read_csv(S.DATA / "feature_sets.csv"))
    edges = len(pd.read_csv(S.CFG / "edges.csv")) + len(pd.read_csv(S.CFG / "simulation_edges.csv"))
    power = pd.read_csv(S.STRUCTURE / "power_by_edge.csv")["edge"].nunique()
    true_edges = len(pd.read_csv(S.STRUCTURE / "difference_dag_edges.csv").query("true_edge"))
    conf = len(pd.read_csv(S.ATTR / "confounding_summary.csv")) * 2
    contrasts = sum(len(f) for f in S.contrast_files())
    return {
        "1": nodes + edges,
        "2": N_DATASETS + power,
        "4": len(S.VARIANTS) * N_DATASETS,
        "5": conf,
        "6": attribution_rows("with_era") * len(S.FEATURE_SETS6),
        "6b": attribution_rows("without_era") * len(S.FEATURE_SETS6),
        "7": reported_contrasts(),
        "8": len(S.LABELS["truth_type"]) * N_DATASETS * 4,
        "9": len(S.METHOD_ORDER9) * N_DATASETS,
    }


def colspec(tex):
    return re.search(r"\\begin\{(?:tabular|longtable)\}\{(.*)\}", tex).group(1)


@pytest.mark.parametrize("k", sorted(S.TABLES))
def test_table_builds_with_expected_rows(k):
    tex, n = S.TABLES[k][1]()
    assert n == expected_counts()[k]
    env = "longtable" if r"\begin{longtable}" in tex else "tabular"
    assert tex.count(rf"\begin{{{env}}}") == 1 and tex.count(rf"\end{{{env}}}") == 1
    assert (env == "longtable" or r"\begin{threeparttable}" in tex) and r"\toprule" in tex and r"\bottomrule" in tex
    assert tex.count(r"\begin{landscape}") == tex.count(r"\end{landscape}") <= 1


@pytest.mark.parametrize("k", sorted(S.TABLES))
def test_cells_per_row_match_column_spec(k):
    tex, _ = S.TABLES[k][1]()
    spec = colspec(tex)
    ncols = len(re.findall(r"[lrc]", re.sub(r"@\{[^}]*\}", "", spec)))
    if r"\begin{longtable}" in tex:
        body = tex.split(r"\endfoot")[1].split(r"\end{longtable}")[0]
    else:
        body = tex.split(r"\toprule")[1].split(r"\bottomrule")[0]
    for line in body.splitlines():
        if line.rstrip().endswith(r"\\") and r"\multicolumn" not in line:
            assert len(re.findall(r"(?<!\\)&", line)) + 1 == ncols, line


def test_attribution_tables_pair_feature_sets_row_by_row():
    s = pd.read_csv(S.ATTR / "summary.csv")
    assert set(s.feature_set) == {"all", "ancestor"}
    for build, scope in ((S.build_t6, "with_era"), (S.build_t6b, "without_era")):
        tex, n = build()
        spec = colspec(tex)
        assert len(re.findall(r"[lrc]", re.sub(r"@\{[^}]*\}", "", spec))) == 3 + len(S.METRICS6)
        truths = len(S.LABELS["truth_type"])
        assert tex.count("Astronaut-set, source population &") == truths
        assert tex.count("Astronaut-set &") == truths
        assert n == truths * 4 * len(S.METHOD_ORDER6) * len(S.FEATURE_SETS6)
        assert "the model was also given the 19 non-ancestors as inputs" in tex
        methods = re.findall(r"& (Standard SHAP|Ordering-only SHAP|Causal interventional SHAP|Causal predictive SHAP) &", tex)
        assert methods[:4] == [S.LABELS["methods"][m] for m in S.METHOD_ORDER6]
        lines = [l for l in tex.splitlines() if " & All & " in l or " & Anc. & " in l]
        assert all(" & All & " in a and " & Anc. & " in b for a, b in zip(lines[::2], lines[1::2]))


def test_attribution_cells_match_summary_and_non_ancestor_share_is_all_features_only():
    tex, _ = S.build_t6()
    s = pd.read_csv(S.ATTR / "summary.csv")
    key = dict(scope="with_era", dataset="full_set", population="source", truth_type="per_unit", method="standard_shap")
    lines = tex.splitlines()
    i = next(i for i, l in enumerate(lines) if l.startswith("Full-set & Standard SHAP & All &"))
    for fset, line in (("all", lines[i]), ("ancestor", lines[i + 1])):
        cells = [c.replace(r" \\", "") for c in line.split(" & ")[3:]]
        r = s[(s.feature_set == fset) & np.logical_and.reduce([s[k] == v for k, v in key.items()])].iloc[0]
        for mi, (m, _) in enumerate(S.METRICS6):
            if m == "non_ancestor_share" and fset == "ancestor":
                assert cells[mi] == "--"
            else:
                assert cells[mi] == S.ms(r[f"{m}_mean"], r[f"{m}_mcse"])


def test_display_labels_from_labels_yaml():
    tex, _ = S.build_t6()
    assert "Intervention-contrast" in tex and "Full-set" in tex
    assert "per_unit" not in tex and "full_set" not in tex


def test_credit_table_lists_methods_in_order_and_marks_zero_by_construction():
    tex, n = S.build_t9()
    names = re.findall(r"^(Standard SHAP|Ordering-only SHAP|Causal predictive SHAP|Causal interventional SHAP) &",
                       tex, re.M)
    assert names == [S.LABELS["methods"][m] for m in S.METHOD_ORDER9]
    credit = pd.read_csv(S.ATTR / "credit_by_role.csv")
    assert n == credit[credit["method"].isin(S.METHOD_ORDER9)].groupby(["method", "dataset"]).ngroups
    pred = [l for l in tex.splitlines() if "0.0 (0.0)" in l and S.DAGGER in l]
    assert len(pred) == N_DATASETS and "Zero by construction" in tex
