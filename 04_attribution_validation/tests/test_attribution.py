from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from attribution_validation import data as D
from attribution_validation.evaluation import (kendall_tau_b, l1_distance, margin_tie_groups, mcse_tie_groups,
                                               non_ancestor_share, proximity_bias_index, tie_adjusted,
                                               tie_adjusted_truth, top_k_recovery, weighted_tau)
from attribution_validation.interventional import (LinearLogisticSCM, NodeSpec, OrderSampler, ancestral_pairs, fit_scm,
                                                   interventional_shap)
from attribution_validation.ng import (CausalSampler, antithetic_orders, causal_weights, edge_strengths,
                                      fitted_ng_inputs, ng_style_shap, permutation_shapley)
from attribution_validation.ordering import ordering_only_shap
from attribution_validation.standard import permutation_shap
from attribution_validation.superlearner import SuperLearner, default_library

CFG = D.load_generation_config()
AG = D.analysis_graph(CFG)
DATA_DIR = D.env_dir("DATA_DIR", "./02_data")


def toy_graph():
    """a -> b -> y, c -> y, d isolated; all continuous except y."""
    return D.AnalysisGraph(nodes=["a", "c", "d", "b", "nephrolithiasis"],
                           edges=[("a", "b"), ("b", "nephrolithiasis"), ("c", "nephrolithiasis")],
                           scale={"a": "continuous", "b": "continuous", "c": "continuous", "d": "continuous",
                                  "nephrolithiasis": "binary"})


def toy_data(n=4000, seed=0):
    rng = np.random.default_rng(seed)
    a, c, d = rng.normal(size=(3, n))
    b = 0.8 * a + rng.normal(scale=0.6, size=n)
    y = (rng.uniform(size=n) < 1 / (1 + np.exp(-(-1 + 1.0 * b + 0.5 * c)))).astype(float)
    return pd.DataFrame({"a": a, "c": c, "d": d, "b": b, "nephrolithiasis": y})


def test_analysis_graph_matches_design_and_shape_validation():
    assert (len(AG.nodes), len(AG.edges)) == (36, 63)
    assert ("individual_factors", "k_citrate") in AG.edges and ("pre_flight_fitness", "bisphosphonates") in AG.edges
    r = D.REPO / "outputs" / "03_shape_validation" / "analysis_graph_edges.csv"
    if r.exists():
        e = pd.read_csv(r)
        assert set(zip(e.parent, e.child)) == set(AG.edges)
    dist = D.distance_to_outcome(AG)
    assert dist["individual_factors"] == 1 and dist["mineralized_renal_material"] == 1


def test_order_sampler_is_uniform_over_consistent_orders():
    rng = np.random.default_rng(1)
    s = OrderSampler(["a", "b", "c"], {("a", "b")}, rng, burn_in=500, thin=20)
    counts = Counter(tuple(s.sample()) for _ in range(6000))
    assert set(counts) == {("a", "b", "c"), ("a", "c", "b"), ("c", "a", "b")}
    for v in counts.values():
        assert abs(v / 6000 - 1 / 3) < 0.03


def test_ancestral_pairs_pass_through_non_features():
    pairs = ancestral_pairs(AG, ["individual_factors", "ureterolithiasis", "pre_flight_fitness"])
    assert ("individual_factors", "ureterolithiasis") in pairs
    assert not any(p[0] == "ureterolithiasis" for p in pairs)


def test_scm_round_trip_and_fit():
    df = toy_data()
    scm = fit_scm(df, toy_graph())
    u = scm.recover_exogenous(df, np.random.default_rng(0))
    sim = scm.simulate(u)
    for col in ["a", "b", "c", "d", "nephrolithiasis"]:
        assert np.allclose(sim[col], df[col])
    assert abs(scm.specs["b"].coefficients[0] - 0.8) < 0.05


def test_interventional_credits_ancestor_and_satisfies_efficiency():
    df = toy_data()
    ag = toy_graph()
    feats = ["a", "b", "c", "d"]
    predict = lambda X: 1 / (1 + np.exp(-(-1 + X[:, 1] + 0.5 * X[:, 2])))  # noqa: E731  uses b and c only
    res = interventional_shap(predict, fit_scm(df, ag), df.iloc[:20], df.iloc[100:164], feats, ag, orders=64,
                              burn_in=100, thin=5, seed=2)
    assert np.max(np.abs(res.efficiency_error)) < 1e-10
    imp = res.values.abs().mean()
    assert imp["a"] > 0.01          # a acts through b, so it receives credit
    assert imp["d"] < 1e-12         # neither in the model nor an ancestor of anything in it


def test_ordering_alone_gives_no_credit_to_an_ancestor_the_model_ignores():
    df = toy_data()
    ag = toy_graph()
    feats = ["a", "b", "c", "d"]
    predict = lambda X: 1 / (1 + np.exp(-(-1 + X[:, 1] + 0.5 * X[:, 2])))  # noqa: E731  uses b and c only
    res = ordering_only_shap(predict, df.iloc[:20], df.iloc[100:164], feats, ag, orders=32, burn_in=100, thin=5, seed=2)
    assert np.max(np.abs(res.efficiency_error)) < 1e-10
    imp = res.values.abs().mean()
    assert imp["a"] < 1e-12         # a is in no coalition's value: the model does not use it and nothing propagates
    assert imp["b"] > 0.01


def test_chain_matches_the_closed_form():
    """x -> m with m = 0.8 x + noise, and a model that uses m only: do-propagation credits x with 0.8 (x - E x)."""
    rng = np.random.default_rng(5)
    x = rng.normal(size=6000)
    df = pd.DataFrame({"x": x, "m": 0.8 * x + rng.normal(scale=0.6, size=6000)})
    ag = D.AnalysisGraph(nodes=["x", "m"], edges=[("x", "m")], scale={"x": "continuous", "m": "continuous"})
    predict = lambda X: X[:, 1]  # noqa: E731
    explain, bg = df.iloc[:30], df.iloc[100:228]
    res = interventional_shap(predict, fit_scm(df, ag), explain, bg, ["x", "m"], ag, orders=4, burn_in=20, thin=2, seed=1)
    expected_x = 0.8 * (explain["x"] - bg["x"].mean())
    assert np.allclose(res.values["x"], expected_x, atol=0.05)
    assert np.allclose(res.values["m"], explain["m"] - bg["m"].mean() - res.values["x"], atol=1e-8)
    ordering = ordering_only_shap(predict, explain, bg, ["x", "m"], ag, orders=4, burn_in=20, thin=2, seed=1)
    assert np.allclose(ordering.values["x"], 0, atol=1e-12)


def test_ng_weights_zero_for_non_ancestors_and_local_accuracy():
    df = toy_data()
    ag = toy_graph()
    feats = ["a", "b", "c", "d"]
    w = causal_weights(ag, edge_strengths(df, ag), feats)
    assert w["d"] == 0 and abs(sum(w.values()) - 1) < 1e-12
    assert w["a"] < w["b"]          # a's path runs through b, product of strengths below 1
    predict = lambda X: 1 / (1 + np.exp(-(-1 + X[:, 1] + 0.5 * X[:, 2])))  # noqa: E731
    expected = float(np.mean(predict(df[feats].to_numpy())))
    res = ng_style_shap(predict, df.iloc[:10], feats, *fitted_ng_inputs(df, ag, feats), samples=32, orders=8, seed=3,
                        expected_value=expected)
    assert np.max(np.abs(res.efficiency_error)) < 1e-10
    assert np.all(res.values["d"] == 0)
    d = res.diagnostics
    assert list(d.columns) == ["shapley_total", "weighted_total", "rescale_factor", *[f"pre_{f}" for f in feats]]
    pre = d[[f"pre_{f}" for f in feats]].to_numpy()
    assert np.allclose(pre.sum(axis=1), d["weighted_total"])
    assert np.allclose(pre * d["rescale_factor"].to_numpy()[:, None], res.values.to_numpy())


def test_permutation_shapley_over_all_orders_is_exact():
    from itertools import combinations, permutations
    from math import factorial
    k = 4
    rng = np.random.default_rng(5)
    table = {S: rng.normal() for r in range(k + 1) for S in combinations(range(k), r)}
    value = lambda C: np.array([table[tuple(np.flatnonzero(c))] for c in C])  # noqa: E731
    exact = np.zeros(k)
    for i in range(k):
        others = [j for j in range(k) if j != i]
        for r in range(k):
            for S in combinations(others, r):
                w = factorial(r) * factorial(k - r - 1) / factorial(k)
                exact[i] += w * (table[tuple(sorted(S + (i,)))] - table[S])
    est = permutation_shapley(value, k, [np.array(p) for p in permutations(range(k))])
    assert np.allclose(est, exact, atol=1e-12)
    orders = antithetic_orders(k, 6, rng)
    assert len(orders) == 6 and np.array_equal(orders[1], orders[0][::-1])


def test_causal_sampler_keeps_binaries_binary_and_shares_noise():
    rng = np.random.default_rng(6)
    n = 3000
    a = rng.normal(size=n)
    b = (rng.uniform(size=n) < 1 / (1 + np.exp(-a))).astype(float)
    df = pd.DataFrame({"a": a, "b": b, "nephrolithiasis": b})
    ag = D.AnalysisGraph(nodes=["a", "b", "nephrolithiasis"], edges=[("a", "b"), ("b", "nephrolithiasis")],
                         scale={"a": "continuous", "b": "binary", "nephrolithiasis": "binary"})
    s = CausalSampler(df, ag, ["a", "b"])
    noise = s.noise(500, rng)
    x = {"a": 2.0, "b": 1.0}
    C = np.array([[False, False], [True, False], [True, True]])
    draws = s.draw(x, C, noise).reshape(3, 500, 2)
    assert set(np.unique(draws[..., 1])) <= {0.0, 1.0}
    assert np.array_equal(draws[0, :, 0], draws[2, :, 0] * 0 + noise["a"])   # same root draws in every coalition
    assert np.all(draws[1, :, 0] == 2.0) and np.all(draws[2] == [2.0, 1.0])
    assert draws[1, :, 1].mean() > draws[0, :, 1].mean()                    # fixing a high raises b downstream


def test_permutation_shap_efficiency():
    df = toy_data()
    feats = ["a", "b", "c", "d"]
    predict = lambda X: 1 / (1 + np.exp(-(-1 + X[:, 1] + 0.5 * X[:, 2])))  # noqa: E731
    res = permutation_shap(predict, df.iloc[:10], df.iloc[100:164], feats, permutations=16, seed=4)
    assert np.max(np.abs(res.efficiency_error)) < 1e-10


def test_super_learner_weights_form_a_convex_combination():
    df = toy_data(2000)
    X, y = df[["a", "b", "c", "d"]].to_numpy(), df["nephrolithiasis"].to_numpy()
    sl = SuperLearner(default_library(0), folds=5, seed=0).fit(X, y)
    w = np.array(list(sl.weights_.values()))
    assert np.all(w >= 0) and abs(w.sum() - 1) < 1e-8
    p = sl.predict(X)
    assert p.min() >= 0 and p.max() <= 1


def test_metrics():
    truth = {"a": 0.3, "b": -0.2, "c": 0.1}
    assert kendall_tau_b({"a": 3, "b": 2, "c": 1}, truth) == pytest.approx(1)
    assert top_k_recovery({"a": 1, "b": 3, "c": 2}, truth, 2) == pytest.approx(0.5)
    dist = {"a": 3, "b": 2, "c": 1}
    assert proximity_bias_index({"a": 0, "b": 0, "c": 1}, truth, dist) > 0     # all credit next to the outcome
    assert proximity_bias_index({"a": 0.3, "b": 0.2, "c": 0.1}, truth, dist) == pytest.approx(0)
    assert non_ancestor_share({"a": 1, "x": 3}, {"a"}) == pytest.approx(0.75)


def test_tied_true_effects_are_not_ordered():
    effect = {"a": 0.30, "b": -0.01040, "c": -0.01035, "d": 0.001}
    se = {"a": 1e-4, "b": 4e-5, "c": 4e-5, "d": 1e-4}
    adj = tie_adjusted_truth(effect, se)
    assert adj["b"] == adj["c"] == pytest.approx(0.010375) and adj["a"] == 0.30 and adj["d"] == 0.001
    assert kendall_tau_b({"a": 4, "b": 3, "c": 2, "d": 1}, adj) == kendall_tau_b({"a": 4, "b": 2, "c": 3, "d": 1}, adj)
    # a tie straddling the top-k boundary: both members count as correct
    assert top_k_recovery({"a": 4, "b": 1, "c": 3, "d": 2}, adj, 2) == pytest.approx(1.0)
    assert top_k_recovery({"a": 4, "b": 3, "c": 1, "d": 2}, adj, 2) == pytest.approx(1.0)
    assert tie_adjusted_truth({"a": 0.3, "b": 0.1}, {"a": 1e-4, "b": 1e-4}) == {"a": 0.3, "b": 0.1}


def test_stale_data_refused(tmp_path):
    (tmp_path / "provenance.txt").write_text("config hash: not-current\nreplicates: 2\n")
    with pytest.raises(RuntimeError, match="make data"):
        D.check_data_current(tmp_path, CFG)


def small_params():
    import yaml
    params = yaml.safe_load((Path(__file__).parents[1] / "config" / "params.yaml").read_text())
    params["super_learner"]["folds"] = 3
    params["explain"].update(records=12, background=16)
    params["standard_shap"]["permutations"] = 4
    params["ng_style_shap"].update(samples=8, orders=4)
    params["oracle"].update(pool=2000, background=16)
    params["ordering_only_shap"].update(orders=4, burn_in=50, thin=5)
    params["interventional_shap"].update(orders=4, burn_in=50, thin=5)
    return params


needs_data = pytest.mark.skipif(not (DATA_DIR / "replicates" / "r001").exists(), reason="generated data not present")


def small_context(params):
    from attribution_validation.pipeline import truth_context
    fs = D.feature_sets(DATA_DIR)
    return {"cfg": CFG, "ag": AG, "fsets": fs, "distance": D.distance_to_outcome(AG),
            "truth": truth_context(D.truth_values(DATA_DIR), fs["ancestor"], params)}


@needs_data
def test_one_run_end_to_end_with_small_budgets_and_oracle():
    from attribution_validation.pipeline import oracle_context, run_one
    params = small_params()
    df = D.load_replicate(DATA_DIR, 1, "reference_subsample")
    ctx = small_context(params)
    anc = ctx["fsets"]["ancestor"]
    res, imp = run_one(df, ctx, anc, params, 7, params["methods"], oracle_context(CFG, AG, anc, params))
    rows = pd.DataFrame(res["rows"])
    assert set(rows["method"]) == set(params["methods"]) | {"ng_style_shap_all_records"}
    assert len(rows) == (len(params["methods"]) + 1) * len(D.TRUTH_TYPES) * len(D.POPULATIONS) * 2
    assert (rows["max_efficiency_error"] < 1e-8).all()
    assert rows["kendall_tau_b"].between(-1, 1).all() and rows["l1_distance"].between(0, 2).all()
    assert res["info"]["explained_events"] == int(res["values"].drop_duplicates("record")["outcome"].sum())
    assert set(imp["feature"]) == set(anc)
    o = pd.DataFrame(res["oracle"])
    assert set(o["method"]) == set(params["methods"]) and len(o) == len(params["methods"]) * 6
    assert o[["tau_oracle_truth", "tau_fitted_truth", "tau_fitted_oracle"]].abs().le(1).all().all()


@needs_data
def test_parallel_cells_match_sequential(tmp_path):
    from attribution_validation.pipeline import run
    params = small_params()
    params["oracle"]["replicates"] = 1
    args = dict(replicates=1, datasets=["reference_subsample", "astronaut_set"], feature_sets=["ancestor"])
    run(params, DATA_DIR, tmp_path / "seq", workers=1, **args)
    run(params, DATA_DIR, tmp_path / "par", workers=2, **args)
    seq, par = (pd.read_csv(tmp_path / d / "per_run.csv").drop(columns=["seconds", "fit_seconds"]) for d in ("seq", "par"))
    assert len(seq) == 2 * (len(params["methods"]) + 1) * len(D.TRUTH_TYPES) * len(D.POPULATIONS) * 2
    for f in ("oracle_decomposition.csv", "tie_groups.csv", "contrasts.csv", "summary.csv"):
        assert (tmp_path / "seq" / f).exists()
    pd.testing.assert_frame_equal(seq, par, rtol=1e-8)
    imp = [pd.read_csv(tmp_path / d / "importance.csv") for d in ("seq", "par")]
    pd.testing.assert_frame_equal(*imp, rtol=1e-8)


def test_cell_directory_is_never_silently_reused_or_cleared(tmp_path):
    from attribution_validation.pipeline import _prepare_cell_dir
    cells = tmp_path / "cells"
    want = {"params": {"seed": 1}, "code_hash": "a"}
    _prepare_cell_dir(cells, want, resume=False, fresh=False)
    (cells / "r001_x_y_runs.csv").write_text("done")
    with pytest.raises(RuntimeError, match="--resume"):
        _prepare_cell_dir(cells, want, resume=False, fresh=False)
    _prepare_cell_dir(cells, want, resume=True, fresh=False)
    assert (cells / "r001_x_y_runs.csv").exists()
    with pytest.raises(RuntimeError, match="different settings"):
        _prepare_cell_dir(cells, {**want, "code_hash": "b"}, resume=True, fresh=False)
    _prepare_cell_dir(cells, {**want, "code_hash": "b"}, resume=False, fresh=True)
    assert not (cells / "r001_x_y_runs.csv").exists()


def test_margin_ties_chain_neighbours_and_mcse_ties_use_combined_error():
    v = {"a": 0.30, "b": 0.295, "c": 0.290, "d": 0.10, "e": -0.099}
    assert margin_tie_groups(v, 0.006) == [["a", "b", "c"], ["d", "e"]]     # chained: a and c differ by 0.01
    assert margin_tie_groups(v, 0.001) == [["a"], ["b"], ["c"], ["d"], ["e"]]
    adj = tie_adjusted(v, margin_tie_groups(v, 0.006))
    assert adj["a"] == adj["c"] == pytest.approx(0.295) and adj["e"] == pytest.approx(0.0995)
    g = mcse_tie_groups({"a": 0.3, "b": 0.299, "c": 0.1}, {"a": 1e-3, "b": 1e-3, "c": 1e-3})
    assert g == [["a", "b"], ["c"]]


def test_weighted_tau_weights_top_disagreements_more():
    truth = {f: v for f, v in zip("abcdef", [6, 5, 4, 3, 2, 1])}
    assert weighted_tau(dict(truth), truth) == pytest.approx(1)
    top_swap = {**truth, "a": 5, "b": 6}
    bottom_swap = {**truth, "e": 1, "f": 2}
    assert weighted_tau(top_swap, truth) < weighted_tau(bottom_swap, truth) < 1
    assert kendall_tau_b(top_swap, truth) == pytest.approx(kendall_tau_b(bottom_swap, truth))


def test_l1_distance_compares_shares():
    truth = {"a": 0.3, "b": -0.1}
    assert l1_distance({"a": 3, "b": 1}, truth) == pytest.approx(0)
    assert l1_distance({"a": 0, "b": 1}, truth) == pytest.approx(1.5)
    assert l1_distance({"a": 1, "b": 0}, {"a": 0, "b": 1}) == pytest.approx(2)


def test_true_scm_reproduces_simcausal_records():
    from attribution_validation.truesim import TrueSCM
    f = DATA_DIR / "scm_check.csv"
    if not f.exists():
        pytest.skip("scm_check.csv not generated")
    s = pd.read_csv(f)
    scm = TrueSCM(CFG)
    v = scm.simulate(scm.recover_exogenous(s))
    for node in scm.order:
        assert np.allclose(v[node], s[node], atol=1e-9), node
    assert np.allclose(scm.outcome_probability(v), D.true_outcome_probability(s, CFG))


def test_true_sampler_fixes_coalition_and_propagates_downstream():
    from attribution_validation.truesim import TrueSampler, TrueSCM
    scm = TrueSCM(CFG)
    feats = ["hydration", "urine_concentration", "water_intake"]
    s = TrueSampler(scm, feats)
    noise = s.noise(400, np.random.default_rng(3))
    C = np.array([[False, False, False], [True, False, False], [True, True, True]])
    d = s.draw({"hydration": 2.0, "urine_concentration": 0.0, "water_intake": -1.0}, C, noise).reshape(3, 400, 3)
    assert np.all(d[1, :, 0] == 2.0) and np.all(d[2] == [2.0, 0.0, -1.0])
    assert np.array_equal(d[0, :, 2], d[1, :, 2])                 # water intake is upstream of hydration: unchanged
    assert d[1, :, 1].mean() < d[0, :, 1].mean()                  # more hydration, lower urine concentration


def test_paired_contrasts_classify_by_both_truths():
    from attribution_validation.pipeline import paired_contrasts
    rows = []
    for r in range(1, 21):
        for tt, vf in (("per_unit", 0.05), ("pop", -0.05), ("rec", 0.0)):
            o = 0.6 + 0.01 * (r % 3)
            for m, t in (("standard_shap", 0.4), ("ordering_only_shap", o),
                         ("interventional_shap", o + vf + 0.01 * (r % 2))):
                rows.append({"dataset": "reference_subsample", "feature_set": "ancestor", "truth_type": tt,
                             "population": "source", "scope": "with_era", "replicate": r, "method": m,
                             "kendall_tau_b": t})
    c = paired_contrasts(pd.DataFrame(rows))
    order = c[(c["contrast"] == "order_effect") & (c["truth_type"] == "per_unit")]
    assert order["mean"].iloc[0] == pytest.approx(0.21, abs=0.01) and order["classification"].iloc[0] == "robust"
    vf = c[c["contrast"] == "value_function_effect"]
    assert (vf["classification"] == "truth_dependent").all()       # opposite signs under the two truths
    flat = pd.DataFrame([{**r, "kendall_tau_b": 0.4 + (0.01 if r["method"] != "standard_shap" else 0) * (-1) ** r["replicate"]}
                         for r in rows])
    assert (paired_contrasts(flat)["classification"] == "inconclusive").all()      # no effect under either truth


def _order_effect_runs(per_unit, pop, rec=None, reps=50, dataset="reference_subsample", population="source"):
    """per-run frame in which the order effect equals the given per-replicate values under each truth type."""
    rows = []
    for r in range(1, reps + 1):
        for tt, f in (("per_unit", per_unit), ("pop", pop), ("rec", rec or pop)):
            for m, t in (("standard_shap", 0.3), ("ordering_only_shap", 0.3 + f(r))):
                rows.append({"dataset": dataset, "feature_set": "ancestor", "truth_type": tt, "population": population,
                             "scope": "with_era", "replicate": r, "method": m, "kendall_tau_b": t})
    return pd.DataFrame(rows)


def _co_primary_class(runs):
    from attribution_validation.pipeline import paired_contrasts
    c = paired_contrasts(runs)
    cls = c[c["contrast"] == "order_effect"]["classification"].unique()
    assert len(cls) == 1
    return cls[0], c


def test_contrast_class_one_interval_excludes_zero_and_null_paired_difference_is_inconclusive():
    sign = lambda r: (-1) ** r                                    # noqa: E731
    runs = _order_effect_runs(per_unit=lambda r: 0.012 + 0.01 * sign(r), pop=lambda r: 0.012 + 0.1 * sign(r))
    cls, c = _co_primary_class(runs)
    pu = c[(c["truth_type"] == "per_unit") & (c["contrast"] == "order_effect")].iloc[0]
    po = c[(c["truth_type"] == "pop") & (c["contrast"] == "order_effect")].iloc[0]
    assert pu["ci_low"] > 0 and po["ci_low"] < 0 < po["ci_high"]       # one truth excludes zero, the other does not
    assert pu["truth_diff_ci_low"] < 0 < pu["truth_diff_ci_high"]      # and the truths do not differ
    assert cls == "inconclusive"


def test_contrast_class_robust_takes_precedence_over_a_significant_truth_difference():
    sign = lambda r: (-1) ** r                                    # noqa: E731
    runs = _order_effect_runs(per_unit=lambda r: 0.3 + 0.01 * sign(r), pop=lambda r: 0.1 + 0.01 * sign(r))
    cls, c = _co_primary_class(runs)
    d = c[(c["truth_type"] == "per_unit") & (c["contrast"] == "order_effect")].iloc[0]
    assert d["truth_diff_mean"] == pytest.approx(0.2) and d["truth_diff_ci_low"] > 0      # significant paired difference
    assert cls == "robust"


def test_contrast_class_significant_truth_difference_without_robustness_is_truth_dependent():
    sign = lambda r: (-1) ** r                                    # noqa: E731
    runs = _order_effect_runs(per_unit=lambda r: 0.05 + 0.01 * sign(r), pop=lambda r: 0.1 * sign(r))
    cls, c = _co_primary_class(runs)
    po = c[(c["truth_type"] == "pop") & (c["contrast"] == "order_effect")].iloc[0]
    d = c[(c["truth_type"] == "per_unit") & (c["contrast"] == "order_effect")].iloc[0]
    assert po["ci_low"] < 0 < po["ci_high"] and d["truth_diff_ci_low"] > 0
    assert cls == "truth_dependent"


def test_contrast_class_opposite_signs_with_both_intervals_excluding_zero_is_truth_dependent():
    sign = lambda r: (-1) ** r                                    # noqa: E731
    runs = _order_effect_runs(per_unit=lambda r: 0.05 + 0.01 * sign(r), pop=lambda r: -0.05 + 0.01 * sign(r))
    assert _co_primary_class(runs)[0] == "truth_dependent"


def test_contrast_class_both_null_is_inconclusive_and_rec_rows_are_labelled_sensitivity():
    sign = lambda r: (-1) ** r                                    # noqa: E731
    runs = _order_effect_runs(per_unit=lambda r: 0.05 * sign(r), pop=lambda r: 0.05 * sign(r),
                              rec=lambda r: 0.3 + 0.01 * sign(r))                   # rec would be robust, and is ignored
    cls, c = _co_primary_class(runs)
    assert cls == "inconclusive"
    assert set(c.loc[c["truth_type"] == "rec", "role"]) == {"sensitivity"}
    assert set(c.loc[c["truth_type"] != "rec", "role"]) == {"co_primary"}
    assert {"truth_diff_mean", "truth_diff_ci_low", "truth_diff_ci_high"} <= set(c.columns)
