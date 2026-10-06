from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from attribution_validation import data as D
from attribution_validation.evaluation import (kendall_tau_b, non_ancestor_share, proximity_bias_index,
                                               tie_adjusted_truth, top_k_recovery)
from attribution_validation.interventional import (LinearLogisticSCM, NodeSpec, OrderSampler, ancestral_pairs, fit_scm,
                                                   interventional_shap)
from attribution_validation.ng import (CausalSampler, antithetic_orders, causal_weights, edge_strengths, ng_causal_shap,
                                      permutation_shapley)
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
    assert (len(AG.nodes), len(AG.edges)) == (36, 59)
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
    res = ng_causal_shap(predict, df, df.iloc[:10], feats, ag, samples=32, orders=8, seed=3, expected_value=expected)
    assert np.max(np.abs(res.efficiency_error)) < 1e-10
    assert np.all(res.values["d"] == 0)
    assert list(res.diagnostics.columns) == ["shapley_total", "weighted_total", "rescale_factor"]


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
    params["ng_causal_shap"].update(samples=8, orders=4)
    params["ordering_only_shap"].update(orders=4, burn_in=50, thin=5)
    params["interventional_shap"].update(orders=4, burn_in=50, thin=5)
    return params


needs_data = pytest.mark.skipif(not (DATA_DIR / "replicates" / "r001").exists(), reason="generated data not present")


@needs_data
def test_one_run_end_to_end_with_small_budgets():
    from attribution_validation.pipeline import run_one
    params = small_params()
    df = D.load_replicate(DATA_DIR, 1, "reference_subsample")
    fs = D.feature_sets(DATA_DIR)
    truth = D.ground_truth(DATA_DIR, "full").loc[fs["ancestor"]]
    truth_se = D.ground_truth_se(DATA_DIR, "full").loc[fs["ancestor"]]
    res, imp = run_one(df, CFG, AG, fs["ancestor"], truth, truth_se, set(fs["ancestor"]), D.distance_to_outcome(AG), params, 7,
                       params["methods"])
    assert {r["method"] for r in res["rows"]} == set(params["methods"])
    assert all(r["max_efficiency_error"] < 1e-8 for r in res["rows"])
    assert res["info"]["explained_events"] >= params["explain"]["min_events"]
    assert set(imp["feature"]) == set(fs["ancestor"])


@needs_data
def test_parallel_cells_match_sequential(tmp_path):
    from attribution_validation.pipeline import run
    params = small_params()
    args = dict(replicates=1, datasets=["reference_subsample", "astronaut_set"], feature_sets=["ancestor"])
    run(params, DATA_DIR, tmp_path / "seq", workers=1, **args)
    run(params, DATA_DIR, tmp_path / "par", workers=2, **args)
    seq, par = (pd.read_csv(tmp_path / d / "per_run.csv").drop(columns=["seconds", "fit_seconds"]) for d in ("seq", "par"))
    assert len(seq) == 2 * len(params["methods"])
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
