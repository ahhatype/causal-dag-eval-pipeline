"""ConfoundingSHAP (Brockschmidt et al., 2026, arXiv:2605.10533) for one binary treatment and the outcome.

The signed global game v(S) = -E[b_S(X_S)] is shapiq 1.7.0's released `GlobalConfoundingXAI`: one S-learner of the
outcome on (X_S, A) per coalition, v(S) = mean(tau_hat) - mean(m_S(1) - m_S(0)), with tau_hat the full-set S-learner
contrast. Shapley credits come from shapiq's exact computer when the budget covers every coalition, otherwise from
RegressionMSR (the paper's approximator) with a decision-tree proxy: shapiq's default xgboost proxy cannot share a
process with torch on this Mac (two OpenMP runtimes segfault), and the proxy only affects the estimator's variance.
The learner is TabPFN with one estimator and the target power transform
disabled, as in the paper's Supplement D; the model version is explicit because tabpfn 9.x defaults to 3.5.

Covariates are the pretreatment ancestors of the outcome: no descendant of the treatment and no other drug.
"""

from __future__ import annotations

import contextlib
import os
import time
from typing import Callable, Sequence

import networkx as nx
import numpy as np
import pandas as pd

from .data import OUTCOME, AnalysisGraph, load_secret_env


def pretreatment_covariates(ag: AnalysisGraph, ancestors: Sequence[str], treatment: str, drugs: Sequence[str]) -> list[str]:
    """Ancestor-set features that are neither the treatment, a descendant of it, nor another drug."""
    desc = nx.descendants(ag.graph, treatment)
    return [f for f in ancestors if f != treatment and f not in desc and f not in set(drugs)]


def common_cause_sets(ag: AnalysisGraph, covariates: Sequence[str], treatment: str) -> dict[str, list[str]]:
    """Role labels for the covariates, from the analysis graph.

    direct: parents of the treatment that reach the outcome without passing through it.
    all: every ancestor of the treatment that also reaches the outcome without passing through it (adds upstream
    common causes, which can earn credit as proxies for the direct ones).
    """
    g = ag.graph.copy()
    g.remove_edges_from([(treatment, c) for c in list(g.successors(treatment))])
    reaches_y = nx.ancestors(g, OUTCOME)
    an_a = nx.ancestors(ag.graph, treatment)
    direct = [c for c in covariates if c in set(ag.graph.predecessors(treatment)) and c in reaches_y]
    every = [c for c in covariates if c in an_a and c in reaches_y]
    return {"direct": direct, "all": every}


def confounder_mass(credits: dict[str, float], confounders: Sequence[str]) -> float:
    """Share of absolute credit on the true confounders (paper, Supplement E, eq. 36)."""
    total = sum(abs(v) for v in credits.values())
    return float(sum(abs(credits[c]) for c in confounders) / total) if total > 0 else float("nan")


def confounder_recovery(credits: dict[str, float], confounders: Sequence[str]) -> float:
    """Share of the true confounders among the |C| covariates with the largest absolute credit (eq. 38)."""
    k = len(confounders)
    if k == 0:
        return float("nan")
    top = sorted(credits, key=lambda f: (-abs(credits[f]), f))[:k]
    return float(len(set(top) & set(confounders)) / k)


def tabpfn_factory(version: str, seed: int) -> Callable:
    """TabPFN regressor factory for shapiq's game: the requested model, one estimator by default, power transform off."""
    import torch
    from shapiq_games.benchmark.causal_xai import base as upstream
    from tabpfn import TabPFNRegressor
    from tabpfn.constants import ModelVersion

    load_secret_env()
    os.environ.setdefault("TABPFN_NO_BROWSER", "1")
    os.environ.setdefault("TABPFN_DISABLE_TELEMETRY", "1")
    os.environ.setdefault("TABPFN_ALLOW_CPU_LARGE_DATASET", "1")
    mv = ModelVersion(version)

    def make(device: str, n_estimators: int = 1):
        return TabPFNRegressor.create_default_for_version(
            mv, device=device, n_estimators=n_estimators, n_jobs=1, random_state=seed,
            inference_precision=torch.float32, inference_config=upstream._TABPFN_INFERENCE_CONFIG)

    return make


@contextlib.contextmanager
def learner(factory: Callable | None):
    """Swap the regressor shapiq's game fits per coalition (None keeps shapiq's own TabPFN default)."""
    from shapiq_games.benchmark.causal_xai import base as upstream
    original = upstream._make_tabpfn
    if factory is not None:
        upstream._make_tabpfn = factory
    try:
        yield
    finally:
        upstream._make_tabpfn = original


def confounding_credits(df: pd.DataFrame, covariates: Sequence[str], treatment: str, factory: Callable | None,
                        budget: int, seed: int, device: str = "cpu", n_estimators: int = 1, proxy: str = "tree") -> dict:
    """Signed global ConfoundingSHAP credits of each covariate for treatment -> outcome on one data set."""
    from shapiq import ExactComputer
    from shapiq.approximator import RegressionMSR
    from shapiq_games.benchmark.causal_xai import base as upstream

    cov = list(covariates)
    p = len(cov)
    X = df[cov].to_numpy(dtype=float)
    A = df[treatment].to_numpy(dtype=float)
    Y = df[OUTCOME].to_numpy(dtype=float)
    t0 = time.time()
    with learner(factory):
        m1, m0 = upstream._fit_s_learner(X, A, Y, device=device, n_estimators=n_estimators)
        tau_hat = m1(X) - m0(X)
        game = upstream.GlobalConfoundingXAI(X, A, Y, tau_hat, mode="signed", device=device, n_estimators=n_estimators)
        if budget >= 2 ** p:
            sv = ExactComputer(n_players=p, game=game)(index="SV", order=1)
            exact = True
        else:
            sv = RegressionMSR(n=p, index="SV", proxy_model=proxy, random_state=seed).approximate(budget=budget, game=game)
            exact = False
        credits = {c: float(sv[(j,)]) for j, c in enumerate(cov)}
        v_empty = float(game(np.zeros((1, p), dtype=bool))[0])
        v_full = float(game(np.ones((1, p), dtype=bool))[0])
    crude = float(Y[A == 1].mean() - Y[A == 0].mean())
    return {"credits": credits, "v_empty": v_empty, "v_full": v_full, "crude": crude, "adjusted_ate": float(tau_hat.mean()),
            "exact": exact, "coalitions_fitted": len(game._cache), "seconds": time.time() - t0,
            "n": len(df), "treated": int(A.sum()), "events": int(Y.sum())}
