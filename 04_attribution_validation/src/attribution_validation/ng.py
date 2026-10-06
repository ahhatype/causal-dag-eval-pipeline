"""Causal SHAP of Ng et al. (arXiv:2509.00846), with the analysis graph supplied in place of PC + IDA.

Attribution phi_i = gamma_i * Shapley_i(v_c), rescaled for local accuracy. The Shapley values of the causal value
function v_c are estimated from sampled feature orders (antithetic pairs), with common random numbers across the
coalitions of each order.

Adapted from causal-shap-spaceflight-renal-stones (python/src/causal_shap_renal/attribution_structural.py).
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from .data import OUTCOME, AnalysisGraph
from .interventional import AttributionResult, sigmoid


def edge_strengths(data: pd.DataFrame, ag: AnalysisGraph) -> dict[tuple[str, str], float]:
    """|coefficient| of each parent in one linear regression of the child on all its parents (standardized data)."""
    z = (data[ag.nodes] - data[ag.nodes].mean()) / data[ag.nodes].std().replace(0, 1)
    out = {}
    for node in ag.nodes:
        pa = ag.parents(node)
        if not pa:
            continue
        coef, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(z)), z[pa].to_numpy()]), z[node].to_numpy(), rcond=None)
        out.update({(p, node): abs(float(c)) for p, c in zip(pa, coef[1:])})
    return out


def causal_weights(ag: AnalysisGraph, strengths: dict, features: Sequence[str]) -> dict[str, float]:
    """gamma_i: normalized sum over directed paths to the outcome of the product of edge strengths."""
    g = ag.graph
    w = {}
    for f in features:
        w[f] = sum(math.prod(strengths[(a, b)] for a, b in zip(p[:-1], p[1:]))
                   for p in nx.all_simple_paths(g, f, OUTCOME)) if nx.has_path(g, f, OUTCOME) else 0.0
    total = sum(abs(v) for v in w.values())
    return {f: abs(v) / total if total > 0 else 0.0 for f, v in w.items()}


class CausalSampler:
    """Draws features outside a coalition in graph order (Ng et al., Eqs. 4-6): roots from their empirical marginal,
    continuous features from a linear regression on their feature parents, binary features from a logistic one."""

    def __init__(self, train: pd.DataFrame, ag: AnalysisGraph, features: Sequence[str], binary_c: float = 1e4):
        self.features = list(features)
        self.col = {f: j for j, f in enumerate(self.features)}
        self.order = list(nx.lexicographical_topological_sort(ag.graph.subgraph(self.features)))
        self.fits = {}
        for f in self.order:
            pa = [p for p in ag.parents(f) if p in self.features]
            y = train[f].to_numpy(dtype=float)
            if not pa or len(np.unique(y)) < 2:
                self.fits[f] = ("marginal", y)
            elif ag.scale[f] == "binary":
                m = LogisticRegression(C=binary_c, max_iter=5000).fit(train[pa].to_numpy(dtype=float), y)
                self.fits[f] = ("logistic", pa, float(m.intercept_[0]), m.coef_[0])
            else:
                X = np.column_stack([np.ones(len(train)), train[pa].to_numpy(dtype=float)])
                beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                self.fits[f] = ("linear", pa, beta, float(np.std(y - X @ beta)))

    def noise(self, m: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        out = {}
        for f in self.order:
            fit = self.fits[f]
            if fit[0] == "marginal":
                out[f] = fit[1][rng.integers(0, len(fit[1]), size=m)]
            elif fit[0] == "logistic":
                out[f] = rng.uniform(size=m)
            else:
                out[f] = rng.normal(0, fit[3], size=m)
        return out

    def draw(self, x: dict[str, float], coalitions: np.ndarray, noise: dict[str, np.ndarray]) -> np.ndarray:
        """coalitions: (c, k) boolean; the same m noise draws serve every coalition. Returns (c * m, k)."""
        c = len(coalitions)
        vals = {}
        for f in self.order:
            fit = self.fits[f]
            e = np.tile(noise[f], c)
            if fit[0] == "marginal":
                drawn = e
            elif fit[0] == "logistic":
                _, pa, b0, b = fit
                drawn = (e < sigmoid(b0 + sum(bj * vals[p] for bj, p in zip(b, pa)))).astype(float)
            else:
                _, pa, beta, _ = fit
                drawn = beta[0] + sum(bj * vals[p] for bj, p in zip(beta[1:], pa)) + e
            vals[f] = np.where(np.repeat(coalitions[:, self.col[f]], len(noise[f])), x[f], drawn)
        return np.column_stack([vals[f] for f in self.features])


def antithetic_orders(k: int, n: int, rng: np.random.Generator) -> list[np.ndarray]:
    """n orders as n // 2 random orders, each followed by its reverse."""
    out = []
    for _ in range(n // 2):
        p = rng.permutation(k)
        out += [p, p[::-1]]
    return out


def permutation_shapley(value: Callable[[np.ndarray], np.ndarray], k: int, orders: Sequence[np.ndarray]) -> np.ndarray:
    """Mean marginal contribution over the given orders; `value` maps (c, k) boolean coalitions to (c,) values."""
    phi = np.zeros(k)
    for order in orders:
        prefixes = np.zeros((k + 1, k), dtype=bool)
        for j, f in enumerate(order):
            prefixes[j + 1:, f] = True
        v = value(prefixes)
        phi[order] += np.diff(v)
    return phi / len(orders)


def ng_causal_shap(predict: Callable[[np.ndarray], np.ndarray], train: pd.DataFrame, explain: pd.DataFrame,
                   features: Sequence[str], ag: AnalysisGraph, samples: int, orders: int, seed: int,
                   expected_value: float) -> AttributionResult:
    features = list(features)
    k = len(features)
    weights = causal_weights(ag, edge_strengths(train, ag), features)
    gamma = np.array([weights[f] for f in features])
    sampler = CausalSampler(train, ag, features)
    rng = np.random.default_rng(seed)
    preds = predict(explain[features].to_numpy(dtype=float))
    out = np.zeros((len(explain), k))
    diag = np.full((len(explain), 3), np.nan)
    if gamma.sum() > 0:
        for r, (_, row) in enumerate(explain[features].iterrows()):
            x = row.to_dict()

            def value(coalitions: np.ndarray) -> np.ndarray:
                draws = sampler.draw(x, coalitions, sampler.noise(samples, rng))
                return predict(draws).reshape(len(coalitions), samples).mean(axis=1)

            phi = permutation_shapley(value, k, antithetic_orders(k, orders, rng))
            weighted = gamma * phi
            total = weighted.sum()
            scale = (preds[r] - expected_value) / total if total != 0 else 0.0
            out[r] = weighted * scale
            diag[r] = (phi.sum(), total, scale)
    eff = out.sum(axis=1) + expected_value - preds if gamma.sum() > 0 else np.zeros(len(explain))
    diagnostics = pd.DataFrame(diag, columns=["shapley_total", "weighted_total", "rescale_factor"], index=explain.index)
    return AttributionResult(pd.DataFrame(out, columns=features, index=explain.index), expected_value, eff, diagnostics)
