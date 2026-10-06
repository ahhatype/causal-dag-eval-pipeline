"""Causal SHAP of Ng et al. (arXiv:2509.00846), Algorithm 1, with the analysis graph supplied in place of PC + IDA.

Adapted from causal-shap-spaceflight-renal-stones (python/src/causal_shap_renal/attribution_structural.py).
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

import networkx as nx
import numpy as np
import pandas as pd

from .data import OUTCOME, AnalysisGraph
from .interventional import AttributionResult


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
    """Draws features outside a coalition in graph order: roots from their empirical marginal,
    others from a linear regression on their parents among the features (Ng et al., Eqs. 4-6)."""

    def __init__(self, train: pd.DataFrame, ag: AnalysisGraph, features: Sequence[str]):
        self.features = list(features)
        sub = ag.graph.subgraph(self.features)
        self.order = list(nx.lexicographical_topological_sort(sub))
        self.fits = {}
        for f in self.order:
            pa = [p for p in ag.parents(f) if p in self.features]
            y = train[f].to_numpy(dtype=float)
            if not pa:
                self.fits[f] = ("marginal", y)
            else:
                X = np.column_stack([np.ones(len(train)), train[pa].to_numpy(dtype=float)])
                beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                self.fits[f] = ("linear", pa, beta, float(np.std(y - X @ beta)))

    def draw(self, x: dict[str, float], coalitions: np.ndarray, m: int, rng: np.random.Generator) -> np.ndarray:
        """coalitions: (c, k) boolean; returns (c * m, k) feature matrix in self.features column order."""
        c = len(coalitions)
        col = {f: j for j, f in enumerate(self.features)}
        vals = {}
        for f in self.order:
            fixed = np.repeat(coalitions[:, col[f]], m)
            fit = self.fits[f]
            if fit[0] == "marginal":
                drawn = rng.choice(fit[1], size=c * m, replace=True)
            else:
                _, pa, beta, sd = fit
                drawn = beta[0] + sum(b * vals[p] for b, p in zip(beta[1:], pa)) + rng.normal(0, sd, c * m)
            vals[f] = np.where(fixed, x[f], drawn)
        return np.column_stack([vals[f] for f in self.features])


def ng_causal_shap(predict: Callable[[np.ndarray], np.ndarray], train: pd.DataFrame, explain: pd.DataFrame,
                   features: Sequence[str], ag: AnalysisGraph, samples: int, iterations: int, seed: int,
                   expected_value: float) -> AttributionResult:
    features = list(features)
    k = len(features)
    weights = causal_weights(ag, edge_strengths(train, ag), features)
    gamma = np.array([weights[f] for f in features])
    sampler = CausalSampler(train, ag, features)
    rng = np.random.default_rng(seed)
    out = np.zeros((len(explain), k))
    preds = predict(explain[features].to_numpy(dtype=float))
    if gamma.sum() > 0:
        for r, (_, row) in enumerate(explain[features].iterrows()):
            x = row.to_dict()
            coalitions, pairs = [], []
            for _ in range(iterations):
                s = rng.integers(0, k)
                S = np.zeros(k, dtype=bool)
                S[rng.choice(k, size=s, replace=False)] = True
                base = len(coalitions)
                coalitions.append(S)
                for i in np.flatnonzero(~S):
                    Si = S.copy()
                    Si[i] = True
                    pairs.append((base, len(coalitions), i, math.factorial(s) * math.factorial(k - s - 1) / math.factorial(k)))
                    coalitions.append(Si)
            v = predict(sampler.draw(x, np.array(coalitions), samples, rng)).reshape(len(coalitions), samples).mean(axis=1)
            phi = np.zeros(k)
            for b, j, i, wk in pairs:
                phi[i] += wk * gamma[i] * (v[j] - v[b])
            total = phi.sum()
            out[r] = phi * (preds[r] - expected_value) / total if total != 0 else 0.0
    eff = out.sum(axis=1) + expected_value - preds
    if gamma.sum() == 0:
        eff = np.zeros(len(explain))
    return AttributionResult(pd.DataFrame(out, columns=features, index=explain.index), expected_value, eff)
