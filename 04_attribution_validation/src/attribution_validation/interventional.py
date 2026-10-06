"""Asymmetric interventional Shapley values on a structural causal model fitted to the data.

Adapted from causal-shap-target-dags (app/causal_shap/structural_value.py).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, LogisticRegression

from .data import AnalysisGraph


def sigmoid(v: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(v, -35.0, 35.0)))


@dataclass(frozen=True)
class NodeSpec:
    name: str
    kind: str
    parents: tuple[str, ...] = ()
    coefficients: tuple[float, ...] = ()
    intercept: float = 0.0
    noise_sd: float = 1.0
    root_probability: float = 0.5


class LinearLogisticSCM:
    def __init__(self, specs: Sequence[NodeSpec]):
        self.specs = {s.name: s for s in specs}
        g = nx.DiGraph()
        g.add_nodes_from(self.specs)
        g.add_edges_from((p, s.name) for s in specs for p in s.parents)
        if not nx.is_directed_acyclic_graph(g):
            raise ValueError("structural model must be acyclic")
        self.order = list(nx.lexicographical_topological_sort(g))

    def _linear(self, spec: NodeSpec, values: Mapping[str, np.ndarray], n: int) -> np.ndarray:
        lin = np.full(n, spec.intercept)
        for p, c in zip(spec.parents, spec.coefficients):
            lin = lin + c * values[p]
        return lin

    def recover_exogenous(self, data: pd.DataFrame, rng: np.random.Generator) -> dict[str, np.ndarray]:
        """Abduction. Binary noise is drawn once, uniformly over the range consistent with the observed value."""
        n = len(data)
        cols = {k: data[k].to_numpy(dtype=float) for k in self.order}
        u = {}
        for node in self.order:
            s = self.specs[node]
            lin = self._linear(s, cols, n)
            if s.kind == "continuous":
                u[node] = (cols[node] - lin) / s.noise_sd
            else:
                p = np.clip(sigmoid(lin) if s.parents else np.full(n, s.root_probability), 1e-9, 1 - 1e-9)
                obs = cols[node]
                u[node] = rng.uniform(np.where(obs == 1, 0.0, p), np.where(obs == 1, p, 1.0))
        return u

    def simulate(self, u: Mapping[str, np.ndarray], interventions: Mapping[str, np.ndarray] | None = None) -> dict:
        interventions = interventions or {}
        n = len(next(iter(u.values())))
        values: dict[str, np.ndarray] = {}
        for node in self.order:
            s = self.specs[node]
            if node in interventions:
                values[node] = np.broadcast_to(np.asarray(interventions[node], dtype=float), (n,))
                continue
            lin = self._linear(s, values, n)
            if s.kind == "continuous":
                values[node] = lin + s.noise_sd * u[node]
            else:
                p = sigmoid(lin) if s.parents else np.full(n, s.root_probability)
                values[node] = (u[node] < p).astype(float)
        return values


def fit_scm(data: pd.DataFrame, ag: AnalysisGraph, binary_c: float = 1e4) -> LinearLogisticSCM:
    """Each node regressed on its analysis-graph parents: linear for continuous, logistic for binary.

    The large-but-finite inverse penalty keeps logistic fits finite under separation (no 1960s drug use).
    """
    specs = []
    for node in ag.nodes:
        pa = tuple(ag.parents(node))
        y = data[node].to_numpy(dtype=float)
        if ag.scale[node] == "continuous":
            if pa:
                m = LinearRegression().fit(data[list(pa)], y)
                resid = y - m.predict(data[list(pa)])
                specs.append(NodeSpec(node, "continuous", pa, tuple(m.coef_), float(m.intercept_), float(resid.std())))
            else:
                specs.append(NodeSpec(node, "continuous", (), (), float(y.mean()), float(y.std())))
        else:
            if pa and 0 < y.mean() < 1:
                m = LogisticRegression(C=binary_c, max_iter=5000).fit(data[list(pa)], y)
                specs.append(NodeSpec(node, "binary", pa, tuple(m.coef_[0]), float(m.intercept_[0])))
            else:
                specs.append(NodeSpec(node, "binary", (), (), 0.0, 1.0, float(np.clip(y.mean(), 1e-6, 1 - 1e-6))))
    return LinearLogisticSCM(specs)


def ancestral_pairs(ag: AnalysisGraph, features: Sequence[str]) -> set[tuple[str, str]]:
    """(a, b) when feature a is an ancestor of feature b in the analysis graph, through any nodes."""
    g = ag.graph
    fs = set(features)
    return {(a, b) for a in features for b in nx.descendants(g, a) if b in fs}


class OrderSampler:
    """Uniform sampling of orders consistent with a partial order, by adjacent transpositions (Karzanov & Khachiyan 1991)."""

    def __init__(self, features: Sequence[str], before: set[tuple[str, str]], rng: np.random.Generator,
                 burn_in: int, thin: int):
        g = nx.DiGraph()
        g.add_nodes_from(features)
        g.add_edges_from(before)
        self.order = list(nx.lexicographical_topological_sort(g))
        self.before = before
        self.rng = rng
        self.thin = thin
        self._step(burn_in)

    def _step(self, steps: int) -> None:
        k = len(self.order)
        if k < 2:
            return
        for i in self.rng.integers(0, k - 1, size=steps):
            a, b = self.order[i], self.order[i + 1]
            if (a, b) not in self.before:
                self.order[i], self.order[i + 1] = b, a

    def sample(self) -> list[str]:
        self._step(self.thin)
        return list(self.order)


@dataclass(frozen=True)
class AttributionResult:
    values: pd.DataFrame
    baseline: float
    efficiency_error: np.ndarray


def interventional_shap(predict: Callable[[np.ndarray], np.ndarray], scm: LinearLogisticSCM, explain: pd.DataFrame,
                        background: pd.DataFrame, features: Sequence[str], ag: AnalysisGraph, orders: int,
                        burn_in: int, thin: int, seed: int) -> AttributionResult:
    """Value of a coalition S: mean prediction with S set to the explained record's values by intervention."""
    features = list(features)
    rng = np.random.default_rng(seed)
    u = scm.recover_exogenous(background, rng)
    n_eval, n_bg = len(explain), len(background)
    u_tiled = {k: np.tile(v, n_eval) for k, v in u.items()}
    base_vals = scm.simulate(u)
    baseline = float(np.mean(predict(np.column_stack([base_vals[f] for f in features]))))
    xv = {f: np.repeat(explain[f].to_numpy(dtype=float), n_bg) for f in features}
    sampler = OrderSampler(features, ancestral_pairs(ag, features), rng, burn_in, thin)
    out = np.zeros((n_eval, len(features)))
    idx = {f: i for i, f in enumerate(features)}
    for _ in range(orders):
        order = sampler.sample()
        fixed: dict[str, np.ndarray] = {}
        prev = np.full(n_eval, baseline)
        for f in order:
            fixed[f] = xv[f]
            sim = scm.simulate(u_tiled, fixed)
            cur = predict(np.column_stack([sim[g] for g in features])).reshape(n_eval, n_bg).mean(axis=1)
            out[:, idx[f]] += cur - prev
            prev = cur
    out /= orders
    eff = out.sum(axis=1) + baseline - predict(explain[features].to_numpy(dtype=float))
    return AttributionResult(pd.DataFrame(out, columns=features, index=explain.index), baseline, eff)
