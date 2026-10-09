"""The data-generating structural equations, mirroring 01_data_generation/R/model_spec.R node by node.

Each node is a function of its parents and its own uniform noise U_<node>, so fixing a node re-evaluates only what
depends on it. `tests` compare this copy against records simulated by simcausal (`scm_check.csv`).
"""

from __future__ import annotations

from typing import Mapping

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import norm

from .data import OUTCOME, GenerationConfig


def _expit(v: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-v))


SPECIAL = {"mission_era": "era_binary", "cumulative_mission_duration": "duration", "ultrasound": "ultrasound",
           "detect_mrm_stone": "detect_mrm", "detect_long_term_health_outcomes": "detect_ltho",
           "loss_of_mission_objectives": "lmo", "loss_of_mission": "lm"}


class TrueSCM:
    def __init__(self, cfg: GenerationConfig):
        n = cfg.nodes
        self.p = cfg.params
        self.cal = cfg.calibration
        sim = set(n.loc[n["simulate"], "id"])
        e = cfg.edges[cfg.edges["parent"].isin(sim) & cfg.edges["child"].isin(sim)]
        self.edges = e
        g = nx.DiGraph()
        g.add_nodes_from(sorted(sim))
        g.add_edges_from(zip(e["parent"], e["child"]))
        self.order = list(nx.topological_sort(g))
        self.scale = dict(zip(n["id"], n["scale"]))
        self.lin = {c: [(p, float(k)) for p, k in zip(grp["parent"], grp["coef"]) if pd.notna(k)]
                    for c, grp in e.groupby("child")}
        self.noise_ids = [*self.order, "duration_mix", "selection"]
        self.rule = {i: self._rule(i) for i in self.order}

    def _rule(self, i: str) -> str:
        if i in SPECIAL:
            return SPECIAL[i]
        if i in self.p["era_values"]:
            return "era_fixed"
        if i in self.p["uptake"]:
            return "uptake_logit" if self.lin.get(i) else "era_prob"
        if not self.edges["child"].eq(i).any():
            return "root_normal"
        if self.scale[i] == "continuous":
            return "cont"
        if i in self.p["base_rates"]:
            return "logit"
        raise ValueError(f"no generation rule for node {i}")

    def noise(self, n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        return {f"U_{i}": rng.uniform(size=n) for i in self.noise_ids}

    def _linear(self, i: str, v: Mapping[str, np.ndarray]) -> np.ndarray:
        return sum(k * v[p] for p, k in self.lin[i])

    def _trunc_norm(self, mean, sd, lo, hi, u):
        pa, pb = norm.cdf((lo - mean) / sd), norm.cdf((hi - mean) / sd)
        return mean + sd * norm.ppf(pa + (pb - pa) * u)

    def _node(self, i: str, v: Mapping[str, np.ndarray], u: Mapping[str, np.ndarray], n: int) -> np.ndarray:
        p, cal, U = self.p, self.cal, u[f"U_{i}"]
        r = self.rule[i]
        if r == "era_binary":
            return (U < p["era"]["share_2000s"]).astype(float)
        if r == "root_normal":
            return norm.ppf(U)
        if r == "era_fixed":
            a, b = p["era_values"][i]
            if a == b:
                return np.zeros(n)
            sh = p["era"]["share_2000s"]
            mu, sd = (1 - sh) * a + sh * b, np.sqrt(sh * (1 - sh)) * abs(b - a)
            return ((a + (b - a) * v["mission_era"]) - mu) / sd
        if r == "duration":
            d = p["duration"]
            e60, s = d["era_1960s"], d["era_2000s"]
            raw = np.where(v["mission_era"] == 0,
                           np.exp(self._trunc_norm(np.log(e60["median"]), e60["log_sd"], np.log(e60["min"]),
                                                   np.log(e60["max"]), U)),
                           np.where(u["U_duration_mix"] < s["shuttle_share"],
                                    self._trunc_norm(s["shuttle"]["mean"], s["shuttle"]["sd"], s["shuttle"]["min"],
                                                     s["shuttle"]["max"], U),
                                    self._trunc_norm(s["iss"]["mean"], s["iss"]["sd"], s["iss"]["min"], s["iss"]["max"], U)))
            k = cal["cont"][i]
            return (raw - k["mu"]) / k["sd"]
        if r == "cont":
            k = cal["cont"][i]
            return ((self._linear(i, v) + k["rs"] * norm.ppf(U)) - k["mu"]) / k["sd"]
        if r == "logit":
            return (U < _expit(cal["logit"][i] + self._linear(i, v))).astype(float)
        if r == "era_prob":
            a, b = p["uptake"][i]
            return (U < a + (b - a) * (v["medical_prevention_capability"] > 0)).astype(float)
        if r == "uptake_logit":
            a60, a00 = (float(x) for x in cal["uptake"][i])
            return (U < _expit(np.where(v["medical_prevention_capability"] > 0, a00, a60) + self._linear(i, v))).astype(float)
        if r == "ultrasound":
            f = p["flip"]
            return (U < np.where(v["medical_monitoring_capability"] > 0, 1 - f, f)).astype(float)
        if r == "detect_mrm":
            d = p["detection"]["mrm_stone"]
            us = v["ultrasound"] == 1
            sens = _expit(np.log(np.where(us, d["with_ultrasound"]["sensitivity"], d["without_ultrasound"]["sensitivity"])
                                 / (1 - np.where(us, d["with_ultrasound"]["sensitivity"],
                                                 d["without_ultrasound"]["sensitivity"])))
                          + d["mrm_slope"] * v["mineralized_renal_material"])
            fpr = np.where(us, 1 - d["with_ultrasound"]["specificity"], 1 - d["without_ultrasound"]["specificity"])
            present = np.maximum(v["nephrolithiasis"], (v["mineralized_renal_material"] > d["mrm_threshold"]).astype(float))
            return (U < np.where(present == 1, sens, fpr)).astype(float)
        if r == "detect_ltho":
            d = p["detection"]["ltho"]
            s_on = v["surveillance"] > 0
            sens = np.where(s_on, d["with_surveillance"]["sensitivity"], d["without_surveillance"]["sensitivity"])
            fpr = np.where(s_on, 1 - d["with_surveillance"]["specificity"], 1 - d["without_surveillance"]["specificity"])
            return (U < np.where(v["long_term_health_outcomes"] > d["threshold"], sens, fpr)).astype(float)
        if r == "lmo":
            k = float(self.edges.loc[(self.edges["parent"] == "task_performance") & (self.edges["child"] == i), "coef"].iloc[0])
            return (U < np.where(v["evacuation"] == 1, 1 - p["flip"], _expit(cal["logit"][i] + k * v["task_performance"]))
                    ).astype(float)
        if r == "lm":
            return (U < np.where(v["loss_of_mission_objectives"] == 1, 1 - p["flip"], p["flip"])).astype(float)
        raise ValueError(r)

    def simulate(self, u: Mapping[str, np.ndarray], interventions: Mapping[str, np.ndarray] | None = None,
                 masks: Mapping[str, np.ndarray] | None = None) -> dict[str, np.ndarray]:
        """All nodes. An intervened node takes the given values; with a mask, only where the mask is True."""
        interventions = interventions or {}
        masks = masks or {}
        n = len(next(iter(u.values())))
        v: dict[str, np.ndarray] = {}
        for i in self.order:
            if i in interventions and i not in masks:
                v[i] = np.broadcast_to(np.asarray(interventions[i], dtype=float), (n,)).copy()
                continue
            x = np.broadcast_to(self._node(i, v, u, n), (n,)).astype(float)
            v[i] = np.where(masks[i], interventions[i], x) if i in masks else x
        return v

    def recover_exogenous(self, data: pd.DataFrame, rng: np.random.Generator | None = None) -> dict[str, np.ndarray]:
        """The true noise, carried as U_<node> columns by records simulated from this model."""
        return {f"U_{i}": data[f"U_{i}"].to_numpy(dtype=float) for i in self.noise_ids}

    def population(self, n: int, rng: np.random.Generator) -> pd.DataFrame:
        u = self.noise(n, rng)
        return pd.DataFrame({**self.simulate(u), **u})

    def outcome_probability(self, v: Mapping[str, np.ndarray]) -> np.ndarray:
        return _expit(self.cal["logit"][OUTCOME] + self._linear(OUTCOME, v))


class TrueSampler:
    """Ng-style sampler using the true equations: features outside a coalition are generated in causal order with the
    coalition's features fixed. Same interface as ng.CausalSampler."""

    def __init__(self, scm: TrueSCM, features: list[str]):
        self.scm = scm
        self.features = list(features)

    def noise(self, m: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        return self.scm.noise(m, rng)

    def draw(self, x: dict[str, float], coalitions: np.ndarray, noise: dict[str, np.ndarray]) -> np.ndarray:
        c, m = len(coalitions), len(next(iter(noise.values())))
        u = {k: np.tile(a, c) for k, a in noise.items()}
        fixed = {f: np.full(c * m, x[f]) for f in self.features}
        masks = {f: np.repeat(coalitions[:, j], m) for j, f in enumerate(self.features)}
        v = self.scm.simulate(u, fixed, masks)
        return np.column_stack([v[f] for f in self.features])
