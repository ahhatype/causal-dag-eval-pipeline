from __future__ import annotations

from typing import Mapping

import numpy as np
from scipy.stats import kendalltau


def _common(a: Mapping[str, float], b: Mapping[str, float]) -> list[str]:
    common = sorted(set(a) & set(b))
    if not common:
        raise ValueError("no shared features")
    return common


def tie_adjusted_truth(effect: Mapping[str, float], se: Mapping[str, float], z: float = 1.96) -> dict[str, float]:
    """|true effect|, with features whose effects differ by less than z paired standard errors set equal.

    Chains of such pairs form one group, whose members take the group's mean |effect|.
    """
    names = sorted(effect)
    parent = {f: f for f in names}

    def root(f: str) -> str:
        while parent[f] != f:
            f = parent[f]
        return f

    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if abs(abs(effect[a]) - abs(effect[b])) < z * np.hypot(se[a], se[b]):
                parent[root(a)] = root(b)
    groups: dict[str, list[str]] = {}
    for f in names:
        groups.setdefault(root(f), []).append(f)
    return {f: float(np.mean([abs(effect[g]) for g in members])) for members in groups.values() for f in members}


def kendall_tau_b(importance: Mapping[str, float], truth: Mapping[str, float]) -> float:
    """Kendall's tau-b between importance and |true total effect|."""
    c = _common(importance, truth)
    return float(kendalltau([importance[f] for f in c], [abs(truth[f]) for f in c]).statistic)


def top_k_recovery(importance: Mapping[str, float], truth: Mapping[str, float], k: int) -> float:
    """Share of the k true top features that are also in the method's top-k; features tied with the k-th all count."""
    c = _common(importance, truth)
    k = min(k, len(c))
    cutoff = sorted((abs(truth[f]) for f in c), reverse=True)[k - 1]
    top_true = {f for f in c if abs(truth[f]) >= cutoff}
    top_method = set(sorted(c, key=lambda f: (-importance[f], f))[:k])
    return min(len(top_true & top_method), k) / k


def proximity_bias_index(importance: Mapping[str, float], truth: Mapping[str, float],
                         distance: Mapping[str, int]) -> float:
    """Importance-weighted mean distance to the outcome under the truth minus under the method.

    Positive values mean the method pools credit closer to the outcome than the truth does.
    """
    c = [f for f in _common(importance, truth) if f in distance]
    d = np.array([distance[f] for f in c], dtype=float)
    wt = np.array([abs(truth[f]) for f in c])
    wm = np.array([importance[f] for f in c])
    if wt.sum() == 0 or wm.sum() == 0:
        return float("nan")
    return float(d @ wt / wt.sum() - d @ wm / wm.sum())


def non_ancestor_share(importance: Mapping[str, float], ancestors: set[str]) -> float:
    """Share of total attribution assigned to features with a true effect of zero."""
    total = sum(importance.values())
    if total == 0:
        return float("nan")
    return float(sum(v for f, v in importance.items() if f not in ancestors) / total)


def auc(y: np.ndarray, score: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, score)) if 0 < np.mean(y) < 1 else float("nan")
