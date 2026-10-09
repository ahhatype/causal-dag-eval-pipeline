from __future__ import annotations

from typing import Mapping

import numpy as np
from scipy.stats import kendalltau, weightedtau


def _common(a: Mapping[str, float], b: Mapping[str, float]) -> list[str]:
    common = sorted(set(a) & set(b))
    if not common:
        raise ValueError("no shared features")
    return common


def margin_tie_groups(values: Mapping[str, float], margin: float) -> list[list[str]]:
    """Features ordered by |value|; neighbours closer than `margin` chain into one tie group."""
    order = sorted(values, key=lambda f: (-abs(values[f]), f))
    groups = [[order[0]]] if order else []
    for a, b in zip(order[:-1], order[1:]):
        if abs(values[a]) - abs(values[b]) < margin:
            groups[-1].append(b)
        else:
            groups.append([b])
    return groups


def mcse_tie_groups(values: Mapping[str, float], mcse: Mapping[str, float], z: float = 1.96) -> list[list[str]]:
    """Any two features whose |values| differ by less than z combined MCSE are tied; ties chain into groups."""
    names = sorted(values)
    parent = {f: f for f in names}

    def root(f: str) -> str:
        while parent[f] != f:
            f = parent[f]
        return f

    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if abs(abs(values[a]) - abs(values[b])) < z * np.hypot(mcse[a], mcse[b]):
                parent[root(a)] = root(b)
    groups: dict[str, list[str]] = {}
    for f in names:
        groups.setdefault(root(f), []).append(f)
    return sorted(groups.values(), key=lambda g: -max(abs(values[f]) for f in g))


def tie_adjusted(values: Mapping[str, float], groups: list[list[str]]) -> dict[str, float]:
    """|value|, with every member of a tie group set to the group's mean |value|."""
    return {f: float(np.mean([abs(values[g]) for g in members])) for members in groups for f in members}


def tie_adjusted_truth(effect: Mapping[str, float], se: Mapping[str, float], z: float = 1.96) -> dict[str, float]:
    return tie_adjusted(effect, mcse_tie_groups(effect, se, z))


def kendall_tau_b(importance: Mapping[str, float], truth: Mapping[str, float]) -> float:
    """Kendall's tau-b between importance and |truth|."""
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


def weighted_tau(importance: Mapping[str, float], truth: Mapping[str, float]) -> float:
    """Vigna's (2015) weighted tau with hyperbolic weights: disagreements among the top-ranked count more."""
    c = _common(importance, truth)
    return float(weightedtau([importance[f] for f in c], [abs(truth[f]) for f in c]).statistic)


def l1_distance(importance: Mapping[str, float], truth: Mapping[str, float]) -> float:
    """Sum of absolute differences between importance and |truth|, each scaled to sum to 1 (range 0 to 2)."""
    c = _common(importance, truth)
    a = np.array([importance[f] for f in c], dtype=float)
    t = np.array([abs(truth[f]) for f in c], dtype=float)
    if a.sum() == 0 or t.sum() == 0:
        return float("nan")
    return float(np.abs(a / a.sum() - t / t.sum()).sum())


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
