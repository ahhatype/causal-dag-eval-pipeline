from __future__ import annotations

from typing import Mapping

import numpy as np
import pandas as pd
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


ROLE_COUNTS = {"ancestor": 16, "descendant": 15, "unrelated": 4}


def roles_from_feature_sets(fsets: pd.DataFrame) -> dict[str, str]:
    """Role of each all-features feature: ancestor (in the ancestor set), descendant or unrelated (neither)."""
    f = fsets[fsets["in_all_features_set"]]
    masks = {"ancestor": f["in_ancestor_set"], "descendant": ~f["in_ancestor_set"] & (f["role"] == "Descendant"),
             "unrelated": ~f["in_ancestor_set"] & (f["role"] == "Neither")}
    assigned = sum(m.astype(int) for m in masks.values())
    if (assigned != 1).any():
        raise ValueError(f"features without exactly one role: {sorted(f.loc[assigned != 1, 'id'])}")
    roles = {i: r for r, m in masks.items() for i in f.loc[m, "id"]}
    counts = pd.Series(roles).value_counts().to_dict()
    if counts != ROLE_COUNTS:
        raise ValueError(f"role counts {counts} differ from {ROLE_COUNTS}")
    return roles


def credit_share_by_role(importance: Mapping[str, float], roles: Mapping[str, str]) -> dict[str, float]:
    """Each role's share of total importance in one all-features run; the shares sum to 1."""
    if set(importance) != set(roles):
        raise ValueError(f"features and roles differ: {sorted(set(importance) ^ set(roles))}")
    total = sum(importance.values())
    out = {r: 0.0 for r in sorted(set(roles.values()))}
    for f, v in importance.items():
        out[roles[f]] += v / total if total else float("nan")
    return out


def credit_by_role(importance: pd.DataFrame, roles: Mapping[str, str], methods: list[str]) -> pd.DataFrame:
    """Credit share by role per data set and method: shares within each replicate, then mean and MC SE over replicates."""
    rows = []
    runs = importance[(importance["feature_set"] == "all") & importance["method"].isin(methods)]
    for (ds, m, rep), g in runs.groupby(["dataset", "method", "replicate"]):
        for r, v in credit_share_by_role(dict(zip(g["feature"], g["importance"])), roles).items():
            rows.append({"dataset": ds, "method": m, "replicate": rep, "role": r, "share": v})
    g = pd.DataFrame(rows).groupby(["dataset", "method", "role"])["share"]
    out = pd.concat([g.mean().rename("share_mean"), (g.std() / np.sqrt(g.count())).rename("share_mcse"),
                     g.count().rename("n_replicates")], axis=1).reset_index()
    return out


def check_non_ancestor_share(credit: pd.DataFrame, summary: pd.DataFrame, own: Mapping[str, str],
                             tol: float = 1e-3) -> pd.DataFrame:
    """Descendant plus unrelated share against the summary's non-ancestor share (all features, with era, own population)."""
    non = credit[credit["role"] != "ancestor"].groupby(["dataset", "method"])["share_mean"].sum().rename("by_role")
    s = summary[(summary["feature_set"] == "all") & (summary["scope"] == "with_era")
                & (summary["population"] == summary["dataset"].map(own)) & (summary["truth_type"] == "per_unit")]
    m = s.set_index(["dataset", "method"])["non_ancestor_share_mean"].rename("summary").to_frame().join(non, how="right")
    m["difference"] = m["by_role"] - m["summary"]
    if not (m["difference"].abs() <= tol).all():
        raise ValueError(f"descendant + unrelated share differs from the summary by more than {tol}:\n{m[m['difference'].abs() > tol]}")
    return m


def auc(y: np.ndarray, score: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, score)) if 0 < np.mean(y) < 1 else float("nan")
