"""Table 1 (ordering of all 35 features by truth and by each method) and the ancestor-displacement check.

Truth columns are |truth| for the per-unit and population-scaled truths, each for the data set's own population,
ranked with ties from the 2% margin tie groups; non-ancestors (truth exactly zero) share the bottom rank. Method
columns are the mean importance over replicates from the all-features run. Every column is also expressed as a share
of its total over the 35 features. A rank is a competition rank: tied features share the best rank in the tie.

Displacement: an ancestor is displaced when its rank among the 16 ancestors differs between the ancestor-set run and
the all-features run (non-ancestors removed before ranking); the per-replicate Kendall tau_b between the two
orderings shows the flag against replicate noise.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import kendalltau

from . import data as D

TRUTH_COLUMNS = ("per_unit", "pop")
TABLE_METHODS = ("standard_shap", "ordering_only_shap", "ng_style_shap", "interventional_shap")
TIE_RULE = "delta2"
SCOPE = "with_era"


def load_labels(path: Path) -> dict:
    return yaml.safe_load(path.read_text())


def competition_rank(values: pd.Series, groups: list[list[str]] | None = None) -> pd.Series:
    """Rank by descending value; features in the same group share the best rank and the next rank skips the group."""
    order = sorted(values.index, key=lambda f: (-values[f], f))
    grouped = {f for g in groups or [] for f in g}
    members = list(groups or []) + [[f] for f in order if f not in grouped]
    rank, seen = {}, 0
    for g in sorted(members, key=lambda g: min(order.index(f) for f in g)):
        for f in g:
            rank[f] = seen + 1
        seen += len(g)
    return pd.Series(rank).reindex(values.index)


def tie_groups_for(ties: pd.DataFrame, truth_type: str, population: str) -> list[list[str]]:
    t = ties[(ties["truth_type"] == truth_type) & (ties["population"] == population) & (ties["scope"] == SCOPE)
             & (ties["rule"] == TIE_RULE)].sort_values("group")
    return [g.split("; ") for g in t["features"]]


def truth_block(truth: pd.DataFrame, ties: pd.DataFrame, fsets: pd.DataFrame, dataset: str, truth_type: str) -> pd.DataFrame:
    """One row per feature: |truth|, rank (ties shared, non-ancestors at the bottom), tie group and share of the total."""
    pop = D.TRUTH_POPULATION[dataset]
    v = D.truth_series(truth, truth_type, pop).abs()
    feats = fsets["id"].tolist()
    v = v.reindex(feats)
    anc = fsets.set_index("id")["in_ancestor_set"]
    groups = tie_groups_for(ties, truth_type, pop)
    ancestors = [f for f in feats if anc[f]]
    rank = competition_rank(v[ancestors], groups)
    bottom = len(ancestors) + 1
    out = pd.DataFrame({"feature": feats, "value": v.to_numpy()})
    out["rank"] = [rank.get(f, bottom) for f in feats]
    gid = {f: i for i, g in enumerate(groups, 1) for f in g}
    out["tie_group"] = [gid.get(f, np.nan) for f in feats]
    out["share_pct"] = 100 * out["value"] / out["value"].sum()
    out["mcse_share_pct"] = np.nan
    out["n_replicates"] = np.nan
    return out


def method_block(importance: pd.DataFrame, fsets: pd.DataFrame, dataset: str, method: str) -> pd.DataFrame:
    """Mean importance over replicates from the all-features run, its rank and share, and the MCSE of the share."""
    feats = fsets["id"].tolist()
    d = importance[(importance["dataset"] == dataset) & (importance["feature_set"] == "all") & (importance["method"] == method)]
    wide = d.pivot_table(index="replicate", columns="feature", values="importance")[feats]
    mean = wide.mean()
    shares = wide.div(wide.sum(axis=1), axis=0) * 100
    n = len(wide)
    out = pd.DataFrame({"feature": feats, "value": mean.to_numpy()})
    out["rank"] = competition_rank(mean).to_numpy()
    out["tie_group"] = np.nan
    out["share_pct"] = 100 * out["value"] / out["value"].sum()
    out["mcse_share_pct"] = (shares.std(ddof=1) / np.sqrt(n)).to_numpy()
    out["n_replicates"] = n
    return out


def table1_long(importance: pd.DataFrame, truth: pd.DataFrame, ties: pd.DataFrame, fsets: pd.DataFrame, labels: dict,
                datasets: list[str]) -> pd.DataFrame:
    names = fsets.set_index("id")["name"]
    anc = fsets.set_index("id")["in_ancestor_set"]
    parts = []
    for ds in datasets:
        cols = [("truth", t, labels["truth_type"][t], truth_block(truth, ties, fsets, ds, t)) for t in TRUTH_COLUMNS]
        cols += [("method", m, labels["methods"][m], method_block(importance, fsets, ds, m)) for m in TABLE_METHODS]
        for kind, key, label, b in cols:
            b = b.assign(dataset=ds, dataset_label=labels["datasets"][ds], column_type=kind, column=key, column_label=label,
                         population=D.TRUTH_POPULATION[ds] if kind == "truth" else "all-features run")
            b["feature_name"] = b["feature"].map(names)
            b["ancestor"] = b["feature"].map(anc).astype(bool)
            parts.append(b)
    cols_out = ["dataset", "dataset_label", "feature", "feature_name", "ancestor", "column_type", "column", "column_label",
                "population", "value", "rank", "tie_group", "share_pct", "mcse_share_pct", "n_replicates"]
    return pd.concat(parts, ignore_index=True)[cols_out]


def displacement(importance: pd.DataFrame, fsets: pd.DataFrame, labels: dict, datasets: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per ancestor, data set and method: ranks among the ancestors from each run and whether they differ; and per
    data set and method: the replicate-wise tau_b between the two ancestor orderings (mean and MCSE)."""
    ancestors = fsets.loc[fsets["in_ancestor_set"], "id"].tolist()
    names = fsets.set_index("id")["name"]
    rows, summ = [], []
    for ds in datasets:
        for m in TABLE_METHODS:
            d = importance[(importance["dataset"] == ds) & (importance["method"] == m)]
            a = d[d["feature_set"] == "ancestor"].pivot_table(index="replicate", columns="feature", values="importance")[ancestors]
            b = d[d["feature_set"] == "all"].pivot_table(index="replicate", columns="feature", values="importance")[ancestors]
            ra, rb = competition_rank(a.mean()), competition_rank(b.mean())
            for f in ancestors:
                rows.append({"dataset": ds, "dataset_label": labels["datasets"][ds], "method": m,
                             "method_label": labels["methods"][m], "feature": f, "feature_name": names[f],
                             "mean_importance_ancestor_run": a[f].mean(), "mean_importance_all_run": b[f].mean(),
                             "rank_ancestor_run": int(ra[f]), "rank_all_run_among_ancestors": int(rb[f]),
                             "displaced": bool(ra[f] != rb[f])})
            taus = np.array([kendalltau(a.loc[r], b.loc[r]).statistic for r in a.index.intersection(b.index)])
            summ.append({"dataset": ds, "dataset_label": labels["datasets"][ds], "method": m, "method_label": labels["methods"][m],
                         "n_ancestors": len(ancestors), "n_displaced": int(sum(r["displaced"] for r in rows[-len(ancestors):])),
                         "n_replicates": len(taus), "tau_b_mean": float(np.nanmean(taus)),
                         "tau_b_mcse": float(np.nanstd(taus, ddof=1) / np.sqrt(np.isfinite(taus).sum()))})
    return pd.DataFrame(rows), pd.DataFrame(summ)


def cell_text(rank: float, share_pct: float, displaced: bool = False) -> str:
    return f"{int(rank)} ({share_pct:.1f}%)" + ("*" if displaced else "")


def wide_table(long: pd.DataFrame, disp: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """One wide table per data set: rows ordered by the first truth column, cells 'rank (share%)', displaced ancestors
    starred in the method columns, ancestors marked with a dagger after the feature name."""
    out = {}
    flag = disp.set_index(["dataset", "method", "feature"])["displaced"]
    for ds, d in long.groupby("dataset", sort=False):
        first = d[d["column_type"] == "truth"]["column"].iloc[0]
        order = d[d["column"] == first].sort_values(["rank", "share_pct", "feature"], ascending=[True, False, True])["feature"]
        table = pd.DataFrame({"feature": order.to_numpy()})
        info = d.drop_duplicates("feature").set_index("feature")
        table.insert(0, "Feature", [info.loc[f, "feature_name"] + (" †" if info.loc[f, "ancestor"] else "") for f in table["feature"]])
        for (kind, key, label), c in d.groupby(["column_type", "column", "column_label"], sort=False):
            c = c.set_index("feature")
            star = [kind == "method" and bool(flag.get((ds, key, f), False)) for f in table["feature"]]
            table[label] = [cell_text(c.loc[f, "rank"], c.loc[f, "share_pct"], s) for f, s in zip(table["feature"], star)]
        out[ds] = table.drop(columns="feature")
    return out


def markdown(wide: dict[str, pd.DataFrame], labels: dict, summary: pd.DataFrame) -> str:
    lines = ["Cells show rank (share of the column total, %). † ancestor of nephrolithiasis. * ancestor displaced: its rank among "
             "the 16 ancestors differs between the ancestor-set run and the all-features run. Truth ranks share a rank within "
             "a tie group (2% margin); non-ancestors share the bottom rank.", ""]
    for ds, t in wide.items():
        lines += [f"### {labels['datasets'][ds]}", "", "| " + " | ".join(t.columns) + " |", "|" + "---|" * len(t.columns)]
        lines += ["| " + " | ".join(str(x) for x in r) + " |" for r in t.itertuples(index=False)]
        s = summary[summary["dataset"] == ds]
        lines += ["", "Ancestor ordering, ancestor-set run vs all-features run (replicate-wise τ_b, mean ± MCSE; displaced ancestors): " +
                  "; ".join(f"{r.method_label} {r.tau_b_mean:.2f} ± {r.tau_b_mcse:.2f} ({r.n_displaced} of {r.n_ancestors})"
                            for r in s.itertuples()), ""]
    return "\n".join(lines)
