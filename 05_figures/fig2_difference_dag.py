"""Figure 2: true analysis graph (a) and edge-level recovery by MGM-PC-Stable in each data set (b-d).

Reads node_layout.csv and outputs/03_shape_validation/difference_dag_edges.csv when it exists. Until the
rerun adds that file, it falls back to edge_recovery.csv, where reversed and unoriented edges cannot be
separated and are drawn together as "orientation unresolved".
"""
import argparse

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

import style

SHAPE = style.ROOT / "outputs" / "03_shape_validation"
VARIANT = "MGM-PC-Stable"



def load_recovery(added_threshold):
    new = SHAPE / "difference_dag_edges.csv"
    if new.exists():
        d = pd.read_csv(new)
        rates = {"correct": "correct_rate", "reversed": "reversed_rate", "unoriented": "unoriented_rate",
                 "overoriented": "overoriented_rate", "absent": "absent_rate"}
        recovery = pd.read_csv(SHAPE / "edge_recovery.csv").query("variant == @VARIANT")
        spurious = recovery[~recovery["true_edge"]].assign(added_rate=lambda x: x["adjacency_rate"])
        d = pd.concat([d[d["true_edge"]], spurious[["dataset", "edge", "true_edge", "adjacency_rate", "added_rate"]]])
        source = "difference_dag_edges.csv"
    else:
        d = pd.read_csv(SHAPE / "edge_recovery.csv").query("variant == @VARIANT").copy()
        d["correct_rate"] = d["cpdag_mark_rate"]
        d["orientation_rate"] = d["adjacency_rate"] - d["cpdag_mark_rate"]
        d["absent_rate"] = 1 - d["adjacency_rate"]
        rates = {"correct": "correct_rate", "orientation": "orientation_rate", "absent": "absent_rate"}
        source = "edge_recovery.csv (fallback: reversed not separated)"
    d["pair"] = d["edge"].str.split(" ")
    true = d[d["true_edge"]].copy()
    true["category"] = true[list(rates.values())].idxmax(axis=1).map({v: k for k, v in rates.items()})
    true["category_rate"] = true[list(rates.values())].max(axis=1)
    added = d[~d["true_edge"] & (d["adjacency_rate"] >= added_threshold)].copy()
    true[list(rates.values())] = true[list(rates.values())].fillna(0)
    added["category"], added["category_rate"] = "added", added["adjacency_rate"]
    return true, added, source


def draw_edge(ax, p, q, color, ls, alpha, arrow=True, lw=0.7, gid=None):
    ax.add_patch(FancyArrowPatch(gid=gid, posA=p, posB=q, arrowstyle="-|>,head_length=2.2,head_width=1.1" if arrow else "-",
                                 color=color, lw=lw, ls=ls, alpha=alpha, mutation_scale=1,
                                 shrinkA=6.2, shrinkB=6.2, zorder=1))


def draw_nodes(ax, nodes, labels):
    face = {"ancestor": "#3B3B3B", "non_ancestor": "white", "outcome": "#0072B2"}
    for _, n in nodes.iterrows():
        ax.scatter(n.x, n.y, s=9 if n.group != "outcome" else 16, c=face[n.group], edgecolors="#3B3B3B",
                   linewidths=0.5, zorder=3, gid=f"node_{n.group}_{n.id}")
        if labels:
            ax.annotate(n.label.replace("\\n", "\n"), (n.x, n.y), xytext=(0, -3.2), textcoords="offset points", ha="center",
                        va="top", fontsize=4.2, zorder=4, wrap=True, gid=f"label_{n.id}")
    ax.axis("off")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--added-threshold", type=float, default=0.10,
                    help="draw added edges found in at least this share of replicates")
    args = ap.parse_args()
    style.apply()

    nodes = pd.read_csv("node_layout.csv")
    xy = {r.id: (r.x, r.y) for r in nodes.itertuples()}
    graph = pd.read_csv(SHAPE / "analysis_graph_edges.csv")
    true, added, source = load_recovery(args.added_threshold)

    fig = plt.figure(figsize=(style.DOUBLE, style.DOUBLE * 0.72))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.12, wspace=0.03)
    ax = fig.add_subplot(grid[0, :])
    small = [fig.add_subplot(grid[1, i]) for i in range(3)]
    for e in graph.itertuples():
        draw_edge(ax, xy[e.parent], xy[e.child], "#8C8C8C", "-", 1, gid=f"true_edge_{e.parent}__{e.child}")
    draw_nodes(ax, nodes, labels=True)
    ax.set_title("a  True analysis graph", loc="left", fontweight="bold")
    ax.set_gid("panel_a_true_graph")

    for ax, ds, letter in zip(small, style.DATASETS, "bcd"):
        t = true[true["dataset"] == ds].drop(columns=["parent", "child"], errors="ignore").merge(graph, how="left", left_on="edge",
                                              right_on=graph["parent"].combine(graph["child"],
                                                                               lambda a, b: " ".join(sorted([a, b]))))
        for e in t.itertuples():
            color, ls, _ = style.EDGE_STYLE[e.category]
            draw_edge(ax, xy[e.parent], xy[e.child], color, ls, 0.35 + 0.65 * e.category_rate,
                      gid=f"{ds}_{e.category}_{e.parent}__{e.child}")
        for e in added[added["dataset"] == ds].itertuples():
            a, b = e.pair
            color, ls, _ = style.EDGE_STYLE["added"]
            draw_edge(ax, xy[a], xy[b], color, ls, 0.35 + 0.65 * e.category_rate, arrow=False, lw=0.9,
                      gid=f"{ds}_added_{a}__{b}")
            if e.category_rate >= 0.3:
                mx, my = (xy[a][0] + xy[b][0]) / 2, (xy[a][1] + xy[b][1]) / 2
                ax.annotate(f"{e.category_rate:.2f}", (mx, my), fontsize=4.5, color=color, ha="center",
                            va="bottom", xytext=(0, 1.5), textcoords="offset points")
        draw_nodes(ax, nodes, labels=False)
        ax.set_gid(f"panel_{letter}_{ds}")
        ax.set_title(f"{letter}  {style.DATASET_LABELS[ds].replace(chr(10), ' ')}", loc="left",
                     fontweight="bold")

    cats = sorted(set(true["category"]) | ({"added"} if len(added) else set()), key=list(style.EDGE_STYLE).index)
    handles = [Line2D([], [], color=style.EDGE_STYLE[c][0], ls=style.EDGE_STYLE[c][1], lw=1,
                      label=style.EDGE_STYLE[c][2]) for c in cats]
    handles += [Line2D([], [], marker="o", ls="", mfc=f, mec="#3B3B3B", mew=0.5, ms=3, label=l)
                for f, l in [("#3B3B3B", "Ancestor of outcome"), ("white", "Non-ancestor"), ("#0072B2", "Outcome")]]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))
    style.save(fig, "fig2_difference_dag")

    summary = true.groupby(["dataset", "category"]).size().unstack(fill_value=0)
    summary["added_drawn"] = added.groupby("dataset").size()
    summary.to_csv(style.OUT / "fig2_category_counts.csv")
    print(f"source: {source}\n{summary}")


if __name__ == "__main__":
    main()
