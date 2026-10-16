"""Figure 2: the true analysis graph on its own, at full double-column width, for layout in a vector editor.

Each node, label and edge is a named SVG group (node_<id>, label_<id>, edge_<parent>__<child>) with live text.
"""
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

import style

SHAPE = style.ROOT / "outputs" / "03_shape_validation"
FACE = {"ancestor": "#3B3B3B", "non_ancestor": "white", "outcome": "#0072B2"}


def main():
    style.apply()
    nodes = pd.read_csv("node_layout.csv")
    xy = {r.id: (r.x, r.y) for r in nodes.itertuples()}
    edges = pd.read_csv(SHAPE / "analysis_graph_edges.csv")

    fig, ax = plt.subplots(figsize=(style.DOUBLE, style.DOUBLE * 0.5))
    for e in edges.itertuples():
        ax.add_patch(FancyArrowPatch(posA=xy[e.parent], posB=xy[e.child], gid=f"edge_{e.parent}__{e.child}",
                                     arrowstyle="-|>,head_length=3,head_width=1.5", mutation_scale=1,
                                     color="#8C8C8C", lw=0.6, shrinkA=4, shrinkB=4, zorder=1))
    for n in nodes.itertuples():
        ax.scatter(n.x, n.y, s=14 if n.group != "outcome" else 24, c=FACE[n.group], edgecolors="#3B3B3B",
                   linewidths=0.5, zorder=3, gid=f"node_{n.id}")
        ax.annotate(n.label.replace("\\n", "\n"), (n.x, n.y), xytext=(0, -4), textcoords="offset points",
                    ha="center", va="top", fontsize=6, zorder=4, gid=f"label_{n.id}")
    ax.axis("off")
    ax.margins(0.04)
    handles = [Line2D([], [], marker="o", ls="", mfc=f, mec="#3B3B3B", mew=0.5, ms=4, label=l)
               for f, l in [("#3B3B3B", "Ancestor of outcome"), ("white", "Non-ancestor"), ("#0072B2", "Outcome")]]
    leg = fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.04))
    leg.set_gid("legend_nodes")
    style.save(fig, "fig2a_true_graph")
    print(f"{len(nodes)} nodes, {len(edges)} edges")


if __name__ == "__main__":
    main()
