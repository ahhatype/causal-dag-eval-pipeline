"""Figure 3: recovery of each true edge by MGM-PC-Stable in each data set, as stacked shares over replicates.

Rows are grouped by how the edge fared across data sets; each cell splits the 50 replicates into outcome categories.
"""
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch, Rectangle

import style

STRUCTURE = style.ROOT / "outputs" / "03_structure_validation"
CATS = ["correct", "reversed", "overoriented", "unoriented", "absent"]
RATE = {c: f"{c}_rate" for c in CATS}
GROUPS = ["Recovered correctly in all data sets", "Lost at n = 900", "Missed in all data sets",
          "Orientation errors", "Added by selection"]


def group(row):
    m = row["modal"]
    if not row["true_edge"]:
        return GROUPS[4]
    if all(m.get(d) == "correct" for d in style.DATASETS):
        return GROUPS[0]
    if m.get("full_set") == "absent":
        return GROUPS[2]
    if any(m.get(d) == "absent" for d in style.DATASETS[1:]):
        return GROUPS[1]
    return GROUPS[3]


def main():
    style.apply()
    labels = pd.read_csv("node_layout.csv").set_index("id")["label"].str.replace("\\n", " ", regex=False)
    d = pd.read_csv(STRUCTURE / "difference_dag_edges.csv")
    cols = list(RATE.values()) + ["added_rate"]
    d[cols] = d[cols].fillna(0)
    rows = []
    for edge, g in d.groupby("edge"):
        g = g.set_index("dataset")
        r = g.iloc[0]
        rows.append({"edge": edge, "parent": r.parent, "child": r.child, "true_edge": bool(r.true_edge),
                     "mark": r.true_mark, "modal": g["modal_category"].to_dict(), "data": g,
                     "full_correct": g["correct_rate"].get("full_set", 0)})
    rows = pd.DataFrame(rows)
    rows["group"] = rows.apply(group, axis=1)
    rows["order"] = rows["group"].map(GROUPS.index)
    rows = rows.sort_values(["order", "full_correct"], ascending=[True, False]).reset_index(drop=True)

    n_lines = len(rows) + rows["group"].nunique()
    fig = plt.figure(figsize=(style.DOUBLE, 0.105 * n_lines + 0.7))
    ax = fig.add_axes([0.36, 0.07, 0.62, 0.88])
    y, ticks, ticklabels = 0, [], []
    width, gap = 1.0, 0.08
    for name in GROUPS:
        sub = rows[rows["group"] == name]
        if sub.empty:
            continue
        ax.text(-0.02, y, f"{name} ({len(sub)})", transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=5.5, fontweight="bold")
        y += 1
        for r in sub.itertuples():
            arrow = "—" if r.mark == "undir" or not r.true_edge else "→"
            a, b = r.parent, r.child
            ticks.append(y)
            ticklabels.append(f"{labels.get(a, a)} {arrow} {labels.get(b, b)}")
            for j, ds in enumerate(style.DATASETS):
                if ds not in r.data.index:
                    continue
                cell, left = r.data.loc[ds], j * (width + gap)
                parts = [("added", cell["added_rate"])] if not r.true_edge else [(c, cell[RATE[c]]) for c in CATS]
                for c, share in parts:
                    if share > 0:
                        ax.add_patch(Rectangle((left, y - 0.38), width * share, 0.76,
                                               color=style.EDGE_STYLE[c][0], lw=0))
                    left += width * share
            y += 1
    ax.set_xlim(0, 3 * width + 2 * gap)
    ax.set_ylim(y - 0.5, -0.5)
    ax.set_yticks(ticks, ticklabels, fontsize=5)
    ax.tick_params(axis="y", length=0)
    ax.set_xticks([j * (width + gap) + width / 2 for j in range(3)], [style.DATASET_SHORT[ds] for ds in style.DATASETS])
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    handles = [Patch(color=style.EDGE_STYLE[c][0], label=style.EDGE_STYLE[c][2]) for c in CATS + ["added"]]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    style.save(fig, "fig2b_edge_strip")
    print(rows.groupby("group").size().reindex(GROUPS))


if __name__ == "__main__":
    main()
