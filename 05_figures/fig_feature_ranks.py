"""Feature-rank dot plot: where each attribution method ranks each feature, against the true rank.

One panel per data set. Rows are the 16 ancestors in order of the full-set intervention-contrast truth,
then the non-ancestors that any method ranks within the top 16 in any data set. The x-axis is rank among
all 35 features by mean |SHAP| over 50 replicates of the all-features run (the ranks in
table1_ordering.csv). Truth is the rank among the 16 ancestors for the data set's own population; tied
truths share the best rank. A black outline marks an ancestor whose order among ancestors differed between
the ancestor-set and all-features runs. Causal predictive SHAP gives non-ancestors no credit, so it has
no marker on those rows.
"""
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

import style

A = style.ROOT / "outputs" / "04_attribution_validation"
TRUTH = "per_unit"
N_ANC = 16
MARKERS = {"standard_shap": "o", "ordering_only_shap": "^", "ng_style_shap": "D", "interventional_shap": "s"}
DODGE = {"standard_shap": -0.27, "ordering_only_shap": -0.09, "ng_style_shap": 0.09, "interventional_shap": 0.27}


def load():
    t = pd.read_csv(A / "table1_ordering.csv")
    d = pd.read_csv(A / "ancestor_displacement.csv")
    disp = set(zip(d.loc[d.displaced, "dataset"], d.loc[d.displaced, "method"], d.loc[d.displaced, "feature"]))
    truth = t[(t.column_type == "truth") & (t.column == TRUTH)]
    meth = t[t.column_type == "method"]
    anc = (truth[(truth.dataset == "full_set") & truth.ancestor]
           .sort_values(["rank", "value"], ascending=[True, False]).feature.tolist())
    intruders = (meth[~meth.ancestor & (meth["rank"] <= N_ANC)].groupby("feature")["rank"].min()
                 .sort_values().index.tolist())
    names = dict(zip(t.feature, t.feature_name))
    return truth, meth, disp, anc, intruders, names


def main():
    style.apply()
    truth, meth, disp, anc, non, names = load()
    rows = anc + non
    ypos = {f: i if i < len(anc) else i + 0.8 for i, f in enumerate(rows)}
    height = (len(rows) + 1.5) * 4.0 * style.MM + 18 * style.MM
    fig, axes = plt.subplots(1, 3, figsize=(style.DOUBLE, height), sharey=True,
                             gridspec_kw={"wspace": 0.08})
    for ax, ds in zip(axes, style.DATASETS):
        ax.axvspan(N_ANC + 0.5, 36, color="0.94", zorder=0, lw=0)
        ax.axhline(len(anc) - 0.1, color="0.6", lw=0.4, ls=(0, (2, 2)))
        for f in rows:
            ax.axhline(ypos[f], color="0.9", lw=0.3, zorder=0)
        tr = truth[truth.dataset == ds].set_index("feature")["rank"]
        for f in anc:
            ax.plot([tr[f], tr[f]], [ypos[f] - 0.42, ypos[f] + 0.42], color="black", lw=1.0, zorder=3,
                    solid_capstyle="butt")
        for m in style.METHODS:
            sub = meth[(meth.dataset == ds) & (meth.column == m)].set_index("feature")
            for f in rows:
                if f not in sub.index or sub.loc[f, "value"] <= 0:
                    continue
                outline = (ds, m, f) in disp
                ax.scatter(sub.loc[f, "rank"], ypos[f] + DODGE[m], marker=MARKERS[m], s=11,
                           facecolor=style.METHOD_COLORS[m], edgecolor="black" if outline else "none",
                           lw=0.6 if outline else 0, zorder=5)
        ax.set_xlim(0.3, 35.7)
        ax.set_xticks([1, 5, 10, 16, 20, 25, 30, 35])
        ax.set_title(style.DATASET_LABELS[ds].replace("\n", " "), fontsize=6.5, pad=3)
        ax.set_xlabel("Rank among 35 features")
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    axes[0].set_yticks([ypos[f] for f in rows])
    axes[0].set_yticklabels([names[f] for f in rows])
    axes[0].set_ylim(max(ypos.values()) + 0.7, -0.7)
    for label, block in (("Ancestors", anc), ("Non-ancestors", non)):
        axes[0].text(-0.80, ypos[block[0]] + (len(block) - 1) / 2, label, rotation=90, ha="center",
                     va="center", transform=axes[0].get_yaxis_transform(), fontsize=6, style="italic")
    handles = [Line2D([], [], marker="|", ls="", color="black", ms=6, mew=1.1, label="Truth (intervention-contrast)")]
    handles += [plt.Rectangle((0, 0), 1, 1, fc="0.94", ec="none", label="Below at least one non-ancestor")]
    handles += [Line2D([], [], marker=MARKERS[m], ls="", ms=3.6, mfc=style.METHOD_COLORS[m], mec="none",
                       label=style.METHOD_LABELS[m]) for m in style.METHODS]
    handles += [Line2D([], [], marker="o", ls="", ms=3.6, mfc="white", mec="black", mew=0.6,
                       label="Ancestor order changed by adding non-ancestors")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0),
               columnspacing=1.2, handletextpad=0.3)
    fig.subplots_adjust(bottom=0.13)
    style.save(fig, "fig_feature_ranks")


if __name__ == "__main__":
    main()
