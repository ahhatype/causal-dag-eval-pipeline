"""Attribution accuracy across the three data conditions, as a 2 x 2 grid.

a, b  Kendall's tau_b against the per-unit and population-scaled truths (ancestor set).
c     Proximity bias index (ancestor set); filled = per-unit truth, open = population-scaled truth.
d     Share of attribution given to non-ancestors (all-features set; truth-independent).
Each data set is scored against its own population (source for full set and reference subsample,
selected for astronaut set), mission era included. Error bars are +/- 1.96 Monte Carlo SE over 50 replicates.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import style

ATTR = style.ROOT / "outputs" / "04_attribution_validation"


def load():
    d = pd.read_csv(ATTR / "summary.csv")
    own = np.where(d["dataset"] == "astronaut_set", "selected", "source")
    d = d[(d["population"] == own) & (d["scope"] == "with_era") & d["method"].isin(style.METHODS)]
    return d


def slope(ax, d, truth, title):
    x = np.arange(3)
    for m in style.METHODS:
        r = d[(d.feature_set == "ancestor") & (d.truth_type == truth) & (d.method == m)].set_index("dataset")
        r = r.loc[style.DATASETS]
        ax.errorbar(x, r.kendall_tau_b_mean, yerr=1.96 * r.kendall_tau_b_mcse, color=style.METHOD_COLORS[m],
                    ls=style.METHOD_LINES[m], lw=1.4, marker="o", ms=4.5, capsize=2, elinewidth=0.8,
                    label=style.METHOD_LABELS[m])
    ax.set_xticks(x, [style.DATASET_SHORT[k] for k in style.DATASETS])
    ax.set_ylim(0, 1)
    ax.set_ylabel(r"Kendall's $\tau_b$ with true effects")
    ax.set_title(title, loc="left", fontweight="bold")


def main():
    style.apply()
    d = load()
    out = d[["dataset", "feature_set", "truth_type", "method", "kendall_tau_b_mean", "kendall_tau_b_mcse",
             "proximity_bias_index_mean", "proximity_bias_index_mcse", "non_ancestor_share_mean",
             "non_ancestor_share_mcse", "top_k_recovery_mean", "top_k_recovery_mcse"]]
    style.OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(style.OUT / "fig3_source_data.csv", index=False)

    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8.5, "xtick.labelsize": 8,
                         "ytick.labelsize": 8, "legend.fontsize": 8, "axes.linewidth": 0.7})
    fig, grid = plt.subplots(2, 2, figsize=(style.DOUBLE, 135 * style.MM),
                             gridspec_kw={"hspace": 0.5, "wspace": 0.3})
    axes = grid.ravel()
    slope(axes[0], d, "per_unit", "a  Intervention-contrast truth")
    slope(axes[1], d, "pop", "b  Population-scaled truth")

    ax = axes[2]
    w = 0.18
    for i, m in enumerate(style.METHODS):
        for truth, filled in (("per_unit", True), ("pop", False)):
            r = d[(d.feature_set == "ancestor") & (d.truth_type == truth) & (d.method == m)].set_index("dataset")
            r = r.loc[style.DATASETS]
            xs = np.arange(3) + (i - 1.5) * w + (0 if filled else 0.06)
            c = style.METHOD_COLORS[m]
            ax.errorbar(xs, r.proximity_bias_index_mean, yerr=1.96 * r.proximity_bias_index_mcse, fmt="o",
                        ms=4.2, color=c, mfc=c if filled else "white", mew=0.9, capsize=1.6, elinewidth=0.8)
    ax.axhline(0, color="#8C8C8C", lw=0.5, zorder=0)
    ax.set_xticks(np.arange(3), [style.DATASET_SHORT[k] for k in style.DATASETS])
    ax.set_ylabel("Proximity bias index (edges)\n+ = credit nearer the outcome")
    ax.set_title("c  Proximity bias", loc="left", fontweight="bold")

    ax = axes[3]
    for i, m in enumerate(style.METHODS):
        r = d[(d.feature_set == "all") & (d.truth_type == "per_unit") & (d.method == m)].set_index("dataset")
        r = r.loc[style.DATASETS]
        ax.bar(np.arange(3) + (i - 1.5) * w, 100 * r.non_ancestor_share_mean, width=w * 0.92,
               yerr=196 * r.non_ancestor_share_mcse, color=style.METHOD_COLORS[m], error_kw={"elinewidth": 0.8,
               "capsize": 1.6})
    ax.set_xticks(np.arange(3), [style.DATASET_SHORT[k] for k in style.DATASETS])
    ax.set_ylim(0, 60)
    ax.text(0.98, 0.98, "Causal predictive SHAP:\n0 by construction", transform=ax.transAxes, fontsize=7.5,
            ha="right", va="top", color=style.METHOD_COLORS["ng_style_shap"])
    ax.set_ylabel("Credit to non-ancestors (%)")
    ax.set_title("d  Credit to non-causes", loc="left", fontweight="bold")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.995))
    style.save(fig, "fig3_attribution")
    print(out[out.feature_set == "ancestor"].query("truth_type == 'per_unit'")[
        ["dataset", "method", "kendall_tau_b_mean", "proximity_bias_index_mean"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
