"""Table 1 as LaTeX (booktabs + threeparttable) from table1_ordering.csv and ancestor_displacement.csv.

Rows per data set: the 16 ancestors in per-unit truth order, then non-ancestors that any method ranked in
its top --top features. Cells are "rank (share%)"; tied truth ranks carry a dagger. An asterisk marks an ancestor
whose order among ancestors differed between the ancestor-set and all-features runs.
Use --full for the supplementary version with all 35 features. The main-text version (default) is one
landscape table with the data sets side by side, ranks only, against the intervention-contrast truth.
"""
import argparse

import pandas as pd

import style

ATTR = style.ROOT / "outputs" / "04_attribution_validation"
TRUTHS = ["per_unit", "pop"]
TIE = r"$^{\dagger}$"
HEAD = {"per_unit": r"\PerUnitLabel", "pop": "Pop.-scaled", "standard_shap": "Standard",
        "ordering_only_shap": "Ordering-only", "ng_style_shap": "Causal pred.", "interventional_shap": "Causal int."}
BLOCK = {"full_set": "Full-set ($n = 10{,}000$)", "reference_subsample": "Reference-set ($n = 900$)",
         "astronaut_set": "Astronaut-set ($n = 900$, selected)"}


def esc(s):
    return s.replace("&", r"\&").replace("%", r"\%").replace("+", "+")


def cell(r, tied, displaced, truth):
    if (truth and not r.ancestor) or (not truth and r.share_pct == 0):
        return "--"
    rank = f"{int(r['rank'])}{TIE if tied else ''}"
    star = r"$^{*}$" if displaced else ""
    return rf"{rank}{star}\,{{\scriptsize({r.share_pct:.1f})}}"


def block(t, disp, ds, top, full):
    d = t[t.dataset == ds]
    truth_rank = d[(d.column == "per_unit")].set_index("feature")
    methods = d[d.column_type == "method"]
    if full:
        keep = set(d.feature)
    else:
        lead = methods[(~methods.ancestor) & (methods["rank"] <= top)].feature
        keep = set(d[d.ancestor].feature) | set(lead)
    ties = d[d.column_type == "truth"].groupby(["column", "rank"]).feature.transform("count") > 1
    d = d.assign(tied=False)
    d.loc[ties.index, "tied"] = ties
    best = methods.groupby("feature")["rank"].min()
    order = sorted(keep, key=lambda f: (not truth_rank.loc[f, "ancestor"], truth_rank.loc[f, "rank"]
                                        if truth_rank.loc[f, "ancestor"] else best[f], f))
    dd = disp[disp.dataset == ds].set_index(["feature", "method"])["displaced"]
    lines = [rf"\multicolumn{{7}}{{l}}{{\textbf{{{BLOCK[ds]}}}}}\\", r"\addlinespace[2pt]"]
    section = None
    for f in order:
        anc = bool(truth_rank.loc[f, "ancestor"])
        if anc != section:
            section = anc
            label = "Ancestors" if anc else (r"Non-ancestors" if full else
                                             rf"Non-ancestors ranked $\le {top}$ by any method")
            lines.append(rf"\multicolumn{{7}}{{l}}{{\textit{{{label}}}}}\\")
        cells = []
        for col in TRUTHS + style.METHODS:
            r = d[(d.feature == f) & (d.column == col)].iloc[0]
            displaced = bool(dd.get((f, col), False)) if col in style.METHODS else False
            cells.append(cell(r, r.tied and col in TRUTHS, displaced, col in TRUTHS))
        name = esc(d[d.feature == f].feature_name.iloc[0])
        lines.append(rf"\quad {name} & " + " & ".join(cells) + r" \\")
    return lines


SHORT = {"per_unit": "Truth", "standard_shap": "Std.", "ordering_only_shap": "Ord.", "ng_style_shap": "Pred.",
         "interventional_shap": "Int."}
COLS = ["per_unit"] + style.METHODS


def wide(t, disp, top):
    keep_anc = set(t[t.ancestor].feature)
    m = t[(t.column_type == "method") & (~t.ancestor) & (t["rank"] <= top)]
    keep_non = set(m.feature)
    ref = t[(t.dataset == "full_set") & (t.column == "per_unit")].set_index("feature")["rank"]
    best = t[t.column_type == "method"].groupby("feature")["rank"].min()
    anc = sorted(keep_anc, key=lambda f: (ref[f], f))
    non = sorted(keep_non, key=lambda f: (best[f], f))
    ties = t[t.column_type == "truth"].groupby(["dataset", "column", "rank"]).feature.transform("count") > 1
    t = t.assign(tied=False)
    t.loc[ties.index, "tied"] = ties
    idx = t.set_index(["dataset", "feature", "column"])
    dd = disp.set_index(["dataset", "feature", "method"])["displaced"]
    names = t.drop_duplicates("feature").set_index("feature")["feature_name"]
    ncol = 1 + len(COLS) * len(style.DATASETS)

    def row(f):
        cells = []
        for ds in style.DATASETS:
            for c in COLS:
                r = idx.loc[(ds, f, c)]
                if (c == "per_unit" and not r.ancestor) or (c != "per_unit" and r.share_pct == 0):
                    cells.append("--")
                    continue
                s = f"{int(r['rank'])}{TIE if (c == 'per_unit' and r.tied) else ''}"
                if c != "per_unit" and bool(dd.get((ds, f, c), False)):
                    s += r"$^{*}$"
                cells.append(s)
        return rf"\quad {esc(names[f])} & " + " & ".join(cells) + r" \\"

    body = [rf"\multicolumn{{{ncol}}}{{l}}{{\textit{{Ancestors}}}}\\", *[row(f) for f in anc],
            r"\addlinespace[2pt]",
            rf"\multicolumn{{{ncol}}}{{l}}{{\textit{{Non-ancestors ranked $\le {top}$ by any method in any data set}}}}\\",
            *[row(f) for f in non]]
    groups = " & ".join(rf"\multicolumn{{{len(COLS)}}}{{c}}{{{style.DATASET_SHORT[ds]}}}" for ds in style.DATASETS)
    rules = "".join(rf"\cmidrule(lr){{{2 + i * len(COLS)}-{1 + (i + 1) * len(COLS)}}}" for i in range(3))
    heads = " & ".join(SHORT[c] for c in COLS)
    spec = "@{}l " + " ".join("r" * len(COLS) for _ in style.DATASETS) + "@{}"
    return [
        r"% Generated by 05_figures/table1_latex.py; edit the script, not this file.",
        r"% Preamble: \usepackage{booktabs,threeparttable,rotating}",
        r"\begin{sidewaystable}", r"\centering", r"\small", r"\begin{threeparttable}",
        r"\caption{Rank of each feature by true effect and by attribution method in each data set.}",
        r"\label{tab:ordering}",
        rf"\begin{{tabular}}{{{spec}}}", r"\toprule",
        f" & {groups} \\\\", rules,
        "Feature & " + " & ".join([heads] * 3) + r" \\", r"\midrule",
        *body, r"\bottomrule", r"\end{tabular}", r"\begin{tablenotes}\footnotesize",
        r"\item Truth, rank of the intervention-contrast true effect among the 16 ancestors, for each data set's own "
        r"population (source for the full-set and reference-set, selected for the astronaut-set); $^{\dagger}$tied rank. "
        r"Std., standard; Ord., ordering-only; Pred., causal predictive; Int., causal interventional SHAP, ranked among "
        r"all 35 features by mean $|\mathrm{SHAP}|$ over 50 replicates of the all-features run. -- marks zero credit "
        r"(causal predictive SHAP assigns none to non-ancestors by construction). $^{*}$Order among ancestors differed "
        r"between the ancestor-set and all-features runs. Shares and the population-scaled truth are given in "
        r"Supplementary Table~\ref{tab:ordering-full}.",
        r"\end{tablenotes}", r"\end{threeparttable}", r"\end{sidewaystable}"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=style.DATASETS)
    a = ap.parse_args()
    t = pd.read_csv(ATTR / "table1_ordering.csv")
    disp = pd.read_csv(ATTR / "ancestor_displacement.csv")
    if not a.full:
        tex = wide(t, disp, a.top)
        (style.OUT / "table1_ordering.tex").write_text("\n".join(tex) + "\n")
        print(f"wrote table1_ordering.tex: {sum(1 for l in tex if l.startswith(chr(92) + 'quad'))} rows")
        return
    body = []
    for i, ds in enumerate(a.datasets):
        if i:
            body.append(r"\midrule")
        body += block(t, disp, ds, a.top, a.full)
    head = [r"\toprule", r" & \multicolumn{2}{c}{True effect} & \multicolumn{4}{c}{SHAP method} \\",
            r"\cmidrule(lr){2-3}\cmidrule(l){4-7}",
            "Feature & " + " & ".join(HEAD[c] for c in TRUTHS + style.METHODS) + r" \\", r"\midrule"]
    notes = (r"Ranks are among all 35 features (methods) or the 16 ancestors (truths); $^{\dagger}$tied truth rank. "
             r"Shares sum to 100\% within each column. True effects are for each data set's own population (source "
             r"for the full-set and reference-set, selected for the astronaut-set). Method ranks use mean "
             r"$|\mathrm{SHAP}|$ over 50 replicates of the all-features run; -- marks zero credit (causal predictive "
             r"SHAP assigns none to non-ancestors by construction). "
             r"$^{*}$Order among ancestors differed between the ancestor-set and all-features runs.")
    tex = [
        r"% Generated by 05_figures/table1_latex.py; edit the script, not this file.",
        r"% Preamble: \usepackage{booktabs,longtable,pdflscape}; define \newcommand{\PerUnitLabel}{Intervention-contrast}",
        r"\begin{landscape}", r"\begingroup", r"\small", r"\setlength{\LTcapwidth}{\linewidth}",
        r"\begin{longtable}{@{}l rr rrrr@{}}",
        r"\caption{Rank of each feature by true effect and by attribution method, with share of total importance (\%).}"
        r"\label{tab:ordering-full}\\",
        *head, r"\endfirsthead", *head, r"\endhead", r"\bottomrule", r"\endfoot",
        *body, r"\end{longtable}", r"\par\footnotesize " + notes, r"\endgroup", r"\end{landscape}"]
    name = "table1_ordering_full.tex" if a.full else "table1_ordering.tex"
    style.OUT.mkdir(parents=True, exist_ok=True)
    (style.OUT / name).write_text("\n".join(tex) + "\n")
    print(f"wrote {name}: {sum(1 for l in body if l.startswith(chr(92) + 'quad'))} rows")


if __name__ == "__main__":
    main()
