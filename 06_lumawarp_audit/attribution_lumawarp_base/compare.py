"""Pairs every LumaWarp-base cell with the pipeline's committed Super Learner result for the same replicate, data set
and method; prints the pre-registered bars."""
import glob, os, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); PIPE = os.path.abspath(os.path.join(H, "..", ".."))
lw = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(H, "cells", "*_runs.csv")))], ignore_index=True)
sl = pd.read_csv(os.path.join(PIPE, "outputs", "04_attribution_validation", "per_run.csv"))
own = {"full_set": "source", "reference_subsample": "source", "astronaut_set": "selected"}
def primary(d): return d[(d.truth_type == "per_unit") & (d.scope == "with_era") & (d.population == d.dataset.map(own))]
lw = primary(lw); lw["base"] = lw.method.str.replace("_lumawarp", ""); sl = primary(sl[sl.feature_set == "ancestor"])
m = lw.merge(sl, left_on=["replicate", "dataset", "base"], right_on=["replicate", "dataset", "method"], suffixes=("_lw", "_sl"))
pd.set_option("display.width", 250); out = []
for (ds, b), g in m.groupby(["dataset", "base"]):
    d = g.kendall_tau_b_lw - g.kendall_tau_b_sl
    out.append(dict(dataset=ds, method=b, reps=len(g), tau_lw=g.kendall_tau_b_lw.mean(), tau_sl=g.kendall_tau_b_sl.mean(), diff=d.mean(), diff_mcse=d.std(ddof=1) / np.sqrt(len(g)) if len(g) > 1 else np.nan,
                    top5_lw=g.top_k_recovery_lw.mean(), top5_sl=g.top_k_recovery_sl.mean(), pbi_lw=g.proximity_bias_index_lw.mean(), pbi_sl=g.proximity_bias_index_sl.mean(),
                    auc_lw=g.auc_lumawarp.mean(), auc_sl=g.auc_superlearner.mean()))
S = pd.DataFrame(out); S.to_csv(os.path.join(H, "summary.csv"), index=False); print(S.round(3).to_string(index=False))
w = lw.pivot_table(index=["dataset", "replicate"], columns="base", values="kendall_tau_b")
for ds in w.index.get_level_values(0).unique():
    x = w.loc[ds]
    if {"interventional_shap", "standard_shap"} <= set(x.columns): print(f"\n{ds}: interventional minus standard (LumaWarp base) = {(x.interventional_shap - x.standard_shap).mean():.3f} over {len(x)} reps")
    if {"ordering_only_shap", "standard_shap"} <= set(x.columns): print(f"{ds}: ordering-only minus standard (LumaWarp base) = {(x.ordering_only_shap - x.standard_shap).mean():.3f}")
f = S[(S.dataset == "full_set") & (S.method == "interventional_shap")]
if len(f):
    f = f.iloc[0]; print(f"\nB1: |tau diff| {abs(f['diff']):.3f} (bar 0.05), |AUC diff| {abs(f.auc_lw - f.auc_sl):.3f} (bar 0.01) -> {'PASS' if abs(f['diff']) <= 0.05 and abs(f.auc_lw - f.auc_sl) <= 0.01 else 'FAIL'}")
    x = w.loc["full_set"]
    if {"interventional_shap", "standard_shap"} <= set(x.columns): c = (x.interventional_shap - x.standard_shap).mean(); print(f"B2: contrast {c:.3f} (bar 0.30) -> {'PASS' if c >= 0.30 else 'FAIL'}")
