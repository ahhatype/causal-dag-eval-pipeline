"""Each importance-file column as a shape reading under folder 14's noise rule, scored against folder 28's truth and centroid flag."""
import os, sys, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, cohen_kappa_score
H = os.path.dirname(os.path.abspath(__file__)); S = os.path.join(H, "..", "28_shape_flag_range")
truth = pd.read_csv(os.path.join(S, "shape_flags.csv")); jobs = pd.read_csv(os.path.join(H, "jobs.csv")); rows = []
READINGS = {"h0_variance": "h0 variance", "h1_variance": "h1 variance", "relvar": None, "importance": "importance"}
for j in jobs.itertuples():
    fp = os.path.join(H, j.dir, "model.lsmdl.feature_importance.csv")
    if not os.path.exists(fp): continue
    fi = pd.read_csv(fp); fi = fi[~fi["feature name"].str.contains(" x ")].set_index("feature name")
    fi["relvar"] = fi["importance"] / fi["theta"].replace(0, np.nan)
    t = truth[(truth.a == j.a) & (truth.s == j.s) & (truth.n == j.n) & (truth.seed == j.seed) & (truth.target == j.target)]
    for name, col in READINGS.items():
        v = fi[col if col else name]; nz = v[[c for c in v.index if c.startswith("null_cont")]].max()
        for c in t.itertuples():
            if c.parent not in v.index: continue
            rows.append(dict(reading=name, a=j.a, s=j.s, n=j.n, seed=j.seed, target=j.target, parent=c.parent, shape=c.shape, nonmonotone=c.nonmonotone, value=float(v[c.parent]), max_noise=float(nz), flagged=bool(v[c.parent] > nz), centroid_flagged=bool(c.flagged)))
R = pd.DataFrame(rows); R.to_csv(os.path.join(H, "flags_by_reading.csv"), index=False); pd.set_option("display.width", 250)
cen_nm = int(truth[truth.nonmonotone].flagged.sum()); out = []
for name, g in R.groupby("reading"):
    nm = g[g.nonmonotone]; mo = g[~g.nonmonotone]
    out.append(dict(reading=name, cases=len(g), nonmono_flagged=int(nm.flagged.sum()), of=len(nm), mono_false=int(mo.flagged.sum()), of_mono=len(mo),
                    auc=roc_auc_score(g.nonmonotone, g.value / g.max_noise), kappa_vs_centroid=cohen_kappa_score(g.flagged, g.centroid_flagged),
                    F1=(nm.flagged.sum() >= 0.8 * cen_nm) and (mo.flagged.sum() <= 4)))
Sx = pd.DataFrame(out); Sx.to_csv(os.path.join(H, "summary.csv"), index=False)
print(f"centroid flag (folder 28): {cen_nm} of {int(truth.nonmonotone.sum())} non-monotone flagged, {int(truth[~truth.nonmonotone].flagged.sum())} of {int((~truth.nonmonotone).sum())} monotone false\n")
print(Sx.round(3).to_string(index=False))
print("\nF1:", ", ".join(f"{r.reading} {'PASS' if r.F1 else 'FAIL'}" for r in Sx.itertuples()))
print("\nBy sample size (non-monotone flagged rate):"); print(R[R.nonmonotone].groupby(["reading", "n"]).flagged.mean().unstack("n").round(2).to_string())
