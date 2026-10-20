"""Edge audit scoring with folder 13's thresholds, against the analysis graph; shape flags from saved models. Verdicts per PREREGISTRATION.md."""
import os, numpy as np, pandas as pd
from graph54 import NODES, PA, TRUE, CONT
from lw_kernel_dump import read_kernel, members
H = os.path.dirname(os.path.abspath(__file__)); THR = {"lumawarp": 5, "glm": 4, "spline": 2, "xgb": 3}; NAMES = {"lumawarp": "LumaWarp", "glm": "logistic/linear", "spline": "spline additive", "xgb": "XGBoost"}
PRED = {"lumawarp": "pred.csv", "glm": "pred_glm.csv", "spline": "pred_spline.csv", "xgb": "pred_xgb.csv"}; jobs = pd.read_csv(os.path.join(H, "jobs.csv")); rows = []; flags = []
for j in jobs.itertuples():
    jd = os.path.join(H, j.dir); m0 = pd.read_csv(os.path.join(jd, "score_meta.csv"))
    for tag, pf in PRED.items():
        pp = os.path.join(jd, pf)
        if not os.path.exists(pp): continue
        p = pd.read_csv(pp); m = m0.copy(); m["p"] = (p["p_1"] if "p_1" in p else p["prediction"]).to_numpy(); s = m.groupby("block").p.agg(lambda x: x.max() - x.min()); nm = s[s.index.str.startswith("null_")]
        for c in s.index:
            isn = c.startswith("null_"); ref = max((nm.drop(c) if isn else nm).max(), 1e-12); rows.append(dict(dataset=j.dataset, rep=j.rep, target=j.target, cand=c, null=isn, learner=tag, score=s[c] / ref, kept=s[c] / ref > THR[tag], true=(not isn) and (c in PA[j.target])))
    mp = os.path.join(jd, "model.lsmdl")
    if os.path.exists(mp) and j.target in CONT and not j.binary:
        acc = {}
        try:
            for mem, secs in members(mp):
                _, axes = read_kernel(secs["KERNEL"]); per = {}
                for ax in axes:
                    if ax["name"].startswith("proj#") or ax["space"] != 0: continue
                    K = ax["A"].T @ ax["A"]; d = np.diag(K); cen = float((np.arange(len(d)) * d).sum() / d.sum()); col = ax["name"].replace("#raw", ""); a_ = per.setdefault(col, [0.0, 0.0]); a_[0] += ax["theta"] * cen; a_[1] += ax["theta"]
                for col, (ws, ts) in per.items(): acc.setdefault(col, []).append(ws / ts)
            cen = {c: float(np.mean(v)) for c, v in acc.items()}; nz = max(v for c, v in cen.items() if c.startswith("null_cont"))
            for pnode in PA[j.target]:
                if pnode in CONT and pnode in cen: flags.append(dict(dataset=j.dataset, rep=j.rep, target=j.target, parent=pnode, centroid=cen[pnode], max_noise=nz, flagged=cen[pnode] > nz))
        except Exception as e: print("slice read failed", j.dir, e)
E = pd.DataFrame(rows); E.to_csv(os.path.join(H, "edges_all.csv"), index=False); out = []
for (ds, rep, tag), g in E.groupby(["dataset", "rep", "learner"]):
    if g.target.nunique() < jobs[(jobs.dataset == ds) & (jobs.rep == rep)].shape[0]: continue
    r = g[~g.null]; called = r[r.kept]; tp = called.true.sum(); out.append(dict(dataset=ds, rep=rep, learner=NAMES[tag], edges_called=len(called), recall=tp / len(TRUE), precision=tp / max(1, len(called)), noise_kept=g[g.null].kept.mean()))
S = pd.DataFrame(out); S.to_csv(os.path.join(H, "scores_by_rep.csv"), index=False)
M = S.groupby(["dataset", "learner"]).agg(reps=("rep", "nunique"), edges_called=("edges_called", "mean"), recall=("recall", "mean"), recall_mcse=("recall", lambda x: x.std() / np.sqrt(len(x))), precision=("precision", "mean"), noise_kept=("noise_kept", "mean")); pd.set_option("display.width", 250); print(M.round(3).to_string()); M.to_csv(os.path.join(H, "summary.csv"))
if flags:
    F = pd.DataFrame(flags); F.to_csv(os.path.join(H, "shape_flags.csv"), index=False); print("\nShape flags (continuous parents of continuous targets):"); print(F.groupby("dataset").flagged.agg(["sum", "size"]).to_string())
lw = M.xs("LumaWarp", level="learner") if "LumaWarp" in M.index.get_level_values(1) else None
if lw is not None and "full_set" in lw.index:
    f = lw.loc["full_set"]; print(f"\nE1 full-set recall {f.recall:.3f} (bar 0.80), precision {f.precision:.3f} (bar 0.90) -> {'PASS' if (f.recall >= .8 and f.precision >= .9) else 'FAIL'}")
    for ds in M.index.get_level_values(0).unique():
        g = M.loc[ds]; best = g.drop(index="LumaWarp").recall.max(); print(f"  E3 {ds}: LumaWarp recall {g.loc['LumaWarp'].recall:.3f} vs best reference {best:.3f} ({g.loc['LumaWarp'].recall - best:+.3f})")
