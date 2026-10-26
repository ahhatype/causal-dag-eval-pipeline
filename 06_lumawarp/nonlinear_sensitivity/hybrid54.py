"""Hybrid discovery on the variant: PC-stable (Fisher-z, alpha 0.01) on the 36 analysis variables, XGBoost gap rule over the 153
continuous pairs (8,000 / 2,000 split, three noise columns), chain veto (folder 29). Replicates 1 to NREP. CPU only."""
import os, sys, warnings, numpy as np, pandas as pd, xgboost as xgb; warnings.filterwarnings("ignore")
from multiprocessing import Pool
from graph54 import *
H = os.path.dirname(os.path.abspath(__file__)); NREP = int(os.environ.get("NREP", 20)); BAND, GAIN = 0.10, 0.02; NTR = 8000
TADJ = {frozenset(e) for e in TRUE}; CN = sorted(CONT); POS = {v: i for i, v in enumerate(TOPO)}
def r2(df, feats, target):
    m = xgb.XGBRegressor(n_estimators=150, max_depth=3, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, n_jobs=1, verbosity=0).fit(df.iloc[:NTR][feats], df.iloc[:NTR][target])
    y = df.iloc[NTR:][target].to_numpy(); yh = m.predict(df.iloc[NTR:][feats]); return 1 - ((y - yh) ** 2).sum() / ((y - y.mean()) ** 2).sum()
def one(rep):
    from causallearn.search.ConstraintBased.PC import pc
    df = full_set(rep); r_ = np.random.default_rng(97000 + rep)
    for k in range(3): df[f"null_cont_{k}"] = r_.standard_normal(len(df))
    nulls = ["null_cont_0", "null_cont_1", "null_cont_2"]
    X = df[NODES].to_numpy(float); keep = [i for i in range(len(NODES)) if X[:, i].std() > 0]; X = X[:, keep] + 1e-9 * r_.standard_normal((len(X), len(keep)))
    G = pc(X, alpha=0.01, indep_test="fisherz", stable=True, show_progress=False).G.graph; pc_adj = set(); pc_dir = set()
    for i in range(len(keep)):
        for j in range(i + 1, len(keep)):
            if G[i, j] == 0 and G[j, i] == 0: continue
            u, v = NODES[keep[i]], NODES[keep[j]]; pc_adj.add(frozenset((u, v)))
            if G[j, i] == 1 and G[i, j] == -1: pc_dir.add((u, v))
            elif G[i, j] == 1 and G[j, i] == -1: pc_dir.add((v, u))
    called = {}; gaps = []
    for i, a in enumerate(CN):
        for b in CN[i + 1:]:
            g = r2(df, [a] + nulls, b) - r2(df, [b] + nulls, a); gaps.append(dict(rep=rep, a=a, b=b, gap=g, true_edge=frozenset((a, b)) in TADJ))
            if g > BAND: called[(a, b)] = g
            elif g < -BAND: called[(b, a)] = -g
    hyb_dir = {e for e in pc_dir if frozenset(e) not in {frozenset(c) for c in called}} | set(called); vetoed = set()
    for (a, c) in list(called):
        for b in [b for (x, b) in hyb_dir if x == a and (b, c) in hyb_dir]:
            if r2(df, [a, b] + nulls, c) - r2(df, [b] + nulls, c) < GAIN: vetoed.add((a, c))
    kept = {e: g for e, g in called.items() if e not in vetoed}
    rows = []
    for name, adj, dire in (("PC-stable alone", pc_adj, pc_dir), ("Hybrid, no veto", pc_adj | {frozenset(e) for e in called}, {e for e in pc_dir if frozenset(e) not in {frozenset(c) for c in called}} | set(called)),
                            ("Hybrid, chain veto", pc_adj | {frozenset(e) for e in kept}, {e for e in pc_dir if frozenset(e) not in {frozenset(c) for c in kept}} | set(kept))):
        tp = len(adj & TADJ); prec = tp / max(1, len(adj)); rec = tp / len(TADJ); ot = [e for e in dire if frozenset(e) in TADJ]
        rows.append(dict(rep=rep, method=name, edges=len(adj), precision=prec, recall=rec, f1=2 * prec * rec / max(1e-9, prec + rec), u_edges_found=sum(frozenset(e) in adj for e in UEDGES), u_edges_oriented=sum(e in dire for e in UEDGES),
                         saturating_found=frozenset(SAT) in adj, false_edges=len(adj - TADJ), added_false=len({frozenset(e) for e in (kept if "veto" in name else called)} - pc_adj - TADJ) if "Hybrid" in name else np.nan, oriented_true=len(ot), oriented_correct=sum(e in TRUE for e in ot)))
    print("rep", rep, "done", flush=True); return rows, gaps
if __name__ == "__main__":
    with Pool(int(os.environ.get("NPROC", 4))) as pool: res = pool.map(one, range(1, NREP + 1))
    S = pd.DataFrame([r for rr, _ in res for r in rr]); Gp = pd.DataFrame([g for _, gg in res for g in gg]); S.to_csv(os.path.join(H, "hybrid_by_rep.csv"), index=False); Gp.to_csv(os.path.join(H, "gaps.csv"), index=False)
    M = S.groupby("method").agg(reps=("rep", "nunique"), edges=("edges", "mean"), precision=("precision", "mean"), recall=("recall", "mean"), f1=("f1", "mean"), u_edges_found=("u_edges_found", "mean"), u_edges_oriented=("u_edges_oriented", "mean"), saturating_found=("saturating_found", "mean"), false_edges=("false_edges", "mean"), added_false=("added_false", "mean"), orient_acc=("oriented_correct", "sum"))
    M["orient_acc"] = S.groupby("method").oriented_correct.sum() / S.groupby("method").oriented_true.sum(); pd.set_option("display.width", 250); print(M.round(3).to_string()); M.to_csv(os.path.join(H, "hybrid_summary.csv"))
    print("\nU-edge gaps (true direction positive):"); print(Gp[Gp.apply(lambda r: (r.a, r.b) in UEDGES or (r.b, r.a) in UEDGES, axis=1)].groupby(["a", "b"]).gap.describe()[["mean", "min", "max"]].round(3).to_string())
