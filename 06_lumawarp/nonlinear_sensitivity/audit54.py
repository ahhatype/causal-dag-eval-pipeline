"""Basis-set audit on the variant: linear LRT vs basis-expanded GCM (folder 30) on (a) the analysis graph, (b) it minus the three
U-shaped edges, (c) the PC skeleton oriented by the ordering. Replicates 1 to NREP."""
import os, warnings, numpy as np, pandas as pd; warnings.filterwarnings("ignore")
from multiprocessing import Pool
from scipy.stats import chi2, norm
from sklearn.preprocessing import SplineTransformer
from sklearn.linear_model import LogisticRegression
from graph54 import *
H = os.path.dirname(os.path.abspath(__file__)); NREP = int(os.environ.get("NREP", 20)); ALPHA = 0.05; POS = {v: i for i, v in enumerate(TOPO)}
def graphs(rep):
    g = {"true": {c: set(p) for c, p in PA.items()}, "u_deleted": {c: {p for p in ps if (p, c) not in UEDGES} for c, ps in PA.items()}}
    pc = pd.read_csv(os.path.join(H, "pc_skeletons.csv")); pc = pc[pc.rep == rep]; gp = {v: set() for v in NODES}
    for r in pc.itertuples():
        a, b = (r.a, r.b) if POS[r.a] < POS[r.b] else (r.b, r.a); gp[b].add(a)
    g["pc_skeleton_ordered"] = gp; return g
def basis(pa):
    adj = {frozenset((p, c)) for c, ps in pa.items() for p in ps}; out = []
    for i, x in enumerate(TOPO):
        for y in TOPO[:i]:
            if frozenset((x, y)) in adj: continue
            out.append((x, y, sorted(pa[x] | pa[y])))
    return out
def design(df, cols, flexible):
    parts = [np.ones((len(df), 1))]
    for c in cols:
        x = df[c].to_numpy(float).reshape(-1, 1); parts.append(SplineTransformer(n_knots=6, degree=3, extrapolation="linear").fit_transform(x) if (flexible and c in CONT and df[c].nunique() > 10) else x)
    return np.hstack(parts)
def resid(df, target, X):
    y = df[target].to_numpy(float)
    if target in CONT and df[target].nunique() > 2: return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    if X.shape[1] == 1 or y.std() == 0: return y - y.mean()
    return y - LogisticRegression(C=1e4, max_iter=2000).fit(X[:, 1:], y).predict_proba(X[:, 1:])[:, 1]
def lrt_p(df, x, y, Z):
    X0 = design(df, Z, False); X1 = design(df, Z + [y], False); t = df[x].to_numpy(float)
    if x in CONT and df[x].nunique() > 2:
        r0 = t - X0 @ np.linalg.lstsq(X0, t, rcond=None)[0]; r1 = t - X1 @ np.linalg.lstsq(X1, t, rcond=None)[0]; stat = len(t) * np.log((r0 ** 2).sum() / (r1 ** 2).sum())
    else:
        if t.std() == 0: return 1.0
        def ll(X):
            p = np.full(len(t), t.mean()) if X.shape[1] == 1 else LogisticRegression(C=1e4, max_iter=2000).fit(X[:, 1:], t).predict_proba(X[:, 1:])[:, 1]
            p = np.clip(p, 1e-9, 1 - 1e-9); return (t * np.log(p) + (1 - t) * np.log(1 - p)).sum()
        stat = 2 * (ll(X1) - ll(X0))
    return chi2.sf(max(stat, 0), 1)
def gcm_basis_p(df, x, y, Z):
    X = design(df, Z, True); rx = resid(df, x, X)
    B = SplineTransformer(n_knots=6, degree=3, extrapolation="linear").fit_transform(df[y].to_numpy(float).reshape(-1, 1))[:, 1:] if (y in CONT and df[y].nunique() > 10) else df[y].to_numpy(float).reshape(-1, 1)
    RB = B - X @ np.linalg.lstsq(X, B, rcond=None)[0]; ps = []
    for k in range(RB.shape[1]):
        prod = rx * RB[:, k]; T = np.sqrt(len(prod)) * prod.mean() / max(prod.std(), 1e-12); ps.append(2 * norm.sf(abs(T)))
    return min(1.0, min(ps) * len(ps))
def bh(p, alpha=ALPHA):
    p = np.asarray(p); n = len(p); o = np.argsort(p); thr = p[o] * n / (np.arange(n) + 1) <= alpha; k = np.max(np.where(thr)[0]) + 1 if thr.any() else 0; rej = np.zeros(n, bool); rej[o[:k]] = True; return rej
def one(rep):
    df = full_set(rep); rows = []
    for gname, pa in graphs(rep).items():
        B = basis(pa)
        for test, fn in (("linear LRT", lrt_p), ("GCM spline basis", gcm_basis_p)):
            ps = [fn(df, x, y, Z) for x, y, Z in B]; rej = bh(ps)
            for (x, y, Z), p, r in zip(B, ps, rej): rows.append(dict(rep=rep, graph=gname, test=test, x=x, y=y, p=p, rejected=bool(r), deletion_pair=(y, x) in UEDGES or (x, y) in UEDGES))
    print("rep", rep, "audit done", flush=True); return rows
if __name__ == "__main__":
    with Pool(int(os.environ.get("NPROC", 4))) as pool: res = pool.map(one, range(1, NREP + 1))
    R = pd.DataFrame([r for rr in res for r in rr]); R.to_csv(os.path.join(H, "implications_all.csv"), index=False)
    S = R.groupby(["graph", "test", "rep"]).agg(implications=("p", "size"), rejected=("rejected", "sum"), deletion_rejected=("rejected", lambda s: int(s[R.loc[s.index, "deletion_pair"]].sum()))).reset_index(); S["rate"] = S.rejected / S.implications
    M = S.groupby(["graph", "test"]).agg(reps=("rep", "nunique"), implications=("implications", "mean"), rejected=("rejected", "mean"), rejection_rate=("rate", "mean"), graph_rejected=("rejected", lambda s: (s > 0).mean()), deletion_rejected=("deletion_rejected", "mean"))
    pd.set_option("display.width", 250); print(M.round(3).to_string()); M.to_csv(os.path.join(H, "audit_summary.csv"))
