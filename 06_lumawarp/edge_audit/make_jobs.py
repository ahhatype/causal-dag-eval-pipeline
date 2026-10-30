import os, numpy as np, pandas as pd
from graph54 import NODES, PA, TOPO, P, DATA_DIR
H = os.path.dirname(os.path.abspath(__file__)); POS = {v: i for i, v in enumerate(TOPO)}; REPS = range(1, 11); DS = ["reference_subsample", "astronaut_set", "full_set"]
def w(df, path):
    df = df.copy()
    for c in df.columns:
        if set(np.unique(df[c])) <= {0.0, 1.0}: df[c] = df[c].astype(int)
    df.to_csv(path, index=False, float_format="%.8g")
jobs = []
for ds in DS:
    for rep in REPS:
        df = pd.read_csv(os.path.join(DATA_DIR, "replicates", f"r{rep:03d}", f"{ds}.csv"))[NODES]; n = len(df); r = np.random.default_rng(99000 + rep); nulls = []
        for k in range(6): df[f"null_cont_{k}"] = r.standard_normal(n); nulls.append(f"null_cont_{k}")
        for k in range(6): df[f"null_bin_{k}"] = (r.random(n) < 0.35).astype(float); nulls.append(f"null_bin_{k}")
        for t in TOPO:
            feats = [v for v in TOPO if POS[v] < POS[t]]
            if not feats or df[t].nunique() < 2: continue
            feats = feats + nulls; jd = os.path.join(H, ds, f"r{rep:03d}", t); os.makedirs(jd, exist_ok=True)
            w(df[feats + [t]], os.path.join(jd, "train.csv")); med = df[feats].median(); blocks = []
            for f in feats:
                grid = np.quantile(df[f], np.linspace(.01, .99, 21)) if df[f].nunique() > 2 else np.array([0., 1.])
                sw = pd.DataFrame(np.tile(med.to_numpy(), (len(grid), 1)), columns=feats); sw[f] = grid; sw["block"] = f; blocks.append(sw)
            sc = pd.concat(blocks, ignore_index=True); w(sc[feats], os.path.join(jd, "score.csv")); sc[["block"]].to_csv(os.path.join(jd, "score_meta.csv"), index=False)
            jobs.append(dict(dataset=ds, rep=rep, target=t, binary=int(df[t].nunique() == 2), n_features=len(feats), dir=os.path.relpath(jd, H)))
J = pd.DataFrame(jobs); J.to_csv(os.path.join(H, "jobs.csv"), index=False); print(J.groupby("dataset").size().to_string())
