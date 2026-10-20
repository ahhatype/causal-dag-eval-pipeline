"""Engine side on the paper's data: for each (dataset, treatment, route, replicate) draw nothing, read the replicate, assign 3
stratified folds, fit g and Q with LumaWarp under the playbook profile, write tmle_inputs.csv. Resumable.
Usage: PLAN=subset|full python run_engine.py"""
import os, sys, json, subprocess, numpy as np, pandas as pd
from graph54 import NODES, PA, TOPO, P, DATA_DIR
TRAIN = os.environ.get("LUMAWARP_TRAIN", "lumawarp_training"); INFER = os.environ.get("LUMAWARP_INFER", "lumawarp_inference"); H = os.path.dirname(os.path.abspath(__file__)); K = 3; Y = "nephrolithiasis"
CFG = open(os.path.join(H, "lumawarp_playbook.cfg")).read(); POS = {v: i for i, v in enumerate(TOPO)}
ROUTES = {"parents": lambda t: sorted(PA[t]), "ordering": lambda t: [v for v in TOPO if POS[v] < POS[t] and v != Y]}
PLANS = {"subset": [("reference_subsample", "k_citrate", r, rt) for rt in ("parents", "ordering") for r in range(1, 11)] + [("astronaut_set", "k_citrate", r, rt) for rt in ("parents", "ordering") for r in range(1, 11)],
         "full": [(ds, t, r, rt) for t in ("k_citrate", "bisphosphonates") for rt in ("parents", "ordering") for ds in ("full_set", "reference_subsample", "astronaut_set") for r in range(1, 51)]}
def wcsv(df, path):
    df = df.copy()
    for c in df.columns:
        if set(np.unique(df[c])) <= {0.0, 1.0}: df[c] = df[c].astype(int)
    df.to_csv(path, index=False, float_format="%.10g")
def folds(a, y, seed):
    r = np.random.default_rng(seed); f = np.zeros(len(a), int); strata = a * 2 + y
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]; f[idx] = r.permutation(np.resize(np.arange(1, K + 1), len(idx)))
    return f
def run(exe, args, log):
    with open(log, "w") as lf: rc = subprocess.call([exe] + args, stdout=lf, stderr=subprocess.STDOUT)
    if rc != 0: raise RuntimeError(f"{exe} exited {rc}, see {log}")
def lw_fit_predict(Xtr, ytr, work, seed, Xtests):
    os.makedirs(work, exist_ok=True); open(os.path.join(work, "effective.cfg"), "w").write(CFG.replace("@TASK@", "binary").replace("@SEED@", str(seed)))
    tr = Xtr.copy(); tr["LW_TARGET"] = ytr.astype(int); wcsv(tr, os.path.join(work, "training.csv")); model = os.path.join(work, "model.lsmdl")
    run(TRAIN, [os.path.join(work, "training.csv"), model, os.path.join(work, "effective.cfg")], os.path.join(work, "training.log")); outs = []
    for tag, Xt in Xtests.items():
        fi = os.path.join(work, f"{tag}_features.csv"); fo = os.path.join(work, f"{tag}_predictions.csv"); wcsv(Xt, fi); run(INFER, [fi, model, fo], os.path.join(work, f"{tag}.log"))
        p = pd.read_csv(fo)["p_1"].to_numpy(); assert len(p) == len(Xt) and np.isfinite(p).all(); outs.append(np.clip(p, 0, 1))
    for f in os.listdir(work):
        if f.startswith("model.lsmdl") or f.startswith("axes_out") or f.endswith(".diag.jsonl") or f.endswith("_features.csv"): os.remove(os.path.join(work, f))
    return outs
def replicate(ds, A, rep, route):
    out = os.path.join(H, "reps", ds, A, route, f"rep{rep:03d}"); fin = os.path.join(out, "tmle_inputs.csv")
    if os.path.exists(fin): return
    df = pd.read_csv(os.path.join(DATA_DIR, "replicates", f"r{rep:03d}", f"{ds}.csv")); Wc = ROUTES[route](A); N = len(df)
    a = df[A].to_numpy(int); y = df[Y].to_numpy(int); seed = 5000 + rep; f = folds(a, y, seed); g = np.full(N, np.nan); q0 = np.full(N, np.nan); q1 = np.full(N, np.nan)
    for k in range(1, K + 1):
        tr, te = f != k, f == k; Wtr, Wte = df.loc[tr, Wc].reset_index(drop=True), df.loc[te, Wc].reset_index(drop=True); fd = os.path.join(out, f"fold_{k:02d}")
        g[te] = lw_fit_predict(Wtr, a[tr], os.path.join(fd, "g"), seed * 10 + 100 * k, {"heldout": Wte})[0]
        x0 = pd.concat([pd.DataFrame({A: np.zeros(te.sum(), int)}), Wte], axis=1); x1 = x0.copy(); x1[A] = 1; xtr = pd.concat([pd.DataFrame({A: a[tr]}), Wtr], axis=1)
        q0[te], q1[te] = lw_fit_predict(xtr, y[tr], os.path.join(fd, "Q"), seed * 10 + 100 * k + 1, {"counterfactual_A0": x0, "counterfactual_A1": x1})
    inputs = pd.concat([pd.DataFrame({"id": np.arange(N), "fold": f}), df[Wc].reset_index(drop=True)], axis=1); inputs[A] = a; inputs[Y] = y; inputs["g_lw"] = g; inputs["Q0_lw"] = q0; inputs["Q1_lw"] = q1
    tmp = fin + ".tmp"; inputs.to_csv(tmp, index=False, float_format="%.10g"); os.replace(tmp, fin)
if __name__ == "__main__":
    plan = PLANS[os.environ.get("PLAN", "subset")]
    for ds, A, rep, route in plan:
        try: replicate(ds, A, rep, route); print(f"{ds} {A} {route} rep{rep:03d} done", flush=True)
        except Exception as e: print(f"{ds} {A} {route} rep{rep:03d} FAILED: {e}", flush=True)
    open(os.path.join(H, "reps", f"ENGINE_DONE_{os.environ.get('PLAN', 'subset')}"), "w").write("done\n"); print("engine side finished")
