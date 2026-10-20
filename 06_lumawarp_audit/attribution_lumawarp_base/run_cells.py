"""The paper's attribution arm with LumaWarp as the explained model. Uses the pipeline package for split, records,
background, orders, SCM, truth and scoring; only the value-function evaluation is batched through the engine CLI.
Env: REPS (10), DATASETS (comma list), METHODS (comma list), BACKGROUND (128), ORDERS (64), NCHUNK (2000000), NPAR (2)."""
import os, sys, time, subprocess, numpy as np, pandas as pd, yaml
from pathlib import Path
H = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(PIPE, "04_attribution_validation", "src"))
from attribution_validation import data as D, pipeline as P
from attribution_validation.interventional import fit_scm, OrderSampler, ancestral_pairs, AttributionResult
from attribution_validation.evaluation import auc
from threadpoolctl import threadpool_limits
threadpool_limits(1)   # as the pipeline's workers: the Super Learner's forests otherwise oversubscribe the cores
TRAIN = os.environ.get("LUMAWARP_TRAIN", "lumawarp_training")
INFER = os.environ.get("LUMAWARP_INFER", "lumawarp_inference")
CFG = os.path.join(H, "cfg_attribution.cfg"); DATA = Path(PIPE) / os.environ.get("DATA_DIR", "02_data").removeprefix("./"); CELLS = os.path.join(H, "cells")
REPS = int(os.environ.get("REPS", 10)); NCHUNK = int(os.environ.get("NCHUNK", 2_000_000)); NPAR = int(os.environ.get("NPAR", 2))
params = yaml.safe_load(open(os.path.join(PIPE, "04_attribution_validation", "config", "params.yaml")))
params["explain"]["background"] = int(os.environ.get("BACKGROUND", 128)); ORDERS = int(os.environ.get("ORDERS", 64))
params["interventional_shap"]["orders"] = params["ordering_only_shap"]["orders"] = params["standard_shap"]["permutations"] = ORDERS
DATASETS = os.environ.get("DATASETS", ",".join(params["datasets"])).split(",")
METHODS = os.environ.get("METHODS", "interventional_shap,standard_shap,ordering_only_shap").split(",")
cfg = D.load_generation_config(); ag = D.analysis_graph(cfg); fsets = D.feature_sets(DATA)
ctx = dict(cfg=cfg, ag=ag, fsets=fsets, distance=D.distance_to_outcome(ag), truth=P.truth_context(D.truth_values(DATA), fsets["ancestor"], params))

class LumaWarp:
    """Fit once per cell; predict any matrix by writing it to disk and calling the engine (NPAR chunks in parallel)."""
    def __init__(self, d, train, features):
        self.d, self.features, self.calls, self.rows, self.seconds = d, list(features), 0, 0, 0.0
        os.makedirs(d, exist_ok=True); train[self.features + [D.OUTCOME]].to_csv(os.path.join(d, "train.csv"), index=False)
        if not os.path.exists(os.path.join(d, "model.lsmdl")):
            t0 = time.time()
            with open(os.path.join(d, "train_log.txt"), "w") as log:
                subprocess.run([TRAIN, "train.csv", "model.lsmdl", CFG], cwd=d, stdout=log, stderr=subprocess.STDOUT, check=True)
            self.fit_seconds = time.time() - t0
            for f in ("model.lsmdl.diag.jsonl", "model.lsmdl.ckpt"): _rm(os.path.join(d, f))
        else: self.fit_seconds = float("nan")
    def predict(self, X):
        X = np.asarray(X, dtype=np.float32); out = np.empty(len(X)); t0 = time.time()
        starts = list(range(0, len(X), NCHUNK))
        for g in range(0, len(starts), NPAR):
            procs = []
            for j, s in enumerate(starts[g:g + NPAR]):
                wd = os.path.join(self.d, f"infer_{j}"); os.makedirs(wd, exist_ok=True)
                pd.DataFrame(X[s:s + NCHUNK], columns=self.features).to_csv(os.path.join(wd, "score.csv"), index=False, float_format="%.6g")
                procs.append((s, wd, subprocess.Popen([INFER, "score.csv", os.path.join("..", "model.lsmdl"), "pred.csv"], cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)))
            for s, wd, p in procs:
                if p.wait() != 0: raise RuntimeError(f"inference failed in {wd}")
                out[s:s + NCHUNK] = pd.read_csv(os.path.join(wd, "pred.csv"), usecols=["p_1"]).p_1.to_numpy()
                for f in ("score.csv", "pred.csv", "pred.csv.diag.jsonl"): _rm(os.path.join(wd, f))
        self.calls += 1; self.rows += len(X); self.seconds += time.time() - t0; return out

def _rm(p):
    try: os.remove(p)
    except FileNotFoundError: pass

def _aggregate(p, order_list, n_eval, n_bg, k, idx, baseline):
    out = np.zeros((n_eval, k)); blk = n_eval * n_bg; b = 0
    for order in order_list:
        prev = np.full(n_eval, baseline)
        for f in order:
            cur = p[b:b + blk].reshape(n_eval, n_bg).mean(axis=1); b += blk; out[:, idx[f]] += cur - prev; prev = cur
    return out / len(order_list)

def interventional_batched(predict, scm, explain, background, features, orders, burn_in, thin, seed):
    """pipeline.interventional.interventional_shap with every value-function row generated first (same rng sequence)."""
    features = list(features); rng = np.random.default_rng(seed); u = scm.recover_exogenous(background, rng)
    n_eval, n_bg, k = len(explain), len(background), len(features); u_tiled = {a: np.tile(v, n_eval) for a, v in u.items()}
    base_vals = scm.simulate(u); base_X = np.column_stack([base_vals[f] for f in features]).astype(np.float32)
    xv = {f: np.repeat(explain[f].to_numpy(dtype=float), n_bg) for f in features}
    sampler = OrderSampler(features, ancestral_pairs(ag, features), rng, burn_in, thin); blocks, order_list = [], []
    for _ in range(orders):
        order = sampler.sample(); order_list.append(order); fixed = {}
        for f in order:
            fixed[f] = xv[f]; sim = scm.simulate(u_tiled, fixed); blocks.append(np.column_stack([sim[g] for g in features]).astype(np.float32))
    E = explain[features].to_numpy(dtype=np.float32); p = predict(np.concatenate([base_X, E] + blocks)); del blocks
    baseline = float(p[:n_bg].mean()); p_explain = p[n_bg:n_bg + n_eval]; p = p[n_bg + n_eval:]
    out = _aggregate(p, order_list, n_eval, n_bg, k, {f: i for i, f in enumerate(features)}, baseline)
    return AttributionResult(pd.DataFrame(out, columns=features, index=explain.index), baseline, out.sum(axis=1) + baseline - p_explain)

def ordering_only_batched(predict, explain, background, features, orders, burn_in, thin, seed):
    """pipeline.ordering.ordering_only_shap, batched (same rng sequence, same orders)."""
    features = list(features); rng = np.random.default_rng(seed); n_eval, n_bg, k = len(explain), len(background), len(features)
    bg = background[features].to_numpy(dtype=float); xv = {f: np.repeat(explain[f].to_numpy(dtype=float), n_bg) for f in features}
    col = {f: j for j, f in enumerate(features)}; sampler = OrderSampler(features, ancestral_pairs(ag, features), rng, burn_in, thin)
    blocks, order_list = [], []
    for _ in range(orders):
        batch = np.tile(bg, (n_eval, 1)); order = sampler.sample(); order_list.append(order)
        for f in order: batch[:, col[f]] = xv[f]; blocks.append(batch.astype(np.float32))
    E = explain[features].to_numpy(dtype=np.float32); p = predict(np.concatenate([bg.astype(np.float32), E] + blocks)); del blocks
    baseline = float(p[:n_bg].mean()); p_explain = p[n_bg:n_bg + n_eval]; p = p[n_bg + n_eval:]
    out = _aggregate(p, order_list, n_eval, n_bg, k, col, baseline)
    return AttributionResult(pd.DataFrame(out, columns=features, index=explain.index), baseline, out.sum(axis=1) + baseline - p_explain)

def standard_batched(predict, explain, background, features, permutations, seed):
    """Permutation SHAP with an independent background masker, antithetic order pairs, same budget as the pipeline's
    shap.PermutationExplainer (which cannot be batched); equal in expectation, not bit-identical."""
    features = list(features); rng = np.random.default_rng(seed); k = len(features); n_eval, n_bg = len(explain), len(background)
    bg = background[features].to_numpy(dtype=float); E = explain[features].to_numpy(dtype=float); orders = []
    for _ in range(permutations // 2):
        p = rng.permutation(k); orders += [list(p), list(p[::-1])]
    blocks = []
    for order in orders:
        batch = np.tile(bg, (n_eval, 1))
        for f in order: batch[:, f] = np.repeat(E[:, f], n_bg); blocks.append(batch.astype(np.float32))
    p = predict(np.concatenate([bg.astype(np.float32), E.astype(np.float32)] + blocks)); del blocks
    baseline = float(p[:n_bg].mean()); p_explain = p[n_bg:n_bg + n_eval]; p = p[n_bg + n_eval:]
    out = _aggregate(p, orders, n_eval, n_bg, k, {i: i for i in range(k)}, baseline)
    return AttributionResult(pd.DataFrame(out, columns=features, index=explain.index), baseline, out.sum(axis=1) + baseline - p_explain)

def run_cell(r, ds):
    name = f"r{r:03d}_{ds}_ancestor"; done = os.path.join(CELLS, name + "_runs.csv")
    if os.path.exists(done): return "done"
    features = fsets["ancestor"]; df = D.load_replicate(DATA, r, ds); seed = P.run_seed(params, r, ds)
    pr = P.prepare_run(df, ctx, features, params, seed)          # the paper's split, Super Learner, records, background
    train, explain, background, info = pr["train"], pr["explain"], pr["background"], pr["info"]
    test = df.drop(train.index)
    lw = LumaWarp(os.path.join(CELLS, name), train, features)
    info.update(lw_fit_seconds=lw.fit_seconds, auc_lumawarp=auc(test[D.OUTCOME].to_numpy(), lw.predict(test[features].to_numpy(dtype=float))),
                background_size=len(background), orders=ORDERS, auc_superlearner=info.pop("auc_model"))
    p = params; rows, imps, vals = [], [], []
    for m in METHODS:
        t0 = time.time(); lw.calls = lw.rows = 0; lw.seconds = 0.0
        if m == "interventional_shap": res = interventional_batched(lw.predict, fit_scm(train, ag), explain, background, features, p[m]["orders"], p[m]["burn_in"], p[m]["thin"], seed)
        elif m == "ordering_only_shap": res = ordering_only_batched(lw.predict, explain, background, features, p[m]["orders"], p[m]["burn_in"], p[m]["thin"], seed)
        elif m == "standard_shap": res = standard_batched(lw.predict, explain, background, features, p[m]["permutations"], seed)
        else: raise ValueError(m)
        imp = res.values.abs().mean().to_dict(); label = m + "_lumawarp"
        minfo = dict(seconds=time.time() - t0, engine_rows=lw.rows, engine_seconds=lw.seconds, max_efficiency_error=float(np.max(np.abs(res.efficiency_error))))
        rows += [{"method": label, **minfo, **sc} for sc in P.score(imp, ctx["truth"], set(features), ctx["distance"], params)]
        imps.append(pd.DataFrame({"method": label, "feature": list(imp), "importance": list(imp.values())}))
        vals.append(res.values.assign(method=label, outcome=explain[D.OUTCOME].to_numpy()).rename_axis("record").reset_index())
        print(f"    {m}: {minfo['engine_rows']:,} rows in {minfo['engine_seconds']/60:.1f} min", flush=True)
    key = {"replicate": r, "dataset": ds, "feature_set": "ancestor"}
    pd.concat(vals, ignore_index=True).assign(**key).to_csv(os.path.join(CELLS, name + "_values.csv.gz"), index=False)
    pd.concat(imps, ignore_index=True).assign(**key).to_csv(os.path.join(CELLS, name + "_importance.csv"), index=False)
    pd.DataFrame([{**key, **info, **row} for row in rows]).to_csv(done + ".tmp", index=False); os.replace(done + ".tmp", done)
    return "ran"

if __name__ == "__main__":
    os.makedirs(CELLS, exist_ok=True); cells = [(r, ds) for r in range(1, REPS + 1) for ds in DATASETS]; T0 = time.time()
    for i, (r, ds) in enumerate(cells, 1):
        t0 = time.time(); status = run_cell(r, ds)
        print(f"[{i}/{len(cells)}] r{r:03d} {ds} {status} ({(time.time()-t0)/60:.1f} min, total {(time.time()-T0)/60:.0f} min)", flush=True)
    print("ALL CELLS COMPLETED")
