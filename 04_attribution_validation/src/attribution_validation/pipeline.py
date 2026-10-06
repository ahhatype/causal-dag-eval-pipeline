from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import time
import traceback
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from . import data as D
from .evaluation import (auc, kendall_tau_b, non_ancestor_share, proximity_bias_index, tie_adjusted_truth,
                         top_k_recovery)
from .interventional import fit_scm, interventional_shap
from .ng import ng_causal_shap
from .ordering import ordering_only_shap
from .standard import permutation_shap
from .superlearner import SuperLearner, default_library


def run_seed(params: dict, replicate: int, dataset: str) -> int:
    """Shared by both feature sets, so they use the same split, explained records and background."""
    return params["seed"] + 1000 * replicate + 10 * params["datasets"].index(dataset)


def choose_explained(test: pd.DataFrame, n: int, min_events: int, rng: np.random.Generator) -> pd.DataFrame:
    idx = rng.choice(len(test), size=min(n, len(test)), replace=False)
    pick = test.iloc[idx]
    events = test.index[test[D.OUTCOME] == 1].difference(pick.index)
    short = min_events - int(pick[D.OUTCOME].sum())
    if short > 0 and len(events):
        add = rng.choice(events, size=min(short, len(events)), replace=False)
        drop = rng.choice(pick.index[pick[D.OUTCOME] == 0], size=len(add), replace=False)
        pick = pd.concat([pick.drop(drop), test.loc[add]])
    return pick


def run_one(df: pd.DataFrame, cfg: D.GenerationConfig, ag: D.AnalysisGraph, features: list[str],
            truth: pd.Series, truth_se: pd.Series, ancestors: set[str], distance: dict, params: dict,
            seed: int, methods: list[str]) -> tuple[dict, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    train, test = train_test_split(df, test_size=params["split"]["test_fraction"], stratify=df[D.OUTCOME],
                                   random_state=seed % (2**32))
    sl = SuperLearner(default_library(seed), params["super_learner"]["folds"], seed)
    t0 = time.time()
    sl.fit(train[features].to_numpy(dtype=float), train[D.OUTCOME].to_numpy())
    info = {"seed": seed, "n_features": len(features), "fit_seconds": time.time() - t0,
            "auc_model": auc(test[D.OUTCOME].to_numpy(), sl.predict(test[features].to_numpy(dtype=float))),
            "auc_true_probability": auc(test[D.OUTCOME].to_numpy(), D.true_outcome_probability(test, cfg)),
            **{f"weight_{k}": v for k, v in sl.weights_.items()},
            **{f"cv_log_loss_{k}": v for k, v in sl.cv_risk_.items()}}
    ex = params["explain"]
    explain = choose_explained(test, ex["records"], ex["min_events"], rng)
    background = train.iloc[rng.choice(len(train), size=min(ex["background"], len(train)), replace=False)]
    expected = float(np.mean(sl.predict(background[features].to_numpy(dtype=float))))
    info["explained_events"] = int(explain[D.OUTCOME].sum())

    ranked_truth = tie_adjusted_truth(truth.to_dict(), truth_se.to_dict())
    runners = {
        "standard_shap": lambda: permutation_shap(sl.predict, explain, background, features,
                                                  params["standard_shap"]["permutations"], seed),
        "ordering_only_shap": lambda: ordering_only_shap(sl.predict, explain, background, features, ag,
                                                         params["ordering_only_shap"]["orders"],
                                                         params["ordering_only_shap"]["burn_in"],
                                                         params["ordering_only_shap"]["thin"], seed),
        "ng_causal_shap": lambda: ng_causal_shap(sl.predict, train, explain, features, ag,
                                                 params["ng_causal_shap"]["samples"],
                                                 params["ng_causal_shap"]["orders"], seed, expected),
        "interventional_shap": lambda: interventional_shap(sl.predict, fit_scm(train, ag), explain, background,
                                                           features, ag, params["interventional_shap"]["orders"],
                                                           params["interventional_shap"]["burn_in"],
                                                           params["interventional_shap"]["thin"], seed),
    }
    rows, imps, vals = [], [], []
    for m in methods:
        t0 = time.time()
        res = runners[m]()
        secs = time.time() - t0
        imp = res.values.abs().mean().to_dict()
        anc = {f: v for f, v in imp.items() if f in ancestors}
        rows.append({"method": m, "seconds": secs,
                     "max_efficiency_error": float(np.max(np.abs(res.efficiency_error))),
                     "kendall_tau_b": kendall_tau_b(anc, ranked_truth),
                     "top_k_recovery": top_k_recovery(anc, ranked_truth, params["evaluation"]["top_k"]),
                     "proximity_bias_index": proximity_bias_index(anc, truth.to_dict(), distance),
                     "non_ancestor_share": non_ancestor_share(imp, ancestors)})
        if res.diagnostics is not None:
            factor = res.diagnostics["rescale_factor"].abs()
            rows[-1].update(rescale_factor_median_abs=float(factor.median()), rescale_factor_max_abs=float(factor.max()))
        imps.append(pd.DataFrame({"method": m, "feature": list(imp), "importance": list(imp.values())}))
        v = res.values.copy()
        if res.diagnostics is not None:
            v = v.join(res.diagnostics)
        vals.append(v.assign(method=m, outcome=explain[D.OUTCOME].to_numpy()).rename_axis("record").reset_index())
    return {"info": info, "rows": rows, "values": pd.concat(vals, ignore_index=True)}, pd.concat(imps, ignore_index=True)


_CONTEXT: dict = {}


def _init_worker(params: dict, data_dir: Path, methods: list[str]) -> None:
    """Every library runs single-threaded, so floating-point results do not depend on `workers`."""
    from threadpoolctl import threadpool_limits
    threadpool_limits(1)
    cfg = D.load_generation_config()
    ag = D.analysis_graph(cfg)
    fsets = D.feature_sets(data_dir)
    _CONTEXT.update(params=params, data_dir=data_dir, methods=methods, cfg=cfg, ag=ag, fsets=fsets,
                    distance=D.distance_to_outcome(ag))


def _cell_name(r: int, ds: str, fs: str) -> str:
    return f"r{r:03d}_{ds}_{fs}"


def _write(df: pd.DataFrame, path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, index=False, compression="gzip" if path.suffix == ".gz" else None)
    os.replace(tmp, path)


def run_cell(cell: tuple[int, str, str], cell_dir: Path) -> tuple[str, str | None]:
    """One replicate, data set and feature set. `_runs.csv` is written last and marks the cell complete."""
    r, ds, fs = cell
    name = _cell_name(r, ds, fs)
    try:
        c = _CONTEXT
        fsets = c["fsets"]
        target = D.GROUND_TRUTH_TARGET[ds]
        truth = D.ground_truth(c["data_dir"], target).loc[fsets["ancestor"]]
        truth_se = D.ground_truth_se(c["data_dir"], target).loc[fsets["ancestor"]]
        df = D.load_replicate(c["data_dir"], r, ds)
        res, imp = run_one(df, c["cfg"], c["ag"], fsets[fs], truth, truth_se, set(fsets["ancestor"]), c["distance"],
                           c["params"], run_seed(c["params"], r, ds), c["methods"])
        key = {"replicate": r, "dataset": ds, "feature_set": fs}
        _write(res["values"].assign(**key), cell_dir / f"{name}_values.csv.gz")
        _write(imp.assign(**key), cell_dir / f"{name}_importance.csv")
        _write(pd.DataFrame([{**key, **res["info"], **row} for row in res["rows"]]), cell_dir / f"{name}_runs.csv")
        return name, None
    except Exception:
        err = traceback.format_exc()
        (cell_dir / f"{name}_error.txt").write_text(err)
        return name, err


def code_hash() -> str:
    h = hashlib.sha256()
    for f in sorted(Path(__file__).parent.glob("*.py")):
        h.update(f.name.encode() + f.read_bytes())
    return h.hexdigest()[:16]


def manifest(params: dict, methods: list[str], cfg: D.GenerationConfig) -> dict:
    """What a cell's results depend on, apart from which cell it is."""
    return {"params": {k: v for k, v in params.items() if k not in ("workers", "methods", "datasets", "feature_sets")},
            "datasets": params["datasets"], "feature_sets": params["feature_sets"], "methods": methods,
            "data_config_hash": cfg.calibration["config_hash"], "code_hash": code_hash()}


def _prepare_cell_dir(cell_dir: Path, want: dict, resume: bool, fresh: bool) -> None:
    mf = cell_dir / "manifest.json"
    existing = list(cell_dir.glob("*.csv*")) if cell_dir.exists() else []
    if fresh:
        for f in existing + list(cell_dir.glob("*_error.txt")):
            f.unlink()
    elif existing:
        if not resume:
            raise RuntimeError(f"{cell_dir} already holds results: pass --resume to continue them or --fresh to discard them")
        if not mf.exists() or json.loads(mf.read_text()) != want:
            raise RuntimeError(f"{cell_dir} was written with different settings or code: pass --fresh to discard it")
    cell_dir.mkdir(parents=True, exist_ok=True)
    for f in cell_dir.glob("*_error.txt"):
        f.unlink()
    mf.write_text(json.dumps(want, indent=1, sort_keys=True))


def write_provenance(out_dir: Path, want: dict, n_cells: int, workers: int) -> None:
    repo = Path(__file__).resolve().parents[3]

    def git(*args: str) -> str:
        try:
            return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True).stdout.strip()
        except OSError:
            return "unavailable"

    pkgs = ["numpy", "pandas", "scikit-learn", "scipy", "shap", "networkx"]
    lines = [f"generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
             f"git commit: {git('rev-parse', 'HEAD')}" + (" (uncommitted changes)" if git("status", "--porcelain") else ""),
             f"code hash: {want['code_hash']}", f"data config hash: {want['data_config_hash']}",
             f"cells: {n_cells}", f"workers: {workers}", f"python: {platform.python_version()}",
             *[f"{p}: {version(p)}" for p in pkgs], "settings:", json.dumps(want["params"], indent=1, sort_keys=True)]
    (out_dir / "provenance.txt").write_text("\n".join(lines) + "\n")


def run(params: dict, data_dir: Path, out_dir: Path, replicates: int | None = None,
        methods: list[str] | None = None, datasets: list[str] | None = None,
        feature_sets: list[str] | None = None, workers: int | None = None, resume: bool = False,
        fresh: bool = False) -> pd.DataFrame:
    """Cells are independent and seeded individually, so results do not depend on `workers`.

    Existing cells in `out_dir/cells` are kept only with `resume`, and only if settings and code are unchanged;
    `fresh` discards them. A failed cell is logged and the rest continue; outputs are combined only when every cell
    is complete.
    """
    cfg = D.load_generation_config()
    available = D.check_data_current(data_dir, cfg)
    reps = range(1, min(available, replicates or available) + 1)
    ag = D.analysis_graph(cfg)
    if (len(ag.nodes), len(ag.edges)) != (36, 59):
        raise RuntimeError(f"analysis graph has {len(ag.nodes)} nodes and {len(ag.edges)} edges, expected 36 and 59")
    methods = methods or params["methods"]
    workers = workers or params.get("workers", 1)
    cell_dir = out_dir / "cells"
    want = manifest(params, methods, cfg)
    _prepare_cell_dir(cell_dir, want, resume, fresh)
    cells = [(r, ds, fs) for r in reps for ds in datasets or params["datasets"] for fs in feature_sets or params["feature_sets"]]
    todo = [c for c in cells if not (cell_dir / f"{_cell_name(*c)}_runs.csv").exists()]
    print(f"{len(cells)} cells, {len(cells) - len(todo)} already done, {workers} worker(s)", flush=True)
    start = time.time()
    failed = []

    def report(i: int, name: str, err: str | None) -> None:
        status = "FAILED, see cells/" + name + "_error.txt" if err else ""
        print(f"[{i}/{len(todo)}] {name} ({(time.time() - start) / 60:.0f} min) {status}".rstrip(), flush=True)
        if err:
            failed.append(name)

    if workers == 1:
        _init_worker(params, data_dir, methods)
        for i, c in enumerate(todo, 1):
            report(i, *run_cell(c, cell_dir))
    else:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(workers, mp_context=mp.get_context("spawn"), initializer=_init_worker,
                                 initargs=(params, data_dir, methods)) as pool:
            futures = [pool.submit(run_cell, c, cell_dir) for c in todo]
            for i, fut in enumerate(as_completed(futures), 1):
                report(i, *fut.result())
    if failed:
        raise RuntimeError(f"{len(failed)} cell(s) failed ({', '.join(failed[:5])}...): fix and re-run with --resume")

    names = [_cell_name(*c) for c in cells]
    per_run = pd.concat([pd.read_csv(cell_dir / f"{n}_runs.csv") for n in names], ignore_index=True)
    _write(per_run, out_dir / "per_run.csv")
    _write(pd.concat([pd.read_csv(cell_dir / f"{n}_importance.csv") for n in names], ignore_index=True),
           out_dir / "importance.csv")
    _write(pd.concat([pd.read_csv(cell_dir / f"{n}_values.csv.gz") for n in names], ignore_index=True),
           out_dir / "values.csv.gz")
    summary = summarize(per_run)
    _write(summary, out_dir / "summary.csv")
    write_provenance(out_dir, want, len(cells), workers)
    return summary


def summarize(per_run: pd.DataFrame) -> pd.DataFrame:
    metrics = ["kendall_tau_b", "top_k_recovery", "proximity_bias_index", "non_ancestor_share", "auc_model",
               "auc_true_probability", "max_efficiency_error", "seconds"]
    g = per_run.groupby(["dataset", "feature_set", "method"])[metrics]
    mean = g.mean().add_suffix("_mean")
    se = (g.std() / np.sqrt(g.count())).add_suffix("_se")
    return pd.concat([g.size().rename("runs"), mean, se], axis=1).reset_index()
