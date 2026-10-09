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
from .evaluation import (auc, kendall_tau_b, l1_distance, margin_tie_groups, mcse_tie_groups, non_ancestor_share,
                         proximity_bias_index, tie_adjusted, top_k_recovery, weighted_tau)
from .interventional import fit_scm, interventional_shap
from .ng import causal_weights, edge_strengths, fitted_ng_inputs, ng_style_shap
from .ordering import ordering_only_shap
from .standard import permutation_shap
from .superlearner import SuperLearner, default_library
from .truesim import TrueSampler, TrueSCM

SCOPES = ("with_era", "without_era")
CONTRASTS = {"order_effect": ("ordering_only_shap", "standard_shap"),
             "value_function_effect": ("interventional_shap", "ordering_only_shap")}


def run_seed(params: dict, replicate: int, dataset: str) -> int:
    """Shared by both feature sets, so they use the same split, explained records and background."""
    return params["seed"] + 1000 * replicate + 10 * params["datasets"].index(dataset)


def choose_explained(test: pd.DataFrame, n: int, rng: np.random.Generator) -> pd.DataFrame:
    return test.iloc[rng.choice(len(test), size=min(n, len(test)), replace=False)]


def rule_names(params: dict) -> list[str]:
    """Tie rules; the first is primary."""
    return [f"delta{round(100 * d)}" for d in params["evaluation"]["tie_margins"]] + ["mcse"]


def truth_context(truth: pd.DataFrame, ancestors: list[str], params: dict) -> dict:
    """Per (truth_type, population, scope): |truth| over the scored ancestors, its tie groups under each rule and
    the tie-adjusted values. The without-era scope drops mission era."""
    ev = params["evaluation"]
    out = {}
    for tt in D.TRUTH_TYPES:
        for pop in D.POPULATIONS:
            v = D.truth_series(truth, tt, pop)
            m = D.truth_series(truth, tt, pop, "mcse")
            for scope in SCOPES:
                fs = [f for f in ancestors if scope == "with_era" or f != "mission_era"]
                vals = {f: abs(float(v[f])) for f in fs}
                top = max(vals.values())
                groups = {name: margin_tie_groups(vals, d * top) for name, d in zip(rule_names(params), ev["tie_margins"])}
                groups["mcse"] = mcse_tie_groups(vals, {f: float(m[f]) for f in fs}, ev["tie_z"])
                out[(tt, pop, scope)] = {"raw": vals, "groups": groups,
                                         "tied": {k: tie_adjusted(vals, g) for k, g in groups.items()}}
    return out


def tie_group_table(tc: dict) -> pd.DataFrame:
    rows = []
    for (tt, pop, scope), t in tc.items():
        for rule, groups in t["groups"].items():
            for i, g in enumerate(groups, 1):
                rows.append({"truth_type": tt, "population": pop, "scope": scope, "rule": rule, "group": i,
                             "size": len(g), "features": "; ".join(g),
                             "max_abs_truth": max(t["raw"][f] for f in g), "min_abs_truth": min(t["raw"][f] for f in g)})
    return pd.DataFrame(rows)


def score(imp: dict[str, float], tc: dict, ancestors: set[str], distance: dict, params: dict) -> list[dict]:
    """Every metric against every truth; rank metrics use tie-adjusted truth, the others raw |truth|."""
    primary, *others = rule_names(params)
    rows = []
    for (tt, pop, scope), t in tc.items():
        anc = {f: imp[f] for f in t["raw"]}
        rows.append({"truth_type": tt, "population": pop, "scope": scope,
                     "kendall_tau_b": kendall_tau_b(anc, t["tied"][primary]),
                     **{f"kendall_tau_b_{k}": kendall_tau_b(anc, t["tied"][k]) for k in others},
                     "weighted_tau": weighted_tau(anc, t["raw"]),
                     "l1_distance": l1_distance(anc, t["raw"]),
                     "top_k_recovery": top_k_recovery(anc, t["tied"][primary], params["evaluation"]["top_k"]),
                     "proximity_bias_index": proximity_bias_index(anc, t["raw"], distance),
                     "non_ancestor_share": non_ancestor_share(imp, ancestors)})
    return rows


def ng_importance(res, threshold: float) -> tuple[dict, dict, float]:
    """Mean |phi| over records whose |pre-rescaling total| reaches the threshold, over all records, and the share
    excluded."""
    keep = res.diagnostics["weighted_total"].abs().to_numpy() >= threshold
    absv = res.values.abs()
    kept = absv[keep].mean() if keep.any() else absv.mean() * np.nan
    return kept.to_dict(), absv.mean().to_dict(), float(1 - keep.mean())


def check_oracle_inputs(cfg: D.GenerationConfig, ancestors: list[str]) -> None:
    """Every parent of the outcome must be an observed member of the ancestor set."""
    pa = cfg.edges.loc[cfg.edges["child"] == D.OUTCOME, "parent"].tolist()
    observed = set(cfg.nodes.loc[cfg.nodes["observed"], "id"])
    bad = [p for p in pa if p not in ancestors or p not in observed]
    if bad:
        raise RuntimeError(f"oracle predictor needs outcome parents {bad} in the ancestor set")


def oracle_context(cfg: D.GenerationConfig, ag: D.AnalysisGraph, ancestors: list[str], params: dict) -> dict:
    o = params["oracle"]
    scm = TrueSCM(cfg)
    pool = scm.population(o["pool"], np.random.default_rng(o["seed"]))
    e = cfg.edges[cfg.edges["child"] == D.OUTCOME]
    idx = [ancestors.index(p) for p in e["parent"]]
    coef = e["coef"].astype(float).to_numpy()
    b0 = float(cfg.calibration["logit"][D.OUTCOME])

    def predict(X: np.ndarray) -> np.ndarray:
        return 1.0 / (1.0 + np.exp(-(b0 + np.asarray(X, dtype=float)[:, idx] @ coef)))

    return {"scm": scm, "pool": pool, "predict": predict,
            "weights": causal_weights(ag, edge_strengths(pool, ag), ancestors),
            "expected": float(np.mean(predict(pool[ancestors].to_numpy(dtype=float))))}


def method_runners(predict, explain, background, train, features, ag, params, seed, expected, scm, ng_inputs):
    p = params
    return {
        "standard_shap": lambda: permutation_shap(predict, explain, background, features,
                                                  p["standard_shap"]["permutations"], seed),
        "ordering_only_shap": lambda: ordering_only_shap(predict, explain, background, features, ag,
                                                         p["ordering_only_shap"]["orders"],
                                                         p["ordering_only_shap"]["burn_in"],
                                                         p["ordering_only_shap"]["thin"], seed),
        "ng_style_shap": lambda: ng_style_shap(predict, explain, features, *ng_inputs(),
                                               p["ng_style_shap"]["samples"], p["ng_style_shap"]["orders"], seed,
                                               expected),
        "interventional_shap": lambda: interventional_shap(predict, scm(), explain, background, features, ag,
                                                           p["interventional_shap"]["orders"],
                                                           p["interventional_shap"]["burn_in"],
                                                           p["interventional_shap"]["thin"], seed),
    }


def importance_of(m: str, res, params: dict) -> tuple[dict, dict]:
    """Primary importance plus method-level diagnostics."""
    info = {"max_efficiency_error": float(np.max(np.abs(res.efficiency_error)))}
    if m != "ng_style_shap":
        return res.values.abs().mean().to_dict(), info
    ng = params["ng_style_shap"]
    imp, imp_all, excluded = ng_importance(res, ng["exclude_below"])
    factor = res.diagnostics["rescale_factor"].abs()
    info.update(ng_excluded_share=excluded, rescale_factor_median_abs=float(factor.median()),
                rescale_factor_max_abs=float(factor.max()), rescale_flag_share=float((factor > ng["rescale_flag"]).mean()),
                importance_all_records=imp_all)
    return imp, info


def prepare_run(df: pd.DataFrame, ctx: dict, features: list[str], params: dict, seed: int) -> dict:
    """Split, super learner fit, explained records and background: shared by every method in a run."""
    cfg = ctx["cfg"]
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
    explain = choose_explained(test, ex["records"], rng)
    background = train.iloc[rng.choice(len(train), size=min(ex["background"], len(train)), replace=False)]
    info["explained_events"] = int(explain[D.OUTCOME].sum())
    info["explained_events_low"] = info["explained_events"] < ex["flag_min_events"]
    return {"sl": sl, "train": train, "explain": explain, "background": background, "info": info,
            "expected": float(np.mean(sl.predict(background[features].to_numpy(dtype=float))))}


def run_one(df: pd.DataFrame, ctx: dict, features: list[str], params: dict, seed: int, methods: list[str],
            oracle: dict | None = None) -> tuple[dict, pd.DataFrame]:
    ag, tc, ancestors = ctx["ag"], ctx["truth"], set(ctx["fsets"]["ancestor"])
    pr = prepare_run(df, ctx, features, params, seed)
    sl, train, explain, background, expected, info = (pr[k] for k in ("sl", "train", "explain", "background",
                                                                         "expected", "info"))
    runners = method_runners(sl.predict, explain, background, train, features, ag, params, seed, expected,
                             lambda: fit_scm(train, ag), lambda: fitted_ng_inputs(train, ag, features))
    rows, imps, vals, fitted = [], [], [], {}
    for m in methods:
        t0 = time.time()
        res = runners[m]()
        secs = time.time() - t0
        imp, minfo = importance_of(m, res, params)
        fitted[m] = imp
        extra = {m + "_all_records": minfo.pop("importance_all_records")} if m == "ng_style_shap" else {}
        for name, im in {m: imp, **extra}.items():
            rows += [{"method": name, "seconds": secs, **minfo, **r} for r in score(im, tc, ancestors, ctx["distance"], params)]
            imps.append(pd.DataFrame({"method": name, "feature": list(im), "importance": list(im.values())}))
        v = res.values.copy()
        if res.diagnostics is not None:
            v = v.join(res.diagnostics)
        vals.append(v.assign(method=m, outcome=explain[D.OUTCOME].to_numpy()).rename_axis("record").reset_index())

    oracle_rows = []
    if oracle is not None:
        o = params["oracle"]
        orng = np.random.default_rng(o["seed"] + seed)
        bg = oracle["pool"].iloc[orng.choice(len(oracle["pool"]), size=o["background"], replace=False)]
        oruns = method_runners(oracle["predict"], explain, bg, train, features, ag, params, seed, oracle["expected"],
                               lambda: oracle["scm"],
                               lambda: (oracle["weights"], TrueSampler(oracle["scm"], features)))
        for m in methods:
            res = oruns[m]()
            imp, _ = importance_of(m, res, params)
            imps.append(pd.DataFrame({"method": m + "_oracle", "feature": list(imp), "importance": list(imp.values())}))
            v = res.values.copy()
            if res.diagnostics is not None:
                v = v.join(res.diagnostics)
            vals.append(v.assign(method=m + "_oracle", outcome=explain[D.OUTCOME].to_numpy()).rename_axis("record").reset_index())
            primary = rule_names(params)[0]
            for (tt, pop, scope), t in tc.items():
                if scope != "with_era":
                    continue
                f_anc = {f: fitted[m][f] for f in t["raw"]}
                o_anc = {f: imp[f] for f in t["raw"]}
                oracle_rows.append({"method": m, "truth_type": tt, "population": pop,
                                    "tau_oracle_truth": kendall_tau_b(o_anc, t["tied"][primary]),
                                    "tau_fitted_truth": kendall_tau_b(f_anc, t["tied"][primary]),
                                    "tau_fitted_oracle": kendall_tau_b(f_anc, o_anc),
                                    "l1_fitted_oracle": l1_distance(f_anc, o_anc)})
    return ({"info": info, "rows": rows, "values": pd.concat(vals, ignore_index=True), "oracle": oracle_rows},
            pd.concat(imps, ignore_index=True))


_CONTEXT: dict = {}


def _init_worker(params: dict, data_dir: Path, methods: list[str]) -> None:
    """Every library runs single-threaded, so floating-point results do not depend on `workers`."""
    from threadpoolctl import threadpool_limits
    threadpool_limits(1)
    cfg = D.load_generation_config()
    ag = D.analysis_graph(cfg)
    fsets = D.feature_sets(data_dir)
    _CONTEXT.update(params=params, data_dir=data_dir, methods=methods, cfg=cfg, ag=ag, fsets=fsets,
                    distance=D.distance_to_outcome(ag),
                    truth=truth_context(D.truth_values(data_dir), fsets["ancestor"], params))


def _oracle() -> dict:
    c = _CONTEXT
    if "oracle" not in c:
        c["oracle"] = oracle_context(c["cfg"], c["ag"], c["fsets"]["ancestor"], c["params"])
    return c["oracle"]


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
        df = D.load_replicate(c["data_dir"], r, ds)
        use_oracle = fs == "ancestor" and r <= c["params"]["oracle"]["replicates"]
        res, imp = run_one(df, c, fsets[fs], c["params"], run_seed(c["params"], r, ds), c["methods"],
                           _oracle() if use_oracle else None)
        key = {"replicate": r, "dataset": ds, "feature_set": fs}
        _write(res["values"].assign(**key), cell_dir / f"{name}_values.csv.gz")
        _write(imp.assign(**key), cell_dir / f"{name}_importance.csv")
        if use_oracle:
            _write(pd.DataFrame([{**key, **row} for row in res["oracle"]]), cell_dir / f"{name}_oracle.csv")
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
    if (len(ag.nodes), len(ag.edges)) != (36, 63):
        raise RuntimeError(f"analysis graph has {len(ag.nodes)} nodes and {len(ag.edges)} edges, expected 36 and 63")
    check_oracle_inputs(cfg, D.feature_sets(data_dir)["ancestor"])
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
    oracle_files = [cell_dir / f"{n}_oracle.csv" for n in names if (cell_dir / f"{n}_oracle.csv").exists()]
    if oracle_files:
        _write(pd.concat([pd.read_csv(f) for f in oracle_files], ignore_index=True), out_dir / "oracle_decomposition.csv")
    tc = truth_context(D.truth_values(data_dir), D.feature_sets(data_dir)["ancestor"], params)
    _write(tie_group_table(tc), out_dir / "tie_groups.csv")
    _write(paired_contrasts(per_run), out_dir / "contrasts.csv")
    summary = summarize(per_run)
    _write(summary, out_dir / "summary.csv")
    write_provenance(out_dir, want, len(cells), workers)
    return summary


KEY = ["dataset", "feature_set", "method", "truth_type", "population", "scope"]


def summarize(per_run: pd.DataFrame) -> pd.DataFrame:
    """Mean and MCSE over replicates; `primary` marks each data set's own population with mission era included."""
    metrics = [c for c in ["kendall_tau_b", "kendall_tau_b_delta5", "kendall_tau_b_mcse", "weighted_tau", "l1_distance",
                           "top_k_recovery", "proximity_bias_index", "non_ancestor_share", "auc_model",
                           "auc_true_probability", "max_efficiency_error", "ng_excluded_share", "rescale_flag_share",
                           "explained_events", "seconds"] if c in per_run]
    g = per_run.groupby(KEY)[metrics]
    mean = g.mean().add_suffix("_mean")
    se = (g.std() / np.sqrt(g.count())).add_suffix("_mcse")
    out = pd.concat([g.size().rename("runs"), mean, se], axis=1).reset_index()
    out["primary"] = (out["population"] == out["dataset"].map(D.TRUTH_POPULATION)) & (out["scope"] == "with_era")
    return out


def paired_contrasts(per_run: pd.DataFrame, z: float = 1.96) -> pd.DataFrame:
    """Order effect tau_b(O) - tau_b(S) and value-function effect tau_b(I) - tau_b(O), paired within replicate.

    Intervals are mean +/- z SD / sqrt(n) over replicates. Classification per data set, feature set and contrast, in
    each population (the data set's own is `primary`), from the two co-primary truths (per_unit and pop):

    - **robust**: both intervals exclude zero, with the same sign (takes precedence);
    - **truth_dependent**: the truths substantively disagree, either opposite signs with both intervals excluding
      zero, or a within-replicate paired difference between the truths' contrasts (per_unit minus pop) whose interval
      excludes zero;
    - **inconclusive**: everything else, including one truth excluding zero and the other not, with no paired
      difference.

    Record-level (`rec`) rows are a sensitivity (`role`) and do not enter the classification. Every row carries the
    paired difference between truths of its group: truth_diff_mean, truth_diff_ci_low, truth_diff_ci_high.
    """
    pr = per_run[per_run["scope"] == "with_era"]
    key = ["dataset", "feature_set", "truth_type", "population"]
    wide = pr.pivot_table(index=key + ["replicate"], columns="method", values="kendall_tau_b")
    grp = ["dataset", "feature_set", "contrast", "population"]

    def interval(x: pd.Series) -> tuple[float, float, float, float, int]:
        n = int(x.notna().sum())
        m = float(x.mean())
        se = float(x.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
        return m, se, m - z * se, m + z * se, n

    rows, diffs = [], []
    for name, (a, b) in CONTRASTS.items():
        if a not in wide or b not in wide:
            continue
        d = (wide[a] - wide[b]).dropna()
        for k, g in d.groupby(level=key):
            m, se, lo, hi, n = interval(g)
            rows.append({**dict(zip(key, k)), "contrast": name, "mean": m, "mcse": se, "n": n, "ci_low": lo, "ci_high": hi})
        co = d.reset_index().rename(columns={0: "d"})
        co = co[co["truth_type"].isin(["per_unit", "pop"])].pivot_table(
            index=["dataset", "feature_set", "population", "replicate"], columns="truth_type", values="d").dropna()
        for k, g in (co["per_unit"] - co["pop"]).groupby(level=["dataset", "feature_set", "population"]):
            m, se, lo, hi, n = interval(g)
            diffs.append({"dataset": k[0], "feature_set": k[1], "population": k[2], "contrast": name,
                          "truth_diff_mean": m, "truth_diff_ci_low": lo, "truth_diff_ci_high": hi})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.merge(pd.DataFrame(diffs), on=grp, how="left")

    def classify(g: pd.DataFrame) -> str:
        co = g[g["truth_type"].isin(["per_unit", "pop"])]
        if len(co) < 2 or pd.isna(co["truth_diff_ci_low"]).any():
            return "incomplete"
        excl = (co["ci_low"] > 0) | (co["ci_high"] < 0)
        same_sign = np.sign(co["mean"]).nunique() == 1
        if excl.all() and same_sign:
            return "robust"
        diff_excl = bool((co["truth_diff_ci_low"].iloc[0] > 0) or (co["truth_diff_ci_high"].iloc[0] < 0))
        if (excl.all() and not same_sign) or diff_excl:
            return "truth_dependent"
        return "inconclusive"

    cls = out.groupby(grp).apply(classify, include_groups=False)
    out = out.merge(cls.rename("classification").reset_index(), on=grp)
    out["role"] = np.where(out["truth_type"].isin(["per_unit", "pop"]), "co_primary", "sensitivity")
    out["primary"] = out["population"] == out["dataset"].map(D.TRUTH_POPULATION)
    return out
