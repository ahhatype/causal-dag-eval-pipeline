"""Usage: uv run python run_ng_convergence.py [--workers N]

Ng-style convergence on phi (replicate 1 of each data set, both feature sets): two independent halves at each order
budget, compared by the correlation of per-feature mean |phi|, its largest absolute difference, and tau-b. If the last
configured budget fails, it is doubled until it passes or a run exceeds `max_budget_factor` times the runtime at the
production budget."""

import argparse
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from attribution_validation import data as D
from attribution_validation.evaluation import kendall_tau_b
from attribution_validation.ng import fitted_ng_inputs, ng_style_shap
from attribution_validation.pipeline import ng_importance, prepare_run, run_seed


def ng_run(pr: dict, features: list[str], ag, params: dict, orders: int, seed: int) -> tuple[dict, float]:
    ng = params["ng_style_shap"]
    t0 = time.time()
    res = ng_style_shap(pr["sl"].predict, pr["explain"], features, *fitted_ng_inputs(pr["train"], ag, features),
                        ng["samples"], orders, seed, pr["expected"])
    return ng_importance(res, ng["exclude_below"])[0], time.time() - t0


def cell(args: tuple) -> list[dict]:
    from threadpoolctl import threadpool_limits
    threadpool_limits(1)
    params, data_dir, ds, fs = args
    cv = params["ng_convergence"]
    cfg = D.load_generation_config()
    ag = D.analysis_graph(cfg)
    features = D.feature_sets(data_dir)[fs]
    seed = run_seed(params, 1, ds)
    pr = prepare_run(D.load_replicate(data_dir, 1, ds), {"cfg": cfg}, features, params, seed)
    _, base_secs = ng_run(pr, features, ag, params, params["ng_style_shap"]["orders"], seed)
    rows, budgets = [], list(cv["orders"])
    while budgets:
        t = budgets.pop(0)
        (a, sa), (b, sb) = (ng_run(pr, features, ag, params, t, cv["seed"] + h * 7919 + t) for h in (1, 2))
        va, vb = np.array([a[f] for f in features]), np.array([b[f] for f in features])
        ok = np.all(np.isfinite(va)) and np.all(np.isfinite(vb))
        r = {"dataset": ds, "feature_set": fs, "orders": t, "seconds_per_half": (sa + sb) / 2,
             "seconds_production_budget": base_secs,
             "correlation": float(np.corrcoef(va, vb)[0, 1]) if ok else np.nan,
             "max_abs_difference": float(np.max(np.abs(va - vb))) if ok else np.nan,
             "tau_b": kendall_tau_b(a, b) if ok else np.nan}
        r["passed"] = bool(r["tau_b"] >= cv["min_tau_b"] and r["correlation"] >= cv["min_correlation"])
        rows.append(r)
        if not budgets and not r["passed"] and r["seconds_per_half"] * 2 <= cv["max_budget_factor"] * base_secs:
            budgets.append(2 * t)
    return rows


def main() -> None:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    params = yaml.safe_load((here / "config" / "params.yaml").read_text())
    data_dir = D.env_dir("DATA_DIR", "./02_data")
    D.check_data_current(data_dir, D.load_generation_config())
    jobs = [(params, data_dir, ds, fs) for ds in params["datasets"] for fs in params["feature_sets"]]
    with ProcessPoolExecutor(a.workers) as pool:
        rows = [r for rs in pool.map(cell, jobs) for r in rs]
    out = D.env_dir("OUTPUT_DIR", "./outputs") / "04_attribution_validation"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "ng_convergence.csv", index=False)
    print(df.to_string())


if __name__ == "__main__":
    main()
