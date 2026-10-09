"""Usage: uv run python run_confounding_shap.py key [--permutations N]
       uv run python run_confounding_shap.py run [--replicates N] [--model-version v2.5] [--budget B] [--resume | --fresh]

`key` computes reference credits from the true model (source and selected populations) for each treatment.
`run` fits ConfoundingSHAP (shapiq's released game, TabPFN S-learners, RegressionMSR) on each replicate of the
configured data sets for each treatment, one cell per (replicate, data set, treatment), and combines the cells.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from attribution_validation import data as D
from attribution_validation.confounding import (common_cause_sets, confounder_mass, confounder_recovery,
                                                confounding_credits, pretreatment_covariates, tabpfn_factory)
from attribution_validation.reference_key import reference_credits
from attribution_validation.truesim import TrueSCM

HERE = Path(__file__).resolve().parent


def settings() -> dict:
    return yaml.safe_load((HERE / "config" / "params.yaml").read_text())


def contexts(cfg_params: dict, data_dir: Path):
    cfg = D.load_generation_config()
    ag = D.analysis_graph(cfg)
    anc = D.feature_sets(data_dir)["ancestor"]
    cs = cfg_params["confounding_shap"]
    drugs = cs["treatments"]
    out = {}
    for t in drugs:
        cov = pretreatment_covariates(ag, anc, t, ["thiazides", *drugs])
        out[t] = {"covariates": cov, "roles": common_cause_sets(ag, cov, t)}
    return cfg, ag, out


def key(args) -> None:
    params = settings()
    cs = params["confounding_shap"]
    data_dir = D.env_dir("DATA_DIR", "./02_data")
    cfg, ag, ctx = contexts(params, data_dir)
    scm = TrueSCM(cfg)
    n = cs["reference"]["records"]
    rng = np.random.default_rng(cs["reference"]["seed"])
    u = scm.noise(n, rng)
    v = scm.simulate(u)
    pop = pd.DataFrame({**v, **u})
    sel = scm.selected(v, u)
    rows, meta = [], {}
    for t, c in ctx.items():
        for popname, mask in (("source", np.ones(n, dtype=bool)), ("selected", sel)):
            t0 = time.time()
            r = reference_credits(scm, pop, u, t, c["covariates"], mask, args.permutations or cs["reference"]["permutations"],
                                  cs["reference"]["seed"])
            for f, val in r["credits"].items():
                rows.append({"treatment": t, "population": popname, "covariate": f, "credit": val,
                             "direct_common_cause": f in c["roles"]["direct"], "any_common_cause": f in c["roles"]["all"]})
            meta[f"{t}/{popname}"] = {k: r[k] for k in ("ate", "crude", "v_empty", "v_full", "n", "subsets")} | {"seconds": time.time() - t0}
            print(t, popname, {k: round(r[k], 5) for k in ("ate", "crude", "v_empty", "v_full")}, f"{time.time() - t0:.0f}s", flush=True)
    out = D.env_dir("OUTPUT_DIR", "./outputs") / "04_attribution_validation"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "confounding_reference_credits.csv", index=False)
    (out / "confounding_reference_meta.json").write_text(json.dumps(meta, indent=1))


def cell_name(r: int, ds: str, t: str) -> str:
    return f"r{r:03d}_{ds}_{t}"


def run(args) -> None:
    params = settings()
    cs = params["confounding_shap"]
    data_dir = D.env_dir("DATA_DIR", "./02_data")
    cfg, ag, ctx = contexts(params, data_dir)
    D.check_data_current(data_dir, cfg)
    out = D.env_dir("OUTPUT_DIR", "./outputs") / "04_attribution_validation"
    cells = out / "confounding_cells"
    want = {"model_version": args.model_version, "budget": args.budget, "n_estimators": cs["n_estimators"],
            "treatments": cs["treatments"], "seed": params["seed"], "proxy_model": cs["proxy_model"]}
    mf = cells / "manifest.json"
    existing = list(cells.glob("*_runs.csv")) if cells.exists() else []
    if args.fresh:
        for f in cells.glob("*"):
            f.unlink()
    elif existing:
        if not args.resume:
            raise SystemExit(f"{cells} holds results: pass --resume to continue or --fresh to discard")
        if json.loads(mf.read_text()) != want:
            raise SystemExit("settings changed since the saved cells were written: pass --fresh")
    cells.mkdir(parents=True, exist_ok=True)
    mf.write_text(json.dumps(want, indent=1))
    reps = range(1, (args.replicates or cs["replicates"]) + 1)
    todo = [(r, ds, t) for r in reps for ds in cs["datasets"] for t in cs["treatments"]
            if not (cells / f"{cell_name(r, ds, t)}_runs.csv").exists()]
    print(f"{len(todo)} cells to run", flush=True)
    start = time.time()
    for i, (r, ds, t) in enumerate(todo, 1):
        seed = params["seed"] + 1000 * r + 10 * params["datasets"].index(ds) + cs["treatments"].index(t)
        df = D.load_replicate(data_dir, r, ds)
        c = ctx[t]
        res = confounding_credits(df, c["covariates"], t, tabpfn_factory(args.model_version, seed), args.budget, seed,
                                  n_estimators=cs["n_estimators"], proxy=cs["proxy_model"])
        cr = res["credits"]
        key = {"replicate": r, "dataset": ds, "treatment": t}
        pd.DataFrame([{**key, "covariate": f, "credit": v, "direct_common_cause": f in c["roles"]["direct"],
                       "any_common_cause": f in c["roles"]["all"]} for f, v in cr.items()]).to_csv(
            cells / f"{cell_name(r, ds, t)}_credits.csv", index=False)
        row = {**key, "model_version": args.model_version, "budget": args.budget, "n_covariates": len(cr),
               "mass_direct": confounder_mass(cr, c["roles"]["direct"]), "recovery_direct": confounder_recovery(cr, c["roles"]["direct"]),
               "mass_any": confounder_mass(cr, c["roles"]["all"]), "recovery_any": confounder_recovery(cr, c["roles"]["all"]),
               **{k: res[k] for k in ("v_empty", "v_full", "crude", "adjusted_ate", "exact", "coalitions_fitted", "seconds", "n", "treated", "events")}}
        pd.DataFrame([row]).to_csv(cells / f"{cell_name(r, ds, t)}_runs.csv", index=False)
        print(f"[{i}/{len(todo)}] {cell_name(r, ds, t)} {res['seconds'] / 60:.1f} min (total {(time.time() - start) / 60:.0f} min)", flush=True)
    names = [cell_name(r, ds, t) for r in reps for ds in cs["datasets"] for t in cs["treatments"]]
    pd.concat([pd.read_csv(cells / f"{n}_credits.csv") for n in names]).to_csv(out / "confounding_credits.csv", index=False)
    runs = pd.concat([pd.read_csv(cells / f"{n}_runs.csv") for n in names])
    runs.to_csv(out / "confounding_runs.csv", index=False)
    g = runs.groupby(["dataset", "treatment"])
    summ = pd.concat([g.size().rename("runs"), g[["mass_direct", "recovery_direct", "mass_any", "recovery_any", "v_empty", "crude",
                                                   "adjusted_ate", "seconds"]].mean().add_suffix("_mean"),
                      (g[["mass_direct", "recovery_direct", "mass_any", "recovery_any"]].std() / np.sqrt(g.size().to_numpy()[:, None])
                       ).add_suffix("_mcse")], axis=1).reset_index()
    summ.to_csv(out / "confounding_summary.csv", index=False)
    print(summ.round(3).to_string())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    k = sub.add_parser("key")
    k.add_argument("--permutations", type=int)
    r = sub.add_parser("run")
    r.add_argument("--replicates", type=int)
    r.add_argument("--model-version", default=None)
    r.add_argument("--budget", type=int, default=None)
    r.add_argument("--resume", action="store_true")
    r.add_argument("--fresh", action="store_true")
    a = ap.parse_args()
    cs = settings()["confounding_shap"]
    if a.cmd == "run":
        a.model_version = a.model_version or cs["model_version"]
        a.budget = a.budget or cs["budget"]
        run(a)
    else:
        key(a)


if __name__ == "__main__":
    main()
