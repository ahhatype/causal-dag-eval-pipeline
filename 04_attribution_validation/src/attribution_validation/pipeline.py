from __future__ import annotations

import time
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


def run_seed(params: dict, replicate: int, dataset: str, feature_set: str) -> int:
    return params["seed"] + 1000 * replicate + 10 * params["datasets"].index(dataset) + params["feature_sets"].index(feature_set)


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
    info = {"fit_seconds": time.time() - t0,
            "auc_model": auc(test[D.OUTCOME].to_numpy(), sl.predict(test[features].to_numpy(dtype=float))),
            "auc_true_probability": auc(test[D.OUTCOME].to_numpy(), D.true_outcome_probability(test, cfg)),
            **{f"weight_{k}": v for k, v in sl.weights_.items()}}
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
                                                 params["ng_causal_shap"]["iterations"], seed, expected),
        "interventional_shap": lambda: interventional_shap(sl.predict, fit_scm(train, ag), explain, background,
                                                           features, ag, params["interventional_shap"]["orders"],
                                                           params["interventional_shap"]["burn_in"],
                                                           params["interventional_shap"]["thin"], seed),
    }
    rows, imps = [], []
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
        imps.append(pd.DataFrame({"method": m, "feature": list(imp), "importance": list(imp.values())}))
    return {"info": info, "rows": rows}, pd.concat(imps, ignore_index=True)


def run(params: dict, data_dir: Path, out_dir: Path, replicates: int | None = None,
        methods: list[str] | None = None, datasets: list[str] | None = None,
        feature_sets: list[str] | None = None) -> pd.DataFrame:
    cfg = D.load_generation_config()
    available = D.check_data_current(data_dir, cfg)
    reps = range(1, min(available, replicates or available) + 1)
    ag = D.analysis_graph(cfg)
    if (len(ag.nodes), len(ag.edges)) != (36, 59):
        raise RuntimeError(f"analysis graph has {len(ag.nodes)} nodes and {len(ag.edges)} edges, expected 36 and 59")
    fsets = D.feature_sets(data_dir)
    ancestors = set(fsets["ancestor"])
    distance = D.distance_to_outcome(ag)
    methods = methods or params["methods"]
    out_dir.mkdir(parents=True, exist_ok=True)
    rows, imps = [], []
    for r in reps:
        for ds in datasets or params["datasets"]:
            df = D.load_replicate(data_dir, r, ds)
            target = D.GROUND_TRUTH_TARGET[ds]
            truth = D.ground_truth(data_dir, target).loc[fsets["ancestor"]]
            truth_se = D.ground_truth_se(data_dir, target).loc[fsets["ancestor"]]
            for fs in feature_sets or params["feature_sets"]:
                seed = run_seed(params, r, ds, fs)
                print(f"[replicate {r}] {ds} / {fs}", flush=True)
                res, imp = run_one(df, cfg, ag, fsets[fs], truth, truth_se, ancestors, distance, params, seed, methods)
                for row in res["rows"]:
                    rows.append({"replicate": r, "dataset": ds, "feature_set": fs, **res["info"], **row})
                imps.append(imp.assign(replicate=r, dataset=ds, feature_set=fs))
                pd.DataFrame(rows).to_csv(out_dir / "per_run.csv", index=False)
    per_run = pd.DataFrame(rows)
    pd.concat(imps, ignore_index=True).to_csv(out_dir / "importance.csv", index=False)
    summary = summarize(per_run)
    summary.to_csv(out_dir / "summary.csv", index=False)
    return summary


def summarize(per_run: pd.DataFrame) -> pd.DataFrame:
    metrics = ["kendall_tau_b", "top_k_recovery", "proximity_bias_index", "non_ancestor_share", "auc_model",
               "auc_true_probability", "max_efficiency_error", "seconds"]
    g = per_run.groupby(["dataset", "feature_set", "method"])[metrics]
    mean = g.mean().add_suffix("_mean")
    se = (g.std() / np.sqrt(g.count())).add_suffix("_se")
    return pd.concat([g.size().rename("runs"), mean, se], axis=1).reset_index()
