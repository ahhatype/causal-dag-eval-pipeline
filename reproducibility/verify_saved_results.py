"""Usage: python3 reproducibility/verify_saved_results.py [--repo .] [--json report.json]

Recomputes the reported summary tables from the committed per-run and per-replicate results and compares them with the
stored tables (nothing is simulated or refit; needs only NumPy and pandas). The exit status is zero only if every
comparison passes. --json writes the full report: input SHA-256 hashes, runtime versions, comparisons, headline values.
"""

import argparse
import ast
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

TOLERANCE = 1e-10
SCOPE = ("Saved-output aggregation only; no simulations, refits, .env reads, fitting-module imports or source edits. "
         "Passing does not verify underlying fitted models, per-record attributions or unavailable nuisance "
         "predictions.")
PIPELINE = "04_attribution_validation/src/attribution_validation/pipeline.py"
AGGREGATORS = {"summarize", "_interval", "_classify", "paired_contrasts", "data_condition_contrasts"}
CONSTANTS = {"KEY", "CONTRASTS", "DATA_CONDITIONS", "CONDITION_METRICS", "CONDITION_METHODS"}
# attribution_validation.data.TRUTH_POPULATION
TRUTH_POPULATION = {"full_set": "source", "reference_subsample": "source", "astronaut_set": "selected"}
HEADLINE_METHODS = ["standard_shap", "ordering_only_shap", "ng_style_shap", "interventional_shap"]
ESTIMATORS = [("lw", "LW-TMLE"), ("sl", "SL-TMLE"), ("gcomp_logistic", "Logistic g-computation"),
              ("naive", "Naive difference")]
# 06_lumawarp/edge_audit/score.py: a candidate edge is kept when its score exceeds the learner's threshold.
EDGE_THRESHOLDS = {"lumawarp": 5, "glm": 4, "spline": 2, "xgb": 3}
EDGE_LEARNERS = {"lumawarp": "Lumawarp", "glm": "logistic/linear", "spline": "spline additive", "xgb": "XGBoost"}


def compare(name: str, expected: pd.DataFrame, actual: pd.DataFrame, key: list[str]) -> dict:
    """One comparison record: row keys and text must match exactly, numbers within TOLERANCE (absolute and relative)."""
    expected = expected.sort_values(key).reset_index(drop=True)
    actual = actual.sort_values(key).reset_index(drop=True)
    record = {"name": name, "expected_rows": len(expected), "actual_rows": len(actual), "key": key,
              "numeric_cells": 0, "max_abs_difference": 0.0, "failures": []}
    if len(expected) != len(actual):
        record["failures"].append("row count")
    if expected[key].astype(str).to_dict("records") != actual[key].astype(str).to_dict("records"):
        record["failures"].append("row keys")
    if not record["failures"]:
        for col in expected:
            if col not in actual:
                record["failures"].append("missing " + col)
                continue
            a, b = expected[col], actual[col]
            if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
                x, y = a.to_numpy(dtype=float), b.to_numpy(dtype=float)
                good = np.isclose(x, y, rtol=TOLERANCE, atol=TOLERANCE, equal_nan=True)
                record["numeric_cells"] += len(x)
                delta = np.abs(x - y)
                finite = delta[np.isfinite(delta)]
                if len(finite):
                    record["max_abs_difference"] = max(record["max_abs_difference"], float(finite.max()))
                if not good.all():
                    examples = [{"row": int(i), "expected": float(x[i]), "recomputed": float(y[i])}
                                for i in np.flatnonzero(~good)[:3]]
                    record["failures"].append({"column": col, "bad_rows": int((~good).sum()), "examples": examples})
            elif not a.fillna("<NA>").astype(str).equals(b.fillna("<NA>").astype(str)):
                record["failures"].append({"column": col, "issue": "text values differ"})
    record["passed"] = not record["failures"]
    return record


def pipeline_aggregators(source: str, filename: str) -> dict:
    """The attribution module's aggregation functions, compiled from its source alone.

    Importing the package would load the fitting libraries; these functions need only pandas and NumPy.
    """
    keep = [n for n in ast.parse(source).body
            if (isinstance(n, ast.FunctionDef) and n.name in AGGREGATORS)
            or (isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id in CONSTANTS for t in n.targets))]
    ns = {"pd": pd, "np": np, "D": SimpleNamespace(TRUTH_POPULATION=TRUTH_POPULATION)}
    exec(compile(ast.Module(body=keep, type_ignores=[]), filename, "exec"), ns)
    return ns


def structure_summary(per_replicate: pd.DataFrame, saved: pd.DataFrame) -> pd.DataFrame:
    """Mean and Monte Carlo SE of each per-replicate column that the saved summary reports."""
    rows = []
    for (dataset, variant), g in per_replicate.groupby(["dataset", "variant"]):
        row = {"dataset": dataset, "variant": variant}
        for col in saved:
            if col.endswith("_mean") and col[:-5] in g:
                row[col] = pd.to_numeric(g[col[:-5]], errors="coerce").mean()
            elif col.endswith("_se") and col[:-3] in g:
                z = pd.to_numeric(g[col[:-3]], errors="coerce")
                row[col] = z.std(ddof=1) / np.sqrt(z.count())
        rows.append(row)
    return pd.DataFrame(rows)


def headlines(structure: pd.DataFrame, summary: pd.DataFrame) -> list[dict]:
    """Per data set: discovery recall, replicates with rejected implications, rank agreement of the four methods."""
    out = []
    for dataset, population in TRUTH_POPULATION.items():
        s = structure[(structure.dataset == dataset) & (structure.variant == "MGM-PC-Stable")]
        a = summary[(summary.dataset == dataset) & (summary.feature_set == "ancestor")
                    & (summary.truth_type == "per_unit") & (summary.population == population)
                    & (summary.scope == "with_era")]
        out.append({
            "dataset": dataset, "replicates": len(s),
            "adjacency_recall_mean": float(s.adjacency_recall.mean()),
            "adjacency_recall_printed": f"{s.adjacency_recall.mean():.2f}",
            "additive_rejected_replicates": int((s.rejected_other > 0).sum()),
            "nonadditive_rejected_replicates": int((s.rejected_nonadditive > 0).sum()),
            "selection_rejected_replicates": int(s.selection_violation_detected.sum()),
            "attribution": {r.method: {"tau": r.kendall_tau_b_mean, "printed": f"{r.kendall_tau_b_mean:.2f}",
                                       "top5_count": 5 * r.top_k_recovery_mean}
                            for r in a.itertuples() if r.method in HEADLINE_METHODS}})
    return out


def estimation_summary(results: pd.DataFrame) -> pd.DataFrame:
    """Bias, empirical SE, RMSE and sign agreement per estimator; mean SE and interval coverage for the two TMLEs."""
    rows = []
    for (dataset, treatment, route), g in results.groupby(["dataset", "treatment", "route"]):
        truth = float(g.truth.iloc[0])
        for est, label in ESTIMATORS:
            x = g[est + "_est"]
            row = {"dataset": dataset, "treatment": treatment, "route": route, "estimator": label, "reps": len(g),
                   "truth": truth, "mean_estimate": x.mean(), "bias": x.mean() - truth, "empirical_se": x.std(),
                   "rmse": np.sqrt(((x - truth) ** 2).mean()), "sign_right": (np.sign(x) == np.sign(truth)).mean()}
            if est in ("lw", "sl"):
                row.update(mean_se=g[est + "_se"].mean(),
                           coverage=((g[est + "_lower"] <= truth) & (truth <= g[est + "_upper"])).mean())
            rows.append(row)
    return pd.DataFrame(rows)


def bounds_summary(results: pd.DataFrame) -> pd.DataFrame:
    """Bias, SEs, coverage and clipped share per propensity bound."""
    key = ["dataset", "treatment", "route", "learner", "bound"]
    rows = []
    for values, g in results.groupby(key):
        rows.append(dict(zip(key, values)) | dict(
            reps=len(g), bias=(g.est - g.truth).mean(), emp_se=g.est.std(), mean_se=g.se.mean(),
            coverage=((g.lower <= g.truth) & (g.truth <= g.upper)).mean(), clipped=g.clipped.mean()))
    return pd.DataFrame(rows)


def primary(runs: pd.DataFrame) -> pd.DataFrame:
    """Runs scored against the per-unit truth of each data set's own population, mission era included."""
    return runs[(runs.truth_type == "per_unit") & (runs.scope == "with_era")
                & (runs.population == runs.dataset.map(TRUTH_POPULATION))]


def explained_predictor_summary(lumawarp_runs: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    """Each method explaining the Lumawarp predictor, paired with the same method on the Super Learner."""
    lw = primary(lumawarp_runs).copy()
    lw["base"] = lw.method.str.replace("_lumawarp", "")
    merged = lw.merge(primary(runs[runs.feature_set == "ancestor"]), left_on=["replicate", "dataset", "base"],
                      right_on=["replicate", "dataset", "method"], suffixes=("_lw", "_sl"))
    rows = []
    for (dataset, base), g in merged.groupby(["dataset", "base"]):
        d = g.kendall_tau_b_lw - g.kendall_tau_b_sl
        rows.append(dict(
            dataset=dataset, method=base, reps=len(g), tau_lw=g.kendall_tau_b_lw.mean(),
            tau_sl=g.kendall_tau_b_sl.mean(), diff=d.mean(), diff_mcse=d.std() / np.sqrt(len(g)),
            top5_lw=g.top_k_recovery_lw.mean(), top5_sl=g.top_k_recovery_sl.mean(),
            pbi_lw=g.proximity_bias_index_lw.mean(), pbi_sl=g.proximity_bias_index_sl.mean(),
            auc_lw=g.auc_lumawarp.mean(), auc_sl=g.auc_superlearner.mean()))
    return pd.DataFrame(rows)


def confounding_summary(runs: pd.DataFrame, saved: pd.DataFrame) -> pd.DataFrame:
    """Mean and Monte Carlo SE of each per-run ConfoundingSHAP metric that the saved summary reports."""
    rows = []
    for (dataset, treatment), g in runs.groupby(["dataset", "treatment"]):
        row = {"dataset": dataset, "treatment": treatment, "runs": len(g)}
        for col in saved:
            if col.endswith("_mean"):
                row[col] = g[col[:-5]].mean()
            elif col.endswith("_mcse"):
                row[col] = g[col[:-5]].std() / np.sqrt(len(g))
        rows.append(row)
    return pd.DataFrame(rows)


def edge_decisions(edge: pd.DataFrame) -> dict:
    """Comparison record: every stored kept flag equals score > threshold for its learner."""
    same = (edge.score > edge.learner.map(EDGE_THRESHOLDS)).equals(edge.kept)
    return {"name": "edge-audit stored scores reproduce retained-edge decisions", "expected_rows": len(edge),
            "actual_rows": len(edge), "key": ["dataset", "rep", "target", "cand", "learner"],
            "numeric_cells": len(edge), "max_abs_difference": 0.0, "failures": [] if same else ["kept flags differ"],
            "passed": same}


def edge_scores(edge: pd.DataFrame, n_edges: int) -> pd.DataFrame:
    """Per replicate and learner: edges called, recall and precision against the analysis graph, noise kept."""
    rows = []
    for (dataset, rep, tag), g in edge.groupby(["dataset", "rep", "learner"]):
        called = g[~g["null"] & g.kept]
        tp = int(called["true"].sum())
        rows.append(dict(dataset=dataset, rep=rep, learner=EDGE_LEARNERS[tag], edges_called=len(called),
                         recall=tp / n_edges, precision=tp / max(1, len(called)),
                         noise_kept=g[g["null"]].kept.mean()))
    return pd.DataFrame(rows)


def edge_summary(scores: pd.DataFrame) -> pd.DataFrame:
    """Means of the per-replicate edge-audit scores, with the Monte Carlo SE of recall."""
    rows = []
    for (dataset, learner), g in scores.groupby(["dataset", "learner"]):
        rows.append(dict(dataset=dataset, learner=learner, reps=g.rep.nunique(), edges_called=g.edges_called.mean(),
                         recall=g.recall.mean(), recall_mcse=g.recall.std() / np.sqrt(len(g)),
                         precision=g.precision.mean(), noise_kept=g.noise_kept.mean()))
    return pd.DataFrame(rows)


def hybrid_summary(by_rep: pd.DataFrame, saved: pd.DataFrame) -> pd.DataFrame:
    """Pooled orientation accuracy and the mean of every other per-replicate metric in the saved summary."""
    rows = []
    for method, g in by_rep.groupby("method"):
        row = {"method": method, "reps": g.rep.nunique(),
               "orient_acc": g.oriented_correct.sum() / g.oriented_true.sum()}
        for col in saved:
            if col not in row:
                row[col] = g[col].mean()
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", type=Path, default=Path.cwd())
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()
    repo = args.repo.resolve()
    inputs, checks = {}, []

    def read_bytes(path: str) -> bytes:
        b = (repo / path).read_bytes()
        inputs[path] = {"bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
        return b

    def read(path: str) -> pd.DataFrame:
        read_bytes(path)
        return pd.read_csv(repo / path)

    def check(name: str, expected: pd.DataFrame, actual: pd.DataFrame, key: list[str]) -> None:
        checks.append(compare(name, expected, actual, key))

    agg = pipeline_aggregators(read_bytes(PIPELINE).decode(), str(repo / PIPELINE))
    attribution = "outputs/04_attribution_validation/"
    runs, summary = read(attribution + "per_run.csv"), read(attribution + "summary.csv")
    check("attribution summary: pure latest-source aggregation", summary, agg["summarize"](runs), agg["KEY"])
    check("paired method contrasts: pure latest-source aggregation", read(attribution + "contrasts.csv"),
          agg["paired_contrasts"](runs), ["dataset", "feature_set", "truth_type", "population", "contrast"])
    check("attribution sample/selection contrasts: pure latest-source aggregation",
          read(attribution + "data_condition_contrasts.csv"), agg["data_condition_contrasts"](runs),
          ["method", "metric", "feature_set", "truth_type", "contrast"])

    structure = read("outputs/03_structure_validation/per_replicate.csv")
    saved = read("outputs/03_structure_validation/summary.csv")
    recomputed = structure_summary(structure, saved)
    check("structure means and MCSE: independent reconstruction", saved[list(recomputed.columns)], recomputed,
          ["dataset", "variant"])

    lumawarp = "outputs/06_lumawarp/"
    results = read(lumawarp + "estimation/results_all_reps.csv")
    check("TMLE and comparator summary from retained estimates", read(lumawarp + "estimation/summary.csv"),
          estimation_summary(results), ["dataset", "treatment", "route", "estimator"])
    results = read(lumawarp + "estimation_bounds/results_all.csv")
    check("propensity-bound summary from retained estimates", read(lumawarp + "estimation_bounds/summary.csv"),
          bounds_summary(results), ["dataset", "treatment", "route", "learner", "bound"])
    lumawarp_runs = read(lumawarp + "attribution_lumawarp_base/per_run.csv")
    check("explained-predictor summary from retained per-run files",
          read(lumawarp + "attribution_lumawarp_base/summary.csv"), explained_predictor_summary(lumawarp_runs, runs),
          ["dataset", "method"])

    confounding_runs = read(attribution + "confounding_runs.csv")
    saved = read(attribution + "confounding_summary.csv")
    check("ConfoundingSHAP summary from saved per-run metrics", saved, confounding_summary(confounding_runs, saved),
          ["dataset", "treatment"])

    edge = read(lumawarp + "edge_audit/edges_all.csv.gz")
    saved_scores = read(lumawarp + "edge_audit/scores_by_rep.csv")
    saved = read(lumawarp + "edge_audit/summary.csv")
    checks.append(edge_decisions(edge))
    scores = edge_scores(edge, len(read("outputs/03_structure_validation/analysis_graph_edges.csv")))
    check("edge-audit per-replicate scores from stored candidate scores", saved_scores, scores,
          ["dataset", "rep", "learner"])
    check("edge-audit summary from reconstructed replicate scores", saved, edge_summary(scores),
          ["dataset", "learner"])

    for arm in ("armA", "armB"):
        by_rep = read(f"{lumawarp}nonlinear_sensitivity/{arm}/hybrid_by_rep.csv")
        saved = read(f"{lumawarp}nonlinear_sensitivity/{arm}/hybrid_summary.csv")
        check(f"nonlinear {arm} hybrid summary from saved replicate metrics", saved, hybrid_summary(by_rep, saved),
              ["method"])

    passed = all(c["passed"] for c in checks)
    report = {"scope": SCOPE, "runtime": {"python": sys.version, "pandas": pd.__version__, "numpy": np.__version__},
              "aggregation_checks": checks, "headlines": headlines(structure, summary), "inputs": inputs,
              "all_aggregation_checks_passed": passed}
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    brief = ["name", "expected_rows", "numeric_cells", "max_abs_difference", "passed", "failures"]
    print(json.dumps({"passed": passed, "checks": [{k: c.get(k) for k in brief} for c in checks], "scope": SCOPE},
                     separators=(",", ":")))
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
