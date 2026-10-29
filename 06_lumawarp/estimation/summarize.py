import glob, os, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); rows = []
for f in glob.glob(os.path.join(H, "results", "*", "*", "*", "rep*.csv")):
    r = pd.read_csv(f); r["dataset"] = f.split("/")[-4]; r["route"] = f.split("/")[-2]; rows.append(r)
R = pd.concat(rows, ignore_index=True)
TV = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", os.environ.get("DATA_DIR", "02_data").removeprefix("./"), "truth_values.csv")); TV = TV[TV.truth_type == "per_unit"]
def tr(ds, t): return float(TV[(TV.feature == t) & (TV.population == ("selected" if ds == "astronaut_set" else "source"))].value.iloc[0])
R["truth"] = [tr(d, t) for d, t in zip(R.dataset, R.treatment)]
for e in ("lw", "sl"): R[f"{e}_covered"] = (R[f"{e}_lower"] <= R.truth) & (R.truth <= R[f"{e}_upper"])
R.to_csv(os.path.join(H, "results_all_reps.csv"), index=False); out = []
for (ds, t, rt), g in R.groupby(["dataset", "treatment", "route"]):
    truth = g.truth.iloc[0]
    for est, lab in (("lw", "LW-TMLE"), ("sl", "SL-TMLE"), ("gcomp_logistic", "Logistic g-computation"), ("naive", "Naive difference")):
        e = g[f"{est}_est"]; r = dict(dataset=ds, treatment=t, route=rt, estimator=lab, reps=len(g), truth=truth, mean_estimate=e.mean(), bias=e.mean() - truth, empirical_se=e.std(), rmse=np.sqrt(((e - truth) ** 2).mean()), sign_right=(np.sign(e) == np.sign(truth)).mean())
        if est in ("lw", "sl"): r.update(mean_se=g[f"{est}_se"].mean(), coverage=g[f"{est}_covered"].mean())
        out.append(r)
S = pd.DataFrame(out); S.to_csv(os.path.join(H, "summary.csv"), index=False); pd.set_option("display.width", 250); print(S.round(4).to_string(index=False))
