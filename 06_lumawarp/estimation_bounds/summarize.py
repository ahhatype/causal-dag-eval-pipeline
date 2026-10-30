import glob, os, numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); rows = []
for f in glob.glob(os.path.join(H, "results", "*", "*", "*", "rep*.csv")):
    r = pd.read_csv(f); r["dataset"] = f.split("/")[-4]; r["route"] = f.split("/")[-2]; rows.append(r)
R = pd.concat(rows, ignore_index=True); R.to_csv(os.path.join(H, "results_all.csv"), index=False)
print("max |LW estimate - folder 34| at the original bound:", R.check_vs_34.max())
S = R.groupby(["dataset", "treatment", "route", "learner", "bound"]).apply(lambda g: pd.Series(dict(reps=len(g), bias=(g.est - g.truth).mean(), emp_se=g.est.std(), mean_se=g.se.mean(), coverage=g.covered.mean(), clipped=g.clipped.mean()))).reset_index()
S.to_csv(os.path.join(H, "summary.csv"), index=False); pd.set_option("display.width", 250)
print(S[S.treatment == "k_citrate"].round(4).to_string(index=False))
k = S[(S.treatment == "k_citrate") & (S.route == "parents") & (S.learner == "lw")].pivot(index="bound", columns="dataset", values="coverage")
print("\nK+ citrate, parents, LW-TMLE coverage by bound:\n", k.round(2))
ok = [b for b in k.index if k.loc[b, "reference_subsample"] >= 0.9 and k.loc[b, "astronaut_set"] >= 0.9 and k.loc[b, "full_set"] >= 0.9]
print("S1:", f"bounds that restore >= 0.90 everywhere: {ok}" if ok else "no bound restores coverage on both 900-row sets -> the sample, not the bounds")
