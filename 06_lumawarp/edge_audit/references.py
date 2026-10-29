"""Folder 13's reference learners on the identical files: logistic/linear, spline additive, XGBoost. Writes pred_<tag>.csv."""
import os, warnings, numpy as np, pandas as pd; warnings.filterwarnings("ignore")
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); TAGS = ("glm", "spline", "xgb")
def model(tag, X, binary):
    import xgboost as xgb
    from sklearn.linear_model import LogisticRegression, LinearRegression
    from sklearn.preprocessing import SplineTransformer
    from sklearn.compose import ColumnTransformer
    from sklearn.pipeline import make_pipeline
    lin = LogisticRegression(C=1e4, max_iter=3000) if binary else LinearRegression()
    if tag == "glm": return lin
    if tag == "spline":
        cont = [c for c in X.columns if X[c].nunique() > 2]
        return make_pipeline(ColumnTransformer([("s", SplineTransformer(n_knots=6, degree=3, extrapolation="linear"), cont)], remainder="passthrough"), LogisticRegression(C=10, max_iter=3000) if binary else LinearRegression())
    return (xgb.XGBClassifier if binary else xgb.XGBRegressor)(n_estimators=150, max_depth=3, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, n_jobs=1, verbosity=0)
def fit_job(jd):
    t = pd.read_csv(os.path.join(jd, "train.csv")); s = pd.read_csv(os.path.join(jd, "score.csv")); X, y = t.iloc[:, :-1], t.iloc[:, -1]; binary = y.nunique() == 2
    for tag in TAGS:
        out = os.path.join(jd, f"pred_{tag}.csv")
        if os.path.exists(out): continue
        try: m = model(tag, X, binary).fit(X, y); p = m.predict_proba(s[X.columns])[:, 1] if binary else m.predict(s[X.columns])
        except Exception as e: print("fail", jd, tag, e); continue
        pd.DataFrame({"p_1": p} if binary else {"prediction": p}).to_csv(out, index=False)
if __name__ == "__main__":
    jobs = pd.read_csv(os.path.join(H, "jobs.csv"))
    with Pool(int(os.environ.get("NPROC", 3))) as p: p.map(fit_job, [os.path.join(H, d) for d in jobs.dir], chunksize=4)
    print("done", len(jobs))
