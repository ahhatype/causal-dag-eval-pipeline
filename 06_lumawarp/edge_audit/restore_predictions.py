"""Writes the archived LumaWarp sweep predictions (outputs/06_lumawarp/edge_audit/lumawarp_predictions.csv.gz) back into the
job directories as pred.csv, so score.py reproduces the LumaWarp rows without the engine. Run after make_jobs.py."""
import os, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
A = pd.read_csv(os.path.join(H, "..", "..", "outputs", "06_lumawarp", "edge_audit", "lumawarp_predictions.csv.gz")); n = 0
for d, g in A.groupby("dir", sort=False):
    jd = os.path.join(H, d)
    if not os.path.isdir(jd): continue
    target = pd.read_csv(os.path.join(jd, "train.csv")).iloc[:, -1]; binary = set(target.unique()) <= {0, 1}
    pd.DataFrame({"p_1" if binary else "prediction": g.sort_values("row").prediction.to_numpy()}).to_csv(os.path.join(jd, "pred.csv"), index=False); n += 1
print(n, "pred.csv files restored")
