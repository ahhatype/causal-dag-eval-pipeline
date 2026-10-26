#!/usr/bin/env bash
# Runs the hybrid and the audit (CPU) once the variant's data exist; the pipeline's own shape validation runs in parallel from run_variant.sh.
cd "$(dirname "$0")"; PY="${PYTHON:-python3}"
until grep -q "data generated" run_log.txt 2>/dev/null; do sleep 120; done
NREP=${NREP:-20} NPROC=4 $PY hybrid54.py > hybrid_log.txt 2>&1
$PY - <<'PY'
# write pc_skeletons.csv for the audit from the hybrid run: re-run PC quickly per replicate (cheap) to keep the two scripts independent
import os, numpy as np, pandas as pd
from graph54 import *
from causallearn.search.ConstraintBased.PC import pc
rows = []
for rep in range(1, int(os.environ.get("NREP", 20)) + 1):
    df = full_set(rep); X = df[NODES].to_numpy(float); keep = [i for i in range(len(NODES)) if X[:, i].std() > 0]; X = X[:, keep] + 1e-9 * np.random.default_rng(97000 + rep).standard_normal((len(X), len(keep)))
    G = pc(X, alpha=0.01, indep_test="fisherz", stable=True, show_progress=False).G.graph
    for i in range(len(keep)):
        for j in range(i + 1, len(keep)):
            if G[i, j] == 0 and G[j, i] == 0: continue
            rows.append(dict(rep=rep, a=NODES[keep[i]], b=NODES[keep[j]]))
pd.DataFrame(rows).to_csv("pc_skeletons.csv", index=False); print("pc skeletons written")
PY
NREP=${NREP:-20} NPROC=4 $PY audit54.py > audit_log.txt 2>&1
echo "HYBRID AND AUDIT FINISHED $(date +%T)" >> run_log.txt
