#!/usr/bin/env bash
# Worker W of NW over folder 34's 600 units. Resumable (skips units with a results file).
set -u
ROOT="$(cd "$(dirname "$0")" && pwd)"; S34="$ROOT/../34_estimation_on_paper_data"; W=${1:-0}; NW=${2:-1}
truth() { "${PYTHON:-python3}" - "$1" "$2" "$ROOT" <<'PY'
import sys, pandas as pd
ds, t = sys.argv[1], sys.argv[2]; pop = "selected" if ds == "astronaut_set" else "source"
import os; v = pd.read_csv(os.path.join(sys.argv[3], "..", "..", os.environ.get("DATA_DIR", "02_data").removeprefix("./"), "truth_values.csv"))
print(v[(v.feature == t) & (v.truth_type == "per_unit") & (v.population == pop)].value.iloc[0])
PY
}
i=0
for F in $(ls "$S34"/reps/*/*/*/rep*/tmle_inputs.csv | sort); do
  i=$((i+1)); [ $((i % NW)) -eq "$W" ] || continue
  R=$(basename "$(dirname "$F")"); RT=$(basename "$(dirname "$(dirname "$F")")"); T=$(basename "$(dirname "$(dirname "$(dirname "$F")")")"); DS=$(basename "$(dirname "$(dirname "$(dirname "$(dirname "$F")")")")")
  O="$ROOT/results/$DS/$T/$RT/$R.csv"; [ -f "$O" ] && continue
  P34="$S34/results/$DS/$T/$RT/$R.csv"; TR=$(truth "$DS" "$T")
  Rscript "$ROOT/rerun_bounds.R" "$F" "$O" "$T" nephrolithiasis "$TR" "$P34" > /dev/null 2>&1 && echo "worker $W: $DS $T $RT $R $(date +%T)" || echo "worker $W: $DS $T $RT $R FAILED"
done
echo "worker $W finished"
