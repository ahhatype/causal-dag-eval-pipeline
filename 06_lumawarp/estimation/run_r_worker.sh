#!/usr/bin/env bash
# R worker W of NW over every reps/<dataset>/<treat>/<route>/repNNN/tmle_inputs.csv; stops when both ENGINE_DONE markers exist... or PLAN marker given.
set -u
ROOT="$(cd "$(dirname "$0")" && pwd)"; W=${1:-0}; NW=${2:-1}; MARK=${3:-ENGINE_DONE_subset}
truth() { "${PYTHON:-python3}" - "$1" "$2" "$ROOT" <<'PY'
import sys, pandas as pd
ds, t = sys.argv[1], sys.argv[2]; pop = "selected" if ds == "astronaut_set" else "source"
import os; v = pd.read_csv(os.path.join(sys.argv[3], "..", "..", os.environ.get("DATA_DIR", "02_data").removeprefix("./"), "truth_values.csv"))
print(v[(v.feature == t) & (v.truth_type == "per_unit") & (v.population == pop)].value.iloc[0])
PY
}
while true; do
  did=0; i=0
  for F in $(ls "$ROOT"/reps/*/*/*/rep*/tmle_inputs.csv 2>/dev/null | sort); do
    i=$((i+1)); [ $((i % NW)) -eq "$W" ] || continue
    R=$(basename "$(dirname "$F")"); RT=$(basename "$(dirname "$(dirname "$F")")"); T=$(basename "$(dirname "$(dirname "$(dirname "$F")")")"); DS=$(basename "$(dirname "$(dirname "$(dirname "$(dirname "$F")")")")")
    O="$ROOT/results/$DS/$T/$RT/$R.csv"; [ -f "$O" ] && continue
    Rscript "$ROOT/tmle_rep.R" "$F" "$O" "$T" nephrolithiasis "$(truth $DS $T)" > "$(dirname "$F")/tmle.log" 2>&1 && echo "worker $W: $DS $T $RT $R $(date +%T)" || echo "worker $W: $DS $T $RT $R FAILED"; did=1
  done
  if [ "$did" -eq 0 ]; then [ -f "$ROOT/reps/$MARK" ] && { echo "worker $W finished"; exit 0; }; sleep 30; fi
done
