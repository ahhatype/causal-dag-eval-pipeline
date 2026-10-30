#!/usr/bin/env bash
# Old build (format-1 models kept for the shape flag), research config. Order: 900-row sets first, then full-set. Resumable.
set -u
TRAIN="${LUMAWARP_TRAIN:-lumawarp_training}"; INFER="${LUMAWARP_INFER:-lumawarp_inference}"; ROOT="$(cd "$(dirname "$0")" && pwd)"
PLAN="${ONLY:-reference_subsample astronaut_set full_set}"; FAILED=""; N=0; T0=$(date +%s)
for PD in $PLAN; do
  for T in $(find "$ROOT/$PD" -name train.csv | sort); do
    JD="$(dirname "$T")"; [ -f "$JD/pred.csv" ] && continue
    cd "$JD"; N=$((N+1)); echo "=== [$N] $(date +%H:%M:%S) +$(( $(date +%s) - T0 ))s ${JD#$ROOT/}"
    if "$TRAIN" train.csv model.lsmdl > train_log.txt 2>&1; then "$INFER" score.csv model.lsmdl pred.csv > infer_log.txt 2>&1 || { rm -f pred.csv; FAILED="$FAILED ${JD#$ROOT/}"; }; else FAILED="$FAILED ${JD#$ROOT/}"; fi
    rm -f model.lsmdl.diag.jsonl pred.csv.diag.jsonl axes_out_*.csv model.lsmdl.ckpt
  done
done
echo "==================================================================="
[ -n "$FAILED" ] && echo "COMPLETED WITH FAILURES:$FAILED" || echo "ALL RUNS COMPLETED ($N new, $(( $(date +%s) - T0 ))s)."
