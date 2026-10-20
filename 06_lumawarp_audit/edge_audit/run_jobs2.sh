#!/usr/bin/env bash
# Two-stream variant of run_jobs.sh for the full-set stage: DIR=fwd|rev; a job is claimed by creating train_log.txt (noclobber), so the two streams never fit the same job.
set -u
TRAIN="${LUMAWARP_TRAIN:-lumawarp_training}"; INFER="${LUMAWARP_INFER:-lumawarp_inference}"; ROOT="$(cd "$(dirname "$0")" && pwd)"
PLAN="${ONLY:-full_set}"; FAILED=""; N=0; T0=$(date +%s); S="${DIR:-fwd}"
for PD in $PLAN; do
  if [ "$S" = rev ]; then L=$(find "$ROOT/$PD" -name train.csv | sort -r); else L=$(find "$ROOT/$PD" -name train.csv | sort); fi
  for T in $L; do
    JD="$(dirname "$T")"; [ -f "$JD/pred.csv" ] && continue
    ( set -o noclobber; : > "$JD/train_log.txt" ) 2>/dev/null || continue
    cd "$JD"; N=$((N+1)); echo "=== [$S $N] $(date +%H:%M:%S) +$(( $(date +%s) - T0 ))s ${JD#$ROOT/}"
    if "$TRAIN" train.csv model.lsmdl > train_log.txt 2>&1; then "$INFER" score.csv model.lsmdl pred.csv > infer_log.txt 2>&1 || { rm -f pred.csv; FAILED="$FAILED ${JD#$ROOT/}"; }; else FAILED="$FAILED ${JD#$ROOT/}"; fi
    rm -f model.lsmdl.diag.jsonl pred.csv.diag.jsonl axes_out_*.csv model.lsmdl.ckpt
  done
done
echo "==================================================================="
[ -n "$FAILED" ] && echo "COMPLETED WITH FAILURES ($S):$FAILED" || echo "ALL RUNS COMPLETED ($S: $N new, $(( $(date +%s) - T0 ))s)."
