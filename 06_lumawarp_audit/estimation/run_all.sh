#!/usr/bin/env bash
# Engine driver plus three R workers, all in the background. Resumable. Logs: engine_log.txt, r_worker_<k>.txt
ROOT="$(cd "$(dirname "$0")" && pwd)"; cd "$ROOT"; mkdir -p reps results; rm -f reps/ENGINE_DONE
PY="${PYTHON:-python3}"
nohup env PLAN=${PLAN:-subset} "$PY" run_engine.py > engine_log.txt 2>&1 &
for k in 0 1 2; do nohup ./run_r_worker.sh $k 3 ENGINE_DONE_${PLAN:-subset} > r_worker_$k.txt 2>&1 & done
echo "launched"
