#!/usr/bin/env bash
# Three R workers over the 600 units (about 4 s per unit; about 15 min wall). Then: python summarize.py
cd "$(dirname "$0")" && for w in 0 1 2; do nohup ./run_worker.sh $w 3 > worker_$w.txt 2>&1 & done; wait; echo ALL WORKERS FINISHED
