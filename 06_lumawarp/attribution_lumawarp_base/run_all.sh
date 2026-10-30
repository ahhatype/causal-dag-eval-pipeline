#!/usr/bin/env bash
# Resumable. Knobs (env): REPS=10 DATASETS=full_set,reference_subsample,astronaut_set METHODS=interventional_shap,standard_shap,ordering_only_shap BACKGROUND=128 ORDERS=64 NPAR=2
# Full budget (the paper's): about 2.3 h per cell -> 30 cells about 70 h.   Reduced: BACKGROUND=32 -> about 18 h; BACKGROUND=32 ORDERS=32 -> about 9 h.
cd "$(dirname "$0")" && export OMP_NUM_THREADS=1 && exec "${PYTHON:-python3}" run_cells.py
