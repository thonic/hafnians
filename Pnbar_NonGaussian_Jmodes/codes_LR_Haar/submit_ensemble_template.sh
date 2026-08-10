#!/bin/bash
# Submit Haar ensemble using N_ENSEMBLE / ENSEMBLE_BATCH_SIZE from config.py.
#
# Usage (from a CLUSTER_*_Haar folder):
#   mkdir -p logs output
#   ./submit_ensemble.sh
#
# Pilot → production: change ONLY config.N_ENSEMBLE (e.g. 10 → 1000), then re-run.
# Do NOT edit this script or the .slurm array range by hand.

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs output

PYTHON="${PYTHON:-/storage/brno12-cerit/home/deeptisharma/.conda/envs/guassian/bin/python}"

# Read array bounds from config.py (single source of truth)
read -r N_ENSEMBLE BATCH_SIZE N_BATCHES LAST_INDEX < <(
  "${PYTHON}" - <<'PY'
import math
import config
n = int(config.N_ENSEMBLE)
b = int(config.ENSEMBLE_BATCH_SIZE)
if b <= 0:
    raise SystemExit("ENSEMBLE_BATCH_SIZE must be positive")
n_batches = max(1, int(math.ceil(n / b)))
last = n_batches - 1
print(n, b, n_batches, last)
PY
)

echo "Submitting Haar ensemble from $(pwd)"
echo "  N_ENSEMBLE=${N_ENSEMBLE}  ENSEMBLE_BATCH_SIZE=${BATCH_SIZE}  PBS array 0-${LAST_INDEX} (${N_BATCHES} tasks)"

qsub -J "0-${LAST_INDEX}" submit_ensemble.slurm
