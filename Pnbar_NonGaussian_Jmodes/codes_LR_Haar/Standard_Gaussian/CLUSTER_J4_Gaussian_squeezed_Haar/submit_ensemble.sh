#!/bin/bash
# Submit Haar ensemble using N_ENSEMBLE / ENSEMBLE_BATCH_SIZE from config.py.
#
# Usage (from a CLUSTER_*_Haar folder):
#   mkdir -p logs output
#   ./submit_ensemble.sh
#
# Pilot → production: change ONLY config.N_ENSEMBLE (e.g. 10 → 1000), then re-run.
# Do NOT edit this script or the .slurm array range by hand.
#
# Note: OpenPBS rejects -J X-Y when X >= Y (e.g. -J 0-0). For a single batch
# we therefore submit without -J; PBS_ARRAY_INDEX defaults to 0 in the .slurm.

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
echo "  N_ENSEMBLE=${N_ENSEMBLE}  ENSEMBLE_BATCH_SIZE=${BATCH_SIZE}  n_batches=${N_BATCHES}"

if [ "${N_BATCHES}" -eq 1 ]; then
  # OpenPBS: -J 0-0 is illegal (start must be < end). One batch → plain job.
  echo "  mode=single-job (no PBS array; batch_index defaults to 0)"
  qsub submit_ensemble.slurm
else
  echo "  mode=PBS array 0-${LAST_INDEX}"
  qsub -J "0-${LAST_INDEX}" submit_ensemble.slurm
fi
