#!/bin/bash
#PBS -N haar_ens_J4_superposition_LR_Haar
#PBS -j oe
# Per-array-index PBS spool (OpenPBS expands ^array_index^ for -J jobs).
# Non-array (pilot) submits: see exec redirect below → ensemble_batch_0.log
#PBS -o logs/ensemble_batch.^array_index^.log
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=8:mem=16GB

# Array range is set by ./submit_ensemble.sh from config.py
# (qsub -J 0-<last>). Do NOT hard-code #PBS -J here.
# Pilot→production: change only N_ENSEMBLE in config.py, then ./submit_ensemble.sh
#
# Logging: each array task writes logs/ensemble_batch_<index>.log
# (diagnostics only; ensemble computation unchanged).

set -euo pipefail
cd "${PBS_O_WORKDIR:-$(pwd)}"
mkdir -p logs output

# Canonical per-task log for both array and single-job (no -J) submits.
ARRAY_IDX="${PBS_ARRAY_INDEX:-0}"
exec >"logs/ensemble_batch_${ARRAY_IDX}.log" 2>&1

NCPUS="${PBS_NCPUS:-${PBS_NP:-8}}"
PYTHON="${PYTHON:-/storage/brno12-cerit/home/deeptisharma/.conda/envs/guassian/bin/python}"
BATCH_INDEX="${ARRAY_IDX}"

echo "HOST=$(hostname)  JOB=${PBS_JOBID:-local}  ARRAY=${BATCH_INDEX}  CPUS=${NCPUS}  ENGINE=haar_ensemble"
"${PYTHON}" run_ensemble.py --batch-index "${BATCH_INDEX}" --workers "${NCPUS}"
