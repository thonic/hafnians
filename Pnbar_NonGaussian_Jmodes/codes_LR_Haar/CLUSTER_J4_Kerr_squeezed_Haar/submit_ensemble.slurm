#!/bin/bash
#PBS -N haar_ens_J4_Kerr_squeezed_Haar
#PBS -j oe
#PBS -o logs/ensemble_batch.log
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=8:mem=16GB

# Array range is set by ./submit_ensemble.sh from config.py
# (qsub -J 0-<last>). Do NOT hard-code #PBS -J here.
# Pilot→production: change only N_ENSEMBLE in config.py, then ./submit_ensemble.sh

set -euo pipefail
cd "${PBS_O_WORKDIR:-$(pwd)}"
mkdir -p logs output

NCPUS="${PBS_NCPUS:-${PBS_NP:-8}}"
PYTHON="${PYTHON:-/storage/brno12-cerit/home/deeptisharma/.conda/envs/guassian/bin/python}"
BATCH_INDEX="${PBS_ARRAY_INDEX:-0}"

echo "HOST=$(hostname)  JOB=${PBS_JOBID:-local}  ARRAY=${BATCH_INDEX}  CPUS=${NCPUS}  ENGINE=haar_ensemble"
"${PYTHON}" run_ensemble.py --batch-index "${BATCH_INDEX}" --workers "${NCPUS}"
