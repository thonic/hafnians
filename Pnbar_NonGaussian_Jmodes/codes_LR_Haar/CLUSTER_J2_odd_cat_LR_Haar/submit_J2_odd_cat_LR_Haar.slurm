#!/bin/bash
#PBS -N pnbar_J2_odd_cat_LR_Haar
#PBS -j oe
#PBS -o logs/j2_odd_cat_LR_Haar.log
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=8:mem=8GB

# Submit: mkdir -p logs output && qsub submit_J2_odd_cat_LR_Haar.slurm

set -euo pipefail
cd "${PBS_O_WORKDIR:-$(pwd)}"
mkdir -p logs output

NCPUS="${PBS_NCPUS:-${PBS_NP:-8}}"
PYTHON="${PYTHON:-/storage/brno12-cerit/home/deeptisharma/.conda/envs/guassian/bin/python}"

echo "HOST=$(hostname)  JOB=${PBS_JOBID:-local}  CPUS=${NCPUS}  ENGINE=low_rank"
"${PYTHON}" run_parallel_J2_odd_cat_LR_Haar.py --workers "${NCPUS}"
