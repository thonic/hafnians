#!/bin/bash
#PBS -N pnbar_J4_Kerr_squeezed_kappa_1_LR
#PBS -j oe
#PBS -o logs/j4_kerr_squeezed_kappa_1_LR.log
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=8:mem=16GB

# Submit: mkdir -p logs output && qsub submit_J4_Kerr_squeezed_kappa_1_LR.slurm

set -euo pipefail
cd "${PBS_O_WORKDIR:-$(pwd)}"
mkdir -p logs output

NCPUS="${PBS_NCPUS:-${PBS_NP:-8}}"
PYTHON="${PYTHON:-/storage/brno12-cerit/home/deeptisharma/.conda/envs/guassian/bin/python}"

echo "HOST=$(hostname)  JOB=${PBS_JOBID:-local}  CPUS=${NCPUS}  ENGINE=low_rank  KAPPA=1.0"
"${PYTHON}" run_parallel_J4_Kerr_squeezed_kappa_1_LR.py --workers "${NCPUS}"
