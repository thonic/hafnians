#!/bin/bash
#PBS -N pnbar_J4_Gaussian_squeezed_LR
#PBS -j oe
#PBS -o logs/j4_gaussian_squeezed_LR.log
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=8:mem=16GB

# Submit: mkdir -p logs output && qsub submit_J4_Gaussian_squeezed_LR.slurm

set -euo pipefail
cd "${PBS_O_WORKDIR:-$(pwd)}"
mkdir -p logs output

NCPUS="${PBS_NCPUS:-${PBS_NP:-8}}"
PYTHON="${PYTHON:-/storage/brno12-cerit/home/deeptisharma/.conda/envs/guassian/bin/python}"

echo "HOST=$(hostname)  JOB=${PBS_JOBID:-local}  CPUS=${NCPUS}  ENGINE=thewalrus_low_rank_hafnian"
"${PYTHON}" run_parallel_J4_Gaussian_squeezed_LR.py --workers "${NCPUS}"
