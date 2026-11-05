#!/bin/bash
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=panel-extend
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=120G --time=06:00:00
#SBATCH --output=${PANEL_ROOT}/slurm-%x-%A_%a.out
: "${PANEL_ROOT:?set PANEL_ROOT to the directory holding these scripts}"
set -euo pipefail
export NSHARDS=8
exec srun --ntasks=1 --cpu-bind=cores timeout --kill-after=60s 21000s \
  python3 -u ${PANEL_ROOT}/panel_extend.py
