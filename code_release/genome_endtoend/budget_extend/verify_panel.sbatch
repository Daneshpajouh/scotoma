#!/bin/bash
#SBATCH --job-name=verify-panel
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --array=0-7
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=2:00:00
#SBATCH --output=${PANEL_ROOT}/slurm-verify-%A_%a.out
: "${PANEL_ROOT:?set PANEL_ROOT to the directory holding these scripts}"
cd ${PANEL_ROOT}
python3 ${PANEL_ROOT}/verify_panel.py ${SLURM_ARRAY_TASK_ID}
