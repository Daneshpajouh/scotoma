#!/bin/bash
#SBATCH --account=def-kwiese_cpu
#SBATCH --job-name=HAPABS
#SBATCH --time=01:40:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=/scratch/amird/umbra_run/hapabs/HAPABS_%A_%a.out
set -uo pipefail
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true
PKG=/scratch/amird/umbra_pkg2
# Append, never assign: scipy-stack supplies numpy through PYTHONPATH.
export PYTHONPATH="$PKG${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
echo "host $(hostname)  task ${SLURM_ARRAY_TASK_ID}  date $(date -u +%FT%TZ)"
python3 "$PKG/cluster/hap_absence.py" \
  --manifest /scratch/amird/umbra_run/hap_manifest.tsv \
  --asm-dir /scratch/amird/hprc_r2/assemblies \
  --universe /scratch/amird/addressability/GS_GRCh38.npz \
  --outdir /scratch/amird/umbra_run/hapabs/results \
  --pam NGG
echo "exit $?"
