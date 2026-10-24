#!/bin/bash
#SBATCH --account=def-kwiese_cpu
#SBATCH --job-name=PANELREF
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/scratch/YOUR_ACCOUNT/umbra_run/PANELREF_%j.out
set -uo pipefail
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true
PKG=${SCRATCH:?set SCRATCH to your own scratch directory}/umbra_pkg2
export PYTHONPATH="$PKG${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
echo "host $(hostname) date $(date -u +%FT%TZ)"
python3 "$PKG/cluster/panel_refcount.py" \
  --universe ${SCRATCH:?set SCRATCH to your own scratch directory}/addressability/GS_GRCh38.npz \
  --panel ${SCRATCH:?set SCRATCH to your own scratch directory}/recomb_mg/data/CHANGEseq_CR_Lazzarotto_2020_dataset.csv --changeseq \
  --out-json ${SCRATCH:?set SCRATCH to your own scratch directory}/umbra_run/PANELREF_GRCh38_NGG.json
echo "PANELREF DONE"
