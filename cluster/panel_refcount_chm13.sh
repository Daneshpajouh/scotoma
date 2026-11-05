#!/bin/bash
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=PANELCHM
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/scratch/YOUR_ACCOUNT/umbra_run/PANELCHM_%j.out
set -uo pipefail
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true
PKG=${SCRATCH:?set SCRATCH to your own scratch directory}/umbra_pkg2
export PYTHONPATH="$PKG${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
echo "host $(hostname) date $(date -u +%FT%TZ)"
python3 "$PKG/cluster/panel_refcount.py" \
  --universe ${SCRATCH:?set SCRATCH to your own scratch directory}/addressability/GS_CHM13.npz \
  --panel ${SCRATCH:?set SCRATCH to your own scratch directory}/recomb_mg/data/CHANGEseq_CR_Lazzarotto_2020_dataset.csv --changeseq \
  --out-json ${SCRATCH:?set SCRATCH to your own scratch directory}/umbra_run/PANELREF_CHM13_NGG.json
echo "PANELCHM DONE"
