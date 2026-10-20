#!/bin/bash
#SBATCH --account=def-kwiese_cpu
#SBATCH --job-name=CHMDECOMP
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --output=/scratch/amird/umbra_run/CHMDECOMP_%j.out
set -uo pipefail
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true
PKG=/scratch/amird/umbra_pkg2
export PYTHONPATH="$PKG${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
echo "host $(hostname) date $(date -u +%FT%TZ)"
python3 "$PKG/cluster/chm13_decompose.py" \
  --chm13 /scratch/amird/addressability/GS_CHM13.npz \
  --grch38 /scratch/amird/addressability/GS_GRCh38.npz \
  --panel /scratch/amird/recomb_mg/data/CHANGEseq_CR_Lazzarotto_2020_dataset.csv --changeseq \
  --out-json /scratch/amird/umbra_run/CHMDECOMP.json
echo "CHMDECOMP DONE"
