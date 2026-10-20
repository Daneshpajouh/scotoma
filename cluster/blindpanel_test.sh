#!/bin/bash
#SBATCH --account=def-kwiese_cpu
#SBATCH --job-name=BLINDPAN
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --output=/scratch/amird/umbra_run/BLINDPANEL_%j.out
set -uo pipefail
HAP="${1:?hap}"
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true
PKG=/scratch/amird/umbra_pkg2
export PYTHONPATH="$PKG${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
echo "host $(hostname) hap $HAP date $(date -u +%FT%TZ)"
python3 -m scotoma.cli scan \
  --assembly /scratch/amird/hprc_r2/assemblies/${HAP}.fa.gz \
  --universe /scratch/amird/addressability/GS_GRCh38.npz \
  --panel /scratch/amird/recomb_mg/data/CHANGEseq_CR_Lazzarotto_2020_dataset.csv --changeseq \
  --keep-present 1.0 \
  --panel-from-blind \
  --seed 20261020 \
  --assume-universe-pam \
  --out-json /scratch/amird/umbra_run/BLINDPANEL_${HAP}.json
echo "BLINDPANEL DONE $HAP"
