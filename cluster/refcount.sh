#!/bin/bash
#SBATCH --account=def-kwiese_cpu
#SBATCH --job-name=REFCOUNT
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --output=/scratch/amird/umbra_run/REFCOUNT_%j.out
set -uo pipefail
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true
PKG=/scratch/amird/umbra_pkg2
export PYTHONPATH="$PKG${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
echo "host $(hostname)  date $(date -u +%FT%TZ)"
echo
echo "### exa-cel / BCL11A +58 enhancer guide, GRCh38, NGG scope"
python3 -m scotoma.cli reference-count \
  --universe /scratch/amird/addressability/GS_GRCh38.npz \
  --spacer CTAACAGTTGCTTTTATCAC \
  --assume-universe-pam \
  --out-json /scratch/amird/umbra_run/REFCOUNT_exacel_NGG.json
echo
echo "### same guide, T2T-CHM13, NGG scope, for the reference-choice comparison"
python3 -m scotoma.cli reference-count \
  --universe /scratch/amird/addressability/GS_CHM13.npz \
  --spacer CTAACAGTTGCTTTTATCAC \
  --assume-universe-pam \
  --out-json /scratch/amird/umbra_run/REFCOUNT_exacel_CHM13_NGG.json
echo "REFCOUNT DONE"
