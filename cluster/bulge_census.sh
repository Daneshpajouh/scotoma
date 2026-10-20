#!/bin/bash
#SBATCH --account=def-kwiese_cpu
#SBATCH --job-name=BULGECEN
#SBATCH --time=11:50:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96000M
#SBATCH --output=/scratch/amird/umbra_run/BULGECEN_%j.out
set -uo pipefail
ENV=/scratch/amird/casoffinder_env
export PATH="$ENV/bin:$PATH"
export LD_LIBRARY_PATH="$ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true
export PYTHONUNBUFFERED=1
W=/scratch/amird/bulge_census
REF=/scratch/amird/cof_h2h/ref           # the 25 primary FASTA files, already prepared and verified
PANEL=/scratch/amird/umbra_run/PANELREF_GRCh38_NGG.json
MM=4                                      # the paper's headline threshold
mkdir -p "$W"
echo "host $(hostname)  date $(date -u +%FT%TZ)"
ls "$REF" | wc -l
for KIND in dna rna; do
  echo "===== $KIND bulge, mm<=$MM ====="
  python3 /scratch/amird/umbra_run/bulge_expand.py "$REF" $KIND "$PANEL" \
          "$W/in_$KIND.txt" $MM "$W/idx_$KIND.json" || exit 1
  date -u +%FT%TZ
  cas-offinder "$W/in_$KIND.txt" C "$W/out_$KIND.txt" < /dev/null
  echo "rc=$?"
  date -u +%FT%TZ
  ls -l "$W/out_$KIND.txt"
done
echo "===== deduplicate and count additional candidates ====="
python3 /scratch/amird/umbra_run/bulge_summarise.py $MM
echo "BULGECEN DONE"
