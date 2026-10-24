#!/bin/bash
#SBATCH --account=def-kwiese_cpu
#SBATCH --job-name=UMBRA_REPRO
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/scratch/YOUR_ACCOUNT/umbra_run/REPRO_%j.out
set -uo pipefail

HAP="${1:?usage: sbatch reproduce.sh <hap_name>}"
echo "host $(hostname)  date $(date -u +%FT%TZ)  hap $HAP"
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true

PKG=${SCRATCH:?set SCRATCH to your own scratch directory}/umbra_pkg
ASM=${SCRATCH:?set SCRATCH to your own scratch directory}/hprc_r2/assemblies/${HAP}.fa.gz
UNI=${SCRATCH:?set SCRATCH to your own scratch directory}/addressability/GS_GRCh38.npz
PANEL=${SCRATCH:?set SCRATCH to your own scratch directory}/recomb_mg/data/CHANGEseq_CR_Lazzarotto_2020_dataset.csv
OUT=${SCRATCH:?set SCRATCH to your own scratch directory}/umbra_run/BLINDSPOT_${HAP}.json

for f in "$ASM" "$UNI" "$PANEL"; do
  if [ ! -s "$f" ]; then echo "GATE FAIL missing input $f" >&2; exit 10; fi
done

echo "PYTHONPATH as the modules left it: ${PYTHONPATH:-<unset>}"
echo "python $(python3 -c 'import sys;print(sys.version.split()[0])')  numpy $(python3 -c 'import numpy;print(numpy.__version__)')"
echo

cd "$PKG" || exit 1
export PYTHONPATH="$PKG${PYTHONPATH:+:$PYTHONPATH}"
echo "PYTHONPATH for the run: $PYTHONPATH"
python3 -c "import numpy; print('numpy still importable under the run PYTHONPATH:', numpy.__version__)" || exit 11
python3 -m scotoma.cli scan \
  --assembly "$ASM" \
  --universe "$UNI" \
  --panel "$PANEL" --changeseq \
  --keep-present 0.12 \
  --seed 20261019 \
  --out-json "$OUT"
status=$?
echo
echo "scan exit status $status"
echo "REPRO DONE $HAP $(date -u +%FT%TZ)"
exit $status
