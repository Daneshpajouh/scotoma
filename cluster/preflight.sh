#!/bin/bash
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=UMBRA_PREFLIGHT
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/scratch/YOUR_ACCOUNT/umbra_run/PREFLIGHT_%j.out
set -uo pipefail

echo "host $(hostname)  date $(date -u +%FT%TZ)"
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true

PKG=${SCRATCH:?set SCRATCH to your own scratch directory}/umbra_pkg
cd "$PKG" || exit 1

echo "=== python and numpy ==="
python3 -c "import sys, numpy; print('python', sys.version.split()[0]); print('numpy', numpy.__version__)" || exit 2

echo
echo "=== unit tests on the cluster ==="
python3 -m pytest tests/ -q 2>&1 | tail -5

echo
echo "=== can umbra read the existing GRCh38 universe? ==="
python3 - <<'PY'
import os, sys
sys.path.insert(0, os.environ["SCRATCH"] + "/umbra_pkg")
from scotoma.universe import load_universe
codes, meta = load_universe(os.environ["SCRATCH"] + "/addressability/GS_GRCh38.npz")
print("codes", f"{codes.size:,}", "dtype", codes.dtype)
print("sorted and unique: enforced by load_universe (no exception raised)")
print("meta", meta)
assert codes.size == 226176633, f"GATE FAIL unexpected universe size {codes.size}"
print("GATE PASS universe size 226,176,633 matches the campaign's assertion")
PY

echo
echo "=== panel: 110 guides from the CHANGE-seq table? ==="
python3 - <<'PY'
import os, sys
sys.path.insert(0, os.environ["SCRATCH"] + "/umbra_pkg")
from scotoma.panel import load_panel
codes, rep = load_panel(os.environ["SCRATCH"] + "/recomb_mg/data/CHANGEseq_CR_Lazzarotto_2020_dataset.csv",
                        changeseq=True)
print(rep)
assert codes.size == 110, f"GATE FAIL expected 110 distinct guides, got {codes.size}"
print("GATE PASS 110 distinct guides")
PY

echo
echo "=== assembly staged? ==="
ls -la ${SCRATCH:?set SCRATCH to your own scratch directory}/hprc_r2/assemblies/HG00235_hap1_hprc_r2_v1.0.1.fa.gz
echo "PREFLIGHT DONE"
