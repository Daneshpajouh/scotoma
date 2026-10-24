#!/bin/bash
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=G38DECOMP
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --output=/scratch/YOUR_ACCOUNT/umbra_run/G38DECOMP_%j.out
set -uo pipefail
# The reverse of job 66001436. chm13_decompose.py is symmetric in its two universe
# arguments: it decomposes the universe given as --chm13 (loaded with multiplicity)
# by membership in the universe given as --grch38. Passing them SWAPPED therefore
# decomposes GRCh38's candidates into
#   shared : code also present in T2T-CHM13
#   novel  : code ABSENT from T2T-CHM13, i.e. the GRCh38-only candidate sites
# which is the term job 66001436 could not separate from copy-number change.
# The argument names in the output are therefore inverted with respect to their
# meaning; the file name and this comment are the record of that.
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true
PKG=${SCRATCH:?set SCRATCH to your own scratch directory}/umbra_pkg2
export PYTHONPATH="$PKG${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
echo "host $(hostname) date $(date -u +%FT%TZ)"
echo "SWAPPED: --chm13 is GS_GRCh38.npz and --grch38 is GS_CHM13.npz, by design"
python3 "$PKG/cluster/chm13_decompose.py" \
  --chm13 ${SCRATCH:?set SCRATCH to your own scratch directory}/addressability/GS_GRCh38.npz \
  --grch38 ${SCRATCH:?set SCRATCH to your own scratch directory}/addressability/GS_CHM13.npz \
  --panel ${SCRATCH:?set SCRATCH to your own scratch directory}/recomb_mg/data/CHANGEseq_CR_Lazzarotto_2020_dataset.csv --changeseq \
  --out-json ${SCRATCH:?set SCRATCH to your own scratch directory}/umbra_run/G38DECOMP_grch38_only.json
echo "G38DECOMP DONE"
