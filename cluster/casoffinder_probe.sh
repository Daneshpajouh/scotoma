#!/bin/bash
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=COFPROBE
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scratch/YOUR_ACCOUNT/umbra_run/COFPROBE_%j.out
set -uo pipefail
ENV=${SCRATCH:?set SCRATCH to your own scratch directory}/casoffinder_env
export PATH="$ENV/bin:$PATH"
export LD_LIBRARY_PATH="$ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OCL_ICD_VENDORS="$ENV/etc/OpenCL/vendors"
echo "host $(hostname)  date $(date -u +%FT%TZ)"
echo "=============== 1. binary identity ==============="
ls -l "$ENV/bin/cas-offinder"
echo "--- conda-meta record ---"
ls "$ENV/conda-meta" | head -20
echo "=============== 2. usage / version ==============="
cas-offinder 2>&1 | head -40 || echo "exit $?"
echo "=============== 3. OpenCL devices ==============="
ls "$ENV/etc/OpenCL/vendors" 2>&1
cas-offinder --help 2>&1 | head -20 || true
echo "=============== 4. candidate GRCh38 FASTA on this filesystem ==============="
for p in ${SCRATCH:?set SCRATCH to your own scratch directory}/addressability ${SCRATCH:?set SCRATCH to your own scratch directory}/hprc_r2 ${SCRATCH:?set SCRATCH to your own scratch directory}/references \
         ${SCRATCH:?set SCRATCH to your own scratch directory}/ot_data ${SCRATCH:?set SCRATCH to your own scratch directory}/validation_data ${SCRATCH:?set SCRATCH to your own scratch directory}/recomb_mg/data; do
  echo "--- $p"
  ls "$p" 2>/dev/null | grep -iE 'fa$|fa\.gz$|fasta|2bit|GRCh38|hg38|chm13' | head -12
done
echo "=============== 5. cvmfs reference grids ==============="
ls /cvmfs/ref.mugqic/genomes/species 2>/dev/null | head -5
ls /cvmfs/soft.computecanada.ca/custom/ref 2>/dev/null | head -5
echo "COFPROBE DONE"
