#!/bin/bash
#SBATCH --account=def-kwiese_cpu
#SBATCH --job-name=COFPROBE2
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --output=/scratch/amird/umbra_run/COFPROBE2_%j.out
set -uo pipefail
ENV=/scratch/amird/casoffinder_env
export PATH="$ENV/bin:$PATH"
export LD_LIBRARY_PATH="$ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
echo "host $(hostname)  date $(date -u +%FT%TZ)"
echo "=== vendors dir ==="
ls -l "$ENV/etc/OpenCL/vendors" || echo "NO VENDORS DIR"
echo "=== pocl present? ==="
ls "$ENV/conda-meta" | grep -i pocl || echo "no pocl conda record"
ls "$ENV/lib" | grep -iE 'pocl|OpenCL' || echo "no pocl/OpenCL libs"
echo "=== usage, stdin closed, hard timeout ==="
timeout 30 cas-offinder < /dev/null 2>&1 | head -30; echo "rc=${PIPESTATUS[0]}"
echo "=== device list (C = CPU) with a 1-guide toy input ==="
WORK=$(mktemp -d)
printf '%s\n' "$WORK/ref" > "$WORK/in.txt"
echo "NNNNNNNNNNNNNNNNNNNNNRG" >> "$WORK/in.txt"
echo "GAACAAGGTGTTCCCACCCGNNN 3" >> "$WORK/in.txt"
mkdir -p "$WORK/ref"
printf '>toy\n' > "$WORK/ref/toy.fa"
printf 'GAACAAGGTGTTCCCACCCGAGGTTTTTTTTTTGAACAAGGTGTTCCCACCCGCGGAAAA\n' >> "$WORK/ref/toy.fa"
cat "$WORK/in.txt"
timeout 120 cas-offinder "$WORK/in.txt" C "$WORK/out.txt" < /dev/null 2>&1 | head -20; echo "rc=${PIPESTATUS[0]}"
echo "--- toy output ---"
cat "$WORK/out.txt" 2>/dev/null || echo "no output file"
echo "=== GRCh38 fasta candidates ==="
ls /scratch/amird/addressability | grep -iE 'fa$|fasta|2bit|GRCh38|hg38' | head
ls /scratch/amird/references 2>/dev/null | head -15
ls /scratch/amird/hprc_r2 2>/dev/null | head -15
echo "COFPROBE2 DONE"
