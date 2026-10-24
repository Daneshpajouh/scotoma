#!/bin/bash
#SBATCH --account=def-kwiese_cpu
#SBATCH --job-name=BULGETOY
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --output=/scratch/YOUR_ACCOUNT/umbra_run/BULGETOY_%j.out
set -uo pipefail
ENV=${SCRATCH:?set SCRATCH to your own scratch directory}/casoffinder_env
export PATH="$ENV/bin:$PATH"
export LD_LIBRARY_PATH="$ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true
export PYTHONUNBUFFERED=1
W=${SCRATCH:?set SCRATCH to your own scratch directory}/bulge_toy
mkdir -p "$W/ref"
GUIDE=GAACAAGGTGTTCCCACCCG
echo "guide $GUIDE"
# planted: exact on-target (23), DNA-bulge site (24, extra A after base 10),
# RNA-bulge site (22, guide base 10 'T' absent), each separated by poly-T spacers
python3 - <<'PYEOF'
g = "GAACAAGGTGTTCCCACCCG"
exact = g + "AGG"
dna   = g[:10] + "A" + g[10:] + "AGG"       # target has one extra base
rna   = g[:10] + g[11:] + "AGG"             # target lacks guide base 10
pad = "T" * 40
seq = pad + exact + pad + dna + pad + rna + pad
with open("${SCRATCH:?set SCRATCH to your own scratch directory}/bulge_toy/ref/toy.fa", "w") as f:
    f.write(">toy\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")
print(f"planted exact={exact} ({len(exact)})")
print(f"planted dna  ={dna} ({len(dna)})  at offset {len(pad)*2+len(exact)}")
print(f"planted rna  ={rna} ({len(rna)})")
PYEOF
echo "$GUIDE" > "$W/panel1.txt"
for KIND in none dna rna; do
  python3 ${SCRATCH:?set SCRATCH to your own scratch directory}/umbra_run/bulge_expand.py "$W/ref" $KIND "$W/panel1.txt" \
          "$W/in_$KIND.txt" 0 "$W/idx_$KIND.json"
  timeout 300 cas-offinder "$W/in_$KIND.txt" C "$W/out_$KIND.txt" < /dev/null > /dev/null 2>&1
  N=$(wc -l < "$W/out_$KIND.txt")
  echo "--- $KIND (mm<=0) rc=$? hits=$N ---"
  cat "$W/out_$KIND.txt"
  case $KIND in
    none) [ "$N" -eq 1 ] && echo "GATE none: PASS (exact on-target only)" || echo "GATE none: FAIL" ;;
    dna)  [ "$N" -eq 1 ] && echo "GATE dna:  PASS (planted DNA bulge only, no end artefact)" || echo "GATE dna:  FAIL" ;;
    rna)  [ "$N" -eq 1 ] && echo "GATE rna:  PASS (planted RNA bulge only, no end artefact)" || echo "GATE rna:  FAIL" ;;
  esac
done
echo "BULGETOY DONE"
