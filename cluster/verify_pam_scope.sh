#!/bin/bash
#SBATCH --account=def-kwiese_cpu
#SBATCH --job-name=UMBRA_PAMSCOPE
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=/scratch/amird/umbra_run/PAMSCOPE_%j.out
set -uo pipefail
echo "host $(hostname)  date $(date -u +%FT%TZ)"
module load StdEnv/2023 2>/dev/null || true
module load scipy-stack 2>/dev/null || true
export PYTHONPATH="/scratch/amird/umbra_pkg2${PYTHONPATH:+:$PYTHONPATH}"
echo "PYTHONPATH $PYTHONPATH"
python3 - <<'PY'
import numpy as np
from scotoma.coding import oriented_sites, canonical
from scotoma.universe import load_universe, contains

# chr20 is small enough to scan quickly and large enough to be decisive.
SEQ = "/scratch/amird/saif_null/restore/r83_saif_probe/ref/chroms/chr20.seq"
raw = open(SEQ, "rb").read()
# .seq files hold bare sequence; strip any whitespace/newlines.
raw = b"".join(raw.split())
bases = np.frombuffer(raw, dtype=np.uint8)
print(f"chr20 bases {bases.size:,}")

uni, meta = load_universe("/scratch/amird/addressability/GS_GRCh38.npz")
print("universe meta pam key:", meta.get("pam", "<ABSENT>"))

narrow = oriented_sites(bases, pam="NGG")
wide   = oriented_sites(bases, pam="NGG+NAG")
print(f"chr20 NGG sites      {narrow.size:,}")
print(f"chr20 NGG+NAG sites  {wide.size:,}   ratio {wide.size/narrow.size:.3f}")

# NAG-only sites: in the wide set, not in the narrow set (compare canonically, since the
# universe is canonical and a site may appear in either orientation).
cn = set(np.unique(canonical(narrow)).tolist())
cw = np.unique(canonical(wide))
nag_only = np.array([c for c in cw.tolist() if c not in cn], dtype=np.uint64)
print(f"chr20 NAG-only distinct canonical codes {nag_only.size:,}")

frac_ngg = float(contains(uni, narrow).mean())
frac_nag = float(contains(uni, nag_only).mean()) if nag_only.size else float("nan")
print(f"\nfraction of chr20 NGG sites present in GS_GRCh38:      {frac_ngg:.6f}")
print(f"fraction of chr20 NAG-ONLY sites present in GS_GRCh38: {frac_nag:.6f}")

print()
if frac_ngg > 0.999 and frac_nag < 0.001:
    print("VERDICT: GS_GRCh38.npz is NGG-ONLY. It contains essentially all NGG sites and")
    print("essentially no NAG-only sites. Scotoma's default NGG scope matches it, and the")
    print("campaign's 8-donor measurement is an NGG-scope result.")
elif frac_ngg > 0.999 and frac_nag > 0.999:
    print("VERDICT: GS_GRCh38.npz includes NAG. Scotoma must be run with --pam NGG+NAG against it.")
else:
    print("VERDICT: INCONCLUSIVE. Neither pattern holds; do not assume a scope.")
PY
echo "PAMSCOPE DONE"
