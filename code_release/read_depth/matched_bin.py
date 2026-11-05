"""Read depth of the 129 exclusively-missed sites against the spared set, INSIDE the
four-mismatch bin. Exact Mann-Whitney U with tie correction; gate on the published set sizes."""
import csv, collections, sys, os
from fractions import Fraction as F
from math import sqrt, erfc

DEP = os.environ["CHANGESEQ_CSV"]
N = 20
csv.field_size_limit(10**7)

def mmd(a, b): return sum(1 for x, y in zip(a, b) if x != y)

at_risk, rescued = [], []
with open(DEP, newline="") as fh:
    for r in csv.DictReader(fh):
        try:
            reads = float((r.get("reads") or "0") or 0)
            nb = int(float((r.get("Align.#Bulges") or "0").strip()))
            nm = int(float((r.get("Align.#Mismatches") or "0").strip()))
        except (TypeError, ValueError):
            continue
        if reads <= 0 or nb != 1 or nm > 4: continue
        g = (r.get("Align.sgRNA") or "").strip().upper()
        t = (r.get("Align.off-target") or "").strip().upper()
        if t.count("-") != 1 or g.count("-") != 0: continue
        if len(g[: t.index("-")].replace("-", "")) != N - 1: continue
        sp = g.replace("-", "")[:N]
        rec = {"reads": reads, "mm": nm, "spacer": sp, "core": t.replace("-", "")[:N - 1]}
        (rescued if sp[N - 1] == sp[N - 2] else at_risk).append(rec)

# exclusively missed: no EMITTED pattern (deletions at positions 1..n-2) within 4 mismatches
def emitted(sp):
    return {sp[:i] + sp[i + 1:] for i in range(1, N - 1)}

excl = [x for x in at_risk if min((mmd(x["core"], p) for p in emitted(x["spacer"])), default=99) > 4]

if (len(at_risk), len(rescued), len(excl)) != (260, 323, 129):
    sys.exit(f"GATE FAILED: {len(at_risk)}, {len(rescued)}, {len(excl)} != 260, 323, 129")
print("GATE PASSED: at_risk 260, rescued 323, exclusive 129")
print("exclusive mm:", dict(sorted(collections.Counter(x['mm'] for x in excl).items())))
print("spared mm:   ", dict(sorted(collections.Counter(x['mm'] for x in rescued).items())))

def mw(a, b):
    """Two-sided Mann-Whitney U, normal approximation with tie and continuity correction."""
    allv = sorted(a + b)
    rank = {}
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1] == allv[i]: j += 1
        r = F(i + j + 2, 2)
        for k in range(i, j + 1): rank[allv[k]] = r
        i = j + 1
    n1, n2 = len(a), len(b)
    R1 = sum(rank[v] for v in a)
    U1 = R1 - F(n1 * (n1 + 1), 2)
    mu = F(n1 * n2, 2)
    ties = collections.Counter(allv)
    n = n1 + n2
    tie = sum(t ** 3 - t for t in ties.values())
    var = F(n1 * n2, 12) * (F(n + 1) - F(tie, n * (n - 1)))
    z = (abs(float(U1 - mu)) - 0.5) / sqrt(float(var))
    return float(U1), erfc(z / sqrt(2))

def med(v):
    s = sorted(v); n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2

ex_r = [x["reads"] for x in excl]
sp_all = [x["reads"] for x in rescued]
sp_4 = [x["reads"] for x in rescued if x["mm"] == 4]

print()
print("AS PUBLISHED, spared set unmatched on mismatch count:")
U, p = mw(ex_r, sp_all)
print(f"  n = {len(ex_r)} vs {len(sp_all)};  median {med(ex_r):g} vs {med(sp_all):g};  U = {U:g};  p = {p:.6f}")
print()
print("MISMATCH-MATCHED, spared set restricted to exactly four mismatches:")
U, p = mw(ex_r, sp_4)
print(f"  n = {len(ex_r)} vs {len(sp_4)};  median {med(ex_r):g} vs {med(sp_4):g};  U = {U:g};  p = {p:.6f}")
print()
print("CONTROL, spared set restricted to FEWER than four mismatches (must differ from the matched row):")
sp_lo = [x["reads"] for x in rescued if x["mm"] < 4]
U, p = mw(ex_r, sp_lo)
print(f"  n = {len(ex_r)} vs {len(sp_lo)};  median {med(ex_r):g} vs {med(sp_lo):g};  U = {U:g};  p = {p:.6f}")
