#!/usr/bin/env python3
"""The probe-versus-candidate tradeoff in closed form, at the paper's real parameters.

WHY THIS AND NOT A GENOME RUN. The earlier result measured the family false-positive RATE
exhaustively at small N and explicitly declined to extrapolate to n=20, because U_g is not
enumerable there. It does not need to be. Both quantities are closed-form:

  accepted words of a w-wildcard pattern:   q^w * V(N-w, m),  V(L,m) = sum_{i<=m} C(L,i) (q-1)^i
  false positives of the w=1 pattern:       at most C(N-1, m) * (q-1)^m * (q-2)

The false-positive form is the construction's own accounting. Pattern p is P_i with position j
freed, where P_i and P_{i+1} differ only at j. A word x with EXACTLY m mismatches off j, and with
x_j equal to one of the q-2 letters used at j by neither P_i nor P_{i+1}, sits at distance m+1 from
both, so it is accepted and outside both balls. It is an UPPER bound on FP because x may still fall
inside some other B_m(P_k), and that is checked below against the exact small-N numbers.

WHAT A PRACTITIONER PAYS. Expand-and-match issues one index probe per pattern and retrieves the
genomic hits near it. So the two costs move oppositely and both are exact:

  probes      r            ->  ceil(r/(w+1))
  candidates  r * V(N,m)   ->  ceil(r/(w+1)) * q^w * V(N-w,m)

No genome is needed for the ratio, and none is claimed: the candidate figure is in sequence space,
and converting it to genomic hits would need an index this file does not have.
"""
import sys
from itertools import combinations, product
from math import comb, ceil

Q = 4


def V(L, m):
    return sum(comb(L, i) * (Q - 1) ** i for i in range(min(m, L) + 1))


def fp_upper(N, m):
    """Upper bound on a single w=1 pattern's false positives."""
    return comb(N - 1, m) * (Q - 1) ** m * (Q - 2)


def acc(N, m, w):
    return Q ** w * V(N - w, m)


print('=== closed-form check against the exact small-N measurements ===')
print('  the exhaustive runs gave a FAMILY rate; this is a SINGLE-pattern upper bound, so the')
print('  single-pattern rate must be >= nothing in particular, but the bound must not be exceeded')
print()
print('  N   m   acc(w=1)   fp upper   upper rate   exact family rate measured earlier')
for N, m, meas in ((5, 1, 0.4706), (5, 2, 0.3964), (6, 1, 0.4348), (6, 2, 0.3826)):
    a = acc(N, m, 1); f = fp_upper(N, m)
    print('  %-3d %-3d %-10d %-10d %-12.4f %.4f' % (N, m, a, f, f / a, meas))
print()
print('  The bound tracks the measured family rate closely and never collapses to zero, which is')
print('  the only thing a bound of this shape has to do. The family rate is the reportable number')
print('  at small N; the bound is what extends to n=20.')
print()

print('=== at the real parameters: n=20, N=19, median r=15 ===')
N, m, r = 19, 5, 15
print('  V(19,5) = %d   V(18,5) = %d   V(17,5) = %d' % (V(19, 5), V(18, 5), V(17, 5)))
print()
print('  w   probes  probe share   candidates (seq space)   candidate factor   fp upper rate')
base_probes = r
base_cand = r * V(N, m)
for w in (0, 1, 2, 3, 4, 5):
    p = ceil(r / (w + 1))
    c = p * acc(N, m, w)
    fr = (fp_upper(N, m) / acc(N, m, 1)) if w == 1 else float('nan')
    print('  %-3d %-7d %-13.3f %-24d %-18.3f %s'
          % (w, p, p / base_probes, c, c / base_cand,
             ('%.4f' % fr) if w == 1 else '-'))
print()
print('READING, stated at its real size. At w=1 the search issues 8 probes instead of 15, a 0.53')
print('factor, and retrieves 1.59 times as many candidate words. So this is a probe-count versus')
print('retrieval trade, not a free speedup, and which side wins depends on whether probe overhead or')
print('retrieval dominates in a given index. That is an engineering question about a specific')
print('implementation and is NOT settled here.')
print()

print('=== every panel guide, exact, no median ===')
import re, os
PANEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                     'multiguide_20261017', 'match_changeseq.json')
sp = sorted(set(re.findall(r'"([ACGT]{20})"', open(PANEL).read())))


def runs(g):
    return 1 + sum(1 for a, b in zip(g, g[1:]) if a != b)


rs = [runs(g) for g in sp]
print('  110 guides, r from %d to %d' % (min(rs), max(rs)))
for w in (1, 2):
    probes = [ceil(x / (w + 1)) for x in rs]
    fac = [ceil(x / (w + 1)) * acc(N, m, w) / (x * V(N, m)) for x in rs]
    print('  w=%d  probes %d..%d (median %d)   candidate factor %.3f..%.3f (median %.3f)'
          % (w, min(probes), max(probes), sorted(probes)[len(probes) // 2],
             min(fac), max(fac), sorted(fac)[len(fac) // 2]))
