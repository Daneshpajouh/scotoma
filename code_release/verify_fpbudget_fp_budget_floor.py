#!/usr/bin/env python3
"""Minimum false-positive count of a wildcard-bearing pattern over a single-deletion library.

Pre-registered in PREREG_2026-10-24_is_the_pattern_floor_robust_to_a_nonzero_fp_budget.md.
Settles, exhaustively at small size, the quantity the zero-budget theorem does not give: it proves
every wildcard pattern has at least ONE false positive, not how many.

Definitions taken from the appendix, not from any summary of it. A pattern fixes F and leaves
wildcards W, accepting every word differing from it in at most m NON-WILDCARD positions.
U_g is the union of radius-m balls about the single-deletion library.

CONTROL. At q=2 there is no spare letter, the escape construction cannot fire, and admissible
patterns (FP = 0) are known to exist. If this code reports a positive minimum at q=2 it is wrong,
and no q=4 number from it means anything.
"""
import sys
from itertools import combinations, product
from math import comb


def library(g):
    """Single-deletion patterns in run order, from run ends."""
    ends = [k for k in range(len(g) - 1) if g[k] != g[k + 1]] + [len(g) - 1]
    return [g[:e] + g[e + 1:] for e in ends]


def ball_union(L, m, alpha):
    """Explicit U_g. Exact, by expanding each ball."""
    N = len(L[0])
    U = set()
    for D in L:
        for k in range(m + 1):
            for pos in combinations(range(N), k):
                for subs in product(*[[c for c in alpha if c != D[i]] for i in pos]):
                    x = list(D)
                    for i, c in zip(pos, subs):
                        x[i] = c
                    U.add(''.join(x))
    return U


def scan(g, m, alpha):
    q = len(alpha)
    L = library(g)
    N = len(L[0])
    U = sorted(ball_union(L, m, alpha))
    Ut = [tuple(u) for u in U]
    best = {}                                   # |W| -> (minFP, minrate, pattern)
    for w in range(1, N + 1):
        nf = N - w
        acc_total = q ** w * sum(comb(nf, i) * (q - 1) ** i for i in range(min(m, nf) + 1))
        for W in combinations(range(N), w):
            Ws = set(W)
            F = [i for i in range(N) if i not in Ws]
            for vals in product(alpha, repeat=nf):
                inside = 0
                for u in Ut:
                    d = 0
                    for j, i in enumerate(F):
                        if u[i] != vals[j]:
                            d += 1
                            if d > m:
                                break
                    if d <= m:
                        inside += 1
                fp = acc_total - inside
                rate = fp / acc_total
                if w not in best or fp < best[w][0]:
                    pat = ''.join('*' if i in Ws else vals[F.index(i)] for i in range(N))
                    best[w] = (fp, rate, pat, acc_total)
    return L, U, best


def report(g, m, alpha, label):
    L, U, best = scan(g, m, alpha)
    print('%s  guide %s  q=%d  m=%d  N=%d  |L|=%d  |U_g|=%d'
          % (label, g, len(alpha), m, len(L[0]), len(L), len(U)))
    print('   |W|   min FP      accepted    min rate   attaining pattern')
    for w in sorted(best):
        fp, rate, pat, tot = best[w]
        print('   %-5d %-11d %-11d %-10.4f %s' % (w, fp, tot, rate, pat))
    allmin = min(v[0] for v in best.values())
    return allmin


print('==== CONTROL, q=2: no spare letter exists, so admissible patterns must appear ====')
c1 = report('ACACA', 1, 'AC', 'CTRL')
c2 = report('ACACAC', 1, 'AC', 'CTRL')
print()
if c1 != 0 or c2 != 0:
    print('REFUSED: the q=2 control found minimum FP %d and %d, not 0. The note records 9 of 65 and'
          ' 5 of 211 admissible patterns there, so this code cannot find FP=0 and no q=4 figure'
          ' from it is reportable.' % (c1, c2))
    sys.exit(1)
print('control passed: FP = 0 is reachable at q=2, so a positive minimum at q=4 is a fact about q=4')
print()

print('==== q=4, exhaustive ====')
for g, m in (('ACGTA', 1), ('ACGTA', 2), ('ACGTAC', 1), ('ACGTAC', 2),
             ('AACGTA', 1), ('AAACCC', 1)):
    report(g, m, 'ACGT', 'q4  ')
    print()
