#!/usr/bin/env python3
"""Is ceil(r/(w+1)) a LOWER bound too? Pre-registered in PREREG_2026-10-24_the_matching_lower_bound.md.

Part A verifies ingredient (i) exhaustively at small N, BOTH DIRECTIONS.
Part B tests the (2m+w)-scattered transversal at real parameters, with the raised-separation control
that must be infeasible.
"""
import os
import re, sys
from itertools import combinations, product
from ortools.sat.python import cp_model

ALPHA = 'ACGT'
LET = ALPHA


def library(g):
    ends = [k for k in range(len(g) - 1) if g[k] != g[k + 1]] + [len(g) - 1]
    return [g[:e] + g[e + 1:] for e in ends]


def ham(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def ball(D, m):
    N = len(D)
    out = set()
    for k in range(m + 1):
        for pos in combinations(range(N), k):
            for subs in product(*[[c for c in ALPHA if c != D[i]] for i in pos]):
                x = list(D)
                for i, c in zip(pos, subs):
                    x[i] = c
                out.add(''.join(x))
    return out


print('=== A. ingredient (i), exhaustive, both directions ===')
print('  a pattern fully covers B_m(P_i) IFF it agrees with P_i on every fixed position')
print('  guide    m   N   patterns  (covers & disagrees)  (agrees & not covers)  verdict')
okA = True
for g, m in (('ACGTA', 1), ('ACGTA', 2), ('ACGTAC', 1), ('ACGTAC', 2), ('AACGTA', 1)):
    L = library(g); N = len(L[0])
    balls = {k: ball(L[k], m) for k in range(len(L))}
    bad1 = bad2 = npat = 0
    for w in range(0, N + 1):
        for W in combinations(range(N), w):
            Ws = set(W); F = [i for i in range(N) if i not in Ws]
            if len(F) <= m:                      # ingredient (i) needs |F| > m
                continue
            for vals in product(ALPHA, repeat=len(F)):
                npat += 1
                pv = {i: vals[j] for j, i in enumerate(F)}
                for k, D in enumerate(L):
                    agrees = all(D[i] == pv[i] for i in F)
                    covers = all(sum(1 for i in F if x[i] != pv[i]) <= m for x in balls[k])
                    if covers and not agrees:
                        bad1 += 1
                    if agrees and not covers:
                        bad2 += 1
    v = 'HOLDS' if (bad1 == 0 and bad2 == 0) else 'FAILS'
    if bad1 or bad2:
        okA = False
    print('  %-8s %-3d %-3d %-9d %-21d %-22d %s' % (g, m, N, npat, bad1, bad2, v))
print()
if not okA:
    print('REFUSED: ingredient (i) failed, so the lower-bound derivation is unsound. Stopping.')
    sys.exit(1)


def transversal(pats, m, sep, limit=60.0):
    K, N = len(pats), len(pats[0])
    mo = cp_model.CpModel()
    x = [[mo.NewIntVar(0, 3, '') for _ in range(N)] for _ in range(K)]
    for i, pat in enumerate(pats):
        ne = []
        for p in range(N):
            b = mo.NewBoolVar('')
            mo.Add(x[i][p] != LET.index(pat[p])).OnlyEnforceIf(b)
            mo.Add(x[i][p] == LET.index(pat[p])).OnlyEnforceIf(b.Not())
            ne.append(b)
        mo.Add(sum(ne) <= m)
    for i, j in combinations(range(K), 2):
        ne = []
        for p in range(N):
            b = mo.NewBoolVar('')
            mo.Add(x[i][p] != x[j][p]).OnlyEnforceIf(b)
            mo.Add(x[i][p] == x[j][p]).OnlyEnforceIf(b.Not())
            ne.append(b)
        mo.Add(sum(ne) >= sep)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = limit
    s.parameters.num_search_workers = 8
    st = s.Solve(mo)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = [''.join(LET[s.Value(x[i][p])] for p in range(N)) for i in range(K)]
        ok = all(ham(sol[i], pats[i]) <= m for i in range(K)) and \
             all(ham(sol[i], sol[j]) >= sep for i, j in combinations(range(K), 2))
        return ('FEASIBLE' if ok else 'WITNESS FAILED RECHECK'), s.WallTime()
    if st == cp_model.INFEASIBLE:
        return 'INFEASIBLE', s.WallTime()
    return 'NO VERDICT', s.WallTime()


print('=== B. the (2m+w)-scattered transversal at real parameters, n=20 ===')
sp = sorted(set(re.findall(r'"([ACGT]{20})"',
            open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'multiguide_20261017', 'match_changeseq.json')).read())))
print('  guide        r   m   w   K=ceil(r/(w+1))  sep=2m+w+1  verdict      control sep+1  wall')
for g in sp[:4]:
    L = library(g); r = len(L)
    for m in (4, 5):
        for w in (1, 2):
            sub = [L[i] for i in range(0, r, w + 1)]
            sep = 2 * m + w + 1
            v, wt = transversal(sub, m, sep)
            vc, _ = transversal(sub, m, sep + 1)
            print('  %-12s %-3d %-3d %-3d %-16d %-11d %-12s %-14s %.1fs'
                  % (g[:12], r, m, w, len(sub), sep, v, vc, wt))
