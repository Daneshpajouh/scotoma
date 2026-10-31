#!/usr/bin/env python3
"""Exact certificates for the FIRST FAILURE of each row, closing the boundary on both sides.

The monotone-deficit lemma makes the provable run counts an INITIAL SEGMENT: d(m,r) = r - tau*(m,r)
never falls, so once the ceiling criterion fails it fails for every larger r. That answers "did you
try r = 20" without computing it. It does NOT say where the segment ends: it gives
d(m,r+1) >= d(m,r), not d(m,r+1) >= 1.

So each row-end needs an exact UPPER bound on tau* at the first failing cell, which means a feasible
PRIMAL fractional cover of value below r - 1 + 1 = r... more precisely, of value whose ceiling is
below r. A float optimum is not that. So: solve the primal in floats, round every variable UP to a
rational, which preserves feasibility and can only raise the objective, then verify feasibility on
every orbit and the objective exactly. If ceil of the rounded value is still below r, the row-end is
proved.
"""
import itertools, json, math, os, sys
from fractions import Fraction as Fr
from math import comb
import numpy as np
from scipy.optimize import linprog

N = 19
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'floor_certificates')
# Smallest first, deliberately. The first version led with (8,11), whose dense LP is
# 10,240 x 10,240 at about 840 MB plus an exact check over 10,240 orbits, and it took the whole
# batch down before any cell was written. A batch ordered worst-first reports nothing when the
# worst cell fails; ordered best-first it reports everything it can.
ROW_ENDS = [(12, 5), (11, 6), (10, 7), (9, 9), (8, 11)]
if '--only' in sys.argv:
    want = sys.argv[sys.argv.index('--only')+1]
    ROW_ENDS = [tuple(int(v) for v in c.split(',')) for c in want.split(';')]


def row_end(M, R, D=10**6):
    AP = R - 1; L = N - AP
    T = [''.join('1'*i + '0'*(AP-i)) for i in range(R)]
    dH = lambda a, b: sum(x != y for x, y in zip(a, b))
    pre = [''.join(b) for b in itertools.product('01', repeat=AP)]
    dmin = {y: min(dH(y, t) for t in T) for y in pre}
    orb = [(y, w) for y in pre for w in range(L+1) if w + dmin[y] <= M]
    shell = {w: comb(L, w)*3**w for w in range(L+1)}
    _f = {}
    def Fc(p, w, s):
        if s < 0: return 0
        k = (p, w, s)
        if k in _f: return _f[k]
        t = 0
        for o in range(0, min(p, w)+1):
            b = comb(p, o)*comb(L-p, w-o)
            if b: t += b*sum(comb(o, e)*2**(o-e) for e in range(o+1) if w+p-o-e <= s)*3**(w-o)
        _f[k] = t; return t
    cen = [(u, p) for u in pre for p in range(L+1)]
    def frac(y, w, u, p):
        s = M - dH(y, u)
        return Fr(Fc(p, w, s), shell[w]) if s >= 0 else Fr(0)
    Af = np.array([[float(frac(y, w, u, p)) for (y, w) in orb] for (u, p) in cen])
    res = linprog(c=np.ones(len(cen)), A_ub=-Af.T, b_ub=-np.ones(len(orb)),
                  bounds=[(0, None)]*len(cen), method='highs')
    if not res.success:
        return {'m': M, 'r': R, 'status': 'lp_failed', 'message': res.message}
    x = [Fr(int(math.ceil(v*D)), D) for v in res.x]        # round UP: feasibility preserved
    total = sum(x)
    nz = [(i, x[i]) for i in range(len(cen)) if x[i]]
    viol = 0
    for (y, w) in orb:
        s = sum(xi*frac(y, w, cen[i][0], cen[i][1]) for i, xi in nz)
        if s < 1: viol += 1
    proved = (viol == 0) and (math.ceil(total) < R)
    return {'N': N, 'm': M, 'r': R, 'free_length': L,
            'status': 'ROW_END_PROVED' if proved else 'not_proved',
            'claim': 'no fractional dual can prove K=r here, because ceil(tau*) < r',
            'float_primal_optimum': float(res.fun),
            'exact_cover_value': str(total), 'exact_cover_value_float': float(total),
            'ceiling': math.ceil(total), 'needed_for_proof': 'ceiling < %d' % R,
            'orbit_constraints_checked': len(orb), 'violations': viol,
            'support_centres': len(nz),
            'rounding': 'primal rounded UP to denominator %d, which preserves feasibility' % D}


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for (M, R) in ROW_ENDS:
        try:
            c = row_end(M, R)
        except (MemoryError, ValueError) as e:
            c = {'m': M, 'r': R, 'status': 'not_attempted',
                 'reason': '%s: %s' % (type(e).__name__, e),
                 'note': 'too large for the dense route; needs row generation'}
            print('m=%2d r=%2d  not attempted: %s' % (M, R, type(e).__name__), flush=True)
        fn = os.path.join(OUT, 'rowend_m%d_r%d.json' % (M, R))
        json.dump(c, open(fn, 'w'), indent=1)
        rows.append(c)
        print('m=%2d r=%2d  %-16s cover %-24s ceil %-3s < r  viol %-3s orbits %-5s -> %s'
              % (M, R, c.get('status'), c.get('exact_cover_value', '-'), c.get('ceiling', '-'),
                 c.get('violations', '-'), c.get('orbit_constraints_checked', '-'),
                 os.path.basename(fn)), flush=True)
    n = sum(1 for r in rows if r.get('status') == 'ROW_END_PROVED')
    print('\n%d of %d row-ends PROVED' % (n, len(rows)))
    return 0 if n == len(rows) else 2


if __name__ == '__main__':
    sys.exit(main())
