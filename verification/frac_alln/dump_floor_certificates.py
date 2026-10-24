#!/usr/bin/env python3
"""Generate and SAVE an exact rational certificate for every proved cell of the boundary map.

Until now these existed only in run output, which is not a release artefact and not something a
referee can check. Each cell writes one JSON holding the exact dual masses as "num/den" strings,
the exact total, the exact maximum ball mass, the number of centre classes verified, and the
verdict. Nothing is read off a float: the LP locates the dual, the masses are rounded DOWN to
rationals (which can only lower what a ball collects, so feasibility is preserved), and both the
total and every ball constraint are then recomputed with Fraction arithmetic.

Criterion: K >= ceil(LP*), and the library attains r, so a feasible dual of exact mass > r-1 proves
K = r. Mass r is NOT required.
"""
import itertools, json, math, os, sys
from fractions import Fraction as Fr
from math import comb

# 2026-10-24: see certify_cell.py. At module level these made a machine without them exit 1 with
# a ModuleNotFoundError traceback, against the appendix's statement that a missing input exits 2
# NOT RUN and that no checker here ends in an uncaught exception.
try:
    import numpy as np
    from scipy.optimize import linprog
except ImportError as _exc:
    sys.stderr.write("NOT RUN: this checker needs numpy and scipy (%s).\n"
                     "install: python3 -m pip install numpy scipy\n"
                     "Every other checker in this folder is standard library only.\n" % _exc)
    raise SystemExit(2)

N = 19
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'floor_certificates')


def certify(R, M, D=10**7, cap=7000):
    AP = R - 1
    L = N - AP
    if L < 0:
        return None
    T = [''.join('1'*i + '0'*(AP-i)) for i in range(R)]
    dH = lambda a, b: sum(x != y for x, y in zip(a, b))
    pre = [''.join(b) for b in itertools.product('01', repeat=AP)]
    dmin = {y: min(dH(y, t) for t in T) for y in pre}
    orb = [(y, w) for y in pre for w in range(L+1) if w + dmin[y] <= M]
    if not orb or len(orb) > cap:
        return {'m': M, 'r': R, 'status': 'skipped', 'orbits': len(orb), 'cap': cap}
    shell = {w: comb(L, w) * 3**w for w in range(L+1)}
    _f = {}
    def Fc(p, w, s):
        if s < 0:
            return 0
        k = (p, w, s)
        if k in _f:
            return _f[k]
        t = 0
        for o in range(0, min(p, w)+1):
            b = comb(p, o) * comb(L-p, w-o)
            if b:
                t += b * sum(comb(o, e) * 2**(o-e) for e in range(o+1) if w+p-o-e <= s) * 3**(w-o)
        _f[k] = t
        return t
    cen = [(u, p) for u in pre for p in range(L+1)]
    A = np.array([[(Fc(p, w, M-dH(y, u)) / shell[w] if M-dH(y, u) >= 0 else 0.0)
                   for (y, w) in orb] for (u, p) in cen])
    res = linprog(c=-np.ones(len(orb)), A_ub=A, b_ub=np.ones(len(cen)),
                  bounds=[(0, None)]*len(orb), method='highs')
    if not res.success:
        return {'m': M, 'r': R, 'status': 'lp_failed', 'message': res.message}
    x = [Fr(int(math.floor(v*D)), D) for v in res.x]
    nz = [(j, x[j]) for j in range(len(orb)) if x[j]]
    total = sum(v for _, v in nz)
    worst, arg = Fr(0), None
    for (u, p) in cen:
        m = Fr(0)
        for j, xi in nz:
            y, w = orb[j]
            s = M - dH(y, u)
            if s >= 0:
                m += xi * Fr(Fc(p, w, s), shell[w])
        if m > worst:
            worst, arg = m, (u, p)
    proved = (worst <= 1) and (total > R - 1)
    return {'N': N, 'm': M, 'r': R, 'free_length': L,
            'status': 'PROVED' if proved else 'not_proved',
            'criterion': 'K >= ceil(LP*); a feasible dual of mass > r-1 proves K = r',
            'float_lp_optimum': float(-res.fun),
            'exact_total_mass': str(total), 'exact_total_mass_float': float(total),
            'ceiling': math.ceil(total),
            'exact_max_ball_mass': str(worst), 'max_ball_at': arg,
            'centre_classes_verified': len(cen), 'support_orbits': len(nz),
            'rounding': 'masses rounded DOWN to denominator %d, which preserves feasibility' % D,
            'dual': [{'prefix': orb[j][0], 'free_weight': orb[j][1], 'mass': str(v)}
                     for j, v in nz]}


CELLS = [(8,7),(8,8),(8,9),(8,10),(9,7),(9,8),(10,3),(10,4),(10,5),(10,6),
         (11,3),(11,4),(11,5),(12,3),(12,4)]

# 2026-10-22: --only lets a slow cell be run on its own rather than taking the batch down with it.
# The first full run stopped after three cells, so the remainder are driven in groups.
if '--only' in sys.argv:
    want = sys.argv[sys.argv.index('--only')+1]
    CELLS = [tuple(int(v) for v in c.split(',')) for c in want.split(';')]

# A cell already on disk is not recomputed, so a restart is cheap and idempotent.
_SKIP_EXISTING = '--force' not in sys.argv

def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for (M, R) in CELLS:
        fn_pre = os.path.join(OUT, 'cert_m%d_r%d.json' % (M, R))
        if _SKIP_EXISTING and os.path.exists(fn_pre):
            c = json.load(open(fn_pre))
            rows.append(c)
            print('m=%2d r=%2d  %-10s  already on disk' % (M, R, c.get('status')), flush=True)
            continue
        c = certify(R, M)
        if c is None:
            continue
        fn = os.path.join(OUT, 'cert_m%d_r%d.json' % (M, R))
        json.dump(c, open(fn, 'w'), indent=1)
        rows.append(c)
        print('m=%2d r=%2d  %-10s  mass %-26s ceil %-3s maxball %-12s centres %-5s -> %s'
              % (M, R, c.get('status'), c.get('exact_total_mass','-'), c.get('ceiling','-'),
                 c.get('exact_max_ball_mass','-'), c.get('centre_classes_verified','-'),
                 os.path.basename(fn)), flush=True)
    idx = os.path.join(OUT, 'INDEX.json')
    json.dump({'criterion': 'K >= ceil(LP*)', 'cells': [
        {k: r.get(k) for k in ('m','r','status','exact_total_mass','ceiling',
                               'exact_max_ball_mass','centre_classes_verified','support_orbits')}
        for r in rows]}, open(idx, 'w'), indent=1)
    nproved = sum(1 for r in rows if r.get('status') == 'PROVED')
    print('\n%d of %d cells PROVED; index at %s' % (nproved, len(rows), idx))
    return 0 if nproved else 2


if __name__ == '__main__':
    sys.exit(main())
