"""Exact UPPER bound on tau*(N,m,r): a feasible fractional cover verified over EVERY point orbit.

This replaces rowend_by_constraint_generation.py, whose five certificates were withdrawn on
2026-10-22. That script built its point orbits from itertools.product('01', repeat=AP), so its
constraint set was the binary-prefix orbits only, while its own metadata recorded "the exact check is
over ALL orbits". Dropping constraints LOWERS a primal optimum, so its values were not upper bounds;
its m=8 r=11 figure 9639109/1000000 = 9.639109 sits BELOW the certified lower bound 10.169959597 for
the same cell, which is the contradiction that exposed it.

WHICH RESTRICTIONS ARE SAFE, since getting this backwards is what cost five certificates:
  - restricting the CENTRES (the variables) removes cover options, so the optimum can only RISE.
    A cover built from a subset of centres is still a cover, so its value is still an upper bound on
    tau*. SAFE, and used here to keep the LP small.
  - restricting the POINT ORBITS (the constraints) lets the LP ignore points it must cover, so the
    optimum FALLS and the result need not be a cover at all. NOT SAFE. Row generation is used, but
    the value is accepted only after an exact check against every orbit in U with zero violations.

Values are rounded UP to a fixed denominator, which preserves coverage and can only raise the
objective, and the whole verification runs in int64 over a common denominator.

TWO-SIDED CONTROL. Where a dual certificate exists for the same cell, the cover value must be at
least the dual mass, because both bracket the same tau*. A cover below the dual is impossible and
means the cover is not a cover. The withdrawn file fails this test by 0.53 at m=8, r=11; this file
refuses to write a certificate that fails it.
"""
import itertools, json, math, os, sys, time
import numpy as np
from fractions import Fraction as Fr
from math import comb, gcd
from scipy.optimize import linprog

N  = int(os.environ.get('CERT_N', 19))
M  = int(os.environ.get('CERT_M', 8))
R  = int(os.environ.get('CERT_R', 11))
NB = int(os.environ.get('COV_NONBIN', 1200))      # non-binary centres offered (safe to cap)
ROUNDS = int(os.environ.get('COV_ROUNDS', 30))
Dy = 10**6
AP, L = R-1, N-(R-1)
t0 = time.time()

LIB   = np.array([[1]*i + [0]*(AP-i) for i in range(R)], dtype=np.int8)
shell = {w: comb(L,w)*3**w for w in range(L+1)}
LCM = 1
for w in shell: LCM = LCM*shell[w]//gcd(LCM, shell[w])
DG = (2**AP) * LCM

_f={}
def F(p,w,s):
    if s<0: return 0
    k=(p,w,s)
    if k in _f: return _f[k]
    tt=0
    for o in range(0,min(p,w)+1):
        b=comb(p,o)*comb(L-p,w-o)
        if b: tt+=b*sum(comb(o,e)*2**(o-e) for e in range(o+1) if w+p-o-e<=s)*3**(w-o)
    _f[k]=tt; return tt

types = np.array(list(itertools.product([0,1,2], repeat=AP)), dtype=np.int8)
nT = len(types); isX = (types == 2); nX = isX.sum(axis=1)
dmin = np.full(nT, 1<<30, dtype=np.int32)
for lib in LIB:
    dmin = np.minimum(dmin, (isX | ((~isX) & (types != lib))).sum(axis=1))

# POINT ORBITS: every (type, free weight) inside U. Types range over {0,1,X}^AP, NOT {0,1}^AP.
_ii = np.repeat(np.arange(nT, dtype=np.int64), L+1)
_ww = np.tile(np.arange(L+1, dtype=np.int64), nT)
_in = (_ww + dmin[_ii]) <= M
pt_i = _ii[_in].copy(); pt_w = _ww[_in].copy()
del _ii, _ww, _in
NPT = int(pt_i.size)
print('point orbits in U: %d over %d types (the withdrawn script used %d binary-prefix orbits)'
      % (NPT, nT, sum(1 for i in range(nT) if nX[i]==0 for w in range(L+1) if w+dmin[i]<=M)),
      flush=True)

_lt={}
def load_table(w,p):
    key=(w,p)
    if key in _lt: return _lt[key]
    T=np.zeros((AP+1,AP+1), dtype=np.int64)
    for c in range(AP+1):
        sc = DG // ((2**c)*shell[w])
        for b in range(AP+1-c):
            tot=0
            for kk in range(c+1): tot += comb(c,kk)*F(p,w,M-b-kk)
            T[b,c]=tot*sc
    _lt[key]=T
    return T

def BC_against_all(j):
    """(B, c) of centre type j against EVERY point type."""
    t=types[j]; tX=isX[j]
    return ((((~tX)&(~isX)&(types!=t)).sum(axis=1) + (tX^isX).sum(axis=1)).astype(np.int32),
            (tX&isX).sum(axis=1).astype(np.int32))

# CENTRES: the variables. Binary types plus a sample of non-binary ones; capping is safe.
rng = np.random.default_rng(0)
bin_types = np.nonzero(nX == 0)[0]
oth_types = np.nonzero(nX >  0)[0]
sel_types = np.concatenate([bin_types,
                            rng.choice(oth_types, size=min(len(oth_types), NB), replace=False)])
cen = [(int(j), p) for j in sel_types for p in range(L+1)]
print('centres offered: %d (%d binary types + %d non-binary types) x %d free weights'
      % (len(cen), len(bin_types), len(sel_types)-len(bin_types), L+1), flush=True)

act = sorted(rng.choice(NPT, size=min(1500, NPT), replace=False).tolist())
best = None
for rnd in range(ROUNDS):
    ai = pt_i[act]; aw = pt_w[act]
    A = np.empty((len(act), len(cen)))
    ta, isXa = types[ai], isX[ai]
    for a,(j,p) in enumerate(cen):
        t=types[j]; tX=isX[j]
        B=(((~tX)&(~isXa)&(ta!=t)).sum(axis=1) + (tX^isXa).sum(axis=1)).astype(np.int32)
        c=((tX)&isXa).sum(axis=1).astype(np.int32)
        colv=np.empty(len(act))
        for w in np.unique(aw):
            sel=(aw==w); T=load_table(int(w), p)
            colv[sel]=T[B[sel], c[sel]]/DG
        A[:,a]=colv
    res = linprog(c=np.ones(len(cen)), A_ub=-A, b_ub=-np.ones(len(act)),
                  bounds=[(0,None)]*len(cen), method='highs')
    if not res.success:
        print('LP failed:', res.message); sys.exit(1)
    ynum = np.ceil(np.asarray(res.x)*Dy).astype(np.int64)     # round UP: preserves coverage
    nz = np.nonzero(ynum)[0]
    total = Fr(int(ynum.sum()), Dy)

    # EXACT coverage sweep over EVERY orbit in U, accumulated per free weight
    cov = [np.zeros(nT, dtype=np.int64) for _ in range(L+1)]
    for a in nz:
        j,p = cen[int(a)]; Bv,cv = BC_against_all(j); ya = int(ynum[a])
        for w in range(L+1):
            cov[w] += ya * load_table(w,p)[Bv,cv]
    need = Dy*DG
    short = cov[0][0]*0
    viol = []
    worst = None
    for w in range(L+1):
        sel = np.nonzero((dmin <= M-w))[0]
        if sel.size == 0: continue
        c_w = cov[w][sel]
        bad = sel[c_w < need]
        if bad.size:
            mn = int(c_w.min())
            if worst is None or mn < worst[0]: worst = (mn, w)
            viol += [(int(i), w) for i in bad[:600]]
    print('round %2d  constraints %6d  cover value %12.7f  support %5d  violations %7d  %5.0fs'
          % (rnd+1, len(act), float(total), len(nz), len(viol), time.time()-t0), flush=True)
    if not viol:
        best = (total, len(nz), rnd+1)
        break
    newact = sorted(set(act) | {int(np.searchsorted(pt_i*(L+1)+pt_w, i*(L+1)+w)) for i,w in viol})
    # searchsorted needs the key array sorted; pt is built in (i,w) lexicographic order so it is
    act = sorted(set(act) | set(np.nonzero(np.isin(pt_i*(L+1)+pt_w,
                   np.array([i*(L+1)+w for i,w in viol], dtype=np.int64)))[0].tolist()))

if best is None:
    print('\nno cover verified in %d rounds; NOTHING is claimed' % ROUNDS)
    sys.exit(0)

total, nsup, rnds = best
print()
print('EXACT cover value %s = %.9f over ALL %d point orbits, zero violations'
      % (total, float(total), NPT))
print('ceil = %d, and r = %d' % (math.ceil(total), R))

# TWO-SIDED CONTROL against the dual certificate for the same cell
here = os.path.dirname(os.path.abspath(__file__))
dual_path = os.path.join(here, 'data', 'floor_certificates', 'cert_m%d_r%d.json' % (M, R))
if os.path.exists(dual_path):
    dmass = Fr(json.load(open(dual_path))['exact_total_mass'])
    print('dual certificate for this cell: %s = %.9f' % (dmass, float(dmass)))
    if total < dmass:
        sys.exit('CONTROL FAILED: cover %s is BELOW the certified dual %s. Both bracket the same\n'
                 '  tau*, so a cover under the dual is not a cover. This is the exact test the\n'
                 '  withdrawn rowend_m8_r11.json fails (9.639109 against 10.169959597).' % (total, dmass))
    print('CONTROL PASSED: cover >= dual, so %.9f <= tau* <= %.9f' % (float(dmass), float(total)))
else:
    print('no dual certificate on disk for this cell, so the two-sided control did not run')

verdict = ('ROW_END: ceil(tau*) = %d < r = %d, so NO fractional dual can prove K=%d here'
           % (math.ceil(total), R, R)) if math.ceil(total) < R else \
          ('consistent with K=%d; this cover does not close the cell' % R)
print(verdict)
