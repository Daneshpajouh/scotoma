"""Certified lower bound on tau*(N,m,r) by a min-max-load LP with an EXACT integer check.

Cell chosen by CERT_N / CERT_M / CERT_R; caps by CERT_BIN / CERT_NONBIN / CERT_MAXROWS / CERT_ROUNDS.

Why not the maximise formulation. Maximising total mass subject to every ball collecting <= 1 is
capped at r = 11, because the r library balls already cover everything and each collects <= 1. The
cutting-plane run hit that cap at round 7 and sat on it for four rounds while the true maximum load
drifted UP (1.129, 1.095, 1.141, 1.163): the LP was walking a degenerate face of mass-11 vertices,
each violating somewhere new. It is not converging and the value it prints is not a bound.

This solves the homogeneous twin instead:

    minimise t   subject to   A x <= t 1  (active centres),  sum x = S,  x >= 0.

Two properties the maximise form does not have.
  1. MONOTONE. Adding centres can only raise t*, so there is no cycling to sit in.
  2. A CERTIFICATE EVERY ROUND. Whatever x comes back, measure its exact maximum load Mex over ALL
     centres; x/Mex is then feasible everywhere, so tau* >= S/Mex. That is a real bound even when
     the separation has not finished, which is the number the appendix has to leave open.

K >= ceil(tau*), so S/Mex > 10 proves K(19,8,11) = 11.

EXACT, not float. Every load is tot / (2^c * shell[w]) with tot an integer, so all of them are
integers over Dglob = 2^AP * lcm_w shell[w]. Rounding x DOWN to a denominator Dx keeps the mass
honest (down) and the loads honest (down), and the whole sweep runs in int64: the largest partial
sum is (sum x_num) * Dglob = 1.9e16 against an int64 ceiling of 9.2e18.
"""
import itertools, sys, time
import numpy as np
from math import comb, gcd
from fractions import Fraction as Fr
from scipy.optimize import linprog

import os
N  = int(os.environ.get('CERT_N', 19))
M  = int(os.environ.get('CERT_M', 8))
R  = int(os.environ.get('CERT_R', 11))
NB = int(os.environ.get('CERT_NONBIN', 1400))     # non-binary columns offered
NBIN = int(os.environ.get('CERT_BIN', 10**9))     # cap on binary columns offered
ROUNDS = int(os.environ.get('CERT_ROUNDS', 40))
AP, L = R-1, N-(R-1)
S  = R           # fixed mass; the ratio S/Mex is homogeneous so the choice is cosmetic
Dx = 10**6       # denominator x is rounded DOWN to
# VALIDITY DOES NOT DEPEND ON OFFERING EVERY COLUMN. A measure supported on ANY subset of point
# orbits, rescaled by its exact maximum load over ALL centres, is a feasible packing, so a capped
# column set can only make the bound WEAKER, never wrong. That is what makes marching up the panel
# safe: at larger r the column set is capped to fit in memory and the bound is still a certificate.
t0 = time.time()

LIB   = np.array([[1]*i + [0]*(AP-i) for i in range(R)], dtype=np.int8)
shell = {w: comb(L,w)*3**w for w in range(L+1)}
LCM = 1
for w in shell: LCM = LCM*shell[w]//gcd(LCM, shell[w])
DG = (2**AP) * LCM
assert all(DG % ((2**c)*shell[w]) == 0 for c in range(AP+1) for w in shell)

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

# BASE-3 DIGITS, NOT itertools.product. The product builds a Python list of nT tuples before numpy
# sees it: 14.3 million tuples at r=16 and 43 million at r=17, several GB of object headers for an
# array that is 688 MB. The i-th tuple of product([0,1,2], repeat=AP) in lexicographic order IS the
# base-3 representation of i, most significant digit first, so the array can be written directly.
# One column at a time, so the int64 temporary is nT and not nT*AP.
# The ORDER MUST MATCH the product exactly, because a certificate's support indices are defined
# against it and verify_ceiling_certificate.py rebuilds them with the product; the assertion below
# checks that on every run, at a size where the product is free.
nT = 3 ** AP
types = np.empty((nT, AP), dtype=np.int8)
_idx = np.arange(nT, dtype=np.int64)
for _j in range(AP):
    types[:, _j] = (_idx // (3 ** (AP-1-_j))) % 3
del _idx
_chk = min(nT, 3**min(AP, 6))
assert np.array_equal(types[:_chk], np.array(list(itertools.islice(
    itertools.product([0,1,2], repeat=AP), _chk)), dtype=np.int8)),     'type ordering diverges from itertools.product; every certificate index would be wrong'
assert nT == len(types)
isX = (types == 2); nX = isX.sum(axis=1)
dmin = np.full(nT, 1<<30, dtype=np.int32)
for lib in LIB:
    dmin = np.minimum(dmin, (isX | ((~isX) & (types != lib))).sum(axis=1))
# NO PYTHON TUPLES FOR THE COLUMN AND CENTRE SETS. At r=14 the support set has 3.6 million columns
# and the centre set 11.2 million classes; held as lists of tuples those are several GB of object
# headers alone, and the r=14 run was killed by the kernel with an empty log before it printed a
# single round. Columns live in two int64 arrays and a centre is the single integer p*nT+i, which is
# the encoding the violation sweep already used.
_ii = np.repeat(np.arange(nT, dtype=np.int64), L+1)
_ww = np.tile(np.arange(L+1, dtype=np.int64), nT)
_keep_col = (_ww + dmin[_ii]) <= M
col_i = _ii[_keep_col].copy(); col_w = _ww[_keep_col].copy()
del _ii, _ww, _keep_col
NCOL = int(col_i.size)
class _Cols:
    __slots__ = ()
    def __getitem__(self, k): return int(col_i[k]), int(col_w[k])
    def __len__(self): return NCOL
cols = _Cols()
NCEN = nT*(L+1)
def centre(c): return int(c) % nT, int(c) // nT          # (active-prefix type index, free weight)
print('types %d  free length %d  support columns %d  centre classes %d  Dglob %d  (%.0fs)'
      % (nT, L, NCOL, NCEN, DG, time.time()-t0), flush=True)

# REMOVED 2026-10-22, and this is why the earlier cap did nothing. There were TWO definitions of
# bc(): an uncapped one here and a capped one further down. Python keeps the LAST, so the cap was
# applied to the shadowed copy while the live one retained a (B, c) pair per distinct column type,
# 115 MB each at r=16, for the whole sweep. RSS climbed 6.4 -> 29.9 GB in one round and the kernel
# killed it twice, while the log reported a 0.7 GB cap that was never in force. The chunked sweep
# needs neither function; both are gone rather than left shadowing each other.
def num_table(w, p, Bv, cv):
    """(B,c) -> INTEGER numerator of the load over Dglob."""
    tab={}
    for bb in np.unique(Bv):
        for cc in np.unique(cv):
            b,c = int(bb), int(cc)
            tot = sum(comb(c,k)*F(p,w,M-b-k) for k in range(c+1))
            tab[(b,c)] = tot * (DG // ((2**c)*shell[w]))
    return tab

_lt = {}
def load_table(w, p):
    """the (AP+1) x (AP+1) integer load table for one (w, p), over Dglob. Tiny, so always cached."""
    key = (w, p)
    if key in _lt: return _lt[key]
    T = np.zeros((AP+1, AP+1), dtype=np.int64)
    for c in range(AP+1):
        sc = DG // ((2**c)*shell[w])
        for b in range(AP+1-c):
            tot = 0
            for kk in range(c+1): tot += comb(c,kk)*F(p,w,M-b-kk)
            T[b,c] = tot * sc
    _lt[key] = T
    return T

CHUNK = int(os.environ.get('CERT_CHUNK', 2_000_000))
print('sweep chunk: %d centre types (%d chunk(s) of %d total)'
      % (CHUNK, -(-nT // CHUNK), nT), flush=True)
rng = np.random.default_rng(0)
_nxc = nX[col_i]
# PRIVACY, AND WHY IT DECIDES THE COLUMN SET. d(z,P_k) = A(k) + X + w with A(k) a +-1 walk in k set
# by the point's type, so the number of library balls containing an orbit is mult = #{k: A(k) <= m-X-w},
# a function of (type, free weight) alone. Summing the library balls' loads for a packing of mass T
# with every load at most 1 gives sum_z x_z*(mult(z)-1) <= r - T, so at T = r-0.73 under one per cent
# of the orbits' worth of excess multiplicity is affordable and over 90 per cent of the mass has to
# sit on PRIVATE orbits, mult = 1. Mean multiplicity over U is 3.3 to 3.9, so a uniform random sample
# of columns is almost entirely the wrong region. Offer the private ones first.
PRIVATE_FIRST = os.environ.get('CERT_PRIVATE', '0') == '1'
if PRIVATE_FIRST:
    # ONLY COMPUTED WHEN IT IS USED. The walk A(k) needs five int arrays of shape (nT, AP+1): at
    # r=16 that is 14.3 million rows and about 4.5 GB at peak, which the r=16 run paid for and then
    # discarded, because private-first seeding is off by default. It was silently OOM-killed after
    # one round with no traceback. int16 holds A, whose range is 0..AP, and the cumulative sums are
    # built one column at a time so the two (nT, AP+1) temporaries never exist at once.
    _A = np.empty((nT, AP+1), dtype=np.int16)
    _zc = np.zeros(nT, dtype=np.int16); _os = (types == 1).sum(axis=1).astype(np.int16)
    for k in range(AP+1):
        _A[:, k] = _zc + _os
        if k < AP:
            _zc = _zc + (types[:, k] == 0)
            _os = _os - (types[:, k] == 1)
    _mult = np.empty(NCOL, dtype=np.int32)
    for w in range(L+1):
        selw = np.nonzero(col_w == w)[0]
        if selw.size:
            _mult[selw] = (_A[col_i[selw]] <= (M - nX[col_i[selw]] - w)[:, None]).sum(axis=1)
    del _A, _zc, _os
    priv = np.nonzero(_mult == 1)[0]
    print('private orbits (mult=1): %d of %d columns (%.1f%%); private-first seeding ON'
          % (priv.size, NCOL, 100.0*priv.size/NCOL), flush=True)
    binary = priv.tolist()
    other  = np.nonzero(_mult > 1)[0].tolist()
    del _mult
else:
    binary = np.nonzero(_nxc == 0)[0].tolist()
    other  = np.nonzero(_nxc >  0)[0].tolist()
del _nxc
if len(binary) > NBIN: binary = sorted(rng.choice(binary, size=NBIN, replace=False))
keep = binary + list(rng.choice(other, size=min(len(other), NB), replace=False))
print('columns: %d binary + %d non-binary' % (len(binary), len(keep)-len(binary)), flush=True)

active  = sorted(rng.choice(NCEN, size=min(2500, NCEN), replace=False))
MAXROWS = int(os.environ.get('CERT_MAXROWS', 11000))

def rows_for(colkeys, act):
    """the LP matrix over the ACTIVE centres only.

    THE ACTIVE SET IS THE POINT. This used to call bc(i), which computes (B, c) against every one of
    the nT centre types, and then index out the few thousand rows it wanted. At r=11 that was
    wasteful; at r=14 it is 1.59 million columns of work per LP column instead of 2,500, about 600
    times too much, and the run never printed a first round before the kernel killed it at 6 GB.
    Comparing against the active slice directly is the same arithmetic on 0.16 percent of the data."""
    ac = np.asarray(act, dtype=np.int64)
    at = (ac % nT); ap = (ac // nT)
    ta, isXa = types[at], isX[at]
    A = np.empty((len(act), len(colkeys)))
    for a,k in enumerate(colkeys):
        i,w = cols[k]; t = types[i]; tX = isX[i]
        Ba = (((~tX)&(~isXa)&(ta!=t)).sum(axis=1) + (tX ^ isXa).sum(axis=1)).astype(np.int32)
        ca = ((tX)&isXa).sum(axis=1).astype(np.int32)
        colv = np.empty(len(act))
        for p in np.unique(ap):
            sel = (ap==p); T = load_table(w, int(p))
            colv[sel] = T[Ba[sel], ca[sel]] / DG
        A[:,a]=colv
    return A

best = Fr(0)
n = len(keep)
for rnd in range(ROUNDS):
    A = rows_for(keep, active)
    Aub = np.hstack([A, -np.ones((len(active),1))])
    res = linprog(c=[0.0]*n+[1.0], A_ub=Aub, b_ub=np.zeros(len(active)),
                  A_eq=[[1.0]*n+[0.0]], b_eq=[float(S)],
                  bounds=[(0,None)]*n+[(0,None)], method='highs')
    if not res.success:
        print('LP failed:', res.message); sys.exit(1)
    t_act = res.x[-1]
    xnum = np.floor(np.asarray(res.x[:n])*Dx).astype(np.int64)
    nz = np.nonzero(xnum)[0]
    mass = Fr(int(xnum.sum()), Dx)

    # one EXACT sweep over every centre: max load, and the worst centres to add as rows
    # SWEEP, ORDERED BY COLUMN AND NOT BY FREE WEIGHT. The (B, c) arrays of a column against all nT
    # centre types cost AP comparisons over nT, which at r=15 is 67 million operations, and they do
    # not depend on the centre's free weight. Looping free weight outermost recomputed them L+1
    # times per column. Accumulating all L+1 totals in one pass over the columns computes them once,
    # at the cost of holding L+1 arrays instead of one: 230 MB at r=15, against an hour of arithmetic.
    worst = 0; cand = []; worst_at = None
    tots = [np.zeros(nT, dtype=np.int64) for _ in range(L+1)]
    # GROUPED BY COLUMN TYPE. (B, c) depends on the column's active-prefix TYPE, not on its free
    # weight, so the several support columns that share a type need it computed once. Ungrouped, the
    # loop recomputed it per column and the cache cannot help: at r=17 it holds three entries against
    # 43 million centre types, and each rebuild is 16 comparisons over all of them. Grouping removes
    # the thrash entirely rather than enlarging the cache.
    bytype = {}
    for a in nz:
        i, w = cols[int(keep[a])]
        bytype.setdefault(i, []).append((w, int(xnum[a])))
    # CHUNKED OVER CENTRE TYPES. Every temporary in this loop is transient and refcounted, yet RSS
    # climbed monotonically from 6.4 GB to 29.9 GB across one r=16 sweep, 3.3 GB per 20 s, and the
    # kernel killed it: BC alone builds several (nT, AP) boolean arrays, 229 MB each at r=16, and
    # allocating hundreds of those does not hand the pages back. Chunking does the same arithmetic
    # on slices, so the temporaries are CH-sized instead of nT-sized and the footprint is flat.
    # types[s0:s1] is a view, not a copy. The per-type (B, c) cache is then pointless and dropped,
    # which also frees its 688 MB.
    for s0 in range(0, nT, CHUNK):
        s1 = min(nT, s0 + CHUNK)
        tsl, isXs = types[s0:s1], isX[s0:s1]
        for i, items in bytype.items():
            t = types[i]; tX = isX[i]
            Bv = (((~tX) & (~isXs) & (tsl != t)).sum(axis=1)
                  + (tX ^ isXs).sum(axis=1)).astype(np.int32)
            cv = ((tX) & isXs).sum(axis=1).astype(np.int32)
            for w, xa in items:
                for p in range(L+1):
                    tots[p][s0:s1] += xa * load_table(w, p)[Bv, cv]
    for p in range(L+1):
        tot = tots[p]
        j = int(tot.argmax())
        if int(tot[j]) > worst: worst, worst_at = int(tot[j]), (j, p)
        # TAKE AT MOST AS MANY CANDIDATES AS THERE ARE CENTRE TYPES. At small active-prefix lengths
        # nT is tiny (243 at m=11 r=6, 81 at m=12 r=5) and a fixed 250 is out of bounds, which killed
        # both of those cells with a ValueError rather than a verdict. The small cells are exactly
        # the cheap ones, so the bug sat where it would be found last.
        kk = min(250, nT)
        top = np.arange(nT) if kk >= nT else np.argpartition(tot, -kk)[-kk:]
        cand += [(int(tot[i]), p*nT+int(i)) for i in top]
    lastload = tots
    Mex = Fr(worst, Dx*DG)                  # exact maximum load of x over ALL centres
    bound = mass / Mex                      # x/Mex is feasible everywhere -> tau* >= this
    best = max(best, bound)
    print('round %2d  rows %5d  t_active %9.6f  exact max load %10.7f  CERTIFIED tau* >= %12.8f  %5.0fs'
          % (rnd+1, len(active), t_act, float(Mex), float(bound), time.time()-t0), flush=True)
    if bound > R - 1:
        print('\n*** CERTIFIED tau*(%d,%d,%d) >= %s = %.9f > %d, so K(%d,%d,%d) = %d.'
              % (N, M, R, bound, float(bound), R-1, N, M, R, R))
        import json
        rec = {'N': N, 'm': M, 'r': R, 'AP': AP, 'L': L, 'Dx': Dx, 'Dglob': DG,
               'mass_numerator_over_Dx': int(xnum.sum()),
               'max_load_numerator_over_Dx_times_Dglob': int(worst),
               'max_load_at': [''.join('01X'[t] for t in types[worst_at[0]]), int(worst_at[1])],
               'bound_num': bound.numerator, 'bound_den': bound.denominator,
               'support': [[int(col_i[int(keep[a])]), int(col_w[int(keep[a])]), int(xnum[a])]
                           for a in nz],
               'note': 'support entries are [active-prefix type index into itertools.product([0,1,2],'
                       ' repeat=AP) in lexicographic order, free weight w, x numerator over Dx].'
                       ' tau* >= (sum numerators)*Dglob/max_load_numerator.'}
        # 2026-10-22: two concurrent r=14 runs both wrote this name and the weaker one landed last,
        # losing a 13.893239388 certificate. The configuration is part of the artefact's identity,
        # so it goes in the name: a run can no longer silently overwrite a different run's proof.
        fn = 'cert_%d_%d_%d_b%d_n%d_r%d.json' % (N, M, R, NBIN if NBIN < 10**9 else 0, NB, MAXROWS)
        open(fn,'w').write(json.dumps(rec))
        canon = 'cert_%d_%d_%d.json' % (N, M, R)
        keep = True
        if os.path.exists(canon):
            try:
                prev = json.load(open(canon))
                keep = (prev.get('r') != R) or (Fr(prev['bound_num'], prev['bound_den']) < bound)
            except Exception:
                keep = True
        if keep:
            open(canon,'w').write(json.dumps(rec))
            print('wrote %s and %s: %d support orbits' % (fn, canon, len(rec['support'])))
        else:
            print('wrote %s: %d support orbits (left %s alone, it holds a STRONGER bound)'
                  % (fn, len(rec['support']), canon))
        break
    cand.sort(reverse=True)
    newrows = [c for _,c in cand if c not in set(active)]
    if not newrows:
        print('\nseparation exhausted; tau* = %.8f on these columns' % float(bound)); break
    active = sorted(set(active) | set(newrows[:1200]))
    if len(active) > MAXROWS:                     # keep the binding ones
        # keep the binding rows. The loads are already in hand from the sweep above, indexed by
        # the same p*nT+i encoding, so this does not re-sweep.
        keepset = sorted(active, key=lambda c: -int(lastload[c // nT][c % nT]))[:MAXROWS]
        active = sorted(keepset)
print('\nbest certified lower bound this run: %s = %.9f  (needs > %d for K=%d)'
      % (best, float(best), R-1, R))
