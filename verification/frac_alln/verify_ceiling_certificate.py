"""Sampled consistency check of a solver ceiling certificate (cert_N_m_r.json). SHARES NO CODE WITH THE
SOLVER. It does NOT establish the maximum load: verify_ceiling_exhaustive.py does that, over every class.

Five row-end certificates were withdrawn on 2026-10-22 because they were verified by the routine
that produced them, and that routine had restricted the constraint set in a way that moved the
bound the wrong way. The defence is a checker built from a different algorithm, so a bug has to
occur twice, independently, in the same direction to survive.

The solver computes a ball's load through F(p,w,s) and a cached (B,c) table over orbit labels.
This file never calls either. It takes EXPLICIT strings over {0,1,2,3}^19, and counts the points of
a support orbit inside an explicit ball by a dynamic program over positions: active positions
contribute a fixed or two-way cost, free positions are a size-w subset with uniform non-zero
letters, and the DP tracks (non-zeros placed, distance so far). All integer arithmetic, so the
fraction it returns is exact.

Three things are checked.
  1. ORBIT INVARIANCE. Two explicit centres in the same (prefix type, free weight) class must have
     EXACTLY equal load. The reduction from 4^N centres to (L+1) x 3^AP classes holds by construction of the blocks;
     this test checks the implementation;
     if it fails the sweep was never a sweep and the bound is void.
  2. THE CEILING HOLDS. No centre sampled, and none of the r library centres, may carry load above
     the certificate's claimed maximum. One counterexample destroys the certificate.
  3. THE ARITHMETIC. The claimed bound must equal (sum of x numerators) * Dglob / (max-load
     numerator), recomputed here from the stored integers.
"""
import json, sys, itertools, random
from fractions import Fraction as Fr

C = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'cert_19_8_11.json'))
N, m, r, AP, L, Dx, DG = C['N'], C['m'], C['r'], C['AP'], C['L'], C['Dx'], C['Dglob']
assert AP + L == N and AP == r - 1
TYPES = list(itertools.product([0, 1, 2], repeat=AP))      # the certificate's stated encoding
LIB = [tuple([1]*k + [0]*(N-k)) for k in range(r)]         # 1^k 0^(N-k), the r deletion patterns

def orbit_counts(typ, w, u):
    """Exact count of orbit-(typ,w) points at each Hamming distance from the explicit centre u.

    Built as two independent distributions and convolved. Nothing here resembles F(p,w,s)."""
    # active prefix, positions 0..AP-1
    act = {0: 1}
    for j in range(AP):
        nxt = {}
        if typ[j] in (0, 1):
            cost = 0 if u[j] == typ[j] else 1
            for d, c in act.items(): nxt[d+cost] = nxt.get(d+cost, 0) + c
        else:                                    # z_j ranges over the two non-library letters
            for letter in (2, 3):
                cost = 0 if u[j] == letter else 1
                for d, c in act.items(): nxt[d+cost] = nxt.get(d+cost, 0) + c
        act = nxt
    # free suffix, positions AP..N-1: exactly w non-zeros, each uniform over {1,2,3}
    free = {(0, 0): 1}
    for j in range(AP, N):
        nxt = {}
        for (k, d), c in free.items():
            cz = 0 if u[j] == 0 else 1                       # this position stays at 0
            nxt[(k, d+cz)] = nxt.get((k, d+cz), 0) + c
            if k < w:                                        # this position carries a non-zero
                for letter in (1, 2, 3):
                    cl = 0 if u[j] == letter else 1
                    nxt[(k+1, d+cl)] = nxt.get((k+1, d+cl), 0) + c
        free = nxt
    out = {}
    for (k, d2), c2 in free.items():
        if k != w: continue
        for d1, c1 in act.items(): out[d1+d2] = out.get(d1+d2, 0) + c1*c2
    return out

def load(u, support):
    """load of the certificate's measure on the explicit ball of radius m about u, exactly."""
    tot = Fr(0)
    for i, w, num in support:
        cnt = orbit_counts(TYPES[i], w, u)
        inside = sum(c for d, c in cnt.items() if d <= m)
        size = sum(cnt.values())      # the DP enumerates the whole orbit, so this IS its size
        tot += Fr(num, Dx) * Fr(inside, size)
    return tot

sup = C['support']
print('certificate: %d support orbits, mass %s, claimed max load %s'
      % (len(sup), Fr(C['mass_numerator_over_Dx'], Dx),
         Fr(C['max_load_numerator_over_Dx_times_Dglob'], Dx*DG)))

# 3. arithmetic, recomputed from the stored integers
mass = Fr(sum(n for _, _, n in sup), Dx)
assert mass == Fr(C['mass_numerator_over_Dx'], Dx), 'stored mass disagrees with the support vector'
Mex = Fr(C['max_load_numerator_over_Dx_times_Dglob'], Dx*DG)
bound = mass / Mex
claimed = Fr(C['bound_num'], C['bound_den'])
print('recomputed bound %s = %.9f   claimed %s   AGREE %s'
      % (bound, float(bound), claimed, bound == claimed))
if bound != claimed: sys.exit('ARITHMETIC MISMATCH')

rng = random.Random(12345)
def centre_in_class(typ, p):
    """an explicit centre whose class is (typ, p), chosen at random inside the class"""
    u = []
    for t in typ: u.append(t if t in (0,1) else rng.choice([2,3]))
    free = [0]*L
    for j in rng.sample(range(L), p): free[j] = rng.choice([1,2,3])
    return tuple(u + free)

# 1. orbit invariance: same class, different representatives, must give EXACTLY equal load
print('\norbit invariance (same class, different explicit centres):')
bad = 0
for trial in range(6):
    typ = tuple(rng.choice([0,1,2]) for _ in range(AP)); p = rng.randrange(L+1)
    a, b = centre_in_class(typ, p), centre_in_class(typ, p)
    while b == a: b = centre_in_class(typ, p)
    la, lb = load(a, sup), load(b, sup)
    ok = (la == lb); bad += (not ok)
    print('  class(nonbinary %d, p=%d)  load %.9f vs %.9f  EQUAL %s'
          % (sum(1 for t in typ if t==2), p, float(la), float(lb), ok))
if bad: sys.exit('ORBIT INVARIANCE FAILS: the reduction to %d classes is not valid, bound is void'
                 % ((L+1)*3**AP))

# 2. the ceiling: nothing may exceed the claimed maximum
print('\nceiling check against the claimed maximum %s:' % float(Mex))
worst = Fr(0); worstu = None
cands = [tuple(x) for x in LIB]
for _ in range(90):
    typ = tuple(rng.choice([0,1,2]) for _ in range(AP))
    cands.append(centre_in_class(typ, rng.randrange(L+1)))
for u in cands:
    lu = load(u, sup)
    if lu > worst: worst, worstu = lu, u
    if lu > Mex:
        sys.exit('COUNTEREXAMPLE: centre %s carries load %s > claimed max %s. CERTIFICATE BROKEN.'
                 % (''.join(map(str,u)), lu, Mex))
print('  %d centres checked (%d library + %d sampled); worst found %.9f, claimed max %.9f'
      % (len(cands), r, len(cands)-r, float(worst), float(Mex)))
print('  library centres carry %s' % [float(load(tuple(x), sup)) for x in LIB[:4]])
# 4. the stored witness must carry the claimed maximum exactly (a claimed maximum below the true one
#    is caught here whenever the witness is the true maximiser; the exhaustive sweep catches the rest)
wp, wf = C['max_load_at']
wu = tuple([{'0': 0, '1': 1, 'X': 2}[ch] for ch in wp] + [1]*int(wf) + [0]*(L-int(wf)))
lw = load(wu, sup)
print('\nstored witness %s, free weight %d: load %s, claimed max %s, EQUAL %s' % (wp, wf, lw, Mex, lw == Mex))
if lw != Mex:
    sys.exit('WITNESS MISMATCH: the stored maximiser carries %s, not the claimed %s.' % (lw, Mex))
print('\nSAMPLED CONSISTENCY CHECKS PASS (arithmetic, invariance, ceiling on %d centres, witness).' % len(cands))
print('This does not establish the maximum load. If verify_ceiling_exhaustive.py confirms it,')
print('tau*(%d,%d,%d) >= %s = %.9f, and ceil = %d.'
      % (N, m, r, bound, float(bound), -(-bound.numerator//bound.denominator)))
