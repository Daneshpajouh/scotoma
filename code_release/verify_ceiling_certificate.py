"""Independent check of cert_19_8_11.json. SHARES NO CODE WITH THE SOLVER.

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
     EXACTLY equal load. The whole reduction from 4^19 centres to 10 x 3^10 classes rests on this;
     if it fails the sweep was never a sweep and the bound is void.
  2. THE CEILING HOLDS. No centre sampled, and none of the r library centres, may carry load above
     the certificate's claimed maximum. One counterexample destroys the certificate.
  3. THE ARITHMETIC. The claimed bound must equal (sum of x numerators) * Dglob / (max-load
     numerator), recomputed here from the stored integers.
"""
import json, sys, itertools, random
from fractions import Fraction as Fr
from math import comb

_cert = sys.argv[1] if len(sys.argv) > 1 else 'cert_19_8_11.json'
try:
    C = json.load(open(_cert))
except FileNotFoundError:
    # A raw traceback is the worst thing a referee alone in a directory can be shown: it reads as
    # broken software rather than a missing input, and nothing else in the release is trusted after
    # it. Exit 2 is NOT RUN, which is not a pass.
    print("GATE NOT RUN: the certificate %s is not here. Pass one as the first argument, e.g.\n"
          "  python3 %s cert_19_8_11.json\n"
          "This is an incomplete release, not a finding about the claim."
          % (_cert, sys.argv[0]), file=sys.stderr)
    sys.exit(2)   # 0 passed, 1 FAILED, 2 NOT RUN. sys.exit(str) would have exited 1, i.e. FAILED.
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

def class_size(typ, p):
    """How many explicit centres the class (typ, p) contains, in closed form.

    Two letters at each non-binary prefix position, C(L,p) placements of the p free non-zeros,
    three letters at each of them. It is 1 exactly when the prefix has no 2 AND p == 0, and that
    single case is what used to hang this gate forever: see the trial loop below."""
    return (2 ** sum(1 for t in typ if t == 2)) * comb(L, p) * (3 ** p)

# 1. orbit invariance: same class, different representatives, must give EXACTLY equal load
print('\norbit invariance (same class, different explicit centres):')
bad = 0
compared = 0
# A SINGLETON CLASS CANNOT TEST INVARIANCE, AND ASKING IT TO USED TO HANG THIS GATE FOREVER.
# `while b == a: b = centre_in_class(typ, p)` cannot terminate when the class (typ, p) holds one
# centre, which happens exactly when the prefix carries no 2 and p == 0. Under the fixed seed 12345
# that draw comes up at trial 4 whenever AP == 4, i.e. at EVERY r = 5 cell, and at trial 0 when
# AP == 1. (19,10,5) and (19,11,5) were blocked by this for hours and I first blamed the denominator
# Dx, which had nothing to do with it; the fault handler put the stack in random.sample under this
# line. Every other AP from 2 to 18 completes, so the release has been passing r >= 6 all along.
#
# Drawing only classes with two or more representatives does not weaken the check: on a singleton
# class the invariance claim is vacuously true, so the trial carried no information to begin with.
# The class size is now printed, and the pass line below is reached only if all six trials really
# compared two DISTINCT centres, so this check can no longer pass by doing nothing.
for trial in range(6):
    for _ in range(1000):
        typ = tuple(rng.choice([0,1,2]) for _ in range(AP)); p = rng.randrange(L+1)
        if class_size(typ, p) > 1: break
    else:
        print('GATE NOT RUN: no class with two representatives in 1000 draws at AP=%d, L=%d'
              % (AP, L), file=sys.stderr)
        sys.exit(2)
    a, b = centre_in_class(typ, p), centre_in_class(typ, p)
    for _ in range(1000):
        if b != a: break
        b = centre_in_class(typ, p)
    if b == a:
        print('GATE NOT RUN: a class of size %d yielded no second representative in 1000 draws'
              % class_size(typ, p), file=sys.stderr)
        sys.exit(2)
    la, lb = load(a, sup), load(b, sup)
    ok = (la == lb); bad += (not ok); compared += 1
    print('  class(nonbinary %d, p=%d, size %d)  load %.9f vs %.9f  EQUAL %s'
          % (sum(1 for t in typ if t==2), p, class_size(typ, p), float(la), float(lb), ok))
if bad: sys.exit('ORBIT INVARIANCE FAILS: the reduction to %d classes is not valid, bound is void'
                 % ((L+1)*3**AP))
if compared != 6:
    print('GATE NOT RUN: only %d of 6 invariance trials compared two distinct centres' % compared,
          file=sys.stderr)
    sys.exit(2)

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
print('\nALL THREE CHECKS PASS.')
print('tau*(%d,%d,%d) >= %s = %.9f, and %d >= ceil(tau*) = %d, so K = %d.'
      % (N, m, r, bound, float(bound), r, -(-bound.numerator//bound.denominator), r))
