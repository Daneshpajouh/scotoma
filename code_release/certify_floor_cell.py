"""Extract the r=8, m=8 dual certificate in EXACT rationals and verify it exactly.

The float LP says the optimum is 8.000000. A float is not a proof. So: take the LP's support,
re-solve the equality system on that support in exact rationals, then verify total mass = 8 and
every one of the 1664 ball constraints <= 1 with Fraction arithmetic.
"""
import numpy as np, itertools, sys
from fractions import Fraction as Fr
from math import comb
from scipy.optimize import linprog

N, R, M = 19, 8, 8
AP, = (R-1,)
L = N - AP
T = [''.join('1'*i + '0'*(AP-i)) for i in range(R)]
def dH(a,b): return sum(x!=y for x,y in zip(a,b))
prefixes = [''.join(b) for b in itertools.product('01', repeat=AP)]
dmin = {y: min(dH(y,t) for t in T) for y in prefixes}
orbits = [(y,w) for y in prefixes for w in range(L+1) if w + dmin[y] <= M]
shell = {w: comb(L,w)*3**w for w in range(L+1)}

_fc = {}
def Fcount(p,w,s):
    if s < 0: return 0
    k=(p,w,s)
    if k in _fc: return _fc[k]
    tot=0
    for o in range(0,min(p,w)+1):
        base=comb(p,o)*comb(L-p,w-o)
        if base:
            tot += base*sum(comb(o,e)*2**(o-e) for e in range(o+1) if w+p-o-e<=s)*3**(w-o)
    _fc[k]=tot
    return tot

centres = [(u,p) for u in prefixes for p in range(L+1)]
A = np.array([[ (Fcount(p,w,M-dH(y,u))/shell[w] if M-dH(y,u)>=0 else 0.0)
                for (y,w) in orbits] for (u,p) in centres])
res = linprog(c=-np.ones(len(orbits)), A_ub=A, b_ub=np.ones(len(centres)),
              bounds=[(0,None)]*len(orbits), method='highs')
sup = [j for j in range(len(orbits)) if res.x[j] > 1e-9]
tight = [i for i in range(len(centres)) if abs(A[i] @ res.x - 1.0) < 1e-7]
print('float optimum %.9f, support %d orbits, %d tight constraints' % (-res.fun, len(sup), len(tight)))

# exact: solve the tight system on the support for mass vector x with A_tight x = 1
import fractions
def Fexact(p,w,s):
    return Fr(Fcount(p,w,s), shell[w])
rows = []
for i in tight:
    u,p = centres[i]
    rows.append([Fexact(p, orbits[j][1], M-dH(orbits[j][0],u)) if M-dH(orbits[j][0],u)>=0 else Fr(0)
                 for j in sup])
# least-squares is not exact; instead pick a square independent subsystem by Gaussian elimination
import copy
Mx = [r[:] + [Fr(1)] for r in rows]
n = len(sup)
piv_rows, piv_cols = [], []
rr = 0
for c in range(n):
    pr = None
    for i in range(rr, len(Mx)):
        if Mx[i][c] != 0: pr = i; break
    if pr is None: continue
    Mx[rr], Mx[pr] = Mx[pr], Mx[rr]
    pivv = Mx[rr][c]
    Mx[rr] = [v/pivv for v in Mx[rr]]
    for i in range(len(Mx)):
        if i != rr and Mx[i][c] != 0:
            f = Mx[i][c]
            Mx[i] = [a - f*b for a,b in zip(Mx[i], Mx[rr])]
    piv_rows.append(rr); piv_cols.append(c); rr += 1
    if rr == n: break
x = [Fr(0)]*n
for k,c in enumerate(piv_cols):
    x[c] = Mx[piv_rows[k]][-1]
# verify: consistency on ALL tight rows, nonnegativity, total, and every ball constraint
consistent = all(sum(a*b for a,b in zip(r, x)) == 1 for r in rows)
nonneg = all(v >= 0 for v in x)
total = sum(x)
worst = Fr(0); arg=None
for i,(u,p) in enumerate(centres):
    m = Fr(0)
    for k,j in enumerate(sup):
        y,w = orbits[j]
        s = M - dH(y,u)
        if s >= 0 and x[k]:
            m += x[k] * Fexact(p,w,s)
    if m > worst: worst, arg = m, (u,p)
print()
print('EXACT CHECK of the reconstructed certificate')
print('  tight rows all satisfied exactly : %s' % consistent)
print('  all masses nonnegative           : %s' % nonneg)
print('  total mass                       : %s = %s' % (total, float(total)))
print('  max ball mass over all %d centres: %s = %.12f at %s' % (len(centres), worst, float(worst), arg))
print()
if consistent and nonneg and total >= R and worst <= 1:
    print('VERDICT: K(%d, %d, %d) = %d. PROVED in exact rationals.' % (N, M, R, R))
    print('  support, mass per orbit:')
    for k,j in enumerate(sup):
        if x[k]: print('    %-8s w=%-2d  %s' % (orbits[j][0], orbits[j][1], x[k]))
else:
    print('VERDICT: reconstruction failed; the float optimum is not yet a proof.')
