"""Walk K(19,8,r) = r up the panel by row generation, with an EXACT proof at each cell.

The full LP is 2^(r-1) x (L+1) in both dimensions, so the dense matrix passes 2 GB around r = 12.
Row generation avoids it: solve on a subset of centre constraints, then verify the resulting dual
against EVERY centre in exact arithmetic. Violated centres are added and it re-solves. The proof at
the end is the exact check over all centres, so the subset never has to be the right one -- it only
has to converge.

Criterion: K >= ceil(LP*), so a feasible dual of exact mass > r - 1 proves K = r. Masses are rounded
DOWN to rationals, which can only lower what a ball collects, so feasibility is preserved.
"""
import numpy as np, itertools, math, sys, time
from fractions import Fraction as Fr
from math import comb
from scipy.optimize import linprog

N, M = 19, 8

def cell(R, D=10**7, max_rounds=24, seed_frac=0.12):
    AP = R-1; L = N-AP
    T=[''.join('1'*i+'0'*(AP-i)) for i in range(R)]
    dH=lambda a,b: sum(x!=y for x,y in zip(a,b))
    pre=[''.join(b) for b in itertools.product('01',repeat=AP)]
    dmin={y:min(dH(y,t) for t in T) for y in pre}
    orb=[(y,w) for y in pre for w in range(L+1) if w+dmin[y]<=M]
    if not orb: return None
    shell={w:comb(L,w)*3**w for w in range(L+1)}
    _f={}
    def Fc(p,w,s):
        if s<0: return 0
        k=(p,w,s)
        if k in _f: return _f[k]
        t=0
        for o in range(0,min(p,w)+1):
            b=comb(p,o)*comb(L-p,w-o)
            if b: t+=b*sum(comb(o,e)*2**(o-e) for e in range(o+1) if w+p-o-e<=s)*3**(w-o)
        _f[k]=t; return t
    # the matrix entry depends only on (p, w, prefix distance), so cache by that
    def col(u, p):
        return np.array([ (Fc(p,w,M-dH(y,u))/shell[w] if M-dH(y,u)>=0 else 0.0) for (y,w) in orb])
    allcen=[(u,p) for u in pre for p in range(L+1)]
    rng=np.random.default_rng(0)
    idx=list(rng.choice(len(allcen), size=max(40,int(len(allcen)*seed_frac)), replace=False))
    active={allcen[i] for i in idx}
    active |= {(t, 0) for t in T} | {(pre[0], 0), (pre[-1], 0)}
    for rnd in range(max_rounds):
        act=sorted(active)
        A=np.array([col(u,p) for (u,p) in act])
        res=linprog(c=-np.ones(len(orb)), A_ub=A, b_ub=np.ones(len(act)),
                    bounds=[(0,None)]*len(orb), method='highs')
        if not res.success: return ('lpfail', res.message)
        x=[Fr(int(math.floor(v*D)), D) for v in res.x]
        nz=[(j,x[j]) for j in range(len(orb)) if x[j]]
        total=sum(v for _,v in nz)
        # EXACT check over every centre; collect violators
        viol=[]; worst=Fr(0)
        for (u,p) in allcen:
            m=Fr(0)
            for j,xi in nz:
                y,w=orb[j]; s=M-dH(y,u)
                if s>=0: m += xi*Fr(Fc(p,w,s), shell[w])
            if m>worst: worst=m
            if m>1: viol.append((u,p))
        if not viol:
            return ('ok', float(-res.fun), total, worst, len(nz), len(act), len(allcen), rnd+1)
        active |= set(viol[:400])
    return ('noconv', float(-res.fun), total, worst, len(viol))

for R in (11, 12, 13):
    t0=time.time()
    out=cell(R)
    if out is None: print('r=%d: no support' % R); continue
    if out[0]=='ok':
        _, fv, total, worst, nz, nact, ncen, rnds = out
        need=R-1
        print('r=%2d  LP*~%.6f  exact dual mass %s = %.7f  max ball %s' %
              (R, fv, total, float(total), '%.9f'%float(worst)))
        print('      %d support orbits, %d of %d centres in the active set, %d rounds, %.0fs' %
              (nz, nact, ncen, rnds, time.time()-t0))
        print('      mass > r-1 = %d: %s   ceil = %d   %s' %
              (need, total>need, math.ceil(total),
               '*** PROVED K(19,8,%d) = %d' % (R,R) if (worst<=1 and total>need) else 'NOT proved'))
    else:
        print('r=%2d  %s %s' % (R, out[0], out[1:]))
    print(flush=True)
