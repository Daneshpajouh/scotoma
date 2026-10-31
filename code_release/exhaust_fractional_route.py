"""Is LP* = r FALSE at (N,m,r) = (19,8,9)? Prove it by an EXACT rational fractional cover.

To prove the floor you need a dual of mass >= r. To prove NO dual can, you need an upper bound on
the LP optimum, i.e. a feasible PRIMAL fractional cover of value < r. Floats cannot prove that, so:
solve the primal in floats, round every variable UP to a rational (which preserves feasibility and
can only increase the objective), then verify feasibility and the objective in exact arithmetic.
If the rounded-up cover still has value < r, LP* < r is PROVED.
"""
import numpy as np, itertools
from fractions import Fraction as Fr
from math import comb
from scipy.optimize import linprog

N = 19
def run(R, M):
    AP = R-1; L = N-AP
    T=[''.join('1'*i+'0'*(AP-i)) for i in range(R)]
    dH=lambda a,b: sum(x!=y for x,y in zip(a,b))
    pre=[''.join(b) for b in itertools.product('01',repeat=AP)]
    dmin={y:min(dH(y,t) for t in T) for y in pre}
    orb=[(y,w) for y in pre for w in range(L+1) if w+dmin[y]<=M]
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
    cen=[(u,p) for u in pre for p in range(L+1)]
    # exact coverage fraction of orbit (y,w) by a centre class (u,p)
    def frac(y,w,u,p):
        s=M-dH(y,u)
        return Fr(Fc(p,w,s), shell[w]) if s>=0 else Fr(0)
    Af=np.array([[float(frac(y,w,u,p)) for (y,w) in orb] for (u,p) in cen])
    # PRIMAL: minimise sum x_c subject to, for each orbit, sum_c frac * x_c >= 1
    res=linprog(c=np.ones(len(cen)), A_ub=-Af.T, b_ub=-np.ones(len(orb)),
                bounds=[(0,None)]*len(cen), method='highs')
    if not res.success: return None
    # round UP to rationals with denominator D: feasibility preserved, objective only grows
    D=10**6
    x=[Fr(int(np.ceil(v*D)), D) for v in res.x]
    total=sum(x)
    # exact feasibility check on every orbit
    viol=0; worst=None
    for (y,w) in orb:
        s=sum(xi*frac(y,w,u,p) for xi,(u,p) in zip(x,cen) if xi)
        if s < 1:
            viol+=1
            if worst is None or s<worst: worst=s
    return float(res.fun), total, viol, worst, len(orb), len(cen)

for (R,M) in [(9,8),(10,8)]:
    out=run(R,M)
    if out is None: print('r=%d m=%d LP failed'%(R,M)); continue
    f,total,viol,worst,no,nc=out
    print('r=%2d m=%d  orbits=%-5d centres=%-5d' % (R,M,no,nc))
    print('   float primal optimum        : %.6f' % f)
    print('   rounded-UP exact cover value: %s = %.6f' % (total, float(total)))
    print('   exact feasibility violations: %d%s' % (viol, '' if not viol else '  (worst %s)'%worst))
    if viol==0 and total < R:
        print('   *** PROVED: LP* <= %s < %d, so NO fractional dual can prove the floor here.' % (total,R))
    elif viol:
        print('   inconclusive: the rounded cover is infeasible, needs a finer denominator')
    else:
        print('   inconclusive: rounded value %s is not below %d' % (total,R))
    print(flush=True)
