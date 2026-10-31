"""The last row-end, (m,r) = (8,11), by CONSTRAINT generation on the primal.

The dense primal is 10,240 orbit constraints by 10,240 centre variables, which killed the batch
route. So: solve on a SUBSET of the orbit constraints, then check the resulting cover against EVERY
orbit in exact arithmetic and add the violated ones. Dropping constraints can only LOWER the primal
optimum, so a subset solution is not automatically feasible; the exact check over all orbits is what
makes it so, and only a cover with zero violations is used.

Rounded UP to rationals, which preserves feasibility and can only raise the objective. Proved if
ceil of the exact value is below r.
"""
import itertools, json, math, os, sys, time
from fractions import Fraction as Fr
from math import comb
import numpy as np
from scipy.optimize import linprog

N, M, R = 19, 8, 11
AP, L = R-1, N-(R-1)
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
def fr(y,w,u,p):
    s=M-dH(y,u)
    return Fr(Fc(p,w,s), shell[w]) if s>=0 else Fr(0)
print('orbits %d, centres %d, free length %d' % (len(orb), len(cen), L), flush=True)
rng=np.random.default_rng(0)
act=set(rng.choice(len(orb), size=600, replace=False).tolist())
act |= {j for j,(y,w) in enumerate(orb) if w in (0, M)}
D=10**6
for rnd in range(40):
    A=np.array([[float(fr(orb[j][0],orb[j][1],u,p)) for (u,p) in cen] for j in sorted(act)])
    res=linprog(c=np.ones(len(cen)), A_ub=-A, b_ub=-np.ones(len(act)),
                bounds=[(0,None)]*len(cen), method='highs')
    if not res.success:
        print('LP failed:', res.message); sys.exit(1)
    x=[Fr(int(math.ceil(v*D)), D) for v in res.x]
    nz=[(i,x[i]) for i in range(len(cen)) if x[i]]
    total=sum(v for _,v in nz)
    viol=[]
    for j,(y,w) in enumerate(orb):
        s=sum(xi*fr(y,w,cen[i][0],cen[i][1]) for i,xi in nz)
        if s<1: viol.append(j)
    print('round %2d  active %5d  cover value %.6f  support %4d  violations %5d'
          % (rnd+1, len(act), float(total), len(nz), len(viol)), flush=True)
    if not viol:
        c={'N':N,'m':M,'r':R,'free_length':L,
           'status':'ROW_END_PROVED' if math.ceil(total)<R else 'not_proved',
           'claim':'no fractional dual can prove K=r here, because ceil(tau*) < r',
           'exact_cover_value':str(total),'exact_cover_value_float':float(total),
           'ceiling':math.ceil(total),'needed_for_proof':'ceiling < %d'%R,
           'orbit_constraints_checked':len(orb),'violations':0,
           'support_centres':len(nz),'generation_rounds':rnd+1,
           'active_constraints_at_end':len(act),
           'rounding':'primal rounded UP to denominator %d, which preserves feasibility'%D,
           'method':'constraint generation on the primal; the exact check is over ALL orbits'}
        # resolve against THIS file, never an absolute home path: a released script carrying
        # one is unrunnable for every reader and names the machine user.
        out=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'data','floor_certificates','rowend_m8_r11.json')
        json.dump(c, open(out,'w'), indent=1)
        print()
        print('EXACT cover value %s = %.7f, ceil %d, needed < %d -> %s'
              % (total, float(total), math.ceil(total), R, c['status']))
        print('written', out)
        break
    act |= set(viol[:500])
else:
    print('did not converge in 40 rounds')
