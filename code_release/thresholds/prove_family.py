#!/usr/bin/env python3
"""Exact infinite-range certificates for constant-prefix packing families.

For each residue and bounded tail distance, certify a degree-p polynomial by
nonnegative forward differences after an explicitly checked finite prefix.
An elementary symmetric mean moment bound covers all remaining distances.
Standard library only; input family weights are checked, never trusted.
"""
import sys
sys.dont_write_bytecode=True
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb,lcm
from pathlib import Path
import hashlib
import json
import time
from verify import Geometry,require

ROOT=Path(__file__).resolve().parent


def choose(n,k):
    return comb(n,k) if 0<=k<=n else 0


def differences(v):
    out=[]
    while v:
        out.append(v[0])
        v=[b-a for a,b in zip(v,v[1:])]
    return out


@lru_cache(None)
def favourable(p,a,k):
    return 3**(p-a)*sum(comb(a,e)*2**(a-e) for e in range(max(0,p+k-a),a+1))


def profiles(family):
    r=family['r']
    G=Geometry(2*r,r,r)
    support=[]
    for term in family['support']:
        a=tuple(2 if x=='X' else int(x) for x in term['prefix'])
        d=term['distance']
        mass=F(term['mass'])
        require(G.dmin[a]==d,'family distance is false')
        require(((a,G.m-d)) in G.regions,'family prefix is not uniquely nearest')
        require(mass>0,'family mass not positive')
        support.append((a,d,mass))
    require(len({a for a,d,v in support})==len(support),'duplicate family prefix')
    sums=[F(0)]*r
    for a,d,v in support:
        sums[G.regions[a,G.m-d]]+=v
    require(sums==[F(1)]*r,'family region masses not one')
    dmax=max(d for a,d,v in support)
    width=r
    scale=lcm(*(v.denominator*2**a.count(2) for a,d,v in support))
    result={}
    for c in G.prefixes:
        h=[0]*((dmax+1)*width)
        for a,d,v in support:
            opts=[(2,3) if x==2 else (x,) for x in a]
            words=tuple(product(*opts))
            mul=v*F(scale,len(words))
            require(mul.denominator==1,'profile scale not integral')
            for b in words:
                delta=sum(x!=y for x,y in zip(c,b))
                h[d*width+delta]+=mul.numerator
        result[c]=tuple(h)
    unique=sorted(set(result.values()))
    cumulative={h:tuple(sum(h[d*width:d*width+j+1])
                         for d in range(dmax+1) for j in range(width)) for h in unique}
    maximal=[h for h in unique if not any(h!=q and all(a<=b for a,b in zip(cumulative[h],cumulative[q]))
                                         for q in unique)]
    domination={''.join(map(str,c)):next(i for i,q in enumerate(maximal)
                                       if all(a<=b for a,b in zip(cumulative[h],cumulative[q])))
                for c,h in result.items()}
    return support,result,maximal,domination,scale,dmax


def moment_proof(actual,scale,r,dmax,beta):
    best=None
    # Choosing a cutoff is exploration only. Each accepted bound is exact.
    for t in (F(7,5),F(3,2),F(4,3),F(5,4),F(8,5)):
        V=max(sum((F(x,scale)*t**(d-delta)
                   for d in range(dmax+1) for delta,x in enumerate(h[d*r:(d+1)*r])),F(0))
              for h in actual.values())
        for cutoff in range(max(20,r),251):
            rho=F(3,5)+F(beta,5*cutoff)
            h=(1-rho)/t+rho*(2+t)/3
            if 0<h<1 and V*h**cutoff < F(99,100):
                candidate=(cutoff,t,V,h,rho)
                if best is None or cutoff<best[0]:
                    best=candidate
                break
    require(best is not None,'no large-distance moment proof found')
    cutoff,t,V,h,rho=best
    return {'cutoff':cutoff,'t':str(t),'prefix_potential':str(V),
            'rho_upper':str(rho),'moment_factor':str(h),
            'bound':str(V*h**cutoff),'bound_less_than_99_over_100':True}


def pvalues(s,p,q,r,beta,scale,dmax,maximal,inverse_C0=None):
    C=5*q+s
    m=3*q+(3*s+beta)//5
    den=choose(C,p)*3**p
    R=[]
    for d in range(dmax+1):
        w=m-d
        require(0<=w<=C,'polynomial evaluation outside tail-weight domain')
        coeff=[choose(w,a)*choose(C-w,p-a) for a in range(p+1)]
        for delta in range(r):
            k=delta-d
            R.append(sum(x*favourable(p,a,k) for a,x in enumerate(coeff)))
    if inverse_C0 is None:
        return [scale*den-sum(x*y for x,y in zip(profile,R)) for profile in maximal]
    width=len(R)
    return [scale*C*den-(C-inverse_C0)*sum(x*y for x,y in zip(profile[:width],R))
            -inverse_C0*sum(x*y for x,y in zip(profile[width:],R)) for profile in maximal]


def inverse_profiles(family):
    maps=[]
    scales=[]
    ds=[]
    for key in ('mass_limit','mass_at_C0'):
        subset=dict(family,support=[dict(t,mass=t[key]) for t in family['support'] if F(t[key])])
        support,actual,maximal,domination,scale,dmax=profiles(subset)
        maps.append(actual)
        scales.append(scale)
        ds.append(dmax)
    scale=lcm(*scales)
    dmax=max(ds)
    r=family['r']
    padded=[]
    for actual,old,d in zip(maps,scales,ds):
        padded.append({c:tuple(x*(scale//old) for x in h)+(0,)*((dmax-d)*r)
                       for c,h in actual.items()})
    joint={c:padded[0][c]+padded[1][c] for c in maps[0]}
    unique=sorted(set(joint.values()))
    cumulative={h:tuple(sum(h[d*r:d*r+j+1])
                         for d in range(2*(dmax+1)) for j in range(r)) for h in unique}
    maximal=[h for h in unique if not any(h!=q and all(a<=b for a,b in zip(cumulative[h],cumulative[q]))
                                         for q in unique)]
    domination={''.join(map(str,c)):next(i for i,q in enumerate(maximal)
                                       if all(a<=b for a,b in zip(cumulative[h],cumulative[q])))
                for c,h in joint.items()}
    moment_actual={(endpoint,c):h for endpoint,actual in enumerate(padded) for c,h in actual.items()}
    return moment_actual,joint,maximal,domination,scale,dmax


def prove(family_path,min_n=None):
    start=time.monotonic()
    family=json.loads(family_path.read_text())
    r,alpha=family['r'],family['alpha']
    min_n=family['min_n_candidate'] if min_n is None else min_n
    beta=3*r-alpha
    require(beta>=0,'moment bound assumes nonnegative beta')
    inverse_C0=family.get('C0') if family.get('mass_type')=='linear_inverse_C' else None
    if inverse_C0 is None:
        support,actual,maximal,domination,scale,dmax=profiles(family)
        moment_actual=actual
    else:
        moment_actual,actual,maximal,domination,scale,dmax=inverse_profiles(family)
        require(min_n-r>=inverse_C0>0,'inverse-C mixture is outside its convex domain')
    require((3*(min_n-r)+beta)//5>=dmax,'initial radius below supported distance')
    moment=moment_proof(moment_actual,scale,r,dmax,beta)
    cutoff=moment['cutoff']
    print(json.dumps({'stage':'profiles','r':r,'unique':len(set(actual.values())),
                      'maximal':len(maximal),'scale':scale,'moment':moment}),flush=True)
    rows=[]
    ncoeff=0
    nfinite=0
    shiftmax=0
    residues=family.get('residues',list(range(5)))
    for s in residues:
        for p in range(cutoff):
            q0=max(0,(min_n-r-s+4)//5,(p-s+4)//5)
            # A finite shift is permissible only after checking each skipped
            # integer q. It does not infer positivity from sampled values.
            cache={}
            def get(q):
                if q not in cache:
                    cache[q]=pvalues(s,p,q,r,beta,scale,dmax,maximal,inverse_C0)
                return cache[q]
            degree=p+(inverse_C0 is not None)
            vals=[get(q0+j) for j in range(degree+1)]
            for i,pr in enumerate(maximal):
                shift=0
                coeff=differences([v[i] for v in vals])
                while any(F(v)<0 for v in coeff):
                    require(F(get(q0+shift)[i])>=0,
                            'actual negative polynomial: residue %d p %d row %d q %d'%(s,p,i,q0+shift))
                    shift+=1
                    require(shift<=100,'Newton positivity not certified within searched shifts')
                    coeff=differences([get(q0+shift+j)[i] for j in range(degree+1)])
                skipped=[get(q0+j)[i] for j in range(shift)]
                require(all(F(v)>=0 for v in skipped),'negative skipped polynomial')
                extra=get(q0+shift+degree+1)[i]
                require(sum(x*comb(degree+1,k) for k,x in enumerate(coeff))==extra,
                        'polynomial degree/reconstruction failed')
                ncoeff+=len(coeff)
                nfinite+=len(skipped)
                shiftmax=max(shiftmax,shift)
                rows.append({'residue':s,'p':p,'profile':i,'q0':q0,'shift':shift,
                             'finite_values':list(map(str,skipped)),
                             'newton_coefficients':list(map(str,coeff))})
            if p%10==0:
                print(json.dumps({'stage':'polynomials','residue':s,'p':p,'rows':len(rows)}),flush=True)
    record={'status':'PROVED: exact computer-assisted infinite-range packing certificate',
            'r':r,'alpha':alpha,'beta':beta,'min_n':min_n,'scale':scale,'dmax':dmax,
            'mass_type':family.get('mass_type','constant'),'C0':inverse_C0,
            'residues':residues,
            'family':family['support'],'prefix_count':3**(r-1),
            'unique_profiles':len(set(actual.values())),'maximal_profiles':maximal,
            'domination_map':domination,'moment':moment,'polynomial_rows':len(rows),
            'checked_newton_coefficients':ncoeff,'checked_finite_values':nfinite,
            'maximum_shift':shiftmax,'polynomials':rows,
            'seconds':round(time.monotonic()-start,3)}
    suffix='' if len(residues)==5 else '_s'+'_'.join(map(str,residues))
    out=ROOT/f'infinite_r{r}{suffix}.certificate.json'
    out.write_text(json.dumps(record,indent=2)+'\n')
    summary={k:v for k,v in record.items() if k not in ('polynomials','domination_map','maximal_profiles')}
    summary['certificate_sha256']=hashlib.sha256(out.read_bytes()).hexdigest()
    (ROOT/f'infinite_r{r}{suffix}.summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)
    return record


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('family',type=Path)
    ap.add_argument('--min-n',type=int)
    args=ap.parse_args()
    prove(args.family,args.min_n)
