#!/usr/bin/env python3
"""Exact checks supporting the analytic all-run packing and period-five theorem."""
import sys
sys.dont_write_bytecode=True
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
import json
from verify import Geometry,require,encode,verify

ROOT=Path(__file__).resolve().parent


def choose(n,k):
    return comb(n,k) if 0<=k<=n else 0


@lru_cache(None)
def conditional_gain_numerator(p,a):
    # Denominator 3^p; A=a occupied positions, E~Binomial(a,1/3).
    return 3**(p-a)*sum(comb(a,e)*2**(a-e)*max(0,2*(a+e-p)+1) for e in range(a+1))


def finite_gain(C,m,p):
    if C==0:
        return F(1)
    numerator=sum(choose(m,a)*choose(C-m,p-a)*conditional_gain_numerator(p,a)
                  for a in range(p+1))
    return F(numerator,comb(C,p)*3**p)


def independent_gain(C,m,p):
    return F(sum(comb(p,a)*m**a*(C-m)**(p-a)*conditional_gain_numerator(p,a)
                 for a in range(p+1)),C**p*3**p)


def main():
    cap=100
    walk={0:1}
    iid=[]
    reflection=0
    recurrence=0
    for p in range(cap+1):
        total=5**p
        expectation=F(sum(v*max(0,2*s+1) for s,v in walk.items()),total)
        require(expectation<=1,'iid gain expectation exceeds one')
        for k in range(1,p+1):
            require(walk.get(-k,0)==2**k*walk.get(k,0),'biased-walk reflection failed')
            reflection+=1
        if p:
            require(expectation-iid[-1]==-F(2,5)*previous_positive_two,
                    'biased-walk drift identity failed')
            recurrence+=1
        iid.append(expectation)
        previous_positive_two=F(sum(v for s,v in walk.items() if s>=2),total)
        nxt={}
        for s,v in walk.items():
            for step,k in ((-1,2),(0,2),(1,1)):
                nxt[s+step]=nxt.get(s+step,0)+v*k
        walk=nxt
    gain_cases=0
    convex_order_cases=0
    for C in range(cap+1):
        m=3*C//5
        for p in range(C+1):
            value=finite_gain(C,m,p)
            require(value<=1,'finite without-replacement gain bound failed')
            gain_cases+=1
            if C:
                independent=independent_gain(C,m,p)
                require(value<=independent<=iid[p],
                        'convex-order or density-monotonicity comparison failed')
                convex_order_cases+=1
    chain_cases=0
    for r in range(1,11):
        L=r-1
        for a in product(range(3),repeat=L):
            ds=[sum(x!=1 for x in a)]
            for x in a:
                ds.append(ds[-1]+(1 if x==1 else -1 if x==0 else 0))
            for k in range(L+1):
                require(sum(d<=k for d in ds)<=2*k+1,'chain count bound failed')
                chain_cases+=1
    cutoff_factor=F(12,5)*F(197,200)**60
    require(cutoff_factor<1,'all-run 60r tail cutoff failed')
    # An above-boundary control retains privacy and all regional masses.
    G=Geometry(7,3,4)
    support=[((a,4),F(1)) for a in G.library]
    witness=sum((v*G.kernel(((0,1),1),q) for q,v in support),F(0))
    require(witness==F(4,3) and witness>1,'general packing control did not overload')
    bad={'N':7,'r':3,'m':4,'kind':'FULL','support':encode(support)}
    (ROOT/'control_general_packing_above_density.invalid.json').write_text(json.dumps(bad,indent=2)+'\n')
    try:
        verify(bad)
    except AssertionError as exc:
        require('packing overload' in str(exc),'wrong general-packing rejection')
        control={'name':'uniform_library_orbits_above_three_fifths_tail_density',
                 'status':'REFUTED','outcome':'REJECTED-AS-REQUIRED',
                 'explicit_witness':'4/3','reason':str(exc)}
    else:
        raise AssertionError('general packing control was accepted')
    record={'status':'VERIFIED-COMPUTATIONALLY','passes':True,'tail_length_cap':cap,
            'finite_gain_cases':gain_cases,'convex_order_comparisons':convex_order_cases,
            'reflection_equalities':reflection,'drift_recurrence_equalities':recurrence,
            'chain_cardinality_cases':chain_cases,'cutoff_factor':str(cutoff_factor),
            'cutoff_factor_strictly_below_one':True,'controls':[control],
            'scope':'finite checks support the separate analytic all-parameter proof; no extrapolation used'}
    (ROOT/'general_verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    main()
