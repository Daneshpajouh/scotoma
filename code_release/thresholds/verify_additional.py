#!/usr/bin/env python3
"""Independent exact checks of the symbolic upper bound and ansatz obstruction.

Standard library only. All polynomials and coefficients are regenerated from
finite, coordinate-wise tail counts, rather than trusting symbolic expressions.
"""
import sys
sys.dont_write_bytecode=True
from fractions import Fraction as F
from math import comb
from functools import lru_cache
from pathlib import Path
import json
from verify import Geometry,require,verify

ROOT=Path(__file__).resolve().parent


def trim(p):
    p=list(map(F,p))
    while len(p)>1 and not p[-1]:
        p.pop()
    return p


def add(a,b):
    return trim([(a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0))
                 for i in range(max(len(a),len(b)))])


def scale(a,t):
    return trim([t*x for x in a])


def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j]+=x*y
    return trim(c)


def shift(a,s):
    return trim([sum((a[j]*comb(j,i)*s**(j-i) for j in range(i,len(a))),F(0))
                 for i in range(len(a))])


def evaluate(a,s):
    out=F(0)
    for x in reversed(a):
        out=out*s+x
    return out


def positive_coefficients(poly,s,name,strict=False):
    shifted=shift(poly,s)
    require(all(x>=0 for x in shifted),name+' coefficient negative')
    if strict:
        require(shifted[0]>0,name+' not strictly positive at the initial integer')
    return len(shifted)


class LimitGeometry(Geometry):
    @staticmethod
    @lru_cache(None)
    def gain_counts(p):
        D={0:1}
        for _ in range(p):
            E={}
            for gain,v in D.items():
                for move,k in ((-1,2),(0,2),(1,1)):
                    E[gain+move]=E.get(gain+move,0)+v*k
            D=E
        require(sum(D.values())==5**p,'limit gain normalization')
        return D

    def probability(self,w,p,forced,both):
        d=self.m-w
        counts=self.gain_counts(p)
        num=sum(comb(both,h)*sum(v for gain,v in counts.items() if gain>=forced+h-d)
                for h in range(both+1))
        return F(num,2**both*5**p)


def parse_prefix(s):
    return tuple(2 if x=='X' else int(x) for x in s)


def constant_obstruction(cert):
    r=cert['r']
    gs=[]
    for item in cert['instances']:
        cls=LimitGeometry if item['kind']=='LIMIT_RHO_3_5' else Geometry
        gs.append(cls(item['N'],r,item['m']))
    prefixes=[a for a,w in gs[0].private]
    require(all([a for a,w in G.private]==prefixes for G in gs),'constant-family prefix sets differ')
    terms=[]
    for term in cert['support']:
        i=term['instance']
        c=parse_prefix(term['prefix']),term['tail_weight']
        v=F(term['mass'])
        require(0<=i<len(gs) and v>=0,'invalid obstruction multiplier')
        terms.append((i,c,v))
    least=min(sum((v*gs[i].kernel(c,(a,gs[i].m-gs[i].dmin[a])) for i,c,v in terms),F(0))
              for a in prefixes)
    mass=sum((v for i,c,v in terms),F(0))
    require(least>=1,'constant-family dual coefficient below one')
    require(mass<r,'constant-family dual is not strictly short')
    return {'status':'PROVED','private_prefixes':len(prefixes),'minimum':str(least),
            'dual_mass':str(mass),'deficit':str(r-mass)}


def upper(cert,check_bridges=True):
    r,n0,m0=cert['r'],cert['n0'],cert['m0']
    z0=cert['positive_from_z']
    require(type(z0) is int and z0>=0,'missing nonnegative positivity shift')
    Q=list(map(F,cert['denominator_coefficients']))
    W=[list(map(F,row)) for row in cert['weight_numerator_coefficients']]
    groups=[[(parse_prefix(t['prefix']),t['tail_weight']) for t in group] for group in cert['groups']]
    require(len(groups)==len(W),'upper support count mismatch')
    require(all(g and len({p for a,p in g})==1 and 0<=g[0][1]<=2 for g in groups),
            'invalid upper centre group')
    positive_coefficients(Q,z0,'upper denominator',True)
    mass=[F(0)]
    ncoeff=0
    for i,(g,w) in enumerate(zip(groups,W)):
        ncoeff+=positive_coefficients(w,z0,'upper weight '+str(i))
        mass=add(mass,scale(w,len(g)))
    gap=add(scale(Q,r),scale(mass,-1))
    ncoeff+=positive_coefficients(gap,z0,'upper strict gap',True)
    C0=n0-r
    C=[F(C0),F(5)]
    Cm1=[F(C0-1),F(5)]
    common=scale(mul(C,Cm1),9)
    quotients=[common,scale(Cm1,3),[F(2)]]
    gs=[Geometry(n0+5*j-1,r,m0+3*j) for j in range(3)]
    prefixes=[a for a,w in gs[0].private]
    require(all([a for a,w in G.private]==prefixes for G in gs),'upper private domain changes')
    alphabet_scale=2**(r-1)
    @lru_cache(None)
    def numerator(a,c):
        b,p=c
        vals=[]
        for j in range(p+1):
            G=gs[j]
            point=a,G.m-G.dmin[a]
            vals.append(G.kernel(c,point)*alphabet_scale*comb(G.C,p)*3**p)
        if p==0:
            return [vals[0]]
        if p==1:
            return [vals[0],vals[1]-vals[0]]
        second=(vals[2]-2*vals[1]+vals[0])/2
        return [vals[0],vals[1]-vals[0]-second,second]
    for a in prefixes:
        total=[F(0)]
        for g,w in zip(groups,W):
            for c in g:
                total=add(total,mul(mul(w,numerator(a,c)),quotients[c[1]]))
        slack=add(total,scale(mul(Q,common),-alphabet_scale))
        ncoeff+=positive_coefficients(slack,z0,'upper coverage '+str(a))
    bridges=cert.get('finite_bridges',[])
    if check_bridges:
        require(len(bridges)==z0,'upper finite bridge missing')
        for j,path in enumerate(bridges):
            bridge=json.loads((ROOT/path).read_text())
            require((bridge['N'],bridge['r'],bridge['m'],bridge['kind'])==
                    (n0+5*j-1,r,m0+3*j,'SHORT'),'wrong upper bridge parameters')
            verify(bridge)
    return {'status':'PROVED','r':r,'line_n0':n0,'line_m0':m0,'shift':z0,
            'private_prefixes':len(prefixes),'positive_coefficients_checked':ncoeff,
            'finite_bridges':len(bridges),'all_z_ge_zero':check_bridges}


def main():
    record={}
    record['constant_mass_obstruction']=constant_obstruction(json.loads((ROOT/'constant_mass_obstruction_r7.json').read_text()))
    record['upper_r7']=upper(json.loads((ROOT/'upper_r7.json').read_text()))
    extra=ROOT/'upper_r7_n46.json'
    if extra.exists():
        record['upper_r7_n46']=upper(json.loads(extra.read_text()))
    extra=ROOT/'upper_r8_n83.json'
    if extra.exists():
        candidate=json.loads(extra.read_text())
        if candidate.get('positive_from_z') is not None:
            record['upper_r8_n83']=upper(candidate)
    (ROOT/'additional_verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    main()
