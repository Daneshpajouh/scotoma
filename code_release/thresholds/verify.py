#!/usr/bin/env python3
"""Standalone exact threshold-certificate verifier; standard library only.

A FULL certificate gives a mass-r private packing and checks every centre.
A SHORT certificate covers every private orbit in the dual of that packing
problem. This rules out mass r, by saturation of the r library constraints.
The latter is an obstruction certificate, not a cover of the whole union.
"""
import argparse
from fractions import Fraction as F
from itertools import product
from math import comb, lcm
from pathlib import Path
import json
import sys
import time
sys.dont_write_bytecode = True


def require(ok, why):
    if not ok:
        raise AssertionError(why)


class Geometry:
    def __init__(self, N, r, m):
        require(all(type(v) is int for v in (N,r,m)), 'noninteger parameter')
        require(1 <= r <= N+1 and 0 <= m <= N, 'inadmissible parameters')
        self.N,self.r,self.m=N,r,m
        self.L,self.C=r-1,N-r+1
        self.prefixes=tuple(product(range(3),repeat=self.L))
        self.library=tuple((0,)*j+(1,)*(self.L-j) for j in range(r))
        self.private=[]
        self.regions={}
        self.dmin={}
        for a in self.prefixes:
            # Consecutive library distances differ by +1,-1,0.
            ds=[sum(x != 1 for x in a)]
            for x in a:
                ds.append(ds[-1]+(1 if x==1 else -1 if x==0 else 0))
            d=min(ds)
            self.dmin[a]=d
            w=m-d
            if ds.count(d)==1 and 0 <= w <= self.C:
                self.private.append((a,w))
                self.regions[(a,w)]=ds.index(d)
        self._tails={}
        self._short_tails={}
        self._prob={}

    def tails(self,p):
        if p not in self._tails:
            # Count actual moving tail words, using per-coordinate polynomials.
            D={(0,0):1}
            for i in range(self.C):
                terms=((0,1,1),(1,0,1),(1,1,2)) if i<p else ((0,0,1),(1,1,3))
                E={}
                for (w,d),v in D.items():
                    for dw,dd,k in terms:
                        key=w+dw,d+dd
                        E[key]=E.get(key,0)+v*k
                D=E
            rows=[]
            for w in range(self.C+1):
                row=[D.get((w,d),0) for d in range(self.C+1)]
                require(sum(row)==comb(self.C,w)*3**w,'tail normalization')
                for d in range(1,self.C+1):
                    row[d]+=row[d-1]
                rows.append(row)
            self._tails[p]=rows
        return self._tails[p]

    def probability(self,w,p,forced,both):
        key=w,p,forced,both
        if key not in self._prob:
            if p<=2:
                # Independently count actual moving tail words. Only p centre
                # positions can overlap, so this count has at most six terms
                # even for an arbitrarily long tail.
                if (w,p) not in self._short_tails:
                    hist={}
                    for a in range(max(0,w-(self.C-p)),min(p,w)+1):
                        base=comb(p,a)*comb(self.C-p,w-a)*3**(w-a)
                        for e in range(a+1):
                            distance=p+w-a-e
                            hist[distance]=hist.get(distance,0)+base*comb(a,e)*2**(a-e)
                    require(sum(hist.values())==comb(self.C,w)*3**w,'short-tail normalization')
                    self._short_tails[w,p]=hist
                hist=self._short_tails[w,p]
                num=sum(comb(both,h)*sum(v for distance,v in hist.items()
                                         if distance<=self.m-forced-h)
                        for h in range(both+1))
                self._prob[key]=F(num,2**both*comb(self.C,w)*3**w)
                return self._prob[key]
            row=self.tails(p)[w]
            num=0
            for h in range(both+1):
                end=self.m-forced-h
                if end>=0:
                    num+=comb(both,h)*row[min(self.C,end)]
            self._prob[key]=F(num,2**both*comb(self.C,w)*3**w)
        return self._prob[key]

    def kernel(self,centre,point):
        a,p=centre
        b,w=point
        forced=sum(x!=y for x,y in zip(a,b))
        both=sum(x==y==2 for x,y in zip(a,b))
        return self.probability(w,p,forced,both)

    def scale_factor(self,w):
        return 2**self.L*comb(self.C,w)*3**w

    def parse(self,terms):
        out=[]
        for t in terms:
            a=tuple(2 if x in ('X','2',2) else int(x) for x in t['prefix'])
            p=t['tail_weight']
            require(len(a)==self.L and all(x in (0,1,2) for x in a),'bad prefix')
            require(type(p) is int and 0<=p<=self.C,'bad tail weight')
            require(not isinstance(t['mass'],float),'float mass forbidden')
            v=F(t['mass'])
            require(v>=0,'negative mass')
            if v:
                out.append(((a,p),v))
        require(len({q for q,v in out})==len(out),'duplicate support')
        return out

    def packing_check(self,support):
        regions=[F(0)]*self.r
        for q,v in support:
            require(q in self.regions,'support not private on the sphere')
            regions[self.regions[q]]+=v
        require(regions==[F(1)]*self.r,'private-region masses must be one')
        scale=lcm(*(v.denominator*self.scale_factor(q[1])
                    for q,v in support))
        width=self.L+1
        tables=[]
        for (b,w),v in support:
            layers=[]
            for p in range(self.C+1):
                tab={}
                for forced in range(width):
                    for both in range(width-forced):
                        x=v*self.probability(w,p,forced,both)
                        quo,rem=divmod(scale,x.denominator)
                        require(rem==0,'inexact common denominator')
                        tab[forced,both]=x.numerator*quo
                layers.append(tab)
            tables.append(layers)
        worst=F(0)
        witness=None
        for a in self.prefixes:
            sig=[(sum(x!=y for x,y in zip(a,b)),sum(x==y==2 for x,y in zip(a,b)))
                 for (b,w),v in support]
            for p in range(self.C+1):
                value=F(sum(t[p][s] for t,s in zip(tables,sig)),scale)
                if value>worst:
                    worst,witness=value,(a,p)
        require(worst<=1,'packing overload: %s at %s'%(worst,witness))
        return {'max_load':str(worst),'witness':encode_orbit(witness),
                'centre_classes':3**self.L*(self.C+1),'region_masses':list(map(str,regions)),
                'mass':str(sum(regions,F(0)))}

    def dual_check(self,support):
        require(self.private,'no private points; use the empty-region certificate')
        worst=None
        witness=None
        for q in self.private:
            value=sum((v*self.kernel(c,q) for c,v in support),F(0))
            if worst is None or value<worst:
                worst,witness=value,q
        mass=sum((v for c,v in support),F(0))
        require(worst>=1,'private dual undercoverage: %s at %s'%(worst,witness))
        require(mass<self.r,'private dual is not strictly short: %s'%mass)
        return {'private_orbits':len(self.private),'minimum':str(worst),
                'mass':str(mass),'deficit':str(self.r-mass),
                'witness':encode_orbit(witness)}


def encode_orbit(q):
    a,p=q
    return {'prefix':''.join('X' if x==2 else str(x) for x in a),'tail_weight':p}


def encode(support):
    return [dict(encode_orbit(q),mass=str(v)) for q,v in support if v]


def verify(cert):
    G=Geometry(cert['N'],cert['r'],cert['m'])
    result={k:cert[k] for k in ('N','r','m','kind')}
    if cert['kind']=='FULL':
        result.update(G.packing_check(G.parse(cert['support'])))
    elif cert['kind']=='SHORT':
        result.update(G.dual_check(G.parse(cert['support'])))
    elif cert['kind']=='EMPTY_REGION':
        missing=set(range(G.r)) - set(G.regions.values())
        require(missing,'claimed empty region has private points')
        result['empty_regions']=sorted(missing)
    else:
        raise AssertionError('unknown certificate kind')
    return dict(result,status='VERIFIED-COMPUTATIONALLY',passes=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths',nargs='+',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    results=[]
    for path in args.paths:
        start=time.monotonic()
        try:
            ans=verify(json.loads(path.read_text()))
        except (AssertionError,ValueError,KeyError,TypeError) as exc:
            print(json.dumps({'passes':False,'path':str(path),'reason':str(exc)}),flush=True)
            return 1
        ans['seconds']=round(time.monotonic()-start,3)
        ans['path']=str(path)
        results.append(ans)
        print(json.dumps(ans),flush=True)
    if args.output:
        args.output.write_text(json.dumps(results,indent=2)+'\n')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
