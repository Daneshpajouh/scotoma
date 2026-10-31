"""Exact support-overlap packing loads, with no discovery dependencies."""
from fractions import Fraction as F
from math import comb, lcm
from exact_tools import prefixes, suffix_counts, library

def parse_terms(terms):
    merged={}
    for t in terms:
        p=tuple(2 if x=='X' else int(x) for x in t['prefix'])
        key=(p,t['suffix_weight'])
        merged[key]=merged.get(key,F(0))+F(t['mass'])
    return [(p,w,mass) for (p,w),mass in sorted(merged.items()) if mass]

def exact_loads(N,r,m,terms):
    L,C=r-1,N-r+1
    lib=[p[:L] for p in library(N,r)]
    for p,w,mass in terms:
        assert mass>=0 and 0<=w<=C and min(sum(x!=y for x,y in zip(p,c)) for c in lib)+w<=m
    scale=lcm(*(mass.denominator*2**L*comb(C,w)*3**w for p,w,mass in terms))
    tables=[]
    for pref,w,mass in terms:
        layers=[]
        for a in range(C+1):
            counts=suffix_counts(C,a)[w];cum=[0]
            for value in counts:cum.append(cum[-1]+value)
            table={}
            for forced in range(L+1):
                for both in range(L+1-forced):
                    num=sum(comb(both,d)*cum[max(0,min(C,m-forced-d)+1)] for d in range(both+1))
                    denom=2**both*comb(C,w)*3**w*mass.denominator
                    assert scale%denom==0
                    table[forced,both]=num*mass.numerator*(scale//denom)
            layers.append(table)
        tables.append(layers)
    loads=[]
    for pref in prefixes(L):
        signatures=[(sum(x!=y for x,y in zip(pref,p)),sum(x==y==2 for x,y in zip(pref,p)))
                    for p,w,mass in terms]
        for a in range(C+1):
            value=sum(table[a][sig] for table,sig in zip(tables,signatures))
            loads.append(value)
    return scale,loads
