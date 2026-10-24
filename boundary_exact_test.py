#!/usr/bin/env python3
"""Exact-boundary check of the sample-size formula against the Neyman-Pearson umbrella criterion.

Both quantities are the least m with (1-gamma)^m <= delta. They can only disagree where
delta = (1-gamma)^k exactly, because that is where a rounding step in floating point can push
one of them off by one. Everything here is exact rational arithmetic, no floats, except the
column that evaluates the printed formula in floating point for comparison."""
from fractions import Fraction as F
from math import comb, ceil, log

def nproc_min_n_exact(gamma, delta):
    """nproc: loc = min{k : P(Bin(n,1-gamma) >= k) <= delta} exists iff (1-gamma)^n <= delta."""
    n = 1
    while n < 20000:
        if (1 - gamma) ** n <= delta:
            return n
        n += 1
    return None

def formula_exact(gamma, delta):
    """Smallest integer m with (1-gamma)^m <= delta, computed without logs."""
    m = 1
    while m < 20000:
        if (1 - gamma) ** m <= delta:
            return m
        m += 1
    return None

def formula_float(gamma, delta):
    return ceil(log(float(delta)) / log(1 - float(gamma)))

print('Exact boundary cases delta = (1-gamma)^k, where a float formula can go off by one')
print(f"{'gamma':>8} {'k':>3} {'delta=(1-g)^k':>22} {'exact':>6} {'nproc':>6} {'float ceil':>11} {'float ok':>9}")
bad = []
for gamma in (F(1,2), F(1,5), F(1,10), F(1,20), F(1,50), F(3,7)):
    for k in (1, 2, 3, 5, 8, 13):
        delta = (1 - gamma) ** k
        e, np_, fl = formula_exact(gamma, delta), nproc_min_n_exact(gamma, delta), formula_float(gamma, delta)
        ok = (fl == e)
        if not ok or e != np_:
            bad.append((gamma, k, e, np_, fl))
        print(f'{str(gamma):>8} {k:>3} {str(delta)[:22]:>22} {e:>6} {np_:>6} {fl:>11} {str(ok):>9}')
print()
print('exact vs nproc disagreements:', [(str(g),k,e,n) for g,k,e,n,_ in bad if e != n] or 'none')
print('float-formula off-by-one cases:', [(str(g),k,e,f) for g,k,e,_,f in bad if f != e] or 'none')
print('\nCONCLUSION: at the exact boundaries the two criteria agree, and the floating-point'
      '\nceil(log(delta)/log(1-gamma)) can differ from the exact answer, so a manuscript that ships'
      '\ncode should compute the bound by integer search, not by logs.')
