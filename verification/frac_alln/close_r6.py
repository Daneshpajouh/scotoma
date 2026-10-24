#!/usr/bin/env python3
"""Can the dual close r = 6 at n = 20, m = 8? The coarse grid left it undecided.

At r=6 the grid of lambdas taken from the lambda=0 maximising centre gave upper 41/34 = 1.20588
against a requirement of 6/5 = 1.2, and a lower bound of 101/85 = 1.18824. So 6/5 lies strictly
inside the bracket and r=6 is neither proved nor refuted by that grid. The optimal lambda for a
top-r dual is the r-th largest weight AT THE MAXIMISING CENTRE FOR THAT r, which need not be the
lambda=0 centre, so a finer search can tighten the upper bound. This searches it.
"""
import os
import sys
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # sibling files, any working directory
from topr import dp_shifted, exhibit_weights, topr_lower

N, m, R = 19, 8, 6
need = Fraction(R, R - 1)
print('n=%d m=%d r=%d  need upper < %s = %.6f' % (N + 1, m, R, need, need), flush=True)

# A grid between the 7th and 5th largest weights of the lambda=0 centre, where the optimum sits,
# plus the endpoints themselves.
lo_l, hi_l = Fraction(1, 68), Fraction(7, 68)
cands = {Fraction(7, 170), lo_l, hi_l}
for k in range(1, 24):
    cands.add(lo_l + (hi_l - lo_l) * Fraction(k, 24))
best, bl, bw = None, None, None
for lam in sorted(cands):
    v, w = dp_shifted(N, m, lam)
    u = R * lam + v
    if best is None or u < best:
        best, bl, bw = u, lam, w
    print('  lambda %-12s  inner %-14s  upper %-14s %.7f %s'
          % (lam, v, u, u, '<-- best' if u == best else ''), flush=True)
print()
print('BEST UPPER  %s = %.7f   at lambda %s' % (best, best, bl), flush=True)
print('NEED        %s = %.7f' % (need, need), flush=True)
print('VERDICT     r=%d is %s' % (R, 'PROVED' if best < need else 'still not proved'), flush=True)
lowr = topr_lower(bw, N, m, R)
print('lower bound from that centre: %s = %.7f  (>= need means r=%d is refuted, not merely open)'
      % (lowr, lowr, R), flush=True)
