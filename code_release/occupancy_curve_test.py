#!/usr/bin/env python3
"""Compare two closed forms with the exact finite-population detection sample size.

Exact target: the smallest m with P(a sample of m from N misses all k bad units) <= delta,
that is C(N-k, m)/C(N, m) <= delta.

Candidate A (plain):     m = ceil(N * (1 - delta^(1/k)))
Candidate B (corrected): m = ceil((1 - delta^(1/k)) * (N - (k-1)/2)), the standard
                         finite-population detection formula with a continuity correction.

B tracks the exact threshold at least as well as A everywhere on the grid below, so A is the
leading-order form of B with the correction dropped. The paper reports the exact integer."""
from math import comb, ceil

def exact_m(N, k, delta):
    lo, hi = 0, N
    while lo < hi:
        mid = (lo + hi) // 2
        p = comb(N - k, mid) / comb(N, mid) if mid <= N - k else 0.0
        if p <= delta:
            hi = mid
        else:
            lo = mid + 1
    return lo

D = 0.025
print(f'delta = {D}')
print(f"{'N':>8} {'k':>5} {'exact':>9} {'A plain':>9} {'B corrected':>12} {'errA':>6} {'errB':>6}")
sa = sb = n = 0
for N in (1000, 5000, 20000):
    for k in (1, 2, 5, 10, 25, 50, 100):
        e = exact_m(N, k, D)
        f = 1 - D ** (1 / k)
        a = ceil(N * f)
        b = ceil(f * (N - (k - 1) / 2))
        sa += abs(a - e); sb += abs(b - e); n += 1
        print(f'{N:>8} {k:>5} {e:>9} {a:>9} {b:>12} {a-e:>6} {b-e:>6}')
print(f'\nmean |error| vs exact:  plain {sa/n:.3f}   corrected {sb/n:.3f}')
print('VERDICT: the corrected (textbook) form is the better approximation to the same exact quantity,'
      '\nso the proposed curve is its leading-order special case, not a new result.' if sb <= sa else
      'VERDICT: the proposed form is closer; re-examine.')
