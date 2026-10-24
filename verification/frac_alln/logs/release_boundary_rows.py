#!/usr/bin/env python3
"""Run step3_dp.scaled_max at the full library r = N+1 at the lengths given, for the finite part of
the weighting-at-every-length check (finite_all_143_rows.log was made with it).
The full-library maximum M bounds every run count and every threshold set, so the floor holds at
every r <= N+1 if (N+1)(M-1) < M when M > 1 (and always when M <= 1); the condition is
sufficient, not necessary.
Usage: python3 release_boundary_rows.py m N [N ...]"""
import sys, time
from fractions import Fraction
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from step3_dp import scaled_max  # released producer, copied unchanged

m = int(sys.argv[1])
for N in map(int, sys.argv[2:]):
    t0 = time.time()
    best, L, word = scaled_max(N, m, N + 1)
    M = Fraction(best, L)
    r = N + 1
    holds_all_r = (M <= 1) or (r * (M - 1) < M)
    print(f"m={m} N={N} n={N+1} r={r} M={M} ({float(M):.6f}) word={''.join('01X'[c] for c in word)} "
          f"every_r_proved={holds_all_r} secs={time.time()-t0:.0f}", flush=True)
