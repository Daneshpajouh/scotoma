#!/usr/bin/env python3
"""Run the run-count-aware bound for the rows the appendix table prints."""
import sys
from fractions import Fraction
sys.path.insert(0, '.')
from topr import dp_shifted, exhibit_weights, topr_lower

CASES = [(19, 8), (20, 8), (22, 8), (19, 7)]

for N, m in CASES:
    n = N + 1
    print('==== n=%d (N=%d), m=%d ====' % (n, N, m), flush=True)
    base, word0 = dp_shifted(N, m, Fraction(0))
    ws = sorted(exhibit_weights(word0, N, m), reverse=True)
    print('  lambda=0 bound (= the r=N+1 maximum): %s = %.6f' % (base, base), flush=True)
    print('  weights of that centre, largest first: %s' % ', '.join(str(w) for w in ws[:8]), flush=True)

    # candidate lambdas: the distinct weights of the exhibited centre, plus 0
    cands = sorted({Fraction(0)} | set(ws))
    cache = {}
    words = {}
    for lam in cands:
        cache[lam], words[lam] = dp_shifted(N, m, lam)
    print('  %3s %12s %14s %14s %14s %s' % ('r', 'need < ', 'best upper', 'at lambda', 'lower bnd', 'verdict'), flush=True)
    largest_proved = 0
    for r in range(1, n + 1):
        need = Fraction(r, r - 1) if r > 1 else None
        best, bl = None, None
        for lam in cands:
            v = r * lam + cache[lam]
            if best is None or v < best:
                best, bl = v, lam
        lo = max(topr_lower(words[lam], N, m, r) for lam in cands)
        if need is None:
            verdict = 'proved' if best <= 1 or True else ''
            ok = True
        else:
            ok = best < need
            verdict = 'PROVED' if ok else ('not proved' + (' (lower bound %s >= need)' % lo if lo >= need else ''))
        if ok:
            largest_proved = r
        print('  %3d %12s %14s %14s %14s %s'
              % (r, str(need) if need else '-', str(best)[:14], str(bl)[:14], str(lo)[:14], verdict), flush=True)
    print('  LARGEST r PROVED: %d of %d' % (largest_proved, n), flush=True)
