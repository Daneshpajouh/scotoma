#!/usr/bin/env python3
"""Exact integer search for the sample sizes used in the paper.

Why this exists. The smallest m0 with (1-gamma)^m0 <= delta1 equals ceil(log(delta1)/log(1-gamma)) as
a real-number statement, but evaluating that ceiling in floating point overshoots by one at exact
boundaries delta1 = (1-gamma)^k (9 of 36 such cases, see boundary_exact_test.py).
Code should therefore compute m0 by integer search with exact comparisons. Every function here
takes rationals (int, fractions.Fraction, or a decimal string such as '0.05', converted exactly).
Floating point is used only to guess a starting bracket, which is verified and widened until valid,
so a bad guess costs time, never correctness.

  m0_iid(gamma, delta)        smallest m >= 0 with (1-gamma)^m <= delta        (Theorem 2(i), i.i.d. negatives)
  m0_lecam(gamma, delta)      smallest m >= 0 with (1-gamma)^m <= 2*delta      (two-point lower bound, Appendix B)
  m0_hypergeom(N, k, delta)   smallest m >= 0 with C(N-k, m)/C(N, m) <= delta  (sampling m of N units
                              without replacement misses all k bad units with probability <= delta)
  m0_float(gamma, delta)      the printed formula evaluated in floating point, for comparison only

How each comparison stays exact and fast.
  * Hypergeometric: C(N-k, m)/C(N, m) = prod_{j<k} (N-m-j)/(N-j) for m <= N-k, so the test is an
    integer comparison of two k-factor products, whatever N is.
  * i.i.d.: (1-gamma)^m <= delta is decided first by m*ln(1-gamma) - ln(delta) in decimal arithmetic at
    50 + len(str(m)) digits. Python's decimal ln, division and product are correctly rounded, so the
    computed value is within a few units of 10^-(P-1) * (m(1 + |ln q|) + |ln delta| + 1) of the true
    one. The sign is accepted only when the value clears that bound by a factor of 10^4; otherwise
    (an exact tie, or a near-tie) the integer comparison q.num^m * d.den <= q.den^m * d.num decides.
    exact_only=True forces the integer comparison everywhere (used by the tests).
"""
from decimal import Decimal, localcontext
from fractions import Fraction
from math import ceil, log, log1p, prod

__all__ = ['as_fraction', 'm0_iid', 'm0_lecam', 'm0_hypergeom', 'm0_float', 'smallest_true']


def as_fraction(x):
    """Exact rational from int, Fraction or str. Floats are refused: 0.05 is not 1/20 in binary."""
    if isinstance(x, float):
        raise TypeError('pass rationals as Fraction, int or a decimal string, not float')
    return Fraction(x)


def smallest_true(pred, guess=0):
    """Smallest integer m >= 0 with pred(m) True, for pred monotone (False...False True...True).
    The bracket starts at the guess and is widened geometrically until pred(lo - 1) is False (or lo
    is 0) and pred(hi) is True, then bisected. Only pred decides; the guess only saves evaluations."""
    lo, hi = max(0, guess - 2), max(0, guess + 2)
    step = 1
    while lo > 0 and pred(lo - 1):          # answer lies below lo
        hi, lo = lo - 1, max(0, lo - 1 - step)
        step *= 2
    step = 1
    while not pred(hi):                     # answer lies above hi
        lo, hi = hi + 1, hi + 1 + step
        step *= 2
    while lo < hi:                          # invariant: answer in [lo, hi], pred(hi) True
        mid = (lo + hi) // 2
        if pred(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def _power_le_exact(q, m, d):
    """q**m <= d for rationals 0 < q < 1, d > 0, by integer arithmetic."""
    return q.numerator ** m * d.denominator <= q.denominator ** m * d.numerator


def _power_le(q, m, d, exact_only=False):
    if exact_only or m == 0:
        return _power_le_exact(q, m, d)
    P = 50 + len(str(m))
    with localcontext() as ctx:
        ctx.prec = P
        lq = (Decimal(q.numerator) / Decimal(q.denominator)).ln()
        ld = (Decimal(d.numerator) / Decimal(d.denominator)).ln()
        x = Decimal(m) * lq - ld
        bound = Decimal(10) ** (5 - P) * (Decimal(m) * (1 + abs(lq)) + abs(ld) + 2)
        if x > bound:
            return False
        if x < -bound:
            return True
    return _power_le_exact(q, m, d)


def m0_iid(gamma, delta, exact_only=False):
    g, d = as_fraction(gamma), as_fraction(delta)
    if not (0 < g <= 1):
        raise ValueError('gamma must be in (0, 1]')
    if d <= 0:
        raise ValueError('delta must be positive')
    if d >= 1:
        return 0
    if g == 1:
        return 1
    q = 1 - g
    try:
        guess = ceil(log(d.numerator / d.denominator) / log1p(-(g.numerator / g.denominator)))
    except (ValueError, OverflowError, ZeroDivisionError):
        guess = 1
    return smallest_true(lambda m: _power_le(q, m, d, exact_only), max(guess, 0))


def m0_lecam(gamma, delta, exact_only=False):
    return m0_iid(gamma, 2 * as_fraction(delta), exact_only)


def m0_hypergeom(N, k, delta):
    if not (isinstance(N, int) and isinstance(k, int) and 1 <= k <= N):
        raise ValueError('need integers 1 <= k <= N')
    d = as_fraction(delta)
    if d <= 0:
        raise ValueError('delta must be positive')
    if d >= 1:
        return 0
    top = prod(N - j for j in range(k))

    def miss_all_le_delta(m):
        if m > N - k:                       # every sample of m units contains a bad unit
            return True
        return prod(N - m - j for j in range(k)) * d.denominator <= top * d.numerator

    f = 1 - (d.numerator / d.denominator) ** (1 / k)
    return smallest_true(miss_all_le_delta, min(N, max(0, ceil(f * N))))


def m0_float(gamma, delta):
    """The printed closed form, evaluated in floating point (natural logs), for comparison only."""
    return ceil(log(float(delta)) / log(1 - float(gamma)))


if __name__ == '__main__':
    for g, d in (('1/5', '64/125'), ('1/10', '729/1000'), ('1/28175', '0.05')):
        print(f'gamma={g:>8} delta={d:>8}  integer search {m0_iid(g, d):>7}  float formula {m0_float(Fraction(g), Fraction(d)):>7}')
