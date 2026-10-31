#!/usr/bin/env python3
"""Tests for m0_search.py. Run: python3 test_m0_search.py  (pytest also works).

Brute force is a plain linear scan with exact rational arithmetic, independent of the bracketing
and bisection in m0_search. Two tests run boundary_exact_test.py and occupancy_curve_test.py as
separate processes and parse their printed tables, so the integer search is checked against what
those scripts actually print."""
import os, re, subprocess, sys
from fractions import Fraction as F
from math import ceil, prod

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from m0_search import m0_iid, m0_lecam, m0_hypergeom, m0_float, smallest_true, as_fraction  # noqa: E402

BOUNDARY = os.path.join(HERE, 'boundary_exact_test.py')
OCCUPANCY = os.path.join(HERE, 'occupancy_curve_test.py')
GAMMAS = [F(1, 2), F(1, 3), F(1, 5), F(1, 10), F(1, 20), F(1, 50), F(3, 7), F(1, 100), F(1, 1000),
          F(2, 3), F(99, 100)]
DELTAS = [F(1, 2), F(1, 10), F(1, 20), F(1, 40), F(1, 100), F(1, 1000), F(7, 9), F(1, 10 ** 6)]
BOUNDARY_CASES = [(g, k) for g in (F(1, 2), F(1, 5), F(1, 10), F(1, 20), F(1, 50), F(3, 7))
                  for k in (1, 2, 3, 5, 8, 13)]


def brute_iid(g, d):
    m = 0
    while (1 - g) ** m > d:
        m += 1
    return m


def brute_hypergeom(N, k, d):
    for m in range(0, N + 1):
        p = F(0) if m > N - k else F(prod(range(N - k - m + 1, N - k + 1)), prod(range(N - m + 1, N + 1)))
        if p <= d:
            return m
    raise AssertionError('unreachable')


def test_iid_equals_brute_force():
    n = 0
    for g in GAMMAS:
        for d in DELTAS:
            b = brute_iid(g, d)
            assert m0_iid(g, d) == b, (g, d)
            assert m0_iid(g, d, exact_only=True) == b, (g, d)
            n += 2
    return n


def test_fast_path_equals_exact_path():
    """The decimal fast path must never change an answer: random rationals plus exact boundaries."""
    import random
    rng = random.Random(20260930)
    n = 0
    cases = [(g, (1 - g) ** k) for g, k in BOUNDARY_CASES]
    cases += [(F(rng.randint(1, 999), rng.randint(1000, 200000)), F(rng.randint(1, 99), 100)) for _ in range(300)]
    for g, d in cases:
        assert m0_iid(g, d) == m0_iid(g, d, exact_only=True), (g, d)
        n += 1
    return n


def test_large_cases():
    """Sizes where exact powers alone would be slow; checked against independent high-precision logs
    (i.i.d.) and against exact binomial coefficients (hypergeometric)."""
    from decimal import Decimal, localcontext
    from math import comb
    n = 0
    for N in (10 ** 6, 10 ** 7):
        for k in (1, 10):
            g, d = F(k, N), F(1, 40)
            m = m0_iid(g, d)
            with localcontext() as ctx:
                ctx.prec = 120
                r = (Decimal(1) / 40).ln() / (Decimal(N - k) / Decimal(N)).ln()
            assert m == int(r.to_integral_value(rounding='ROUND_CEILING')), (N, k, m, r)
            n += 1
    for N, k in ((200000, 10), (200000, 1), (50000, 50)):
        m = m0_hypergeom(N, k, F(1, 40))
        assert comb(N - k, m) * 40 <= comb(N, m)
        assert comb(N - k, m - 1) * 40 > comb(N, m - 1)
        n += 1
    return n


def test_exact_boundaries():
    """delta = (1-gamma)^k exactly: the answer is k, and the float formula overshoots in 9 of 36."""
    over = []
    for g, k in BOUNDARY_CASES:
        d = (1 - g) ** k
        assert m0_iid(g, d) == k, (g, k)
        fl = m0_float(g, d)
        assert fl in (k, k + 1), (g, k, fl)
        if fl != k:
            over.append((str(g), k))
    assert len(BOUNDARY_CASES) == 36
    assert over == [('1/5', 2), ('1/5', 3), ('1/5', 5), ('1/5', 8), ('1/5', 13),
                    ('1/10', 3), ('1/10', 5), ('1/10', 8), ('1/10', 13)], over
    return f'{len(BOUNDARY_CASES)} boundary cases, float formula overshoots in {len(over)}'


def test_lecam_is_iid_at_twice_delta():
    n = 0
    for g in GAMMAS:
        for d in DELTAS:
            if 2 * d < 1:
                assert m0_lecam(g, d) == brute_iid(g, 2 * d)
                assert m0_lecam(g, d) <= m0_iid(g, d)
                n += 1
    return n


def test_hypergeom_equals_brute_force():
    n = 0
    for N in range(1, 41):
        for k in range(1, N + 1):
            for d in (F(1, 2), F(1, 10), F(1, 40), F(1, 100), F(3, 4)):
                assert m0_hypergeom(N, k, d) == brute_hypergeom(N, k, d), (N, k, d)
                n += 1
    return n


def test_hypergeom_never_exceeds_iid_or_census():
    n = 0
    for N in (10, 37, 100, 1000, 5000):
        for k in (1, 2, 3, 5, 10):
            if k > N:
                continue
            for d in (F(1, 20), F(1, 40), F(1, 100)):
                h = m0_hypergeom(N, k, d)
                assert h <= N - k + 1
                assert h <= m0_iid(F(k, N), d)
                n += 1
    return n


def test_monotone():
    for d in DELTAS:
        vals = [m0_iid(g, d) for g in sorted(GAMMAS)]
        assert vals == sorted(vals, reverse=True), d
    for g in GAMMAS:
        vals = [m0_iid(g, d) for d in sorted(DELTAS)]
        assert vals == sorted(vals, reverse=True), g
    return 2


def test_smallest_true_any_guess():
    n = 0
    for target in (0, 1, 2, 7, 100, 12345):
        for guess in (0, 1, target, target + 50, 3 * target + 999, max(0, target - 40)):
            assert smallest_true(lambda m: m >= target, guess) == target
            n += 1
    return n


def test_rejects_float_input():
    try:
        m0_iid(0.2, F(1, 20))
    except TypeError:
        return 1
    raise AssertionError('float input accepted')


def test_against_boundary_exact_test_output():
    out = subprocess.run([sys.executable, BOUNDARY], capture_output=True, text=True, check=True).stdout
    rows = [l.split() for l in out.splitlines() if re.match(r'^\s*\d+/\d+\s+\d+\s', l)]
    assert len(rows) == 36, len(rows)
    for r in rows:
        g, k, exact, flt = F(r[0]), int(r[1]), int(r[3]), int(r[5])
        assert m0_iid(g, (1 - g) ** k) == exact == k
        assert m0_float(g, (1 - g) ** k) == flt
    return len(rows)


def test_against_occupancy_curve_test_output():
    out = subprocess.run([sys.executable, OCCUPANCY], capture_output=True, text=True, check=True).stdout
    assert 'delta = 0.025' in out
    rows = [list(map(int, l.split())) for l in out.splitlines() if re.match(r'^\s*\d+\s+\d+\s+\d+\s+\d+\s+\d+\s', l)]
    assert len(rows) == 21, len(rows)
    for N, k, exact, a, b, ea, eb in rows:
        assert m0_hypergeom(N, k, F(1, 40)) == exact, (N, k)
    return len(rows)


TESTS = [test_iid_equals_brute_force, test_fast_path_equals_exact_path, test_large_cases, test_exact_boundaries, test_lecam_is_iid_at_twice_delta,
         test_hypergeom_equals_brute_force, test_hypergeom_never_exceeds_iid_or_census, test_monotone,
         test_smallest_true_any_guess, test_rejects_float_input, test_against_boundary_exact_test_output,
         test_against_occupancy_curve_test_output]

if __name__ == '__main__':
    failed = 0
    for t in TESTS:
        try:
            n = t()
            print(f'PASS  {t.__name__}  ({n}{" checks" if isinstance(n, int) else ""})')
        except Exception as e:  # report every failure, then exit non-zero
            failed += 1
            print(f'FAIL  {t.__name__}: {type(e).__name__}: {e}')
    print(f'\n{len(TESTS) - failed} of {len(TESTS)} tests passed')
    sys.exit(1 if failed else 0)
