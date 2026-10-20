"""Budget tests, checked against brute force and against the closed-form edge cases.

The sample size is the one number here a reader might quote at a confidence level, so it is
checked by an independent search rather than by re-running the same binary search.
"""
from fractions import Fraction
from math import comb

import pytest

from scotoma.budget import (as_fraction, budget_rows, census_sensitivity, m_star, miss_probability)


def _miss_by_binomials(n, k, m):
    """miss probability from the definition, C(n-k, m) / C(n, m), as an exact Fraction."""
    if m > n - k:
        return Fraction(0)
    return Fraction(comb(n - k, m), comb(n, m))


def _m_star_by_scan(n, k, delta):
    for m in range(0, n - k + 2):
        if miss_probability(n, k, m) <= delta:
            return m
    raise AssertionError("no m satisfied the bound, which cannot happen")


def test_miss_probability_matches_the_binomial_definition():
    for n in (10, 37, 120):
        for k in (1, 2, 5):
            for m in (0, 1, 3, n // 2, n - k):
                assert miss_probability(n, k, m) == _miss_by_binomials(n, k, m)


def test_miss_probability_edge_cases():
    # Sampling nothing misses everything.
    assert miss_probability(100, 3, 0) == 1
    # Sampling all but k-1 cannot miss k of them.
    assert miss_probability(100, 3, 98) == 0
    assert miss_probability(100, 3, 100) == 0
    # With k = 1 the miss probability is exactly the unsampled fraction.
    assert miss_probability(100, 1, 40) == Fraction(60, 100)


def test_miss_probability_is_monotone_decreasing_in_m():
    previous = Fraction(2)
    for m in range(0, 60):
        current = miss_probability(500, 4, m)
        assert current <= previous
        previous = current


def test_m_star_matches_an_independent_scan():
    for n in (50, 200, 999):
        for k in (1, 2, 7):
            for delta in ("0.05", "0.01", "1/3"):
                d = as_fraction(delta)
                assert m_star(n, k, d) == _m_star_by_scan(n, k, d)


def test_m_star_achieves_the_bound_and_is_minimal():
    n, k, delta = 28175, 10, as_fraction("0.05")
    m = m_star(n, k, delta)
    assert miss_probability(n, k, m) <= delta
    assert miss_probability(n, k, m - 1) > delta


def test_k_equals_one_is_the_closed_form():
    """For k = 1 the requirement (N - m) / N <= delta gives m = ceil(N (1 - delta))."""
    for n in (100, 28175, 10 ** 6):
        for delta in ("0.05", "0.1"):
            d = as_fraction(delta)
            expected = -((-(n * (1 - d))) // 1)  # ceiling of an exact Fraction
            assert m_star(n, 1, d) == int(expected)


def test_larger_k_needs_a_smaller_sample():
    """More cleaved sites are easier to catch, so the budget falls as k rises."""
    n, delta = 28175, as_fraction("0.05")
    sizes = [m_star(n, k, delta) for k in (1, 2, 5, 10, 25)]
    assert sizes == sorted(sizes, reverse=True)


def test_tighter_confidence_needs_a_larger_sample():
    n, k = 28175, 5
    assert m_star(n, k, as_fraction("0.01")) > m_star(n, k, as_fraction("0.05"))


def test_float_confidence_is_refused():
    """0.05 is not 1/20 in binary, and a sample size off by one at a published level is a defect."""
    with pytest.raises(TypeError):
        as_fraction(0.05)
    with pytest.raises(TypeError):
        m_star(1000, 5, 0.05)


def test_string_and_fraction_agree():
    assert as_fraction("0.05") == Fraction(1, 20)
    assert m_star(5000, 3, "0.05") == m_star(5000, 3, Fraction(1, 20))


def test_invalid_arguments_are_rejected():
    with pytest.raises(ValueError):
        m_star(100, 0, "0.05")          # k must be at least 1
    with pytest.raises(ValueError):
        m_star(100, 101, "0.05")        # k cannot exceed N
    with pytest.raises(ValueError):
        m_star(100, 5, "1.5")           # delta must be a probability
    with pytest.raises(ValueError):
        miss_probability(100, 5, -1)


def test_budget_rows_are_self_consistent():
    rows = budget_rows(28175, k_values=(1, 5, 10, 25), delta1="0.05")
    assert [r["k"] for r in rows] == [1, 5, 10, 25]
    for r in rows:
        assert r["miss_probability"] <= Fraction(1, 20)
        assert r["fraction"] == Fraction(r["m_star"], 28175)
        assert 0 < r["fraction_float"] <= 1.0


def test_census_sensitivity_is_the_root_not_the_level():
    """The census threshold is 1 - delta^(1/k), which at k = 10 is far below 1 - delta.

    Conflating the two was a real error in an earlier draft: the two models differ, one being
    exact hypergeometric sampling with perfect detection and the other a full census with
    independent per-site detection.
    """
    s10 = float(census_sensitivity(10, "0.05"))
    assert 0.25 < s10 < 0.27
    s1 = float(census_sensitivity(1, "0.05"))
    assert abs(s1 - 0.95) < 1e-9
    assert s10 < s1


def test_census_sensitivity_rejects_bad_k():
    with pytest.raises(ValueError):
        census_sensitivity(0, "0.05")


def test_cross_check_against_the_reference_implementation_when_available():
    """If the companion m0_search module is importable, the two must agree exactly.

    Skipped rather than failed when it is absent, so this package stands alone.
    """
    m0 = pytest.importorskip("m0_search")
    for n in (500, 28175):
        for k in (1, 5, 10):
            assert m_star(n, k, "0.05") == m0.m0_hypergeom(n, k, Fraction(1, 20))
