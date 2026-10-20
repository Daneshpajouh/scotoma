"""From a candidate count to a validation budget, in exact arithmetic.

A screen that assays a subset of a guide's candidate off-target sites and finds nothing has
bounded, not eliminated, the risk. The bound is a sampling statement: if at least ``k`` of the
``N`` enumerated candidates were genuinely cleaved, how many sites ``m`` must be assayed before
missing all of them has probability at most ``delta1``?

Sampling ``m`` of ``N`` without replacement, the chance of missing a fixed set of ``k`` is
hypergeometric, and

    C(N - k, m) / C(N, m) = prod_{i=0}^{k-1} (N - m - i) / (N - i),

which is strictly decreasing in ``m``. So the required budget is

    m*(N, k, delta1) = min { m : prod_{i=0}^{k-1} (N - m - i) / (N - i) <= delta1 }.

Everything here is integer and :class:`~fractions.Fraction` arithmetic. A float ``delta1`` is
refused, because 0.05 is not 1/20 in binary and a sample size that is off by one at a published
confidence level is a defect, not a rounding detail.

**What a row licenses, and what it does not.** ``m*(N, k, delta1)`` certifies "fewer than ``k``
cleaved sites among the enumerated candidates" at confidence ``1 - delta1``. It licenses "none"
only under the additional assumption that an unclean guide cleaves at least ``k`` enumerated
sites. That assumption is about the enzyme and the enumeration, not about the sampling, and this
module cannot check it. Two further limits are structural: the budget is defined over the
candidate set you enumerated, so a candidate your reference could not represent is outside both
``N`` and the guarantee; and the whole statement is about sites, never about whether a sampled
site was correctly called.
"""
from __future__ import annotations

from fractions import Fraction


def as_fraction(x) -> Fraction:
    """Coerce to an exact :class:`~fractions.Fraction`, refusing floats.

    Accepts ``Fraction``, ``int``, or a string such as ``"0.05"`` or ``"1/20"``, both of which
    parse exactly.
    """
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, str):
        return Fraction(x)
    raise TypeError(
        f"refusing a {type(x).__name__} for an exact confidence level; pass a Fraction, an int, "
        'or a string such as "0.05" or "1/20"'
    )


def miss_probability(n_sites: int, k: int, m: int) -> Fraction:
    """Exact probability that a uniform sample of ``m`` of ``n_sites`` misses all of a fixed ``k``.

    Returns 0 once the sample is large enough that missing ``k`` is impossible.
    """
    n_sites, k, m = int(n_sites), int(k), int(m)
    if k <= 0:
        return Fraction(1)
    if m < 0 or m > n_sites:
        raise ValueError("m must lie in 0..n_sites")
    if k > n_sites:
        raise ValueError("k cannot exceed n_sites")
    if m > n_sites - k:
        return Fraction(0)
    num, den = 1, 1
    for i in range(k):
        num *= (n_sites - m - i)
        den *= (n_sites - i)
    return Fraction(num, den)


def m_star(n_sites: int, k: int, delta1) -> int:
    """Smallest sample size whose miss probability is at most ``delta1``.

    Found by binary search on a monotone predicate, so the cost is logarithmic in ``n_sites``
    and every comparison is exact.
    """
    n_sites, k = int(n_sites), int(k)
    delta = as_fraction(delta1)
    if not 0 < delta < 1:
        raise ValueError("delta1 must lie strictly between 0 and 1")
    if k <= 0:
        raise ValueError("k must be at least 1")
    if k > n_sites:
        raise ValueError("k cannot exceed n_sites")
    lo, hi = 0, n_sites - k + 1  # at m = n_sites - k + 1 the miss probability is exactly 0
    while lo < hi:
        mid = (lo + hi) // 2
        if miss_probability(n_sites, k, mid) <= delta:
            hi = mid
        else:
            lo = mid + 1
    return lo


def budget_rows(n_sites: int, k_values=(1, 5, 10, 25), delta1="0.05"):
    """A budget table: one row per ``k``, with the sample size and the fraction of ``N`` it is.

    Returns
    -------
    list of dict
        Each with ``k``, ``m_star``, ``fraction`` (an exact ``Fraction``), ``fraction_float``
        for display, and the exact ``miss_probability`` achieved at ``m_star``.
    """
    delta = as_fraction(delta1)
    rows = []
    for k in k_values:
        m = m_star(n_sites, k, delta)
        frac = Fraction(m, n_sites) if n_sites else Fraction(0)
        rows.append({
            "k": int(k),
            "m_star": m,
            "fraction": frac,
            "fraction_float": float(frac),
            "miss_probability": miss_probability(n_sites, k, m),
        })
    return rows


def census_sensitivity(k: int, delta1) -> Fraction:
    """Per-site sensitivity a full census must reach to match ``delta1``, if detection were independent.

    Assaying every one of the ``N`` sites with per-site detection probability ``s`` misses all
    ``k`` cleaved sites with probability ``(1 - s)^k``, so the requirement is
    ``s >= 1 - delta1 ** (1 / k)``.

    This is a **different model** from :func:`m_star`, which assumes perfect detection on a
    sampled subset. The two must not be presented as two readings of one column: one is exact
    hypergeometric sampling, the other assumes per-site detection events are independent, which
    is an assumption about the assay that nothing here verifies. Returned as a ``Fraction``
    bound only when ``k`` divides cleanly; otherwise the caller should use the float form, since
    a real ``k``-th root is not rational in general.
    """
    delta = as_fraction(delta1)
    k = int(k)
    if k <= 0:
        raise ValueError("k must be at least 1")
    return 1 - Fraction(float(delta) ** (1.0 / k)).limit_denominator(10 ** 12)
