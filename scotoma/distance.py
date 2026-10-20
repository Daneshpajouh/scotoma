"""Spacer mismatch distance between packed windows and a guide panel.

Distance is counted over the 20-base spacer only. The NGG PAM is excluded: every window this
package extracts has a valid PAM by construction, so including it would add a constant zero and
make the metric depend on the N base, which is not part of recognition.

The kernel is branch-free and vectorised over sites. For a packed window ``w`` and guide ``g``,
``w ^ g`` has a nonzero 2-bit group exactly where the bases differ. Shifting right by
``2 * PAM_LEN`` discards the PAM, after which collapsing each 2-bit group to a single bit and
population-counting gives the number of mismatched spacer positions.

Nothing here counts cleavage. A "candidate" is a site within a mismatch threshold of some guide
in the panel; whether it is cut is a separate question this package does not address.
"""
from __future__ import annotations

import numpy as np

from .coding import PAM_LEN, SPACER_LEN

_LOW = np.uint64(0x5555555555555555)
_M2 = np.uint64(0x3333333333333333)
_M4 = np.uint64(0x0F0F0F0F0F0F0F0F)
_ONES = np.uint64(0x0101010101010101)
_SHIFT_PAM = np.uint64(2 * PAM_LEN)
_1, _2, _4, _56 = (np.uint64(x) for x in (1, 2, 4, 56))

#: Returned for sites no guide comes near; larger than any achievable distance.
NO_MATCH = np.int8(99)


def spacer_mismatches(codes, guide_code):
    """Mismatch count between every code and one guide, over the spacer.

    Parameters
    ----------
    codes : numpy.ndarray
        uint64 **oriented** window codes. Passing canonical codes here is the error this
        package exists to prevent; see :mod:`scotoma.coding`.
    guide_code : int
        A packed 23-mer for the guide. Only its spacer bits are used.

    Returns
    -------
    numpy.ndarray
        int8, in ``0 .. SPACER_LEN``.
    """
    codes = np.asarray(codes, dtype=np.uint64)
    x = (codes ^ np.uint64(guide_code)) >> _SHIFT_PAM
    # Collapse each 2-bit group to its low bit if either bit is set, then popcount the result.
    y = (x | (x >> _1)) & _LOW
    c = y - ((y >> _1) & _LOW)
    c = (c & _M2) + ((c >> _2) & _M2)
    c = (c + (c >> _4)) & _M4
    return ((c * _ONES) >> _56).astype(np.int8)


def min_mismatches(codes, guide_codes, chunk: int = 4_000_000):
    """For every code, the smallest spacer mismatch count over the whole panel.

    Processed in chunks so peak memory stays near ``chunk`` elements regardless of input size;
    a haplotype contributes tens of millions of sites.

    Returns
    -------
    numpy.ndarray
        int8 of the same length as ``codes``. :data:`NO_MATCH` where the panel is empty.
    """
    codes = np.asarray(codes, dtype=np.uint64)
    guides = [int(g) for g in np.asarray(guide_codes, dtype=np.uint64).tolist()]
    best = np.full(codes.size, NO_MATCH, np.int8)
    if not guides or codes.size == 0:
        return best
    for start in range(0, codes.size, chunk):
        stop = min(start + chunk, codes.size)
        block = codes[start:stop]
        acc = np.full(block.size, NO_MATCH, np.int8)
        for g in guides:
            np.minimum(acc, spacer_mismatches(block, g), out=acc)
        best[start:stop] = acc
    return best


def threshold_counts(best, thresholds=(0, 1, 2, 3, 4, 5, 6)):
    """Cumulative counts of sites at or below each mismatch threshold.

    Parameters
    ----------
    best : numpy.ndarray
        Output of :func:`min_mismatches`.

    Returns
    -------
    dict
        ``{threshold: count}``, counts being cumulative, so each is at most the next.
    """
    best = np.asarray(best)
    return {int(t): int((best <= t).sum()) for t in thresholds}


def density_per_million(counts, n_sites: int):
    """Convert threshold counts to candidates per million sites.

    This is the quantity that makes two arms of different size comparable. Returns NaN for an
    empty arm rather than dividing by zero.
    """
    if n_sites <= 0:
        return {int(t): float("nan") for t in counts}
    return {int(t): counts[t] / (n_sites / 1e6) for t in counts}
