# K(19, 8, 7) = 7, proved by a fractional dual certificate

2026-10-22. The first cell of the previously open region (n = 20, m = 8 to 12, run count r >= 7)
to be closed. Six radius-8 balls cannot cover the union of the seven radius-8 balls around the
canonical staircase at N = 19 over four letters.

## The certificate

Twelve support classes, each an orbit of (binary active prefix y, passive suffix weight w), given
mass a/5 spread uniformly over the orbit. Numerators sum to 35, so total mass is exactly 7.

  000100 7 4   001000 7 1   010111 6 1   011000 7 4
  011110 7 1   100001 7 1   100011 6 4   110000 8 5
  111000 8 1   111100 8 5   111110 8 4   111111 8 4

Every triple satisfies w + min_i d(y, 1^i 0^(6-i)) = 8, so the support lies inside U.

The ball constraint reduces to 896 cases: an unrestricted centre may be replaced by one with a
binary active prefix without increasing any distance to a supported word, and the suffix
intersection depends on the centre only through its suffix weight.

## How it was verified, and the part that matters

`verify_m8r7.py` checks the 896 inequalities in exact rationals: total mass 7, maximum ball mass
exactly 1, exit 0.

**That verifier cannot check its own closed-form intersection F.** F claims: among suffix words of
weight w, the fraction within distance h of a centre suffix of weight k. The verifier uses F in
every inequality, so if F were wrong the whole certificate would pass and mean nothing.

So F was brute-forced directly against enumeration at L = 3, 4, 5 over every (k, w, h):
**482 cells, zero mismatches.** And the comparison can fail: a deliberately shifted F disagrees on
90 cells. That control is the reason this is recorded as proved rather than as reported.

## What it closes, and what it does not

A positive result propagates DOWN the down-set, so K = r now holds at every cell with m <= 8 and
r <= 7. It does **not** reach the far corner at m = 12, r = 20, and nothing here should be read as
closing the interval.

## Provenance

Constructed from a problem statement that supplied the staircase, the passive-coordinate
decomposition and the warning about the large-radius regime. Verified and controlled independently,
against the certificate alone, before it was recorded.
