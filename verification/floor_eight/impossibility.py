import sys

if not __debug__:
    print("NOT RUN: the checks are assertions and python -O disables them; run without -O")
    sys.exit(2)
from math import comb
from fractions import Fraction as Q

def binary_bound(N, m, r):
    h = r - 1
    L = N - h
    V = sum(
        comb(h, i) * comb(L, j) * 3**j
        for i in range(h + 1)
        for j in range(L + 1)
        if i + j <= m
    )
    return Q(2**h * 4**L, V)

def full_bound(N, m):
    V = sum(comb(N, j) * 3**j
            for j in range(min(N, m) + 1))
    return Q(4**N, V)

for r in (9, 10, 11, 12, 13, 15, 20):
    b = binary_bound(19, 8, r)
    print(r, b, "barrier:", b <= r - 1)

# First failure certified by this bound, not necessarily
# the first failure of the actual binary-support LP.
first = next(
    r for r in range(7, 21)
    if binary_bound(19, 8, r) <= r - 1
)
assert first == 13
assert all(binary_bound(19, 8, r) < r - 1
           for r in range(13, 21))

b = full_bound(19, 12)
assert b == Q(68719476736, 12017379241)
assert b < 6
print("Unrestricted far-corner upper bound:", b)
