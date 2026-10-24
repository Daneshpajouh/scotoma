import sys

if not __debug__:
    print("NOT RUN: the checks are assertions and python -O disables them; run without -O")
    sys.exit(2)
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import comb

C = [
    ("000100", 7, 4), ("001000", 7, 1),
    ("010111", 6, 1), ("011000", 7, 4),
    ("011110", 7, 1), ("100001", 7, 1),
    ("100011", 6, 4), ("110000", 8, 5),
    ("111000", 8, 1), ("111100", 8, 5),
    ("111110", 8, 4), ("111111", 8, 4),
]

def d(x, y):
    return sum(a != b for a, b in zip(x, y))

T = ["1"*i + "0"*(6-i) for i in range(7)]

# Verify that the support belongs to U and has total mass seven.
assert all(w + min(d(y, t) for t in T) == 8
           for y, w, a in C)
assert sum(Q(a, 5) for y, w, a in C) == 7

@cache
def F(k, w, h):
    if h < 0:
        return Q(0)
    count = 0
    for j in range(max(0, w-(13-k)), min(k, w)+1):
        for b in range(j+1):
            if k+w-j-b <= h:
                count += (
                    comb(k, j) * comb(13-k, w-j)
                    * comb(j, b) * 2**(j-b) * 3**(w-j)
                )
    return Q(count, comb(13, w) * 3**w)

def mass(p, k):
    return sum(
        (Q(a, 5) * F(k, w, 8-d(p, y)) for y, w, a in C),
        Q(0),
    )

maximum = Q(0)
checked = 0
for bits in product("01", repeat=6):
    p = "".join(bits)
    for k in range(14):
        value = mass(p, k)
        assert value <= 1, (p, k, value)
        maximum = max(maximum, value)
        checked += 1

assert all(mass(t, 0) == 1 for t in T)
assert maximum == 1
assert checked == 896
print("total mass = 7; maximum ball mass = 1; checks = 896")
