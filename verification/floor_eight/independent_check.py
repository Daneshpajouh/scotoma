"""
INDEPENDENT check of the 564-template dual, by direct summation.

Shares NOTHING with the supplied verifier except the coefficient payload itself:
no transform, no eigenvalues, no modular arithmetic, no worker threads. Plain
expansion of the templates into actual words of Z_4^19, exact Fractions, and
brute-force Hamming distance.

Checks:
  1. every support word lies in U = union over i of B_8(p_i), p_i = 1^i 0^(19-i),
     computed by direct distance, NOT by the supplied nonzero-minus-peak formula
  2. total mass is exactly 30348/1573
  3. the ball at the reported attaining centre collects exactly 1
  4. a sample of other centres collects at most 1
  5. CONTROL: the same check on a deliberately corrupted certificate must FAIL

Check 4 is a SAMPLE: the staircase words and 120 random centres, 140 in all. It is not a sweep and
this script does not establish the maximum. The exhaustive maximisation over all 3^19 centre
classes is dual20.js, which shares nothing with this script but the payload.

Usage: python3 independent_check.py   (reads payload.bin from this script's own folder)
Exit 0 when every check passes; 1 when one fails; 2 if it could not run.
"""
import os
import random
import sys
from fractions import Fraction as Q

if not __debug__:
    print("NOT RUN: the checks are assertions and python -O disables them; run without -O")
    sys.exit(2)

N = 19
DEN = 1573
try:
    raw = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "payload.bin"), "rb").read()
except OSError as exc:
    print(f"NOT RUN: {exc}")
    sys.exit(2)
assert len(raw) == 564 * 5, len(raw)

# staircase: p_i = 1^i 0^(N-i), i = 0..N
STAIR = [tuple([1] * i + [0] * (N - i)) for i in range(N + 1)]

templates = []
for off in range(0, len(raw), 5):
    index = int.from_bytes(raw[off:off + 4], "little")
    b = raw[off + 4]
    trits, p = [], index
    for _ in range(N):
        trits.append(p % 3)
        p //= 3
    assert p == 0
    templates.append((tuple(trits), b))
assert len(templates) == 564
assert len({t for t, _ in templates}) == 564, "templates must be distinct"


def expand(trits):
    """star (trit 2) -> letter 2 or 3; trits 0,1 -> letters 0,1."""
    words = [[]]
    for c in trits:
        if c == 2:
            words = [w + [L] for w in words for L in (2, 3)]
        else:
            words = [w + [c] for w in words]
    return [tuple(w) for w in words]


def d(a, b):
    return sum(x != y for x, y in zip(a, b))


def dist_to_U(x):
    return min(d(x, p) for p in STAIR)


omega = {}
for trits, b in templates:
    ws = expand(trits)
    share = Q(b, DEN * len(ws))
    for w in ws:
        assert w not in omega, "support words must be distinct"
        omega[w] = share

print(f"support words: {len(omega):,}")

# 1. support lies in U, by direct distance
far = [x for x in omega if dist_to_U(x) > 8]
print(f"support words with distance to the staircase above 8: {len(far)}")
assert not far, "support escapes U"
dists = {dist_to_U(x) for x in omega}
print(f"distances from support to the staircase: {sorted(dists)}")

# 2. total mass
total = sum(omega.values(), Q(0))
print(f"total mass = {total} = {float(total):.10f}")
assert total == Q(30348, 1573)
assert total > 19


def ball_mass(c, table):
    return sum((v for x, v in table.items() if d(x, c) <= 8), Q(0))


CENTRE = tuple(int(ch) for ch in "1001010000001000000")
m = ball_mass(CENTRE, omega)
print(f"mass in the ball at the reported attaining centre = {m}")
assert m == 1, m

# 4. sample of other centres, including the staircase itself
random.seed(3)
probes = list(STAIR) + [tuple(random.randrange(4) for _ in range(N)) for _ in range(120)]
worst = max(ball_mass(c, omega) for c in probes)
print(f"worst ball mass over {len(probes)} probe centres = {worst} (must be <= 1)")
assert worst <= 1

# 5. CONTROL: corrupt one coefficient and confirm the attaining ball now exceeds 1.
# The template must have a word INSIDE that ball, or raising it cannot move the sum.
# (A first version corrupted templates[0], which lies outside the ball, so it could
# not fail; the control was wrong, not the certificate.)
inside = [(t, b) for t, b in templates if any(d(w, CENTRE) <= 8 for w in expand(t))]
outside = [(t, b) for t, b in templates if not any(d(w, CENTRE) <= 8 for w in expand(t))]
print(f"templates touching the attaining ball: {len(inside)}; not touching: {len(outside)}")


def corrupt(tb, by):
    trits, b = tb
    bad = dict(omega)
    ws = expand(trits)
    for w in ws:
        bad[w] = Q(b + by, DEN * len(ws))
    return bad


mb = ball_mass(CENTRE, corrupt(inside[0], 1))
print(f"CONTROL A, a touching template raised by 1: ball mass = {mb} -> exceeds 1: {mb > 1}")
assert mb > 1, "the control did not fail, so this check cannot detect a bad certificate"

mo = ball_mass(CENTRE, corrupt(outside[0], 40))
print(f"CONTROL B, a non-touching template raised by 40: ball mass = {mo} (must stay 1)")
assert mo == 1

print()
print("INDEPENDENT CHECK PASSED: support inside U, mass 30348/1573 > 19, "
      "attaining ball exactly 1, no probe above 1, and the control fails as it must.")
