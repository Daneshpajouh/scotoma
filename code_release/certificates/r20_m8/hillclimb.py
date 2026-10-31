"""
A requested load check, independent of the orbit reduction and the transform:
explicit words over all four letters, exact integer weights, plain summation.
Weights: omega(x) = b / (1573 * 2^k); scaled by 1573 * 256 to integers W(x) = b * 2^(8-k), so a
ball is feasible iff its integer load <= 402,688 (= 256 * 1573). k <= 8 for every template.
Checks: (1) all 20 staircase words; (2) every centre in their radius-2 neighbourhoods (20 x 1,597);
(3) steepest-ascent hill climbing from many random starts over Z_4^19, each step trying all 57
neighbours, plus random restarts seeded near the known maximiser. Reports the largest load found.
CONTROL: the same search on a certificate with one coefficient raised must find a load > 402,688.
"""
import itertools, random, sys
import numpy as np
N = 19; CAP = 256 * 1573
raw = open(sys.argv[1], "rb").read()
words, wts = [], []
for off in range(0, len(raw), 5):
    idx = int.from_bytes(raw[off:off+4], "little"); b = raw[off+4]
    t, p = [], idx
    for _ in range(N): t.append(p % 3); p //= 3
    k = t.count(2)
    for fill in itertools.product((2, 3), repeat=k):
        it = iter(fill); words.append([next(it) if c == 2 else c for c in t]); wts.append(b * 2 ** (8 - k))
S = np.array(words, dtype=np.int8); W = np.array(wts, dtype=np.int64)
print("support words", len(S), "total integer mass", int(W.sum()), "= 30348*256 ->", int(W.sum()) == 30348 * 256)

def load(c, Wv):
    return int(Wv[(S != c).sum(1) <= 8].sum())

def neighbours(c):
    for j in range(N):
        for a in range(4):
            if a != c[j]:
                d = c.copy(); d[j] = a; yield d

def search(Wv, starts, seed):
    rng = random.Random(seed); best = (0, None)
    stair = [np.array([1]*i + [0]*(N-i), dtype=np.int8) for i in range(20)]
    for c in stair:
        best = max(best, (load(c, Wv), tuple(c)))
    for c in stair:                                   # radius-2 neighbourhoods, exhaustive
        for d in neighbours(c):
            best = max(best, (load(d, Wv), tuple(d)))
            for e in neighbours(d):
                best = max(best, (load(e, Wv), tuple(e)))
    known = np.array([int(x) for x in "1001010000001000000"], dtype=np.int8)
    for s in range(starts):
        if s % 2 == 0:
            c = np.array([rng.randrange(4) for _ in range(N)], dtype=np.int8)
        else:                                          # perturb the known maximiser
            c = known.copy()
            for _ in range(rng.randrange(1, 7)):
                c[rng.randrange(N)] = rng.randrange(4)
        cur = load(c, Wv)
        while True:
            nb = max(((load(d, Wv), tuple(d)) for d in neighbours(c)), key=lambda x: x[0])
            if nb[0] <= cur: break
            cur, c = nb[0], np.array(nb[1], dtype=np.int8)
        best = max(best, (cur, tuple(c)))
    return best

starts = int(sys.argv[2]) if len(sys.argv) > 2 else 200
best = search(W, starts, 1)
print(f"largest load found: {best[0]} vs cap {CAP} -> {'OK, never above' if best[0] <= CAP else 'VIOLATION'}; "
      f"at {''.join(map(str, best[1]))}")
# CONTROL: raise one coefficient that touches the known maximiser's ball
known = np.array([int(x) for x in "1001010000001000000"], dtype=np.int8)
inball = np.where((S != known).sum(1) <= 8)[0]
W2 = W.copy(); W2[inball[0]] += 1
b2 = search(W2, 20, 2)
print(f"CONTROL one weight raised by 1 inside that ball: largest load {b2[0]} -> {'caught (> cap)' if b2[0] > CAP else 'MISSED'}")
