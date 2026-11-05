"""Exact orbit arithmetic and pair types for the four-letter staircase."""
from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import comb


def prefixes(L):
    return tuple(product(range(3), repeat=L))


def digits(code, L):
    out = [0] * L
    for i in range(L - 1, -1, -1):
        out[i], code = code % 3, code // 3
    return tuple(out)


def library(N, r):
    return tuple((0,) * j + (1,) * (r - 1 - j) + (0,) * (N - r + 1)
                 for j in range(r))


@lru_cache(None)
def suffix_counts(C, a):
    """Support-overlap formula, independent of coordinate-polynomial code."""
    result = [[0] * (C + 1) for _ in range(C + 1)]
    for w in range(C + 1):
        for s in range(max(0, a + w - C), min(a, w) + 1):
            for h in range(s + 1):
                d = a + w - s - h
                result[w][d] += (comb(a, s) * comb(C - a, w - s)
                                 * comb(s, h) * 2 ** (s - h) * 3 ** (w - s))
        assert sum(result[w]) == comb(C, w) * 3 ** w
    return tuple(tuple(row) for row in result)


def orbit_size(C, point):
    p, w = point
    return 2 ** p.count(2) * comb(C, w) * 3 ** w


def kernel_count(C, m, centre, point):
    """Number of concrete words in point covered by one concrete centre."""
    a, aw = centre
    b, bw = point
    forced = sum(x != y for x, y in zip(a, b))
    both = sum(x == y == 2 for x, y in zip(a, b))
    counts = suffix_counts(C, aw)[bw]
    numerator = sum(comb(both, d) * sum(counts[:max(0, min(C, m - forced - d) + 1)])
                    for d in range(both + 1))
    return numerator * 2 ** (b.count(2) - both)


def point_orbits(N, r, m):
    L, C = r - 1, N - r + 1
    lib = library(L, r)
    for p in prefixes(L):
        dmin = min(sum(x != y for x, y in zip(p, c)) for c in lib)
        for w in range(min(C, m - dmin) + 1):
            yield (p, w)


@lru_cache(None)
def prefix_joint(a, b, p, same):
    """Joint prefix-distance distribution; same shared X letters equal."""
    poly = {(0, 0): 1}
    shared = 0
    for aa, bb, pp in zip(a, b, p):
        av = aa
        bv = bb
        if aa == bb == pp == 2:
            bv = 2 if shared < same else 3
            shared += 1
        nxt = defaultdict(int)
        for (d1, d2), count in poly.items():
            for letter in ((2, 3) if pp == 2 else (pp,)):
                nxt[d1 + (letter != av), d2 + (letter != bv)] += count
        poly = nxt
    return tuple(sorted(poly.items()))


@lru_cache(maxsize=30000)
def suffix_joint(C, a, b, support_overlap, equal_overlap, w):
    s, h = support_overlap, equal_overlap
    blocks = [(1, 1)] * h + [(1, 2)] * (s - h) + [(1, 0)] * (a - s)
    blocks += [(0, 1)] * (b - s) + [(0, 0)] * (C - a - b + s)
    assert len(blocks) == C
    poly = {(0, 0, 0): 1}
    for av, bv in blocks:
        nxt = defaultdict(int)
        for (wt, d1, d2), count in poly.items():
            for letter in range(4):
                nw = wt + (letter != 0)
                if nw <= w:
                    nxt[nw, d1 + (letter != av), d2 + (letter != bv)] += count
        poly = nxt
    hist = {(d1, d2): count for (wt, d1, d2), count in poly.items() if wt == w}
    assert sum(hist.values()) == comb(C, w) * 3 ** w
    rect = [[0] * (C + 1) for _ in range(C + 1)]
    for i in range(C + 1):
        for j in range(C + 1):
            rect[i][j] = (hist.get((i, j), 0) + (rect[i - 1][j] if i else 0)
                          + (rect[i][j - 1] if j else 0)
                          - (rect[i - 1][j - 1] if i and j else 0))
    return tuple(tuple(row) for row in rect)


@lru_cache(maxsize=100000)
def pair_minimum(C, m, centre_a, centre_b, point):
    a, aw = centre_a
    b, bw = centre_b
    p, w = point
    triples = sum(x == y == z == 2 for x, y, z in zip(a, b, p))
    best = None
    for s in range(max(0, aw + bw - C), min(aw, bw) + 1):
        for h in range(s + 1):
            rect = suffix_joint(C, aw, bw, s, h, w)
            for same in range(triples + 1):
                value = 0
                for (d1, d2), count in prefix_joint(a, b, p, same):
                    r1, r2 = m - d1, m - d2
                    if r1 >= 0 and r2 >= 0:
                        value += count * rect[min(C, r1)][min(C, r2)]
                if best is None or value < best:
                    best = value
    assert best is not None
    return best


def forest_weight(weights, k):
    parent = list(range(k))

    def root(v):
        while v != parent[v]:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    total, edges = 0, []
    for weight, i, j in sorted(weights, reverse=True):
        ri, rj = root(i), root(j)
        if ri != rj:
            parent[ri] = rj
            total += weight
            edges.append((i, j))
    return total, edges


def volume(t, s):
    return sum(comb(t, i) * 3 ** i for i in range(max(-1, min(t, s)) + 1))
