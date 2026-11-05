#!/usr/bin/env python3
"""LEMMA. The single-deletion library is isometric to the integer line: d(P_i, P_j) = |i - j|.

PROOF, as derived, and checked below rather than trusted.

Let g have runs R_1..R_r, with run i ending at position e_i (0-indexed). Deleting any character of
run i yields the same word, namely g with position e_i removed, so D_1(g) = {P_1,..,P_r} and
P_i = g[:e_i] + g[e_i+1:].

Take e_i < e_j. The two patterns agree on 0..e_i-1. For e_i <= k <= e_j-1, P_i[k] = g[k+1] while
P_j[k] = g[k], so they differ exactly when g[k] != g[k+1], that is exactly when k is the end of a
run. For k >= e_j both equal g[k+1] and they agree. The run ends in [e_i, e_j-1] are precisely
e_i, e_{i+1}, .., e_{j-1}, which is j-i of them. Hence

    d(P_i, P_j) = j - i.

COROLLARY. Patterns pairwise within distance t correspond to indices pairwise within t, so the
largest such set is an interval of t+1 runs. At t=2 the maximum is exactly 3, for every word with
r >= 3, attained by any three consecutive runs. That is the clique-3 ceiling, now with a reason.

WHY IT MATTERS. The library is a geodesic path, maximally spread for its size, which is the property
the 2m-scattered transversal exploits. For b=2 the index set is two-dimensional and the ceiling
rises, which is why the same route dies there.
"""
import os, random, re, sys
from itertools import combinations, product


def dele(g, b):
    """Distinct words from deleting exactly b characters. Inlined so this file stands alone in the
    release: it previously imported a sibling that was never published, and failed on a clean
    clone with ModuleNotFoundError."""
    out = set()
    n = len(g)
    for idx in combinations(range(n), b):
        s = set(idx)
        out.add(''.join(c for i, c in enumerate(g) if i not in s))
    return sorted(out)


def ham(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)



def runs_and_patterns(g):
    """Patterns IN RUN ORDER, built from run ends, so the index i is the run index."""
    ends = [k for k in range(len(g) - 1) if g[k] != g[k + 1]] + [len(g) - 1]
    return [g[:e] + g[e + 1:] for e in ends]


def check(g):
    P = runs_and_patterns(g)
    if sorted(P) != sorted(dele(g, 1)):
        return 'PATTERN SET MISMATCH'
    if len(set(P)) != len(P):
        return 'PATTERNS NOT DISTINCT'
    for i, j in combinations(range(len(P)), 2):
        if ham(P[i], P[j]) != j - i:
            return 'ISOMETRY FAILS at (%d,%d): %d != %d' % (i, j, ham(P[i], P[j]), j - i)
    return None


bad, tested = [], 0
for g in sorted(set(re.findall(r'"([ACGT]{20})"', open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'panel_match.json')).read()))):
    tested += 1
    e = check(g)
    if e: bad.append((g, e))
print('panel: %d guides, failures %d' % (tested, len(bad)))

for alpha in ('AC', 'ACG', 'ACGT'):
    n0 = tested
    for L in range(2, 11):
        if len(alpha) ** L > 300000:
            continue
        for t in product(alpha, repeat=L):
            tested += 1
            e = check(''.join(t))
            if e: bad.append((''.join(t), e))
    print('alphabet %-4s exhaustive to the budget: cumulative tested %d, failures %d'
          % (alpha, tested, len(bad)))

rng = random.Random(20261024)
for L in (20, 23, 30, 47):
    for _ in range(3000):
        g = ''.join(rng.choice('ACGT') for _ in range(L))
        tested += 1
        e = check(g)
        if e: bad.append((g, e))
print('random lengths 20/23/30/47: cumulative tested %d, failures %d' % (tested, len(bad)))

print()
if bad:
    print('LEMMA FALSIFIED on %d words. First: %s -> %s' % (len(bad), bad[0][0], bad[0][1]))
    sys.exit(1)
print('d(P_i, P_j) = |i - j| held on every one of %d words, with the pattern set matching the' % tested)
print('independent deletion enumeration each time. The clique-3 ceiling follows as a corollary,')
print('so it is a LEMMA with a proof, not a panel accident.')
