#!/usr/bin/env python3
"""The non-chain centres that also attain unit load, and why the theorem still holds.

WHY THIS EXISTS. The body once said the private-sphere weighting attains Lambda = 1 "exactly at
the chain words". That converse is false. At the threshold length N*(m) there are non-chain
centres with Lambda exactly 1, and this script produces them so the boundary is on the record
instead of being left to a reader to discover. The theorem is unaffected: it asserts Lambda < 1
for centres with u(c) >= 2, and every centre below carries u(c) < 2.

Exact rational arithmetic throughout. No search, no sampling: the two cells are enumerated in
full at m = 1 and the single named centre is evaluated exactly at m = 6.
"""
import json
import math
import os
from fractions import Fraction
from itertools import combinations, product
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))


def nstar(m):
    return math.ceil((3 * m + 3 + math.isqrt(5 * m * m - 2 * m + 1)) / 2) \
        if (5 * m * m - 2 * m + 1) == math.isqrt(5 * m * m - 2 * m + 1) ** 2 \
        else math.ceil((3 * m + 3 + math.sqrt(5 * m * m - 2 * m + 1)) / 2)


def support(N, j):
    """A_j = [N] \\ {j, j+1}, read 1-based as the appendix writes it."""
    drop = {k for k in (j - 1, j) if 0 <= k < N}
    return [i for i in range(N) if i not in drop]


def load(N, m, c):
    total = Fraction(0)
    for j in range(N + 1):
        A = support(N, j)
        if len(A) < m:
            continue
        pat = [1] * j + [0] * (N - j)
        hits = 0
        for pos in combinations(A, m):
            fixed = sum(1 for i in range(N) if i not in pos and pat[i] != c[i])
            if fixed > m:
                continue
            for letters in product((2, 3), repeat=m):
                if fixed + sum(1 for i, L in zip(pos, letters) if L != c[i]) <= m:
                    hits += 1
        total += Fraction(hits, comb(len(A), m) * 2 ** m)
    return total


def u(c):
    return sum(1 for x in c if x in (2, 3))


def chain(N):
    return {tuple([1] * j + [0] * (N - j)) for j in range(N + 1)}


def main():
    out = {'what': 'non-chain centres with unit private-sphere load at the threshold length',
           'theorem_scope': 'the theorem asserts Lambda(c) < 1 for u(c) >= 2; every centre here '
                            'has u(c) < 2, so none of them contradicts it',
           'cells': []}
    # m = 1: enumerate the whole cube at the threshold length and the two above it
    for m, widths in ((1, (4, 5, 6)),):
        for N in widths:
            ch = chain(N)
            ties = []
            worst_u2 = None
            for c in product((0, 1, 2, 3), repeat=N):
                v = load(N, m, c)
                if u(c) >= 2 and (worst_u2 is None or v > worst_u2):
                    worst_u2 = v
                if v == 1 and c not in ch:
                    ties.append({'centre': ''.join(map(str, c)), 'u': u(c), 'load': '1'})
            assert all(load(N, m, c) == 1 for c in ch), 'a chain word lost unit load'
            out['cells'].append({'m': m, 'N': N, 'nstar': nstar(m), 'exhaustive': True,
                                 'non_chain_unit_load': ties,
                                 'max_load_at_u_at_least_2': str(worst_u2)})
    # m = 6 at its threshold length: the cube is far too large to enumerate, so one named centre
    # is evaluated exactly. This is a witness, not a census, and is labelled as such.
    m, N = 6, nstar(6)
    c = tuple([1] * 7 + [0] * (N - 8) + [2])
    v = load(N, m, c)
    out['cells'].append({'m': m, 'N': N, 'nstar': N, 'exhaustive': False,
                         'witness_only': True,
                         'non_chain_unit_load': [{'centre': ''.join(map(str, c)),
                                                  'u': u(c), 'load': str(v)}],
                         'max_load_at_u_at_least_2': None})
    assert v == 1, 'the named m=6 witness no longer attains unit load'
    path = os.path.join(HERE, 'boundary_profiles.json')
    open(path, 'w', encoding='utf-8').write(json.dumps(out, indent=2) + '\n')
    for cell in out['cells']:
        print('m=%d N=%d  non-chain unit-load centres: %d  %s  max load at u>=2: %s'
              % (cell['m'], cell['N'], len(cell['non_chain_unit_load']),
                 'exhaustive' if cell['exhaustive'] else 'one named witness',
                 cell['max_load_at_u_at_least_2']))
    print('written: boundary_profiles.json')


if __name__ == '__main__':
    main()
