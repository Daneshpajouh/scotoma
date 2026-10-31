#!/usr/bin/env python3
"""Exact standalone verification of the five-pattern full-mass construction.

Standard library only. The infinite proof uses a proved degree bound and exact
Newton coefficients, followed by a sampling-without-replacement moment bound.
Finite checks count actual tail words by a separate intersection formula.
"""
import sys
sys.dont_write_bytecode = True
import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
PREFIXES = tuple(product(range(3), repeat=4))
LIBRARY = tuple((0,)*j + (1,)*(4-j) for j in range(5))
# (prefix, difference from the closest library prefix, total orbit mass).
SUPPORT = tuple((t, 0, F(1, 5)) for t in LIBRARY) + tuple(
    (tuple(map(int, t)), 1, F(a, 5))
    for t, a in [('0012', 2), ('2011', 2), ('0112', 4),
                 ('2001', 4), ('0200', 4), ('1121', 4)])
MAXIMAL = (
    ((1, 2, 2, 0, 0), (0, 4, 8, 8, 0)),
    ((0, 1, 3, 1, 0), (1, 3, 4, 8, 4)),
    ((0, 2, 3, 0, 0), (0, 0, 16, 4, 0)),
    ((0, 1, 2, 2, 0), (2, 3, 1, 10, 4)),
)
DEFAULT_N = list(range(10, 23)) + [27, 32, 52, 102, 103, 137, 202]
CUTOFF = 46


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def distance(a, b):
    return sum(x != y for x, y in zip(a, b))


@lru_cache(None)
def prefix_histogram(fixed, counted):
    """Literal enumeration of the at most two free-letter assignments."""
    h = [F(0)] * 5
    options = [(2, 3) if x == 2 else (x,) for x in counted]
    words = tuple(product(*options))
    for word in words:
        h[distance(fixed, word)] += F(1, len(words))
    return tuple(h)


def profile(c):
    h = [[F(0)] * 5 for _ in range(2)]
    for t, d, mass in SUPPORT:
        for j, prob in enumerate(prefix_histogram(c, t)):
            h[d][j] += 5 * mass * prob
    require(all(v.denominator == 1 for row in h for v in row),
            'nonintegral prefix profile')
    return tuple(tuple(int(x) for x in row) for row in h)


def dominates(upper, lower):
    return all(sum(a[:j]) >= sum(b[:j])
               for a, b in zip(upper, lower) for j in range(1, 6))


def structural_check():
    per_region = [F(0)] * 5
    for t, d, mass in SUPPORT:
        ds = [distance(t, p) for p in LIBRARY]
        require(min(ds) == d and ds.count(d) == 1,
                'prefix is not uniquely private: %s' % (t,))
        require(mass > 0, 'nonpositive mass')
        per_region[ds.index(d)] += mass
    require(per_region == [F(1)] * 5, 'private-region masses are not all one')
    actual = {c: profile(c) for c in PREFIXES}
    require(len(set(actual.values())) == 44, 'unexpected prefix-profile count')
    for upper in MAXIMAL:
        require(upper in actual.values(), 'proposed maximal profile is absent')
    domination_map = {}
    for c, pr in actual.items():
        choices = [i for i, upper in enumerate(MAXIMAL) if dominates(upper, pr)]
        require(choices, 'prefix profile has no dominating row: %s' % (c,))
        domination_map[''.join(map(str, c))] = choices[0]
    return actual, domination_map


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(None)
def favorable(p, a, k):
    """3^p times conditional gain-tail probability with a occupied positions."""
    return 3**(p-a) * sum(comb(a, b) * 2**(a-b)
                          for b in range(max(0, p+k-a), a+1))


def gain_numerator(C, w, p, k):
    # Denominator is choose(C,p)*3^p. This is a polynomial of degree <= p
    # in q after C=5q+s and w=3q+b: each term has degrees a and p-a.
    return sum(choose(w, a) * choose(C-w, p-a) * favorable(p, a, k)
               for a in range(p+1))


def polynomial_values(residue, p, q):
    C = 5*q + residue
    m = 3*q + (3*residue+4)//5
    den = comb(C, p) * 3**p
    x = [gain_numerator(C, m, p, j) for j in range(5)]
    y = [gain_numerator(C, m-1, p, j-1) for j in range(5)]
    return [5*den - sum(a*u+b*v for a, b, u, v in zip(A, B, x, y))
            for A, B in MAXIMAL]


def forward_differences(values):
    out = []
    while values:
        out.append(values[0])
        values = [b-a for a, b in zip(values, values[1:])]
    return out


def polynomial_proof():
    rows = []
    checked_coefficients = 0
    for residue in range(5):
        for p in range(CUTOFF):
            q0 = max(1, (p-residue+4)//5)
            values = [polynomial_values(residue, p, q0+i)
                      for i in range(p+1)]
            # An extra evaluation checks the declared degree implementation.
            extra = polynomial_values(residue, p, q0+p+1)
            for i in range(4):
                coeff = forward_differences([v[i] for v in values])
                require(all(F(v) >= 0 for v in coeff),
                        'Newton coefficient FAILED: residue=%d p=%d profile=%d'
                        % (residue, p, i))
                require(sum(v*comb(p+1, k) for k, v in enumerate(coeff)) == extra[i],
                        'degree/reconstruction check FAILED')
                checked_coefficients += len(coeff)
                rows.append({'residue': residue, 'p': p, 'profile': i,
                             'q0': q0, 'newton_coefficients': list(map(str, coeff))})
    require(len(rows) == 920, 'wrong polynomial-row count')
    return rows, checked_coefficients


def large_distance_proof(actual_profiles):
    t = F(7, 5)
    potentials = [sum((F(a, 5)*t**(-j) + F(b, 5)*t**(1-j)
                      for j, (a, b) in enumerate(zip(A, B))), F(0))
                  for A, B in actual_profiles.values()]
    potential = max(potentials)
    require(potential == F(169, 49), 'prefix exponential potential FAILED')
    rho = F(3, 5) + F(4, 5*CUTOFF)
    h = (1-rho)/t + rho*(2+t)/3
    require(h == F(11749, 12075) and 0 < h < 1, 'moment factor FAILED')
    bound = potential * h**CUTOFF
    require(bound < F(99, 100), 'large-distance bound FAILED')
    return {'cutoff': CUTOFF, 't': str(t), 'rho_upper': str(rho),
            'moment_factor': str(h), 'prefix_potential': str(potential),
            'uniform_load_bound': str(bound), 'below_99_over_100': True,
            'strictly_below_one': True}


@lru_cache(None)
def tail_distance_cumulative(C, p, w):
    """Count actual weight-w tail words, rather than sampling p positions.

    Select the overlap with the p nonzero centre positions, then its matches.
    The normalization is choose(C,w)*3^w, different from gain_numerator.
    """
    h = [0] * (C+1)
    for overlap in range(max(0, w-(C-p)), min(p, w)+1):
        base = comb(p, overlap) * comb(C-p, w-overlap) * 3**(w-overlap)
        for matches in range(overlap+1):
            d = p+w-overlap-matches
            h[d] += base * comb(overlap, matches) * 2**(overlap-matches)
    total = comb(C, w) * 3**w
    require(sum(h) == total, 'tail histogram normalization FAILED')
    for j in range(1, len(h)):
        h[j] += h[j-1]
    return tuple(h), total


def tail_share(C, p, w, budget):
    if budget < 0:
        return F(0)
    counts, total = tail_distance_cumulative(C, p, w)
    return F(counts[min(C, budget)], total)


def load(C, m, c, p, support=SUPPORT, wrong_tail=False):
    value = F(0)
    for t, d, mass in support:
        w = m if wrong_tail else m-d
        value += mass * sum((prob*tail_share(C, p, w, m-j)
                             for j, prob in enumerate(prefix_histogram(c, t))), F(0))
    return value


def finite_check(n, actual_profiles, radius=None, wrong_tail=False):
    require(n >= 10, 'n must be at least ten')
    C = n-5
    m = (3*n-11)//5 if radius is None else radius
    per_region = [F(0)] * 5
    for t, d, mass in SUPPORT:
        w = m if wrong_tail else m-d
        ds = [distance(t, lib)+w for lib in LIBRARY]
        require(0 <= w <= C and min(ds) == m and ds.count(m) == 1,
                'privacy FAILED at prefix %s, tail weight %d, radius %d'
                % (''.join(map(str, t)), w, m))
        per_region[ds.index(m)] += mass
    require(per_region == [F(1)]*5, 'finite region mass FAILED')
    worst = F(0)
    witness = None
    for p in range(C+1):
        x = [tail_share(C, p, m, m-j) for j in range(5)]
        y = [tail_share(C, p, m-1, m-j) for j in range(5)]
        for c, (A, B) in actual_profiles.items():
            value = sum((F(a, 5)*u+F(b, 5)*v
                         for a, b, u, v in zip(A, B, x, y)), F(0))
            if value > worst:
                worst, witness = value, (''.join(map(str, c)), p)
    require(worst <= 1, 'load FAILED: %s at %s' % (worst, witness))
    classes = 81*(C+1)
    tail_distance_cumulative.cache_clear()
    return {'n': n, 'N': n-1, 'C': C, 'm': m, 'mass': '5',
            'private_region_masses': ['1']*5, 'max_load': str(worst),
            'argmax': witness, 'centre_classes': classes,
            'status': 'VERIFIED-COMPUTATIONALLY'}


def negative_controls(actual_profiles):
    controls = []
    # Keep mass five and privacy, but move the construction beyond the bound.
    witness_load = load(12, 9, (0, 0, 1, 1), 1)
    require(witness_load == F(41, 36) and witness_load > 1,
            'above-boundary witness unexpectedly passed')
    try:
        finite_check(17, actual_profiles, radius=9)
    except AssertionError as exc:
        require('load FAILED' in str(exc), 'wrong above-boundary rejection')
        controls.append({'control': 'n=17, radius 9 instead of 8',
                         'outcome': 'REJECTED-AS-REQUIRED', 'reason': str(exc),
                         'explicit_witness': {'prefix': '0011', 'tail_weight': 1,
                                              'load': str(witness_load)}})
    else:
        raise AssertionError('above-boundary control was accepted')
    try:
        finite_check(17, actual_profiles, wrong_tail=True)
    except AssertionError as exc:
        require('privacy FAILED' in str(exc), 'wrong wrong-tail rejection')
        controls.append({'control': 'give free-letter orbits tail weight m',
                         'outcome': 'REJECTED-AS-REQUIRED', 'reason': str(exc)})
    else:
        raise AssertionError('wrong-tail control was accepted')
    tail_distance_cumulative.cache_clear()
    return controls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, nargs='+', default=DEFAULT_N)
    parser.add_argument('--write-results', action='store_true',
                        help='write proof_certificate.json and verification_results.json beside this script')
    args = parser.parse_args()
    actual_profiles, domination_map = structural_check()
    coefficients, count = polynomial_proof()
    large = large_distance_proof(actual_profiles)
    certificate = {'status': 'VERIFIED-COMPUTATIONALLY',
                   'degree_bound': 'p', 'basis': 'binomial(q-q0,k)',
                   'prefix_domination': domination_map, 'polynomials': coefficients}
    canonical = json.dumps(certificate, sort_keys=True, separators=(',', ':')).encode()
    result = {'theorem_status': 'PROVED',
              'range': 'integers n >= 10, r=5, 0 <= m <= floor((3n-11)/5)',
              'construction_orbits': len(SUPPORT),
              'prefix_classes': len(PREFIXES), 'distinct_prefix_profiles': 44,
              'maximal_prefix_profiles': 4, 'polynomial_rows': len(coefficients),
              'nonnegative_newton_coefficients_checked': count,
              'certificate_sha256': sha256(canonical).hexdigest(),
              'large_distance_bound': large,
              'finite_checks': [finite_check(n, actual_profiles) for n in args.n],
              'controls': negative_controls(actual_profiles)}
    if args.write_results:
        (ROOT/'proof_certificate.json').write_text(json.dumps(certificate, indent=2)+'\n')
        (ROOT/'verification_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    try:
        main()
    except AssertionError as exc:
        print('FAIL:', exc, file=sys.stderr)
        sys.exit(1)
