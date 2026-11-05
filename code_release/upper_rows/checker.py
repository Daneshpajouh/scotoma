#!/usr/bin/env python3
"""Independent, standard-library-only orbit coverage and load checker.

Masses are TOTAL orbit masses, spread uniformly. Prefix X (or 2) denotes
either unused letter. No privacy or prescribed mass condition is imposed.
All mathematical decisions use fractions.Fraction and arbitrary-size integers.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, lcm
from pathlib import Path
import argparse
import json
import time


class Cell:
    def __init__(self, N, r, m):
        if not all(type(v) is int for v in (N, r, m)):
            raise ValueError('N, r, m must be integers')
        if not (0 <= m <= N and 1 <= r <= N + 1):
            raise ValueError('invalid cell')
        self.N, self.r, self.m = N, r, m
        self.L, self.C = r - 1, N - r + 1
        self.prefixes = tuple(product(range(3), repeat=self.L))
        self.library = tuple((0,) * j + (1,) * (self.L-j) for j in range(r))
        self.dmin = {a: min(sum(x != y for x, y in zip(a, b))
                            for b in self.library) for a in self.prefixes}
        self.centres = tuple((a, w) for a in self.prefixes for w in range(self.C+1))
        self.points = tuple((a, w) for a, w in self.centres if self.in_union((a, w)))
        self._suffix_counts = {}
        self._prefix_signatures = {}

    def in_union(self, orbit):
        a, w = orbit
        return self.dmin[a] + w <= self.m

    def size(self, orbit):
        a, w = orbit
        return 2**a.count(2) * comb(self.C, w) * 3**w

    def suffix_counts(self, p):
        """Coefficients [z^w t^d] of
        (1+3*z*t)^(C-p) * (t+z+2*z*t)^p.
        Constructed by coordinate-wise polynomial multiplication, no overlap
        counting formula or numerical special function.
        """
        if p in self._suffix_counts:
            return self._suffix_counts[p]
        poly = {(0, 0): 1}
        for coord in range(self.C):
            terms = ((0, 1, 1), (1, 0, 1), (1, 1, 2)) if coord < p else (
                (0, 0, 1), (1, 1, 3))
            nxt = {}
            for (w, d), n in poly.items():
                for dw, dd, multiplicity in terms:
                    key = (w+dw, d+dd)
                    nxt[key] = nxt.get(key, 0) + n*multiplicity
            poly = nxt
        counts = tuple(tuple(poly.get((w, d), 0) for d in range(self.C+1))
                       for w in range(self.C+1))
        for w in range(self.C+1):
            if sum(counts[w]) != comb(self.C, w)*3**w:
                raise ArithmeticError('suffix orbit count mismatch')
        self._suffix_counts[p] = counts
        return counts

    def signature(self, a, b):
        key = (a, b)
        if key not in self._prefix_signatures:
            both = sum(x == y == 2 for x, y in zip(a, b))
            forced = sum(x != y for x, y in zip(a, b))
            self._prefix_signatures[key] = forced, both
        return self._prefix_signatures[key]

    @lru_cache(maxsize=None)
    def probability(self, moving_weight, fixed_weight, forced, both):
        """Fraction of the moving orbit inside a ball around a fixed word.

        At a both-X coordinate, match and mismatch each have probability 1/2.
        forced counts every deterministically mismatching active coordinate.
        """
        budget = self.m - forced
        if budget < 0:
            return F(0)
        counts = self.suffix_counts(fixed_weight)[moving_weight]
        numerator = 0
        for mismatch in range(both+1):
            end = min(self.C, budget-mismatch)
            if end >= 0:
                numerator += comb(both, mismatch) * sum(counts[:end+1])
        denominator = 2**both * comb(self.C, moving_weight) * 3**moving_weight
        return F(numerator, denominator)

    def kernel(self, fixed, moving):
        a, p = fixed
        b, w = moving
        forced, both = self.signature(a, b)
        return self.probability(w, p, forced, both)

    def parse(self, terms, packing=False):
        merged = {}
        for term in terms:
            raw = term['prefix']
            a = tuple(2 if x in ('X', '2', 2) else int(x) for x in raw)
            w = term['suffix_weight']
            mass = term['mass']
            if isinstance(mass, float):
                raise ValueError('floating-point masses are forbidden')
            mass = F(mass)
            if len(a) != self.L or any(x not in (0, 1, 2) for x in a):
                raise ValueError('bad active prefix')
            if type(w) is not int or not 0 <= w <= self.C or mass < 0:
                raise ValueError('bad suffix weight or negative mass')
            orbit = (a, w)
            if packing and not self.in_union(orbit):
                raise ValueError('packing mass outside U')
            merged[orbit] = merged.get(orbit, F(0)) + mass
        return tuple((a, mass) for a, mass in merged.items() if mass)

    def sweep(self, weights, packing=False, collect=False, method='auto'):
        orbits = self.centres if packing else self.points
        if method == 'scaled' or (method == 'auto' and len(orbits)*len(weights) > 500000):
            return self._sweep_scaled(weights, packing, collect)
        return self._sweep_direct(weights, packing, collect)

    def _sweep_direct(self, weights, packing=False, collect=False):
        """Enumerate EVERY relevant orbit, never sample or early accept."""
        orbits = self.centres if packing else self.points
        worst = None
        at = None
        violations = []
        nbad = 0
        for orbit in orbits:
            value = sum((mass*self.kernel(orbit, support)
                         for support, mass in weights), F(0))
            if worst is None or (value > worst if packing else value < worst):
                worst, at = value, orbit
            bad = value > 1 if packing else value < 1
            if bad:
                nbad += 1
                if collect:
                    violations.append((value, orbit))
        return {'passes': nbad == 0, 'extremum': str(worst),
                'at': orbit_json(at), 'violations': nbad,
                'orbits_checked': len(orbits),
                'total_mass': str(sum((mass for _, mass in weights), F(0))),
                'support_orbits': len(weights)}, violations

    def _sweep_scaled(self, weights, packing=False, collect=False):
        """Same Fraction kernel, with exact common-denominator accumulation.

        Every weighted kernel denominator divides mass.denominator *
        2^L * binomial(C,w) * 3^w. Taking their LCM therefore gives an
        exact integer accumulation scale. No bounded machine integer is used.
        Every extremum and feasibility decision still compares Fractions.
        """
        orbits = self.centres if packing else self.points
        scale = lcm(*(mass.denominator * 2**self.L * comb(self.C, w) * 3**w
                      for (_, w), mass in weights)) if weights else 1
        width = self.L+1
        tables = []
        for (_, w), mass in weights:
            layers = []
            for p in range(self.C+1):
                table = [0]*(width*width)
                for forced in range(width):
                    for both in range(width-forced):
                        probability = self.probability(w, p, forced, both)
                        denominator = mass.denominator*probability.denominator
                        quotient, remainder = divmod(scale, denominator)
                        if remainder:
                            raise ArithmeticError('common denominator was not exact')
                        table[forced*width+both] = mass.numerator*probability.numerator*quotient
                layers.append(table)
            tables.append(layers)
        worst, at, last_prefix = None, None, None
        signatures = []
        violations = []
        nbad = 0
        for orbit in orbits:
            a, p = orbit
            if a != last_prefix:
                signatures = []
                for (b, _), _mass in weights:
                    forced = sum(x != y for x, y in zip(a, b))
                    both = sum(x == y == 2 for x, y in zip(a, b))
                    signatures.append(forced*width+both)
                last_prefix = a
            integer_total = sum(table[p][signature] for table, signature in zip(tables, signatures))
            value = F(integer_total, scale)
            if worst is None or (value > worst if packing else value < worst):
                worst, at = value, orbit
            bad = value > 1 if packing else value < 1
            if bad:
                nbad += 1
                if collect:
                    violations.append((value, orbit))
        return {'passes': nbad == 0, 'extremum': str(worst),
                'at': orbit_json(at), 'violations': nbad,
                'orbits_checked': len(orbits),
                'total_mass': str(sum((mass for _, mass in weights), F(0))),
                'support_orbits': len(weights)}, violations


def orbit_json(orbit):
    a, w = orbit
    return {'prefix': ''.join('X' if x == 2 else str(x) for x in a),
            'suffix_weight': w}


def terms_json(weights):
    return [dict(orbit_json(orbit), mass=str(mass)) for orbit, mass in weights]


def verify(cert):
    start = time.monotonic()
    cell = Cell(cert['N'], cert['r'], cert['m'])
    answer = {k: cert[k] for k in ('N', 'r', 'm')}
    answer['status'] = 'VERIFIED-COMPUTATIONALLY'
    answer['total_centre_orbits'] = len(cell.centres)
    answer['total_point_orbits'] = len(cell.points)
    answer['centre_word_count'] = sum(map(cell.size, cell.centres))
    answer['union_word_count'] = sum(map(cell.size, cell.points))
    if answer['centre_word_count'] != 4**cell.N:
        raise ArithmeticError('orbit enumeration omitted words')
    for key, packing in (('cover', False), ('packing', True)):
        if key in cert:
            weights = cell.parse(cert[key], packing)
            answer[key], _ = cell.sweep(weights, packing)
    answer['passes'] = all(answer[key]['passes'] for key in ('cover', 'packing')
                           if key in answer)
    if 'cover' in answer and 'packing' in answer:
        upper = F(answer['cover']['total_mass'])
        lower = F(answer['packing']['total_mass'])
        answer['matched_exactly'] = answer['passes'] and upper == lower
        if answer['matched_exactly']:
            answer['tau_star'] = str(upper)
            answer['deficit'] = str(cell.r-upper)
            # K <= r by the explicit library cover. Only assert equality
            # when the exact packing lower bound forces it.
            if cell.r-1 < lower <= cell.r:
                answer['K'] = cell.r
                answer['K_reason'] = 'ceil(exact packing mass) = r; library cover has r balls'
    answer['elapsed_seconds'] = round(time.monotonic()-start, 3)
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        result = verify(json.loads(args.certificate.read_text()))
    except (ValueError, KeyError, ArithmeticError, TypeError) as exc:
        print(json.dumps({'passes': False, 'error': str(exc)}))
        return 1
    output = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(output)
    print(output, end='')
    return 0 if result['passes'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
