#!/usr/bin/env python3
"""Replay the exact staircase seeds, saturation caps, and guard arithmetic.

Run: python3 -I certificates/infinite_families/verify.py verify
All inputs resolve beside this file. Temporary builds are removed on exit.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from fractions import Fraction as F
from math import comb

sys.dont_write_bytecode = True
BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

import numpy as np
from exact_tools import digits, kernel_count, orbit_size, point_orbits, prefixes, suffix_counts, volume
from orbit_checker import Cell
from packing_tools import exact_loads, parse_terms
from overlap_check import dense_minimum, literal_calibration

SEEDS = [(19, 4, 13), (19, 5, 12), (20, 6, 12), (19, 11, 10),
         (18, 7, 10), (19, 7, 11), (21, 6, 13)]


def read(name):
    return json.loads((BASE / name).read_text())


def check_hashes():
    count = 0
    for line in (BASE / 'SHA256SUMS').read_text().splitlines():
        expected, name = line.split('  ', 1)
        path = BASE / name
        if path.parent != BASE or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise AssertionError('checksum mismatch: ' + name)
        count += 1
    print('SHA256_VERIFIED', count, flush=True)


def run(command, cwd):
    result = subprocess.run([str(x) for x in command], cwd=cwd,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise AssertionError('command failed: ' + ' '.join(map(str, command)) + '\n' + result.stdout)
    return result.stdout


def triples(path):
    return {tuple(map(int, fields[1:4])) for line in path.read_text().splitlines()
            if (fields := line.split()) and fields[0] == '3'}


def kernel_check(N, r, m, work):
    """Compare every native entry with the independent support-overlap formula."""
    L, C = r-1, N-r+1
    pref = np.array(prefixes(L), dtype=np.int64)
    E = len(pref)*(C+1)
    ids = work / f'classes_{N}_{r}_{m}.txt'
    ids.write_text(''.join(str(i)+'\n' for i in range(E)))
    matrix_path = work / f'kernel_{N}_{r}_{m}.bin'
    run([work/'dump_selected', N, r, m, ids, matrix_path], work)
    with matrix_path.open('rb') as handle:
        header = tuple(map(int, np.fromfile(handle, np.int32, 5)))
        assert header[:4] == (N, r, m, E)
        P = header[4]
        actual_ids = np.fromfile(handle, np.int32, E)
        pp = np.fromfile(handle, np.int32, P)
        ww = np.fromfile(handle, np.int32, P)
        rhs = np.fromfile(handle, np.int64, P)
        matrix = np.fromfile(handle, np.int64, E*P).reshape(E, P)
        assert not handle.read(1)
    assert np.array_equal(actual_ids, np.arange(E))
    points = list(point_orbits(N, r, m))
    assert len(points) == P
    for i, (p, w) in enumerate(points):
        assert tuple(pref[pp[i]]) == p and ww[i] == w
        assert int(rhs[i]) == 2**L*comb(C, w)*3**w
    # Every entry, including a sum of r copies, fits in signed int64.
    assert r*2**L*4**C < 2**63
    table = np.zeros((C+1, C+1, L+1, L+1), dtype=np.int64)
    for a in range(C+1):
        for w in range(C+1):
            counts = suffix_counts(C, a)[w]
            for forced in range(L+1):
                for both in range(L-forced+1):
                    value = sum(comb(both, j)*sum(counts[:max(0, min(C, m-forced-j)+1)])
                                for j in range(both+1))
                    table[a, w, forced, both] = value*2**(L-both)
    target = pref[pp]
    for i, p in enumerate(pref):
        forced = np.count_nonzero(target != p, axis=1)
        both = np.count_nonzero((target == 2) & (p == 2), axis=1)
        for a in range(C+1):
            assert np.array_equal(table[a, ww, forced, both], matrix[i*(C+1)+a])
    print('KERNEL_VERIFIED', [N, r, m], 'entries', E*P, flush=True)
    return {'cell': [N, r, m], 'classes': E, 'target_orbits': P, 'entries': E*P}


def packing_check(filename, beta=None, cap=None, lower=None):
    data = read(filename)
    N, r, m = (data[key] for key in ('N', 'r', 'm'))
    C = N-r+1
    terms = parse_terms(data['packing'])
    v = sum((term[2] for term in terms), F(0))
    scale, values = exact_loads(N, r, m, terms)
    maximum = F(max(values), scale)
    assert maximum <= 1
    # A second kernel uses coordinate-polynomial multiplication.
    cell = Cell(N, r, m)
    independent, _ = cell.sweep(cell.parse(data['packing'], packing=True), packing=True)
    assert independent['passes'] and F(independent['extremum']) == maximum
    result = {'cell': [N, r, m], 'mass': str(v), 'max_load': str(maximum),
              'classes': len(values), 'support': len(terms)}
    if beta is None:
        assert v > (r-1 if filename == 'packing_18_7_10.json' else r-2)
    else:
        adjusted = [F(value, scale)+beta*(i % (C+1)) for i, value in enumerate(values)]
        actual_cap = max(adjusted)
        assert actual_cap <= cap and v >= lower
        margin = lower+3*C*beta-(r-1)*cap
        assert margin > 0
        result.update(beta=str(beta), adjusted_max=str(actual_cap),
                      adjusted_cap=str(cap), mass_lower_bound=str(lower),
                      saturation_weight=3*C, margin_lower_bound=str(margin))
        if (N, r, m) == (21, 6, 13):
            assert actual_cap > cap-F(1, 50000)
            print('NEGATIVE_CONTROL_REJECTED adjusted cap reduced by 1/50000', flush=True)
    assert maximum*2 > 1
    print('PACKING_VERIFIED', json.dumps(result), flush=True)
    return result, values, scale


def cover_check(filename):
    data = read(filename)
    N, r, m = (data[key] for key in ('N', 'r', 'm'))
    cell = Cell(N, r, m)
    result, _ = cell.sweep(cell.parse(data['cover']), packing=False)
    assert result['passes'] and F(result['total_mass']) < r-1
    assert F(result['extremum'])*F(9, 10) < 1
    print('FRACTIONAL_COVER_VERIFIED', [N, r, m], json.dumps(result), flush=True)
    return {'cell': [N, r, m], **result}


def overlap_check(work):
    fresh = work / 'triples.txt'
    output = run([work/'enum_triples', 19, 4, 13, 3, fresh], work)
    assert 'SOLUTIONS 3325' in output
    assert 'done 1=0 2=0 3=3325' in output
    accepted = triples(fresh)
    assert len(accepted) == 3325 and accepted == triples(BASE/'class_triples_19_4_13.txt')
    classes = [(p, w) for p in prefixes(3) for w in range(17)]
    points = set(point_orbits(19, 4, 13))
    seen = set()
    smallest = None
    for row in csv.DictReader((BASE/'overlap_witnesses.tsv').open(), delimiter='\t'):
        triple = tuple(int(row[x]) for x in ('a', 'b', 'c'))
        assert triple in accepted and triple not in seen
        seen.add(triple)
        point = (digits(int(row['point_prefix']), 3), int(row['weight']))
        assert point in points
        i, j = int(row['pair_i']), int(row['pair_j'])
        assert 0 <= i < j < 3
        covered = sum(kernel_count(16, 13, classes[v], point) for v in triple)
        lower = dense_minimum(16, 13, classes[triple[i]], classes[triple[j]], point)
        size = orbit_size(16, point)
        gap = size-covered+lower
        assert (covered, lower, size, gap) == tuple(int(row[x]) for x in
                                                   ('sum', 'intersection', 'orbit_size', 'gap'))
        assert gap > 0
        smallest = gap if smallest is None else min(smallest, gap)
    assert seen == accepted
    result = {'triples_checked': len(seen), 'minimum_gap': smallest}
    print('OVERLAP_SEED_VERIFIED', json.dumps(result), flush=True)
    return result


def family_check():
    guards = []
    for r in range(2, 12):
        d = (r+1)//3
        assert volume(d, 1) > (r-1)*volume(d, 0)
        guards.append({'r': r, 'd': d, 'guard_words': volume(d, 1), 'slots': r-1})

    def certified(N, r, m):
        return any(2 <= r <= r0 and N >= N0 and 0 <= m <= m0+(N-N0)//((r+1)//3)
                   for N0, r0, m0 in SEEDS)

    for h in (0, 1, 2, 1000):
        for cell in [(19+h, 4, 13+h), (19+2*h, 5, 12+h), (20+2*h, 6, 12+h),
                     (21+2*h, 6, 13+h), (19+2*h, 7, 11+h)]:
            assert certified(*cell)
    assert not certified(19, 2, 19)
    assert not certified(25, 5, 19)
    assert all(not certified(19, r, 13) for r in range(17, 21))
    # Integer transport does not transport a strict fractional deficit.
    import threshold_control
    profiles, _ = threshold_control.structural_check()
    control = threshold_control.finite_check(42, profiles)
    assert control['N'] == 41 and control['m'] == 23
    assert F(control['mass']) == 5 and F(control['max_load']) == 1
    assert 12+11 == control['m'] and 19+2*11 == control['N']
    print('FALSE_FRACTIONAL_RAY_CLAIM_REJECTED tau*(41,5,23)=5 at h=11', flush=True)
    print('GUARDS_VERIFIED', json.dumps(guards), flush=True)
    return {'seeds': SEEDS, 'guards': guards, 'fractional_ray_counterexample': control,
            'strict_fractional_gap_ray': 'K(19+h,4,13+h)=4 with tau*<3 for every h>=0'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['verify'])
    parser.parse_args()
    start = time.monotonic()
    check_hashes()
    compiler = shutil.which('clang++') or shutil.which('c++')
    assert compiler, 'C++17 compiler unavailable'
    with tempfile.TemporaryDirectory(prefix='staircase_families_') as directory:
        work = Path(directory)
        for name in ('enum_k', 'enum_triples', 'dump_selected'):
            run([compiler, '-O3', '-std=c++17', '-pthread', BASE/(name+'.cpp'), '-o', work/name], work)
        literal_tests = literal_calibration()
        print('LITERAL_PAIR_TESTS', literal_tests, 'mismatches 0', flush=True)
        kernels = [kernel_check(*cell, work) for cell in [(19, 4, 13), (19, 5, 12), (20, 6, 12)]]
        overlap = overlap_check(work)
        enumeration = []
        for cell, filename, expected in [((19, 5, 12), 'matched_19_5_12.json', 426),
                                         ((20, 6, 12), 'packing_20_6_12.json', 371)]:
            N, r, m = cell
            result, values, scale = packing_check(filename)
            v = F(result['mass'])
            eligible = [i for i, value in enumerate(values) if F(value, scale) >= v-(r-2)]
            assert len(eligible) == expected
            path = work/f'eligible_{N}_{r}_{m}.txt'
            path.write_text(''.join(str(i)+'\n' for i in eligible))
            output = run([work/'enum_k', N, r, m, r-1, path], work)
            assert 'solutions 0' in output and 'SOL ' not in output
            item = {'cell': cell, 'eligible': len(eligible), 'solutions': 0}
            enumeration.append(item)
            print('INTEGER_SEED_VERIFIED', json.dumps(item), flush=True)
        saturation = [packing_check('packing_18_7_10.json')[0],
                      packing_check('packing_19_7_11.json', F(23, 600), F(59061, 50000), F(2883, 500))[0],
                      packing_check('packing_21_6_13.json', F(31, 1000), F(60379, 50000), F(9101, 2000))[0]]
        covers = [cover_check(name) for name in ['cover_19_4_13.json', 'matched_19_5_12.json',
                                                'cover_20_6_12.json', 'cover_19_7_11.json',
                                                'matched_21_6_13.json']]
        families = family_check()
    print(json.dumps({'kernels': kernels, 'overlap': overlap, 'integer_seeds': enumeration,
                      'saturation': saturation, 'fractional_covers': covers, 'families': families,
                      'literal_pair_tests': literal_tests, 'seconds': time.monotonic()-start}, indent=2))
    print('INFINITE_FAMILIES_VERIFIED five integer rays; four-run strict fractional-gap ray', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (AssertionError, ValueError, ArithmeticError) as exc:
        print('FAIL:', exc, file=sys.stderr)
        sys.exit(1)
