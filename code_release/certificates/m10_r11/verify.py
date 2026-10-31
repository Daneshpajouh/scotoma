#!/usr/bin/env python3
"""Independent exact audit of the length-19 staircase covering certificate."""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'
import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from math import comb, gcd, lcm
from pathlib import Path
import sys
import time
import sys
sys.dont_write_bytecode = True

import numpy as np

sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
SOURCE = HERE
PACKING = HERE / 'cell_19_11_10.certificate.json'
L, C, M = 10, 9, 10
SCALE = 10**8


def say(*items):
    print(*items, flush=True)


def encode(prefix):
    value = 0
    for letter in prefix:
        value = 3 * value + (2 if letter in ('X', '2', 2) else int(letter))
    return value


def decode(code):
    out = [0] * L
    for i in range(L - 1, -1, -1):
        out[i], code = code % 3, code // 3
    return out


def near_stair(prefix):
    return min(sum(a != (0 if i < j else 1) for i, a in enumerate(prefix)) for j in range(L + 1))


def suffix_distribution(n, centre_weight, word_weight):
    """Count by support intersection s and equal nonzero letters h."""
    out = [0] * (n + 1)
    for s in range(max(0, word_weight + centre_weight - n), min(word_weight, centre_weight) + 1):
        supports = comb(centre_weight, s) * comb(n - centre_weight, word_weight - s)
        outside_labels = 3 ** (word_weight - s)
        for h in range(s + 1):
            distance = centre_weight + word_weight - s - h
            out[distance] += supports * outside_labels * comb(s, h) * 2 ** (s - h)
    assert sum(out) == comb(n, word_weight) * 3**word_weight
    return out


def make_table(m=M, n=C, length=L):
    probabilities = {}
    denominator = 1
    for word_w in range(n + 1):
        for centre_w in range(n + 1):
            hist = suffix_distribution(n, centre_w, word_w)
            for forced in range(length + 1):
                for shared_x in range(length + 1 - forced):
                    count = sum(comb(shared_x, t) * sum(hist[:max(0, min(n + 1, m - forced - t + 1))])
                                for t in range(shared_x + 1))
                    probability = Fraction(count, 2**shared_x * comb(n, word_w) * 3**word_w)
                    probabilities[word_w, centre_w, forced, shared_x] = probability
                    denominator = lcm(denominator, probability.denominator)
    table = np.zeros((n + 1, n + 1, length + 1, length + 1), dtype=np.int64)
    for index, q in probabilities.items():
        table[index] = q.numerator * (denominator // q.denominator)
    return table, denominator


def terms():
    certificate = json.loads(PACKING.read_text())
    assert (certificate['N'], certificate['r'], certificate['m']) == (19, 11, 10)
    assert certificate['normalization'] == 'total mass, uniform within each orbit'
    result = []
    seen = set()
    for row in certificate['packing']:
        assert len(row['prefix']) == L and set(row['prefix']) <= set('01X2')
        code, weight, mass = encode(row['prefix']), row['suffix_weight'], Fraction(row['mass'])
        assert type(weight) is int and 0 <= weight <= C and mass > 0
        assert (code, weight) not in seen
        seen.add((code, weight))
        assert near_stair(decode(code)) + weight <= M, ('packing orbit outside U', row)
        result.append((code, weight, mass))
    return result


def small_checks():
    # Exhaust actual four-letter words, independently of the support-intersection formula.
    checks = 0
    for n in range(5):
        for a in range(n + 1):
            centre = (1,) * a + (0,) * (n - a)
            counts = [[0] * (n + 1) for _ in range(n + 1)]
            for word in itertools.product(range(4), repeat=n):
                b = sum(x != 0 for x in word)
                d = sum(x != y for x, y in zip(word, centre))
                counts[b][d] += 1
            for b in range(n + 1):
                assert counts[b] == suffix_distribution(n, a, b)
                checks += 1
    # Include shared X coordinates as actual independent choices 2 and 3.
    kernel_checks = 0
    length, n = 3, 3
    for prefix in itertools.product(range(3), repeat=length):
        options = [(2, 3) if x == 2 else (x,) for x in prefix]
        for pw in range(n + 1):
            orbit = [p + s for p in itertools.product(*options)
                     for s in itertools.product(range(4), repeat=n) if sum(x != 0 for x in s) == pw]
            for cw in range(n + 1):
                for cp in ((0, 1, 2), (3, 2, 1), (2, 2, 2)):
                    cc = tuple(2 if x == 3 else x for x in cp)
                    f = sum(x != y for x, y in zip(cc, prefix))
                    q = sum(x == y == 2 for x, y in zip(cc, prefix))
                    suffix = tuple((1, 2, 3)[i % 3] if i >= n - cw else 0 for i in range(n))
                    centre = cp + suffix
                    distances = [sum(x != y for x, y in zip(word, centre)) for word in orbit]
                    hist = suffix_distribution(n, cw, pw)
                    for budget in range(length + n + 1):
                        direct = Fraction(sum(d <= budget for d in distances), len(orbit))
                        calculated = Fraction(sum(comb(q, t) * sum(hist[:max(0, min(n+1, budget-f-t+1))])
                                                  for t in range(q + 1)), 2**q * comb(n, pw) * 3**pw)
                        assert direct == calculated, (prefix, pw, cw, cp, budget)
                        kernel_checks += 1
    say('SMALL_EXHAUSTIVE_OK suffix_histograms', checks, 'kernel_probabilities', kernel_checks)
    return {'suffix_histograms': checks, 'kernel_probabilities': kernel_checks}


def check_dual(A, b, z, banned, subtracted, multiple, D, expected=None):
    assert len(z) == len(b) and all(type(q) is int and q >= 0 for q in z), 'negative or malformed dual'
    # A sufficient nonnegative overflow bound licenses fast signed-int64 arithmetic.
    bound = sum(z) * int(A.max())
    assert bound < 2**63, 'dual accumulation requires arbitrary precision'
    products = np.asarray(z, dtype=np.int64) @ A.astype(np.int64)
    allowed = [j for j in range(A.shape[1]) if j not in banned]
    worst_col = max(allowed, key=lambda j: int(products[j]))
    worst = int(products[worst_col])
    assert worst <= D, f'column cap violated at {worst_col}: {worst} > {D}'
    residual = [int(b[i]) - sum(int(A[i, j]) for j in subtracted) for i in range(len(b))]
    dot = sum(zr * br for zr, br in zip(z, residual))
    if expected is not None:
        assert dot == expected, f'stored dot mismatch: {dot} != {expected}'
    assert dot > multiple * D, f'insufficient residual dot {dot}'
    return {'worst_product': worst, 'worst_column': worst_col, 'signed_dot': dot,
            'strict_margin': dot - multiple * D, 'min_z': min(z), 'nnz': sum(q != 0 for q in z),
            'safe_int64_bound': bound, 'negative_residual_rows': sum(q < 0 for q in residual)}


def build():
    start = time.monotonic()
    report = {'small_exhaustive': small_checks()}
    packing = terms()
    mass = sum((q for _, _, q in packing), Fraction())
    threshold = mass - 9
    assert 9 < mass < 10 and 0 < threshold < 1
    table, D = make_table()
    say('PACKING', len(packing), 'all_orbits_in_U', True)
    say('PACKING_MASS', str(mass))
    say('THRESHOLD', str(threshold))
    say('KERNEL_DENOMINATOR', D)
    reported = json.loads((SOURCE / 'threshold.json').read_text())
    assert Fraction(reported['packing_mass']) == mass and Fraction(reported['threshold']) == threshold
    assert reported['dprob'] == D
    # Read data only. No external source script is loaded or executed.
    tokens = (SOURCE / 'filter_instance.txt').read_text().split()
    cursor = 0
    def take():
        nonlocal cursor
        q = int(tokens[cursor]); cursor += 1
        return q
    assert [take() for _ in range(5)] == [L, C, M, 10**18, D]
    orig_scale = 10**18
    scaled_threshold = threshold * orig_scale * D
    assert [take(), take()] == [-(-scaled_threshold.numerator // scaled_threshold.denominator),
                               scaled_threshold.numerator // scaled_threshold.denominator]
    assert take() == len(packing)
    for code, w, q in packing:
        assert [take(), take(), take(), take()] == [code, w, (q * orig_scale).__floor__(), (q * orig_scale).__ceil__()]
    covers = take()
    stored_cover_ids = [(take(), take()) for _ in range(covers)]
    assert take() == table.size
    stored_table = np.asarray([take() for _ in range(table.size)], dtype=np.int64).reshape(table.shape)
    assert cursor == len(tokens)
    mismatches = int(np.count_nonzero(stored_table != table))
    assert mismatches == 0
    say('FILTER_DATA_OK table_entries', table.size, 'table_mismatches', mismatches)

    allprefix = np.asarray([decode(i) for i in range(3**L)], dtype=np.uint8)
    ppre = np.asarray([decode(code) for code, _, _ in packing], dtype=np.uint8)
    pw = np.asarray([w for _, w, _ in packing], dtype=np.intp)
    forced = np.zeros((len(packing), len(allprefix)), dtype=np.uint8)
    shared = np.zeros_like(forced)
    for coordinate in range(L):
        a, b = ppre[:, coordinate, None], allprefix[None, :, coordinate]
        forced += a != b
        shared += (a == 2) & (b == 2)
    low = np.asarray([(q * SCALE).__floor__() for _, _, q in packing], dtype=np.int64)
    high = np.asarray([(q * SCALE).__ceil__() for _, _, q in packing], dtype=np.int64)
    cap = SCALE * D
    assert sum(int(q) for q in high) * D < 2**63
    th_num, th_den = threshold.numerator, threshold.denominator
    th_scaled = threshold * cap
    th_floor, th_ceil = th_scaled.__floor__(), th_scaled.__ceil__()
    chosen = []
    interval_eligible = interval_ineligible = elig_uncertain = load_uncertain = 0
    exact_eligible = exact_ineligible = exact_tight = 0
    exact_tight_ids = []
    max_excluded_upper = -1
    min_included_lower = cap * 10
    nearest_excluded = nearest_included = None
    max_strict_upper = -1
    exact_settlements = []
    for cw in range(C + 1):
        K = table[pw[:, None], cw, forced, shared]
        lower, upper = low @ K, high @ K
        eq_or_excess = np.flatnonzero(upper > cap)
        load_uncertain += len(eq_or_excess)
        exact_loads = {}
        for code0 in eq_or_excess:
            code = int(code0)
            load = sum((q * int(K[i, code]) for i, (_, _, q) in enumerate(packing)), Fraction()) / D
            exact_loads[code] = load
            assert load <= 1, ('packing infeasible', code, cw, str(load))
            if load == 1:
                exact_tight += 1
                exact_tight_ids.append([code, cw])
            exact_settlements.append({'prefix': code, 'weight': cw, 'load': str(load), 'reason': 'load_cap'})
        strict = upper[upper <= cap]
        if strict.size:
            max_strict_upper = max(max_strict_upper, int(strict.max()))
        accept = lower >= th_ceil
        reject = upper <= th_floor
        # Strict exclusion is essential. Threshold is nonintegral at this scale.
        assert th_floor < th_scaled < th_ceil
        unresolved = ~(accept | reject)
        interval_eligible += int(accept.sum())
        interval_ineligible += int(reject.sum())
        elig_uncertain += int(unresolved.sum())
        for code0 in np.flatnonzero(unresolved):
            code = int(code0)
            if code not in exact_loads:
                exact_loads[code] = sum((q * int(K[i, code]) for i, (_, _, q) in enumerate(packing)), Fraction()) / D
            load = exact_loads[code]
            accept[code] = load >= threshold
            reject[code] = load < threshold
            exact_eligible += int(accept[code])
            exact_ineligible += int(reject[code])
            exact_settlements.append({'prefix': code, 'weight': cw, 'load': str(load), 'reason': 'eligibility'})
        for code0 in np.flatnonzero(accept):
            chosen.append((int(code0), cw))
        if np.any(reject):
            idx = np.flatnonzero(reject)
            code = int(idx[np.argmax(upper[idx])])
            if int(upper[code]) > max_excluded_upper:
                max_excluded_upper = int(upper[code]); nearest_excluded = (code, cw)
        if np.any(accept):
            idx = np.flatnonzero(accept)
            code = int(idx[np.argmin(lower[idx])])
            if int(lower[code]) < min_included_lower:
                min_included_lower = int(lower[code]); nearest_included = (code, cw)
        say('SCAN suffix_weight', cw, 'eligible', int(accept.sum()), 'ineligible', int(reject.sum()),
            'threshold_exact_fallback', int(unresolved.sum()), 'load_exact_fallback', len(eq_or_excess))
    chosen.sort()
    independent_index = np.asarray(chosen, dtype=np.int32)
    reported_index = np.load(SOURCE / 'elig_index.npy', allow_pickle=False)
    filter_ids = [tuple(map(int, line.split()[:2])) for line in (SOURCE / 'filter_out_candidates.txt').read_text().splitlines() if line.strip()]
    assert len(filter_ids) == len(set(filter_ids)) == len(chosen)
    assert set(chosen) == set(filter_ids), 'eligible set mismatch'
    assert np.array_equal(reported_index, np.asarray(filter_ids)), 'supplied column order mismatch'
    independent_position = {pair: j for j, pair in enumerate(chosen)}
    source_to_independent = np.asarray([independent_position[pair] for pair in filter_ids], dtype=np.int32)
    say('FULL_SCAN_OK', (3**L) * (C + 1), 'eligible', len(chosen), 'ineligible', 3**L * (C + 1) - len(chosen),
        'threshold_fallbacks', elig_uncertain, 'load_fallbacks', load_uncertain, 'tight', exact_tight)
    assert len(set(chosen)) == len(chosen)
    assert set(stored_cover_ids).issubset(set(chosen))
    assert set(stored_cover_ids) == set(map(tuple, exact_tight_ids)), 'tight load orbit set differs from supplied cover support'

    indices = independent_index[:, 0]
    cweights = independent_index[:, 1]
    K = table[pw[:, None], cweights[None, :], forced[:, indices], shared[:, indices]]
    oldK = np.load(SOURCE / 'pack_matrix.npy', allow_pickle=False)
    assert oldK.shape == K.shape and np.array_equal(oldK, K[:, source_to_independent]), 'kernel matrix mismatch'
    row_gcd = [gcd(D, *(int(v) for v in row)) for row in K]
    A = (K // np.asarray(row_gcd)[:, None]).astype(np.int32)
    b = np.asarray([D // g for g in row_gcd], dtype=np.int32)
    assert np.array_equal(A[:, source_to_independent], np.load(SOURCE / 'A_reduced.npy', allow_pickle=False)), 'reduced matrix mismatch'
    assert np.array_equal(b, np.load(SOURCE / 'rhs_reduced.npy', allow_pickle=False)), 'right hand side mismatch'
    say('MATRIX_OK rows', A.shape[0], 'columns', A.shape[1], 'exact_entries', A.size, 'mismatches', 0)
    # Canonical order and all five case definitions are fixed here, independently of the metadata.
    idx = {pair: j for j, pair in enumerate(chosen)}
    endpoints = [idx[(encode('0000000000'), 0)], idx[(encode('1111111111'), 0)]]
    neighbours = [idx[(encode('0000000001'), 0)], idx[(encode('0111111111'), 0)]]
    cases = [
        ('omit_both_ends', endpoints, [], 10),
        ('end0_omit_1513', [neighbours[0]], [endpoints[0]], 9),
        ('end0_and_1513', [], [endpoints[0], neighbours[0]], 8),
        ('end7074_omit_2065', [neighbours[1]], [endpoints[1]], 9),
        ('end7074_and_2065', [], [endpoints[1], neighbours[1]], 8),
    ]
    meta = json.loads((SOURCE / 'cert_index.json').read_text())
    assert meta['D'] == 10**12
    assert [int(source_to_independent[j]) for j in meta['columns'].values()] == [endpoints[0], endpoints[1], neighbours[0], neighbours[1]]
    assert len(meta['certs']) == len(cases)
    certificates = []
    for (name, banned, subtracted, multiple), stored in zip(cases, meta['certs']):
        assert (name, banned, subtracted, multiple) == (stored['name'], [int(source_to_independent[j]) for j in stored['banned']],
                                                       [int(source_to_independent[j]) for j in stored['subtracted']], stored['multiple'])
        zz = np.load(SOURCE / ('cert_' + name + '_z.npy'), allow_pickle=False)
        assert zz.dtype.kind in 'iu' and zz.ndim == 1
        z = [int(q) for q in zz]
        checked = check_dual(A, b, z, set(banned), subtracted, multiple, meta['D'], int(stored['signed_dot']))
        assert checked['nnz'] == stored['z_nnz']
        certificates.append({'name': name, 'banned': banned, 'subtracted': subtracted, 'multiple': multiple, **checked})
        say('DUAL_OK', name, 'worst', checked['worst_product'], 'dot', checked['signed_dot'], 'margin', checked['strict_margin'])
    # Enumerate the truncated integer case split, allowing multiplicity > 1.
    branch_counts = [0] * 5
    checked_tuples = 0
    for e0, e1, n0, n1 in itertools.product(range(11), repeat=4):
        if e0 + e1 + n0 + n1 > 10:
            continue
        memberships = [e0 == 0 and e1 == 0,
                       e0 >= 1 and n0 == 0,
                       e0 >= 1 and n0 >= 1,
                       e0 == 0 and e1 >= 1 and n1 == 0,
                       e0 == 0 and e1 >= 1 and n1 >= 1]
        assert sum(memberships) == 1
        branch_counts[memberships.index(True)] += 1
        checked_tuples += 1
    say('CASE_SPLIT_OK integer_tuples', checked_tuples, 'branch_counts', branch_counts)
    report.update({'packing_terms': len(packing), 'packing_mass': str(mass), 'threshold': str(threshold),
                   'DPROB': D, 'fixed_point_scale': SCALE, 'fixed_point_int64_bound': sum(int(q) for q in high) * D,
                   'classes_scanned': 3**L * (C + 1), 'eligible': len(chosen),
                   'ineligible': 3**L * (C + 1) - len(chosen),
                   'interval_eligible': interval_eligible, 'interval_ineligible': interval_ineligible,
                   'eligibility_exact_fallbacks': elig_uncertain, 'load_exact_fallbacks': load_uncertain,
                   'exact_load_one': exact_tight, 'exact_fallback_eligible': exact_eligible,
                   'exact_fallback_ineligible': exact_ineligible,
                   'excluded_gap_lower_bound': str(threshold - Fraction(max_excluded_upper, cap)),
                   'included_gap_lower_bound': str(Fraction(min_included_lower, cap) - threshold),
                   'nearest_excluded_by_interval': nearest_excluded, 'nearest_included_by_interval': nearest_included,
                   'max_strict_load_upper_bound': str(Fraction(max_strict_upper, cap)),
                   'matrix_shape': list(A.shape), 'matrix_entries': A.size, 'kernel_mismatches': 0,
                   'certificates': certificates, 'case_split_tuples': checked_tuples,
                   'case_split_branch_counts': branch_counts,
                   'elapsed_seconds': time.monotonic() - start})
    verify_fractional(table, D, mass)
    wrong, wrongD = make_table(m=9)
    unequal = np.argwhere(wrong.astype(object) * D != table.astype(object) * wrongD)
    try:
        assert len(unequal) == 0, 'wrong-budget kernel probability mismatch'
    except AssertionError:
        say('WRONG_BUDGET_REJECTED changed_probabilities', len(unequal))
    else:
        raise AssertionError('wrong-budget control unexpectedly passed')
    z = [int(q) for q in np.load(SOURCE / 'cert_omit_both_ends_z.npy', allow_pickle=False)]
    row = next(i for i, q in enumerate(z) if q > 0)
    z[row] += 10**9
    try:
        check_dual(A, b, z, set(endpoints), [], 10, 10**12)
    except AssertionError as exc:
        assert 'column cap violated' in str(exc), str(exc)
        say('INFLATED_DUAL_REJECTED', str(exc))
    else:
        raise AssertionError('inflated-dual control unexpectedly passed')
    say('M10_R11_VERIFIED K(19,11,10)=11 ceil(tau_star)=10 seconds', round(time.monotonic() - start, 3))


def verify_fractional(table, D, packing_mass):
    data = json.loads(PACKING.read_text())
    cover = []
    seen = set()
    for row in data['cover']:
        code, w, mass = encode(row['prefix']), row['suffix_weight'], Fraction(row['mass'])
        assert len(row['prefix']) == L and set(row['prefix']) <= set('01X2')
        assert type(w) is int and 0 <= w <= C and mass > 0
        assert (code, w) not in seen
        seen.add((code, w))
        cover.append((code, w, mass))
    v = sum((q for _, _, q in cover), Fraction())
    assert v == packing_mass and 9 < v < 10
    cap = SCALE * D
    allprefix = np.asarray([decode(code) for code in range(3**L)], dtype=np.uint8)
    cprefix = np.asarray([decode(code) for code, _, _ in cover], dtype=np.uint8)
    cw = np.asarray([w for _, w, _ in cover], dtype=np.intp)
    forced = np.zeros((len(cover), len(allprefix)), dtype=np.uint8)
    shared = np.zeros_like(forced)
    for i in range(L):
        a, b = cprefix[:, i, None], allprefix[None, :, i]
        forced += a != b
        shared += (a == 2) & (b == 2)
    point_distance = np.asarray([near_stair(prefix) for prefix in allprefix], dtype=np.int8)
    low = np.asarray([(q*SCALE).__floor__() for _, _, q in cover], dtype=np.int64)
    high = np.asarray([(q*SCALE).__ceil__() for _, _, q in cover], dtype=np.int64)
    assert sum(int(q) for q in high) * D < 2**63
    rows = fallbacks = tight = 0
    settled = []
    for pw in range(C+1):
        ids = np.flatnonzero(point_distance + pw <= M)
        K = table[pw, cw[:, None], forced[:, ids], shared[:, ids]]
        lower = low @ K
        uncertain = np.flatnonzero(lower < cap)
        for j0 in uncertain:
            j = int(j0)
            value = sum((q*int(K[i,j]) for i, (_, _, q) in enumerate(cover)), Fraction()) / D
            assert value >= 1, ('fractional cover fails', int(ids[j]), pw, str(value))
            tight += value == 1
            settled.append({'prefix': int(ids[j]), 'weight': pw, 'coverage': str(value)})
        rows += len(ids)
        fallbacks += len(uncertain)
        print('FRACTIONAL_COVER suffix_weight', pw, 'U_orbits', len(ids), 'exact_fallbacks', len(uncertain), flush=True)
    result = {'label': 'VERIFIED-COMPUTATIONALLY', 'cover_classes': len(cover), 'U_orbits': rows,
              'exact_fallbacks': fallbacks, 'tight_coverage_orbits': tight,
              'mass': str(v), 'strictly_below_ten': v < 10, 'matched_packing_mass': True}
    print('FRACTIONAL_UPPER_OK', json.dumps(result), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['verify'], nargs='?', default='verify')
    parser.parse_args()
    build()
