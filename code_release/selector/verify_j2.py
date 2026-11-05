#!/usr/bin/env python3
"""Exact, standalone verification of the J2 counterexamples.

No repository imports; no files are written by this program. All comparisons use
integers or fractions.Fraction. A DP exhaustively covers centre classes over
{0,1,X}; exchanging 2 and 3 independently at each coordinate preserves loads.
"""
import argparse
import itertools
import json
import math
import sys
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache


def choose(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def allowed(N, j):
    return tuple(p for p in range(N) if p not in (j - 1, j))


def sphere_size(N, m, j):
    return choose(len(allowed(N, j)), m) * 2**m


def pattern(N, j):
    return (1,) * j + (0,) * (N - j)


def sphere(N, m, j):
    for positions in itertools.combinations(allowed(N, j), m):
        for letters in itertools.product((2, 3), repeat=m):
            x = list(pattern(N, j))
            for p, letter in zip(positions, letters):
                x[p] = letter
            yield tuple(x)


def distance(a, b):
    return sum(x != y for x, y in zip(a, b))


@lru_cache(None)
def count_profile(M, m, u, e, a, binary_out, double_wild=False):
    """Count accepted sphere words, not their probability.

    u: global number of free coordinates in the centre;
    e: global number of binary mismatches against this pattern;
    a: free coordinates inside the allowed set;
    binary_out: binary mismatches outside the allowed set.
    For selected A free and B binary-mismatch coordinates, acceptance is
    s <= 2*A+B-u-e. There are C(A,s)*2**(m-A) letter assignments.
    """
    b = e - binary_out
    same = M - a - b
    if min(a, b, same) < 0:
        return 0
    total = 0
    # Negative control: count excluded free coordinates a second time.
    extra = u - a if double_wild else 0
    for A in range(min(a, m) + 1):
        for B in range(min(b, m - A) + 1):
            S = m - A - B
            cap = 2 * A + B - u - e - extra
            if S > same or cap < 0:
                continue
            letters = sum(choose(A, s) for s in range(min(A, cap) + 1))
            total += choose(a, A) * choose(b, B) * choose(same, S) * letters * 2**(m - A)
    return total


def counts(word, m, double_wild=False):
    N = len(word)
    u = sum(c >= 2 for c in word)
    result = []
    for j in range(N + 1):
        p = pattern(N, j)
        aset = set(allowed(N, j))
        e = sum(c < 2 and c != v for c, v in zip(word, p))
        a = sum(word[k] >= 2 for k in aset)
        binary_out = sum(word[k] < 2 and word[k] != p[k] for k in range(N) if k not in aset)
        result.append(count_profile(len(aset), m, u, e, a, binary_out, double_wild))
    return result


def load(word, m, double_wild=False):
    return sum((F(v, sphere_size(len(word), m, j))
                for j, v in enumerate(counts(word, m, double_wild))), F(0))


def L(m, N):
    return F(m, N - 2) + F(m, N - m - 1)


def nstar(m):
    N = m + 2
    while L(m, N) > 1:
        N += 1
    return N


def word_string(word):
    return ''.join('01XX'[x] for x in word)


def middle_word(N, a, u=4):
    assert 0 <= a <= N - u
    return (1,) * a + (2,) * u + (0,) * (N - a - u)


def middle4_formula(N, m):
    return F(m * (m - 1) * (9*N*N + 12*N*m - 87*N - 5*m*m - 23*m + 174),
             4 * (N - 2) * (N - 3) * (N - 4) * (N - m - 1))


def H(N, m):
    """Sign polynomial: middle4_formula - L = m*H / positive denominator."""
    return (-8*N**3 + (13*m + 59)*N*N + (12*m*m - 127*m - 93)*N
            - 5*m**3 - 18*m*m + 245*m - 30)


def scaled_dp(N, m, u_values):
    """Exact maximum for every specified u by a max-plus DP.

    Fix global totals Z,O,u. State (z,o,last) after p coordinates records
    prefix totals. Adding the next letter completes pattern j=p. Hence every
    centre is a path and every path is a centre. Equal states have the same
    possible future contributions, so retaining the greatest partial sum is
    exact. Integer scaling is common to endpoint and interior spheres.
    """
    D = math.lcm(sphere_size(N, m, 0), sphere_size(N, m, 1))
    scale_i = D // sphere_size(N, m, 1)
    scale_e = D // sphere_size(N, m, 0)
    output = {}
    for u in u_values:
        if u > 2*m:
            # Every private word has only m free letters; d(word,centre)>=u-m.
            output[u] = dict(value='0', word='X'*u + '0'*(N-u),
                             equivalent_classes=choose(N, u)*2**(N-u), states=0,
                             reason='PROVED: u>2m gives zero load')
            continue
        M = N - 2
        interior = {}
        for e in range(N - u + 1):
            for left in range(3):
                for right in range(3):
                    a = u - (left == 2) - (right == 2)
                    bout = (left == 0) + (right == 1)
                    interior[e, left, right] = count_profile(M, m, u, e, a, bout) * scale_i
        best = -1
        best_word = None
        state_count = 0
        for O in range(N - u + 1):
            Z = N - u - O
            states = {}
            for c in range(3):
                z, o, x = int(c == 0), int(c == 1), int(c == 2)
                if z <= Z and o <= O and x <= u:
                    val = count_profile(N - 1, m, u, O, u - x, int(c == 1)) * scale_e
                    states[z, o, c] = (val, (c,))
            for p in range(1, N):
                state_count += len(states)
                nxt = {}
                for (z, o, last), (acc, word) in states.items():
                    e = z + O - o
                    for c in range(3):
                        z2, o2 = z + (c == 0), o + (c == 1)
                        x2 = p + 1 - z2 - o2
                        if z2 > Z or o2 > O or x2 > u:
                            continue
                        val = acc + interior[e, last, c]
                        key = (z2, o2, c)
                        old = nxt.get(key)
                        if old is None or val > old[0]:
                            nxt[key] = (val, word + (c,))
                states = nxt
            state_count += len(states)
            for (z, o, last), (acc, word) in states.items():
                assert z == Z and o == O
                val = acc + count_profile(N - 1, m, u, Z, u - (last == 2), int(last == 0)) * scale_e
                if val > best:
                    best, best_word = val, word
        exact = F(best, D)
        assert exact == load(best_word, m)
        output[u] = dict(value=str(exact), word=word_string(best_word),
                         equivalent_classes=choose(N, u)*2**(N-u), states=state_count)
    return output


def controls(mutant=False):
    # Actual word enumeration checks endpoints and excluded free coordinates.
    N, m = 7, 2
    spheres = [list(sphere(N, m, j)) for j in range(N + 1)]
    maxima = defaultdict(lambda: F(-1))
    comparisons = mutant_mismatches = 0
    first_mutant = None
    for word in itertools.product(range(3), repeat=N):
        brute_counts = [sum(distance(word, x) <= m for x in sj) for sj in spheres]
        got = counts(word, m, mutant)
        if mutant:
            for j, (a, b) in enumerate(zip(got, brute_counts)):
                if a != b:
                    mutant_mismatches += 1
                    if first_mutant is None:
                        first_mutant = dict(word=word_string(word), j=j, mutant=a, correct=b)
        else:
            assert got == brute_counts, (word, got, brute_counts)
            v = sum((F(x, len(sj)) for x, sj in zip(brute_counts, spheres)), F(0))
            u = word.count(2)
            maxima[u] = max(maxima[u], v)
        comparisons += N + 1
    if mutant:
        print(json.dumps(dict(status='EXPECTED-FAIL', control='double-count excluded free letters',
                              mismatches=mutant_mismatches, first=first_mutant)), flush=True)
        assert mutant_mismatches > 0
        return 1
    dp = scaled_dp(N, m, range(N + 1))
    assert all(F(dp[u]['value']) == maxima[u] for u in maxima)
    # Letter-2/3 symmetry checked on all actual four-letter centres in a second cell.
    symmetry_checks = 0
    for word in itertools.product(range(4), repeat=5):
        collapsed = tuple(min(x, 2) for x in word)
        assert counts(word, 1) == counts(collapsed, 1)
        symmetry_checks += 6
    return dict(status='VERIFIED-COMPUTATIONALLY', actual_sphere_comparisons=comparisons,
                small_cell=[N, m], brute_classes=3**N, all_u_dp_equals_brute=True,
                letter_symmetry_checks=symmetry_checks)


def witness_record(N, m, a):
    word = middle_word(N, a)
    v = load(word, m)
    ell = L(m, N)
    return dict(N=N, m=m, a=a, u=4, word=word_string(word), load=str(v), L=str(ell),
                load_minus_L=str(v-ell), S_holds_for_this_word=(v < max(F(1), ell)),
                Nstar=nstar(m), per_pattern=[str(F(x, sphere_size(N,m,j)))
                                          for j,x in enumerate(counts(word,m))])


def verify():
    result = dict(controls=controls())
    records = [witness_record(41,19,19), witness_record(39,18,18),
               witness_record(76,35,36), witness_record(74,34,35)]
    assert records[0]['load'] == '2571333812653/1837953717600'
    assert not records[0]['S_holds_for_this_word']
    assert records[1]['S_holds_for_this_word']
    assert not records[2]['S_holds_for_this_word']
    assert records[3]['S_holds_for_this_word']
    result['witnesses_and_adjacent_controls'] = records
    # Independent per-pattern counting checks the derived hockey-stick formula.
    checks = 0
    for m in range(4, 46):
        for N in range(2*m+6, 2*m+13):
            for a in (m+1, N-m-5):
                assert load(middle_word(N,a),m) == middle4_formula(N,m)
                den = 4*(N-2)*(N-3)*(N-4)*(N-m-1)
                assert middle4_formula(N,m)-L(m,N) == F(m*H(N,m),den)
                checks += 1
    result['middle4_formula_checks'] = checks
    # Exact sign characterisation of the fully interior family.
    table = []
    for m in (19,22,25,30,34,35,40,50,75,100,200):
        ns = nstar(m)
        failing = [N for N in range(2*m+6, ns+1)
                   if middle4_formula(N,m) >= max(F(1),L(m,N))]
        table.append(dict(m=m, Nstar=ns, first_failing=min(failing) if failing else None,
                          last_failing=max(failing) if failing else None,
                          count=len(failing)))
    result['interior_u4_failure_band'] = table
    # Boundary onset for this particular interior family: m=35, not a global onset claim.
    for m in range(4, 1001):
        assert (H(2*m+6,m)>0) == (m>=35)
        ns = nstar(m)
        assert F(ns-1)>F(5*m,2)
        assert H(ns-1,m)<0
    result['polynomial_sign_checks_m4_through_1000'] = True
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='command', required=True)
    sub.add_parser('verify')
    sub.add_parser('mutant-control')
    cell = sub.add_parser('cell')
    cell.add_argument('N', type=int)
    cell.add_argument('m', type=int)
    cell.add_argument('--u', type=int, nargs='+')
    args = ap.parse_args()
    if args.command == 'mutant-control':
        return controls(mutant=True)
    if args.command == 'verify':
        print(json.dumps(verify(), indent=2), flush=True)
    else:
        u_values = args.u if args.u is not None else list(range(args.N + 1))
        results = {}
        for u in u_values:
            data = scaled_dp(args.N, args.m, [u])[u]
            results[u] = data
            print(json.dumps(dict(N=args.N,m=args.m,u=u,**data)),flush=True)
        best_u = max(results, key=lambda u: F(results[u]['value']))
        print(json.dumps(dict(summary=True,N=args.N,m=args.m,best_u=best_u,
                              maximum=results[best_u]['value'],
                              word=results[best_u]['word'],
                              L=str(L(args.m,args.N)),
                              classes=sum(d['equivalent_classes'] for d in results.values()))),flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
