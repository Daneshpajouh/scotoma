#!/usr/bin/env python3
"""Second independent implementation of the private-sphere fractional weighting at N=19, m=6.

Written from the DESCRIPTION, sharing no code with the originating seat, because that is the whole
value of a second check: two copies of one implementation agree whenever that implementation is
wrong.

THE OBJECT, as described and as I read it.
  Patterns P_j = 1^j 0^(19-j) for j = 0..19, so r = 20, over the alphabet {0,1,2,3}. The patterns
  use only {0,1}; letters 2 and 3 are the two no pattern uses.
  P_j's boundary sits between coordinates j and j+1 (1-indexed), so the two coordinates next to it
  are j and j+1. ALLOWED_j = {1..19} minus those, which is 17 coordinates for j = 1..18 and 18 for
  j = 0 and j = 19.
  PRIVATE SPHERE S_j = the words obtained from P_j by overwriting exactly m = 6 coordinates of
  ALLOWED_j with letters from {2,3}. Every such word is at distance exactly 6 from P_j.
  WEIGHT: one unit spread evenly over each S_j, so each point of S_j carries 1 / |S_j|.

PRIVACY, which I re-derived rather than took, because it explains the exclusion.
  x in S_j overwrites a set T of 6 allowed coordinates. P_i differs from P_j exactly on the
  interval B between them, |B| = |i - j|, and that interval is adjacent to P_j's boundary.
  d(x, P_i) = 6 + |i - j| - |T intersect B|, because the overwritten coordinates carry letters no
  pattern uses and so differ from P_i too. Privacy needs d(x, P_i) > 6, i.e. |T intersect B| <
  |i - j|. B always CONTAINS the excluded coordinate adjacent to the boundary, so |T intersect B|
  <= |i - j| - 1. Hence privacy holds for every i != j, with no computation.
  That is also exactly why the whole-sphere version leaks: with j+1 allowed, a T containing it puts
  x at distance 6 from P_{j+1} as well, so two units land in one ball.

THE BALL COUNT IN CLOSED FORM, derived here independently.
  Fix a centre c. Split ALLOWED_j by what c carries:
    a        coordinates where c is in {2,3}  (all of these disagree with P_j)
    b_diff   coordinates where c is in {0,1} and disagrees with P_j
    b_same   coordinates where c is in {0,1} and agrees with P_j
  and let e = how many of the EXCLUDED coordinates c disagrees with P_j on.
  Choose T with t_a, t_d, t_s coordinates from those three groups, t_a + t_d + t_s = 6, and let k of
  the t_a coordinates receive a letter different from c's. Then
    d(c, x) = a + b_diff + e - t_a + t_s + k
  and the number of such x is C(a,t_a) C(b_diff,t_d) C(b_same,t_s) C(t_a,k) 2^(t_d + t_s).
  So the ball count depends on c only through (a, b_diff, b_same, e), which is the reduction the
  originating seat describes. Everything below is exact rational arithmetic.

CONTROLS, run before any number is reported:
  C1 |S_j| must come to 16,633,344 in total, and every point of S_j must be at distance exactly 6
     from P_j, checked by sampling actual words rather than trusting the formula.
  C2 the closed-form ball count must agree with a BRUTE-FORCE count over an explicitly enumerated
     S_j at a smaller N and m, where enumeration is possible. If the formula is wrong this fires.
  C3 privacy must hold: no point of S_j may lie within 6 of any other pattern. Checked by sampling.
  C4 the weighting must FAIL on the paper's own ACGTA counterexample, where the floor is false. A
     construction that succeeded there would be wrong.
"""
import itertools
import random
import sys
from fractions import Fraction
from math import comb

random.seed(20261022)


def allowed(j, N):
    excl = {p for p in (j, j + 1) if 1 <= p <= N}
    return [p for p in range(1, N + 1) if p not in excl]


def pattern(j, N):
    return [1] * j + [0] * (N - j)          # 1-indexed coordinate p is pattern[p-1]


def sphere_size(j, N, m):
    return comb(len(allowed(j, N)), m) * 2 ** m


def summaries(c, j, N):
    """(a, b_diff, b_same, e) for centre c against pattern j."""
    P = pattern(j, N)
    A = set(allowed(j, N))
    a = b_diff = b_same = e = 0
    for p in range(1, N + 1):
        agree = (c[p - 1] == P[p - 1])
        if p in A:
            if c[p - 1] >= 2:
                a += 1
            elif agree:
                b_same += 1
            else:
                b_diff += 1
        else:
            if not agree:
                e += 1
    return a, b_diff, b_same, e


def ball_count(c, j, N, m):
    """#{x in S_j : d(c,x) <= m}, closed form."""
    a, b_diff, b_same, e = summaries(c, j, N)
    tot = 0
    for t_a in range(0, min(a, m) + 1):
        for t_d in range(0, min(b_diff, m - t_a) + 1):
            t_s = m - t_a - t_d
            if t_s < 0 or t_s > b_same:
                continue
            base = a + b_diff + e - t_a + t_s
            for k in range(0, t_a + 1):
                if base + k <= m:
                    tot += (comb(a, t_a) * comb(b_diff, t_d) * comb(b_same, t_s)
                            * comb(t_a, k) * 2 ** (t_d + t_s))
    return tot


def ball_weight(c, N, m, r):
    return sum(Fraction(ball_count(c, j, N, m), sphere_size(j, N, m)) for j in range(r))


def enumerate_sphere(j, N, m):
    P = pattern(j, N)
    A = allowed(j, N)
    out = []
    for T in itertools.combinations(A, m):
        for letters in itertools.product((2, 3), repeat=m):
            x = P[:]
            for p, L in zip(T, letters):
                x[p - 1] = L
            out.append(tuple(x))
    return out


def d(x, y):
    return sum(p != q for p, q in zip(x, y))


def controls():
    ok = True
    N, m, r = 19, 6, 20

    # C1 total
    total = sum(sphere_size(j, N, m) for j in range(r))
    if total != 16633344:
        print(f'[C1 FAIL] total private points {total:,}, described as 16,633,344')
        ok = False
    else:
        n18 = sum(1 for j in range(r) if len(allowed(j, N)) == 17)
        n17 = sum(1 for j in range(r) if len(allowed(j, N)) == 18)
        print(f'controls: C1 total private points {total:,} as described; '
              f'{n18} patterns exclude two coordinates, {n17} exclude one')

    # C1b every point at distance exactly m, by sampling real words
    for j in (0, 7, 19):
        A = allowed(j, N)
        P = pattern(j, N)
        for _ in range(200):
            T = random.sample(A, m)
            x = P[:]
            for p in T:
                x[p - 1] = random.choice((2, 3))
            if d(x, P) != m:
                print(f'[C1b FAIL] a sampled point of S_{j} is at distance {d(x, P)}, not {m}')
                return False
    print(f'controls: C1b sampled points of S_0, S_7, S_19 are all at distance exactly {m}')

    # C2 closed form against brute force at a size where enumeration is possible
    for (n, mm) in ((7, 2), (8, 3)):
        rr = n + 1
        for j in range(rr):
            S = enumerate_sphere(j, n, mm)
            if len(S) != sphere_size(j, n, mm):
                print(f'[C2 FAIL] |S_{j}| enumerated {len(S)} vs formula {sphere_size(j,n,mm)}')
                return False
            for _ in range(6):
                c = [random.randrange(4) for _ in range(n)]
                brute = sum(1 for x in S if d(c, x) <= mm)
                closed = ball_count(c, j, n, mm)
                if brute != closed:
                    print(f'[C2 FAIL] N={n} m={mm} j={j} centre {c}: brute {brute} '
                          f'vs closed form {closed}')
                    return False
    print('controls: C2 closed-form ball count agrees with brute force at (N,m) = (7,2) and (8,3), '
          'every pattern, random centres')

    # C3 privacy by sampling
    for j in (0, 5, 12, 19):
        A = allowed(j, N)
        P = pattern(j, N)
        for _ in range(300):
            T = random.sample(A, m)
            x = P[:]
            for p in T:
                x[p - 1] = random.choice((2, 3))
            for i in range(r):
                if i != j and d(x, pattern(i, N)) <= m:
                    print(f'[C3 FAIL] a point of S_{j} lies within {m} of P_{i}')
                    return False
    print('controls: C3 privacy holds on every sampled point against every other pattern')
    return ok


def main():
    if not controls():
        print('\nCONTROLS FAILED, nothing below is reportable.')
        return 1
    N, m, r = 19, 6, 20

    print(f'\nball weight at the patterns (each must be exactly 1):')
    bad = False
    for j in (0, 7, 19):
        w = ball_weight(pattern(j, N), N, m, r)
        print(f'  P_{j:<2d} -> {w}  ({float(w):.6f})')
        if w != 1:
            bad = True
    for j in range(r):
        if ball_weight(pattern(j, N), N, m, r) != 1:
            print(f'  [fail] P_{j} does not carry exactly 1')
            bad = True
    if bad:
        return 1
    print(f'  all {r} patterns carry exactly 1')

    # The claimed best non-pattern centre.
    for tail in (2, 3):
        c = [1] * 7 + [0] * 11 + [tail]
        w = ball_weight(c, N, m, r)
        print(f'\nclaimed best non-pattern centre 1^7 0^11 {tail}: {w} = {float(w):.6f}')
        print(f'  described as 29/34 = {float(Fraction(29,34)):.6f}  '
              f'-> {"MATCHES" if w == Fraction(29,34) else "DOES NOT MATCH"}')

    # C4 the weighting must fail where the floor is false: ACGTA, N = 4, m = 2, r = 5.
    print('\nC4, the counterexample it must fail on (N=4, m=2, r=5):')
    n4, m4, r4 = 4, 2, 5
    worst = max(ball_weight([*cc], n4, m4, r4)
                for cc in itertools.product(range(4), repeat=n4))
    print(f'  maximum ball weight over all {4**n4} centres: {worst} = {float(worst):.4f}')
    if worst <= 1:
        print('  [C4 FAIL] the weighting is feasible here, but the floor is FALSE at these '
              'parameters, so a feasible weighting would prove something untrue')
        return 1
    print(f'  [ok  ] infeasible, as it must be: some ball carries {worst} > 1')
    return 0


if __name__ == '__main__':
    sys.exit(main())
