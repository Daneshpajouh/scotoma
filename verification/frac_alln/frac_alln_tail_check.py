#!/usr/bin/env python3
"""The tail of the private-sphere weighting, checked: the length lemma, the neighbour lemma
and the tail bound B(N), as in the appendix subsection on the weighting at every length.

The tail, lemma and centre modes use only the standard library; --finite also needs the released
step3_dp.py beside this file. For release beside step3_dp.py.

Setting. Centres are words over {0, 1, X} of length N (X = either letter no pattern uses). The
chain patterns are P_j = 1^j 0^(N-j), j = 0 .. N. For a centre c, s = number of X coordinates and
e_j = number of 0/1 coordinates where c differs from P_j, so d(c, P_j) = s + e_j.

Share. The share of P_j's private sphere inside B_m(c) depends on the pool size p = |A_j| (written
n in the code below; it is not the spacer length), the X coordinates sigma and 0/1 differences eps
of c inside A_j, and D = d(c, P_j). Excluded coordinates still count in D but are never overwritten.

f(N, s, e): the largest share over pools n in {N-2, N-1, N} and over every placement of the stars
and differences on the h = N - n excluded coordinates (s_out + d_out <= h); infeasible cases count
as share 0. Admitting n = N can only raise f; for N >= 2m + 1 it changes nothing, by the length lemma.
g(E) = max over e >= E of f(N, s, e).

Tail bound. For a centre with s X coordinates whose nearest chain pattern is at 0/1 distance E0,
    sum_j w_j(c) <= B(N, s, E0) = (2 E0 + s + 1) g(E0) + 2 sum_{E > E0} g(E)
(neighbour lemma, then summation by parts). B(N) = max over (s, E0) != (0, 0) with s + E0 <= 2m.
By the length lemma, B(N) is nonincreasing in N for N >= 2m + 2 (where every profile with
s + e <= 2m fits the smaller pool), so B(N0) < 1 settles every N >= N0 >= 2m + 2: chain patterns
collect exactly 1 and every other centre less.

Usage: frac_alln_tail_check.py            N0 table for m = 6..12, monotonicity from N = 2m + 2 to 100
       frac_alln_tail_check.py --lemmas   also the length lemma sweep and the neighbour lemma brute force
       frac_alln_tail_check.py --centres m N k
                                          per-centre test: every non-pattern centre over {0,1,X} within
                                          Hamming distance k of a chain pattern, against its own
                                          B(N, s, E0); control: the count E + 1 in place of 2E + s + 1
                                          at every E, whose violations are counted (positive at the three
                                          (m, N) the appendix reports)
       frac_alln_tail_check.py --finite m N [N ...]
                                          the full-chain maximum M at r = n = N + 1 from the released
                                          step3_dp.scaled_max(N, m, N + 1), imported from beside this file, which returns
                                          (best, L, word) with M = best / L; and R, the
                                          largest run count r <= n with M <= 1 or r (M - 1) < M

Exit codes. 0: every check passed and every printed value equals its published value. 1: a lemma
test found a violation or increase, a centre exceeded its bound, the tail table or a finite row
differs from the published one (PUBLISHED_TAIL and PUBLISHED_FINITE below, generated from the
released logs; a finite row is compared on M and every_r), or a table row is inconsistent
(B(N0) >= 1, B(N0 - 1) < 1, a rise, or a value >= 1 past N0). 2: a usage error. A pass is decided
by this program, not by reading its output.
"""
import itertools
import sys
from fractions import Fraction
from math import comb


def share(n, m, sigma, eps, D):
    """Numerator of the share of the private sphere inside the ball, over C(n, m) 2^m."""
    total = 0
    rest = n - sigma - eps
    for x in range(0, min(sigma, m) + 1):
        for y in range(0, min(eps, m - x) + 1):
            k = m - x - y
            if k > rest:
                continue
            ways = comb(sigma, x) * comb(eps, y) * comb(rest, k) * 2 ** (m - x)
            need = D - x - y                      # z must be at least this
            good = sum(comb(x, z) for z in range(max(need, 0), x + 1)) if need <= x else 0
            total += ways * good
    return total


def F(n, m, sigma, eps, D):
    if sigma + eps > n or n < m:
        return Fraction(0)
    return Fraction(share(n, m, sigma, eps, D), comb(n, m) * 2 ** m)


def f(N, m, s, e):
    best = Fraction(0)
    for n in (N - 2, N - 1, N):
        out = N - n
        for so in range(0, min(s, out) + 1):          # stars on excluded positions
            for eo in range(0, min(e, out - so) + 1):  # differences on excluded positions
                best = max(best, F(n, m, s - so, e - eo, s + e))
    return best


def B(N, m, pattern_count=lambda E, s: 2 * E + s + 1):
    best, arg = Fraction(0), None
    for s in range(0, 2 * m + 1):
        fs = [f(N, m, s, e) for e in range(0, 2 * m - s + 2)]
        g = [max(fs[e:]) for e in range(len(fs))]
        gg = g + [Fraction(0)]
        for E0 in range(0, 2 * m - s + 1):
            if s == 0 and E0 == 0:
                continue
            # summation by parts with a general count c(E): sum c(E) (g(E) - g(E+1))
            val = sum(pattern_count(E, s) * (gg[E] - gg[E + 1]) for E in range(E0, len(g)))
            if val > best:
                best, arg = val, (s, E0)
    return best, arg


def length_lemma(m, nmax):
    bad = tested = 0
    for n in range(m, nmax):
        for sigma in range(0, 2 * m + 1):
            for eps in range(0, 2 * m + 1 - sigma):
                if sigma + eps > n:
                    continue
                for D in range(sigma + eps, 2 * m + 1):   # D = d(c, P_j) >= sigma + eps
                    tested += 1
                    if F(n + 1, m, sigma, eps, D) > F(n, m, sigma, eps, D):
                        bad += 1
    return bad, tested


def neighbour_lemma(N):
    """Over all centres in {0,1,X}^N: for every E, patterns with e_j <= E number at most 2E+s+1."""
    bad = tested = 0
    for c in itertools.product((0, 1, 2), repeat=N):
        s = sum(1 for x in c if x == 2)
        es = []
        for j in range(N + 1):
            P = [1] * j + [0] * (N - j)
            es.append(sum(1 for a, b in zip(c, P) if a != 2 and a != b))
        for E in range(0, N + 1):
            tested += 1
            if sum(1 for x in es if x <= E) > 2 * E + s + 1:
                bad += 1
    return bad, tested


def ball_sum(c, m):
    """Full-chain load W(c) at a centre over {0,1,2}, 2 = X, using the step 2 share formula; with s and E0."""
    N = len(c)
    s = sum(1 for x in c if x == 2)
    tot, es = Fraction(0), []
    for j in range(N + 1):
        A = [p for p in range(1, N + 1) if p not in (j, j + 1)]
        e = sigma = eps = 0
        for pos in range(1, N + 1):
            a, b = c[pos - 1], 1 if pos <= j else 0
            if a == 2:
                sigma += pos in A
            elif a != b:
                e += 1
                eps += pos in A
        es.append(e)
        tot += F(len(A), m, sigma, eps, s + e)
    return tot, s, min(es)


def centres(m, N, k):
    gs = {}

    def bound(s, E0, count):
        if s not in gs:
            fs = [f(N, m, s, e) for e in range(0, 2 * m - s + 2)]
            gs[s] = [max(fs[e:]) for e in range(len(fs))] + [Fraction(0)]
        g = gs[s]
        return sum(count(E, s) * (g[E] - g[E + 1]) for E in range(E0, len(g) - 1))

    seen, tested, viol, ctrl = set(), 0, 0, 0
    for base in range(N + 1):
        P = [1] * base + [0] * (N - base)
        for t in range(1, k + 1):
            for pos in itertools.combinations(range(N), t):
                for kinds in itertools.product((2, "flip"), repeat=t):
                    c = list(P)
                    for p, kd in zip(pos, kinds):
                        c[p] = 2 if kd == 2 else 1 - c[p]
                    key = tuple(c)
                    if key in seen:
                        continue
                    seen.add(key)
                    v, s, E0 = ball_sum(key, m)
                    if s == 0 and E0 == 0:
                        continue                      # a chain pattern itself
                    tested += 1
                    if v > bound(s, E0, lambda E, s_: 2 * E + s_ + 1):
                        viol += 1
                    if v > bound(s, E0, lambda E, s_: E + 1):
                        ctrl += 1
    print(f"per-centre m={m} N={N} k={k}: {tested} non-pattern centres; {viol} above B(N,s,E0); "
          f"control (count E+1): {ctrl} above", flush=True)
    return 0 if tested > 0 and viol == 0 else 1


def finite(m, Ns):
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from step3_dp import scaled_max                # returns (best, L, word); M = best / L
    bad = 0
    for N in Ns:
        best, L, _ = scaled_max(N, m, N + 1)
        M, n = Fraction(best, L), N + 1
        R = n if M <= 1 else max([r for r in range(1, n + 1) if r * (M - 1) < M], default=0)
        want = PUBLISHED_FINITE.get((m, N))
        if want is None:
            note = 'no published row'
        elif (Fraction(want[0]), want[1]) == (M, R == n):
            note = 'equals the published row'
        else:
            note = 'DIFFERS from the published row M=%s every_r=%s' % want
            bad += 1
        print(f"m={m} n={n} M={M} R={R} every_r={R == n}  {note}", flush=True)
    return 1 if bad else 0


PUBLISHED_TAIL = {  # m: (N0, B(N0), B(N0-1), worst (s, E0)); generated from logs/tail_and_lemmas.log
    6: (27, '24/25', '77/76', '(2, 0)'),
    7: (31, '651/667', '45/44', '(2, 0)'),
    8: (35, '424/429', '103/100', '(2, 0)'),
    9: (39, '1071/1073', '29/28', '(2, 0)'),
    10: (44, '75/77', '165/164', '(2, 0)'),
    11: (48, '407/414', '319/315', '(2, 0)'),
    12: (52, '322/325', '948/931', '(2, 0)'),
}

PUBLISHED_FINITE = {  # (m, N): (M, every_r); generated from logs/finite_all_143_rows.log, 143 rows
    (6, 19): ('1', True),
    (6, 20): ('1', True),
    (6, 21): ('1', True),
    (6, 22): ('1', True),
    (6, 23): ('1', True),
    (6, 24): ('1', True),
    (6, 25): ('1', True),
    (6, 26): ('1', True),
    (7, 19): ('196/187', True),
    (7, 20): ('1', True),
    (7, 21): ('1', True),
    (7, 22): ('1', True),
    (7, 23): ('1', True),
    (7, 24): ('1', True),
    (7, 25): ('1', True),
    (7, 26): ('1', True),
    (7, 27): ('1', True),
    (7, 28): ('1', True),
    (7, 29): ('1', True),
    (7, 30): ('1', True),
    (8, 19): ('108/85', False),
    (8, 20): ('116/99', False),
    (8, 21): ('62/57', False),
    (8, 22): ('66/65', True),
    (8, 23): ('1', True),
    (8, 24): ('1', True),
    (8, 25): ('1', True),
    (8, 26): ('1', True),
    (8, 27): ('1', True),
    (8, 28): ('1', True),
    (8, 29): ('1', True),
    (8, 30): ('1', True),
    (8, 31): ('1', True),
    (8, 32): ('1', True),
    (8, 33): ('1', True),
    (8, 34): ('1', True),
    (9, 19): ('46987/29920', False),
    (9, 20): ('517315/369512', False),
    (9, 21): ('270/209', False),
    (9, 22): ('6/5', False),
    (9, 23): ('102/91', False),
    (9, 24): ('81/77', False),
    (9, 25): ('1', True),
    (9, 26): ('1', True),
    (9, 27): ('1', True),
    (9, 28): ('1', True),
    (9, 29): ('1', True),
    (9, 30): ('1', True),
    (9, 31): ('1', True),
    (9, 32): ('1', True),
    (9, 33): ('1', True),
    (9, 34): ('1', True),
    (9, 35): ('1', True),
    (9, 36): ('1', True),
    (9, 37): ('1', True),
    (9, 38): ('1', True),
    (10, 19): ('34551/14144', False),
    (10, 20): ('17072293/8868288', False),
    (10, 21): ('8205/5168', False),
    (10, 22): ('5467093/3879876', False),
    (10, 23): ('55/42', False),
    (10, 24): ('175/143', False),
    (10, 25): ('185/161', False),
    (10, 26): ('13/12', False),
    (10, 27): ('41/40', True),
    (10, 28): ('1', True),
    (10, 29): ('1', True),
    (10, 30): ('1', True),
    (10, 31): ('1', True),
    (10, 32): ('1', True),
    (10, 33): ('1', True),
    (10, 34): ('1', True),
    (10, 35): ('1', True),
    (10, 36): ('1', True),
    (10, 37): ('1', True),
    (10, 38): ('1', True),
    (10, 39): ('1', True),
    (10, 40): ('1', True),
    (10, 41): ('1', True),
    (10, 42): ('1', True),
    (10, 43): ('1', True),
    (11, 19): ('6694579/1584128', False),
    (11, 20): ('53651089/17199104', False),
    (11, 21): ('46315249/19348992', False),
    (11, 22): ('36173889/18811520', False),
    (11, 23): ('188225/117572', False),
    (11, 24): ('15323547/10816624', False),
    (11, 25): ('396/299', False),
    (11, 26): ('209/168', False),
    (11, 27): ('88/75', False),
    (11, 28): ('231/208', False),
    (11, 29): ('484/459', False),
    (11, 30): ('253/252', True),
    (11, 31): ('1', True),
    (11, 32): ('1', True),
    (11, 33): ('1', True),
    (11, 34): ('1', True),
    (11, 35): ('1', True),
    (11, 36): ('1', True),
    (11, 37): ('1', True),
    (11, 38): ('1', True),
    (11, 39): ('1', True),
    (11, 40): ('1', True),
    (11, 41): ('1', True),
    (11, 42): ('1', True),
    (11, 43): ('1', True),
    (11, 44): ('1', True),
    (11, 45): ('1', True),
    (11, 46): ('1', True),
    (11, 47): ('1', True),
    (12, 19): ('67603/8704', False),
    (12, 20): ('187909627/34398208', False),
    (12, 21): ('17180561/4299776', False),
    (12, 22): ('226403179/75246080', False),
    (12, 23): ('1487209/632320', False),
    (12, 24): ('57083797/29745716', False),
    (12, 25): ('4980497/3090464', False),
    (12, 26): ('24051386/16900975', False),
    (12, 27): ('234/175', False),
    (12, 28): ('82/65', False),
    (12, 29): ('43/36', False),
    (12, 30): ('135/119', False),
    (12, 31): ('94/87', False),
    (12, 32): ('98/95', False),
    (12, 33): ('1', True),
    (12, 34): ('1', True),
    (12, 35): ('1', True),
    (12, 36): ('1', True),
    (12, 37): ('1', True),
    (12, 38): ('1', True),
    (12, 39): ('1', True),
    (12, 40): ('1', True),
    (12, 41): ('1', True),
    (12, 42): ('1', True),
    (12, 43): ('1', True),
    (12, 44): ('1', True),
    (12, 45): ('1', True),
    (12, 46): ('1', True),
    (12, 47): ('1', True),
    (12, 48): ('1', True),
    (12, 49): ('1', True),
    (12, 50): ('1', True),
    (12, 51): ('1', True),
}


def main():
    if sys.argv[1:2] == ["--finite"]:
        if len(sys.argv) < 4:
            print("usage: frac_alln_tail_check.py --finite m N [N ...]   (N = n - 1)")
            return 2
        return finite(int(sys.argv[2]), map(int, sys.argv[3:]))
    if sys.argv[1:2] == ["--centres"]:
        if len(sys.argv) != 5:
            print("usage: frac_alln_tail_check.py --centres m N k")
            return 2
        return centres(*map(int, sys.argv[2:5]))
    if sys.argv[1:] not in ([], ["--lemmas"]):
        print(__doc__)
        return 2
    lemmas = "--lemmas" in sys.argv[1:]
    fails = 0
    print("m  N0  n0  B(N0)  B(N0-1)  worst(s,E0)  B nonincreasing 2m+2..100  B(N)<1 for N0..100")
    for m in range(6, 13):
        N = 2 * m + 2                             # N0 is defined as the least N >= 2m + 2
        while B(N, m)[0] >= 1:                    # first N with B(N) < 1
            N += 1
        N0 = N
        b0, arg = B(N0, m)
        b1, _ = B(N0 - 1, m)
        prev, rises, all_below = None, 0, True
        for NN in range(2 * m + 2, 101):
            val, _ = B(NN, m)
            if prev is not None and val > prev:
                rises += 1
            if NN >= N0 and val >= 1:
                all_below = False
            prev = val
        print(f"{m}  {N0}  {N0 + 1}  {b0}  {b1}  {arg}  rises={rises}  {all_below}", flush=True)
        if (N0, str(b0), str(b1), str(arg)) != PUBLISHED_TAIL[m] or rises or not all_below or b0 >= 1 or b1 < 1:
            print(f"  FAIL m={m}: published {PUBLISHED_TAIL[m]}, rises={rises}, all_below={all_below}", flush=True)
            fails += 1
    if lemmas:
        for m in range(6, 13):
            bad, tested = length_lemma(m, 62)
            print(f"length lemma m={m}: {bad} increases in {tested} (n, sigma, eps, D) tests, pools n up to 61 against n + 1",
                  flush=True)
            fails += bad > 0 or tested == 0
        for N in range(2, 9):
            bad, tested = neighbour_lemma(N)
            print(f"neighbour lemma N={N}: {bad} violations in {tested} (centre, E) tests", flush=True)
            fails += bad > 0 or tested == 0
    print("ALL PASS" if not fails else f"FAILED: {fails} check(s)", flush=True)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
