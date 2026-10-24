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


def finite(m, Ns):
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from step3_dp import scaled_max                # returns (best, L, word); M = best / L
    for N in Ns:
        best, L, _ = scaled_max(N, m, N + 1)
        M, n = Fraction(best, L), N + 1
        R = n if M <= 1 else max([r for r in range(1, n + 1) if r * (M - 1) < M], default=0)
        print(f"m={m} n={n} M={M} R={R} every_r={R == n}", flush=True)


def main():
    failures = []
    if sys.argv[1:2] == ["--finite"]:
        if len(sys.argv) < 4:
            print("usage: frac_alln_tail_check.py --finite m N [N ...]   (N = n - 1)")
            return
        finite(int(sys.argv[2]), map(int, sys.argv[3:]))
        return
    if sys.argv[1:2] == ["--centres"]:
        if len(sys.argv) != 5:
            print("usage: frac_alln_tail_check.py --centres m N k")
            return
        centres(*map(int, sys.argv[2:5]))
        return
    if sys.argv[1:] not in ([], ["--lemmas"]):
        print(__doc__)
        return
    lemmas = "--lemmas" in sys.argv[1:]
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
        if rises or not all_below:
            failures.append(f"B(N) at m={m}: rises={rises}, all_below={all_below}")
    if lemmas:
        for m in range(6, 13):
            bad, tested = length_lemma(m, 62)
            print(f"length lemma m={m}: {bad} increases in {tested} (n, sigma, eps, D) tests, pools n up to 61 against n + 1",
                  flush=True)
            if bad:
                failures.append(f"length lemma m={m}: {bad} increases in {tested} tests")
        for N in range(2, 9):
            bad, tested = neighbour_lemma(N)
            print(f"neighbour lemma N={N}: {bad} violations in {tested} (centre, E) tests", flush=True)
            if bad:
                failures.append(f"neighbour lemma N={N}: {bad} violations in {tested} tests")

    # A verdict and an exit code. This script counted its violations, printed them and exited 0,
    # so a run with a real violation was indistinguishable from a clean one to anything but a human
    # reading the lines. Found on 2026-10-23 by a sweep for this shape across the whole repository,
    # after two seats hit it independently in two other folders.
    print("\nVERDICT:", "PASS" if not failures else "FAIL")
    for f in failures:
        print("   ", f)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main() or 0)
