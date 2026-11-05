#!/usr/bin/env python3
"""The sound-pattern theorem, every part that needs no input file.

THE STATEMENT. Fix a spacer, let L be its complete single-deletion library (words of length N = n - 1),
U the union of the radius-m balls around L, and m <= N - 2. A pattern over the four letters plus a
wildcard ACCEPTS a word when the word differs from it at no more than m fixed (non-wildcard)
coordinates, and is SOUND when every word it accepts lies in U. The theorem says that at four letters
no sound pattern carries a wildcard, and the only sound wildcard-free patterns (centres) are the r
members of L. On the 110 real guides, steps (iii) and (iv) of the proof are checked by
sound_panel.py, which reads the CHANGE-seq table.

WHAT THIS SCRIPT CHECKS, ALL EXHAUSTIVELY.
  1. Every guide of length n = 4, 5, 6 over ACGT (N = 3, 4, 5), every m from 0 to N - 1, every
     pattern of length N over ACGT plus a wildcard, and every word each pattern accepts. Inside
     m <= N - 2: no guide has a sound wildcard pattern, and every guide's sound centres are exactly
     L. At m = N - 1 the wildcard half still holds, and the centre half fails: sound centres outside L
     appear for exactly the guides with run count 3 or more.
  2. The same enumeration over a two-letter alphabet, N = 3 to 6, as the control that the
     enumeration is not vacuous: there, sound wildcard patterns exist inside m <= N - 2. The proof
     needs a letter no library word uses at each coordinate, which two letters cannot supply.
  3. ACGTA at m = 2, the paper's smallest case where the covering number is below the library size:
     its sound centres are exactly the five members of L, no four of them cover U, and the printed
     four-ball cover {ACGA, ACTT, AGGT, CGTA} covers U using unsound centres (ACTT accepts 4 words
     outside U, AGGT 8). That the covering number is 4 is floor_counterexample.py's check.
  4. CAGAG at m = 3 = N - 1: eleven sound centres outside L.
  5. Step (iii)'s disagreement claim: for every unsorted binary word b of block length k = 1 to 19,
     with u the first 1 and v the last 0 (so u < v), every sorted word disagrees with b at u or at v.
     That is sum over k of (2^k - k - 1)(k + 1) = 19,920,075 cases.
  6. Step (iv)'s private points, built as the proof builds them and checked by distance, for every
     guide of length 4 to 7 that is not a single repeated letter and every 1 <= m <= N - 2: 427,800.

CONTROLS, RUN EVERY TIME. The two-letter enumeration of part 2 must find sound wildcard patterns.
Part 5 is rerun with v moved to a coordinate that carries the same bit as u, and must fail. Part 6 is
rerun with one coordinate fewer changed, and must fail.

Usage:  sound_patterns.py
Python 3.9 or later, standard library only, about 20 s.
Exit 0 when every part holds and every control fails as it should; 1 otherwise.
"""
import itertools
import sys
import time

ALPHA = "ACGT"


def ham(a, b):
    return sum(x != y for x, y in zip(a, b))


def library(g):
    return sorted({g[:i] + g[i + 1:] for i in range(len(g))})


def runs(g):
    return 1 + sum(g[p] != g[p + 1] for p in range(len(g) - 1))


def enumerate_patterns(q, N):
    """Per m: (guides, guides with a sound wildcard pattern, guides whose sound centres are not L,
    guides with sound centres outside L, guides with r >= 3, the guides with sound centres outside L
    that have r >= 3)."""
    wild = q
    words = list(itertools.product(range(q), repeat=N))
    pats = list(itertools.product(range(q + 1), repeat=N))
    word_ix = {w: i for i, w in enumerate(words)}
    wild_mask = sum(1 << i for i, p in enumerate(pats) if wild in p)
    centre_of_word = [pats.index(w) for w in words]
    centre_mask = sum(1 << i for i in centre_of_word)
    # accepts[m][o]: patterns accepting word o at budget m; ball[m][w]: words within m of w
    accepts = [[0] * len(words) for _ in range(N)]
    ball = [[0] * len(words) for _ in range(N)]
    for o, w in enumerate(words):
        for i, p in enumerate(pats):
            d = sum(1 for a, b in zip(p, w) if a != wild and a != b)
            for m in range(d, N):
                accepts[m][o] |= 1 << i
        for j, x in enumerate(words):
            d = sum(1 for a, b in zip(w, x) if a != b)
            for m in range(d, N):
                ball[m][o] |= 1 << j
    every_word = (1 << len(words)) - 1
    tally = {m: [0, 0, 0, 0, 0, 0] for m in range(N)}
    for g in itertools.product(range(q), repeat=N + 1):
        lib = {g[:i] + g[i + 1:] for i in range(N + 1)}
        lib_ix = [word_ix[w] for w in lib]
        lib_centres = sum(1 << centre_of_word[i] for i in lib_ix)
        r = 1 + sum(g[p] != g[p + 1] for p in range(N))
        for m in range(N):
            union = 0
            for i in lib_ix:
                union |= ball[m][i]
            unsound = 0
            outside = every_word & ~union
            while outside:
                low = outside & -outside
                unsound |= accepts[m][low.bit_length() - 1]
                outside ^= low
            centres = centre_mask & ~unsound
            sound_wild = bool(wild_mask & ~unsound)
            t = tally[m]
            t[0] += 1
            t[1] += sound_wild
            t[2] += centres != lib_centres
            t[3] += bool(centres & ~lib_centres)
            t[4] += r >= 3
            t[5] += bool(centres & ~lib_centres) and r >= 3
    return tally


def sound_centres(g, m):
    N = len(g) - 1
    L = library(g)
    words = ["".join(t) for t in itertools.product(ALPHA, repeat=N)]
    U = {w for w in words if min(ham(w, l) for l in L) <= m}
    ball = {c: frozenset(w for w in words if ham(w, c) <= m) for c in words}
    return L, U, ball, [c for c in words if ball[c] <= U]


def disagreement(K, broken=False):
    """Step (iii) over every unsorted binary word to block length K: cases and failures."""
    cases = fails = 0
    for k in range(1, K + 1):
        sorted_words = [((1 << i) - 1) << (k - i) for i in range(k + 1)]  # zeros, then i ones
        for b in range(1 << k):
            ones = [p for p in range(k) if (b >> p) & 1]
            zeros = [p for p in range(k) if not (b >> p) & 1]
            if not ones or not zeros or max(zeros) < min(ones):
                continue                                          # sorted: zeros then ones
            u, v = min(ones), max(zeros)
            if broken:
                same = [p for p in range(k) if p != u and ((b >> p) & 1) == 1]
                v = same[0] if same else u
            for s in sorted_words:
                cases += 1
                if ((s >> u) & 1) == ((b >> u) & 1) and ((s >> v) & 1) == ((b >> v) & 1):
                    fails += 1
    return cases, fails


def private_points(changed_less=0):
    """Step (iv) for every guide of length 4 to 7 with r > 1 and every 1 <= m <= N - 2."""
    built = fails = 0
    for n in range(4, 8):
        N = n - 1
        for t in itertools.product(ALPHA, repeat=n):
            g = "".join(t)
            if runs(g) == 1:
                continue
            L = library(g)
            free = [next(x for x in ALPHA if all(w[p] != x for w in L)) for p in range(N)]
            for m in range(1, N - 1):
                for j, D in enumerate(L):
                    nbrs = [w for w in L if ham(w, D) == 1]
                    uv = {next(p for p in range(N) if w[p] != D[p]) for w in nbrs}
                    pos = [p for p in range(N) if p not in uv][:m - changed_less]
                    x = "".join(free[p] if p in pos else c for p, c in enumerate(D))
                    built += 1
                    if ham(x, D) > m or any(ham(x, w) <= m for w in L if w != D):
                        fails += 1
    return built, fails


def main():
    t0 = time.time()
    problems = []

    print("1. four letters, every guide, every pattern with or without a wildcard")
    for N in (3, 4, 5):
        tally = enumerate_patterns(4, N)
        for m, (gs, wild, not_l, outside, r3, outside_r3) in sorted(tally.items()):
            zone = "m <= N-2" if m <= N - 2 else "m = N-1 "
            print(f"   N {N} m {m} [{zone}] guides {gs:5d}  with a sound wildcard pattern {wild}"
                  f"  sound centres not L {not_l}  with centres outside L {outside}  (r >= 3: {r3})")
            if wild:
                problems.append(f"N {N} m {m}: a sound wildcard pattern at four letters")
            if m <= N - 2 and not_l:
                problems.append(f"N {N} m {m}: sound centres differ from L inside m <= N-2")
            if m == N - 1 and not (outside == r3 == outside_r3):
                problems.append(f"N {N} m {m}: centres outside L are not exactly the r >= 3 guides")

    print("2. control: two letters, where the free-letter argument has nothing to use")
    found = 0
    for N in range(3, 7):
        for m, row in sorted(enumerate_patterns(2, N).items()):
            if m <= N - 2:
                found += row[1]
    print(f"   guides with a sound wildcard pattern inside m <= N-2, N = 3 to 6: {found}")
    if not found:
        problems.append("control: the two-letter enumeration found no sound wildcard pattern")

    print("3. ACGTA at m = 2")
    L, U, ball, sound = sound_centres("ACGTA", 2)
    four_cover = any(set().union(*(ball[c] for c in combo)) >= U for combo in itertools.combinations(sound, 4))
    cover = ["ACGA", "ACTT", "AGGT", "CGTA"]
    covers = set().union(*(ball[c] for c in cover)) >= U
    outside = {c: len(ball[c] - U) for c in cover}
    print(f"   sound centres {sorted(sound)} equal L {sorted(sound) == L} | four sound centres cover U"
          f" {four_cover} | printed cover covers U {covers} | accepts outside U {outside}")
    if sorted(sound) != L or four_cover or not covers or outside != {"ACGA": 0, "ACTT": 4, "AGGT": 8, "CGTA": 0}:
        problems.append("ACGTA does not read as printed")

    print("4. CAGAG at m = 3 = N - 1")
    L, U, ball, sound = sound_centres("CAGAG", 3)
    extra = len(set(sound) - set(L))
    print(f"   sound centres outside L: {extra}")
    if extra != 11:
        problems.append(f"CAGAG: {extra} sound centres outside L, not 11")

    print("5. step (iii), every unsorted binary word to block length 19")
    cases, fails = disagreement(19)
    expected = sum((2 ** k - k - 1) * (k + 1) for k in range(1, 20))
    print(f"   cases {cases:,} (expected {expected:,}) failures {fails}")
    if cases != expected or cases != 19920075 or fails:
        problems.append("step (iii) does not hold as printed")
    c_cases, c_fails = disagreement(8, broken=True)
    print(f"   control, v moved to a coordinate carrying u's bit, to length 8: failures {c_fails} of {c_cases}")
    if not c_fails:
        problems.append("control: step (iii) with a wrong v did not fail")

    print("6. step (iv) private points, guides of length 4 to 7, 1 <= m <= N-2")
    built, fails = private_points()
    print(f"   built {built:,} failures {fails}")
    if built != 427800 or fails:
        problems.append("private points do not hold as printed")
    c_built, c_fails = private_points(changed_less=1)
    print(f"   control, one coordinate fewer changed: failures {c_fails:,} of {c_built:,}")
    if not c_fails:
        problems.append("control: private points with one change fewer did not fail")

    print(f"seconds {time.time() - t0:.0f}")
    for p in problems:
        print(f"FAIL: {p}")
    print("VERDICT:", "FAIL" if problems else "PASS")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
