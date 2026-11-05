#!/usr/bin/env python3
"""The column rule, checked over every column, and the counting bound that rests on it.

The appendix paragraph "A counting bound that settles m = 6 without any search" rests on one
inequality about a single coordinate:

    Take w words and one coordinate. The w pattern entries there are j ones and w - j zeros, and t
    of the w entries are overwritten with other letters. Then D, the number of disagreeing pairs in
    that column, is at most

        j(w - j) + (w - 2)t,

    plus 1 when the column is constant (j = 0 or j = w) and t >= 1.

This file checks that by enumeration rather than by argument: every w from 2 to 6, every j, every
choice of which entries are overwritten, and every letter each overwritten entry can take over the
four-letter alphabet. It then checks the two structural facts the appendix count uses, and the
arithmetic that turns them into N >= 5m - 6.

WHY A SEPARATE CHECK. The inequality is small and the paper's conclusion is not: it rules out every
scattered transversal, cyclic or otherwise, at N = 4m - 1 for m >= 6, which had previously needed a
SAT run and an exhaustive search. A reader should be able to test the one step the argument turns
on without taking either on trust.

Usage:  check_column_rule.py [--quiet]
Exit 0 when every part passes, 1 when one does not.
"""
import sys
from itertools import combinations, product

ALPHABET = (0, 1, 2, 3)
WMAX = 6


def bound(w, j, t):
    b = j * (w - j) + (w - 2) * t
    if (j == 0 or j == w) and t >= 1:
        b += 1
    return b


def columns(w, j, t):
    """Every column with j ones, w-j zeros, and t entries overwritten by another letter."""
    pattern = [1] * j + [0] * (w - j)
    for changed in combinations(range(w), t):
        choices = [tuple(a for a in ALPHABET if a != pattern[i]) for i in changed]
        for newvals in product(*choices):
            col = list(pattern)
            for i, v in zip(changed, newvals):
                col[i] = v
            yield col


def disagreeing(col):
    w = len(col)
    return sum(1 for a in range(w) for b in range(a + 1, w) if col[a] != col[b])


def part_a(quiet):
    """The rule itself, over every column."""
    seen = tight = bad = 0
    for w in range(2, WMAX + 1):
        for j in range(0, w + 1):
            for t in range(0, w + 1):
                b = bound(w, j, t)
                best = -1
                for col in columns(w, j, t):
                    d = disagreeing(col)
                    seen += 1
                    best = max(best, d)
                    if d > b:
                        bad += 1
                        if not quiet:
                            print(f"  [A] VIOLATION w={w} j={j} t={t} column={col} D={d} > {b}")
                if best == b:
                    tight += 1
    print(f"  [A] {'PASS' if not bad else 'FAIL'}: {seen:,} columns at w = 2 to {WMAX}, "
          f"{bad} violations, {tight} (w,j,t) cases attain the bound")
    return bad == 0


def part_b(quiet):
    """The control. A slope of w-3 must be violated, or the rule is not tight and proves less."""
    bad = 0
    for w in range(2, WMAX + 1):
        for j in range(0, w + 1):
            for t in range(0, w + 1):
                b = j * (w - j) + (w - 3) * t + (1 if (j in (0, w) and t >= 1) else 0)
                for col in columns(w, j, t):
                    if disagreeing(col) > b:
                        bad += 1
    print(f"  [B] {'PASS' if bad else 'FAIL'}: slope w-3 instead of w-2 gives {bad} violations "
          f"(must be positive, else the rule is generous and the bound weaker than claimed)")
    return bad > 0


def part_c(quiet):
    """The two structural facts: four split columns summing to 20, and N-4 constant ones."""
    bad = 0
    for N in (19, 23, 24, 31):
        for k in range(0, N - 4 + 1):
            pats = [[1] * (k + i) + [0] * (N - k - i) for i in range(5)]
            splits, const = [], 0
            for p in range(N):
                j = sum(q[p] for q in pats)
                if 0 < j < 5:
                    splits.append(j * (5 - j))
                else:
                    const += 1
            if sorted(splits) != [4, 4, 6, 6] or const != N - 4:
                bad += 1
                if not quiet:
                    print(f"  [C] VIOLATION N={N} k={k}: splits={splits} const={const}")
    print(f"  [C] {'PASS' if not bad else 'FAIL'}: five consecutive staircase words always split in "
          f"exactly four columns contributing 4+6+6+4 = 20, with N-4 constant")
    return bad == 0


def part_d(quiet):
    """The arithmetic: 10(2m+1) <= 20 + 15m + (N-4) is N >= 5m-6, and it binds from m = 6."""
    bad = 0
    for m in range(2, 21):
        implied = 10 * (2 * m + 1) - 20 - 15 * m + 4
        if implied != 5 * m - 6:
            bad += 1
        binds = (5 * m - 6) > (4 * m - 1)
        if binds != (m >= 6):
            bad += 1
    print(f"  [D] {'PASS' if not bad else 'FAIL'}: 10(2m+1) <= 20 + 15m + (N-4) is N >= 5m-6 for "
          f"every m, and 5m-6 exceeds 4m-1 exactly from m = 6")
    return bad == 0


def main(argv):
    quiet = "--quiet" in argv
    ok = all([part_a(quiet), part_b(quiet), part_c(quiet), part_d(quiet)])
    print()
    print("ALL PARTS PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
