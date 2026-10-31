#!/usr/bin/env python3
"""Recompute the rescaling table's maxima, at the published lengths or beyond them.

Four of ten blind reviewers in one round gave the same limiting sentence, the one reading "at every
length checked, n = 20 through 25, and 30, and every m <= 7". One asked in terms for independently
executable exact verification records for those cells. This is that record: it recomputes each
maximum from the released dynamic program in exact rational arithmetic, and it reproduces the
published values first so a reader can see the harness agree with the paper before it is asked to
say anything the paper does not already contain.

WHAT THE NUMBER IS. M is the step-3 maximum at run count r = n, the largest the length admits.
The convention is worth stating because it is easy to get wrong and the number does not carry it:
dp_shifted takes N = n-1 and returns the maximum at r = N+1, which is n, matching the appendix row
that prints 196/187 "at r=20" for n=20. A maximum at the largest run count bounds every smaller
one, since the best r-subset of weights is a subset of the best n-subset, which is what lets a
single row speak for every run count. M is:
the greatest total weight any radius-m ball collects, where each pattern spreads one unit over its
private sphere. The weighting is feasible unrescaled exactly when M = 1; when M > 1 the rescaled
bound still proves the floor for every r < M/(M-1).

WHY M = 1 IS THE FLOOR, WHICH IS WHAT MAKES A ROW OF ONES MEANINGFUL RATHER THAN SUSPICIOUS. Put
the centre on a pattern. Every word of that pattern's private sphere lies at distance exactly m from
it, so the whole sphere is inside the radius-m ball and that pattern alone contributes a full unit.
So M >= 1 everywhere, with no computation, and a cell reading exactly 1 is at the trivial minimum
rather than merely small. This script treats any M < 1 as a fault in itself, not as a better bound.

Usage:
  verify_rescaling_grid.py                 controls only, a few seconds
  verify_rescaling_grid.py --paper         the lengths and budgets the paper claims
  verify_rescaling_grid.py --to 33         contiguous lengths 20..33 at m = 7, minutes
Exit 0 when every control reproduces and nothing is below 1, 1 on a disagreement, 2 if topr.py is
not beside this script.
"""
import sys
from fractions import Fraction

sys.path.insert(0, ".")
try:
    from topr import dp_shifted
except ImportError:
    sys.stderr.write("NOT RUN: topr.py is not in the current directory. It ships beside this\n"
                     "script in code_release; run this from there.\n")
    raise SystemExit(2)

# Fixed by the paper's own rescaling table. If these do not reproduce, nothing below is reportable.
CONTROLS = {
    (19, 6): Fraction(1),          # n=20, m<=6: "1, at every r"
    (19, 7): Fraction(196, 187),   # n=20, m=7
    (20, 7): Fraction(1),          # n=21, m<=7: "1, at every r"
    (19, 8): Fraction(108, 85),    # n=20, m=8: above 1, so a harness that only prints 1 fails here
}
PAPER_LENGTHS = [20, 21, 22, 23, 24, 25, 30]


def maximum(n, m):
    M, _ = dp_shifted(n - 1, m, Fraction(0))
    return M


def run_controls():
    print("CONTROLS, against the published table:")
    bad = 0
    for (N, m), want in sorted(CONTROLS.items()):
        got = maximum(N + 1, m)
        ok = got == want
        bad += not ok
        print(f"  n={N+1:<3} m={m}   got {str(got):>10}   published {str(want):>10}   "
              f"{'ok' if ok else 'MISMATCH'}")
    print("  The m=8 control is above 1 on purpose: a harness that answers 1 everywhere fails it.")
    return bad


def report(cells):
    print(f"\nMaxima at r = n, the largest run count each length admits.")
    print(f"{'n':>4} {'m':>3} {'maximum M':>14}   note")
    above = []
    for n, m in cells:
        M = maximum(n, m)
        if M < 1:
            print(f"{n:>4} {m:>3} {str(M):>14}   IMPOSSIBLE")
            print("\nM >= 1 holds by putting the centre on a pattern, so a value below 1 means this\n"
                  "computation is wrong rather than the bound improved.")
            return 1
        note = "feasible unrescaled" if M == 1 else \
               f"rescaled bound proves every r < {float(M / (M - 1)):.3f}"
        if M > 1:
            above.append((n, m, M))
        print(f"{n:>4} {m:>3} {str(M):>14}   {note}", flush=True)
    print(f"\n{len(above)} cell(s) above 1 out of {len(cells)}:")
    for n, m, M in above:
        print(f"  n={n}, m={m}: {M} = {float(M):.6f}")
    if not above:
        print("  none: every cell sits at the floor, so no rescaling is used anywhere here.")
    return 0


def main(argv):
    if run_controls():
        print("\nA published value did not reproduce, so no grid is reported.")
        return 1
    if "--paper" in argv:
        cells = [(n, m) for n in PAPER_LENGTHS for m in (5, 6, 7)]
    elif "--to" in argv:
        top = int(argv[argv.index("--to") + 1])
        cells = [(n, 7) for n in range(20, top + 1)]
    else:
        print("\nControls only. Pass --paper for the claimed cells, or --to N for lengths 20..N "
              "at m=7.")
        return 0
    return report(cells)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
