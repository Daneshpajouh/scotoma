#!/usr/bin/env python3
"""Rebuild the paper's frontier table from released files alone, in one command.

Four reviewers asked for this in four wordings: a single documented command that regenerates the
length-and-budget table from the released certificates, so the all-length theorem is checkable
without trusting any path that was not published.

Everything it uses is in this folder. It reads no cluster path, no working-tree file and no
precomputed log; it recomputes.

TWO BLOCKS, BECAUSE THE TABLE HAS TWO KINDS OF ROW AND THEY REST ON DIFFERENT OBJECTS.

  1. The all-length rows. For each budget m, the threshold is the smallest library length N whose
     step-3 maximum M satisfies M < 1 + 1/N. Centre the weighting on a pattern and it collects its
     own private sphere, so M >= 1 always and the question is whether M is close enough to 1 that
     the rescaled bound beats r-1 at every run count the length admits. The largest run count at
     length N is r = N+1, which is why the criterion carries 1/N. Computed here with topr.dp_shifted.

     The spacer length printed in the paper is n = N+1. Mixing the two is the single easiest
     mistake in this table: the closed form ceil((3m+3+sqrt(5m^2-2m+1))/2) is stated in N and is
     SUFFICIENT rather than tight, so it sits one above the measured threshold at m = 5, 7, 8, 10
     and 11. This script prints the measured value and says which is which.

  2. The n = 20 certificate rows, from the four released duals. Each is swept exhaustively by
     verify_ceiling_exhaustive.py, which shares no code with the solver.

SHARPNESS IS REPORTED, NOT ASSUMED. For every budget the script also evaluates the criterion one
length BELOW the threshold and prints whether it fails there. A threshold whose predecessor also
passed would not be a threshold, and the table would be quoting a weaker length than the paper has.

Usage:
    python3 rebuild_frontier_table.py            rebuild, print, and compare against the paper
    python3 rebuild_frontier_table.py --check    the same, plus a closing confirmation line

Both modes compare against the published rows and exit 1 on any disagreement, so one run is
enough; --check only adds the last line. The whole run takes about 10 to 12 minutes on a busy
machine, most of it in the last three rows (m = 10, 11 and 12 take roughly 1, 2.5 and 5 to 6
minutes), and each row announces itself before its search starts.

Exit 0 if the rebuild succeeds and matches every row the paper prints; 1 on a mismatch, including a
published row that was not compared; 2 if it could not run.
"""
import math
import os
import subprocess
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

try:
    from topr import dp_shifted
except ImportError:
    sys.stderr.write("NOT RUN: topr.py must sit beside this script.\n")
    raise SystemExit(2)

# What the paper prints, in spacer length n. Written here so a drift in either direction is caught
# rather than silently accepted.
PUBLISHED_N = {8: 23, 9: 26, 10: 28, 11: 31, 12: 34}
# The budgets the paper prints a threshold for. Checked against PUBLISHED_N's keys and against the
# number of rows actually compared, because a comparison only covers the keys it is given: with a
# key deleted the old gate still reported that the rebuilt rows equal the published rows. Found by
# a clean-clone negative control (nc_drop) on 2026-10-23.
PUBLISHED_BUDGETS = (8, 9, 10, 11, 12)
CERTS = ("cert_m9_r9.json", "cert_m10_r7.json", "cert_m11_r6.json", "cert_m12_r4.json")


def closed_form_N(m):
    """The sufficient condition, in N. Not the threshold; printed to show the gap."""
    return math.ceil((3 * m + 3 + math.sqrt(5 * m * m - 2 * m + 1)) / 2)


def threshold(m, span=40):
    """Smallest N with M < 1 + 1/N, with M there and M one length below."""
    for N in range(2 * m + 3, 2 * m + 3 + span):
        M, _ = dp_shifted(N, m, Fraction(0))
        if M < 1 + Fraction(1, N):
            prev, _ = dp_shifted(N - 1, m, Fraction(0))
            return N, M, prev
    return None, None, None


def sweep_certificates():
    """Run the exhaustive sweep on the four released duals. Returns (ok, lines)."""
    verifier = os.path.join(HERE, "verify_ceiling_exhaustive.py")
    paths = [os.path.join(HERE, "data", "floor_certificates", c) for c in CERTS]
    missing = [p for p in paths if not os.path.exists(p)]
    if not os.path.exists(verifier) or missing:
        return None, ["  certificates or verifier not present beside this script"]
    r = subprocess.run([sys.executable, "-I", verifier] + paths,
                       capture_output=True, text=True)
    lines = [l for l in r.stdout.strip().splitlines() if l.strip()]
    # Keep the verifier's own errors. Without them a verifier that died before printing (a Python
    # version error, say) would surface only as "a certificate did not pass", hiding the cause.
    lines += ["stderr: " + l for l in r.stderr.strip().splitlines() if l.strip()]
    return r.returncode == 0, lines


def main(argv):
    check = "--check" in argv
    print("Rebuilding the frontier table from released files in this folder.\n", flush=True)

    print("BLOCK 1  all-length rows, recomputed from topr.dp_shifted", flush=True)
    print(f"  {'m':>3} {'N':>4} {'n=N+1':>6} {'M at N':>12} {'1+1/N':>10} "
          f"{'fails at N-1':>13} {'N* sufficient':>14} {'paper n':>8}", flush=True)
    bad = 0
    could_not_run = False
    rows = {}
    compared = 0
    if tuple(sorted(PUBLISHED_N)) != PUBLISHED_BUDGETS:
        print(f"  PUBLISHED_N names budgets {sorted(PUBLISHED_N)}, the paper prints "
              f"{list(PUBLISHED_BUDGETS)}. A missing key is a row nobody compares.", flush=True)
        bad += 1
    for m in range(5, 13):
        print(f"  m = {m}: searching from N = {2 * m + 3} ...", flush=True)
        N, M, prev = threshold(m)
        if N is None:
            print(f"  {m:>3}  no threshold in range")
            bad += 1
            continue
        sharp = not (prev < 1 + Fraction(1, N - 1))
        pub = PUBLISHED_N.get(m)
        rows[m] = N + 1
        flag = ""
        if N + 1 < 20:
            # Outside the paper's scope (n >= 20). Printed for completeness, not claimed.
            flag += "  (n < 20: below the paper's scope, not a claim it makes)"
        if pub is not None:
            compared += 1
            if pub != N + 1:
                flag = f"  <- paper says {pub}"
                bad += 1
        print(f"  {m:>3} {N:>4} {N+1:>6} {float(M):>12.6f} {float(1+Fraction(1,N)):>10.6f} "
              f"{str(sharp):>13} {closed_form_N(m):>14} {'' if pub is None else pub:>8}{flag}", flush=True)
        if not sharp:
            print(f"       NOT SHARP at m={m}: the criterion also holds at N-1, so this is not "
                  "the first length that works")
            bad += 1

    print(f"\n  {compared} of {len(PUBLISHED_BUDGETS)} published rows compared.", flush=True)
    if compared != len(PUBLISHED_BUDGETS):
        print("  A published row was not compared, so this run cannot say the table matches.")
        bad += 1

    print("\nBLOCK 2  n = 20 certificate rows, swept exhaustively", flush=True)
    ok, lines = sweep_certificates()
    for l in lines:
        print("  " + l)
    if ok is None:
        print("  BLOCK 2 NOT RUN: the released certificates are not beside this script.")
        could_not_run = True
    elif not ok:
        print("  BLOCK 2 FAILED: a certificate did not pass its exhaustive sweep.")
        bad += 1

    print()
    if bad:
        print(f"{bad} problem(s). The rebuilt table does not match what the paper prints, or a "
              "check did not run. Neither is a pass.")
        if could_not_run:
            print("A released file is also missing, so the run was not complete either.")
        return 1
    if could_not_run:
        print("Could not run in full: a released file is missing. That is not a pass.")
        return 2
    print("Rebuilt. Each threshold is the first length at which THIS WEIGHTING certifies every "
          "run count; that is not the onset of the floor, which a different argument can "
          "reach sooner (Theorem 4 settles m = 8 at n = 20). Every published threshold "
          "reproduces, and all four duals pass an exhaustive sweep.")
    if check:
        print("--check: the rebuilt rows equal the published rows.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
