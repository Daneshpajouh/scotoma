#!/usr/bin/env python3
"""Run step3_dp.scaled_max at the full library r = N+1 at the lengths given, for the finite part of
the weighting-at-every-length check (finite_all_143_rows.log was made with it).
The full-library maximum M bounds every run count and every threshold set, so the floor holds at
every r <= N+1 if (N+1)(M-1) < M when M > 1 (and always when M <= 1); the condition is
sufficient, not necessary.

Usage: python3 release_boundary_rows.py m N [N ...]
       python3 release_boundary_rows.py --all [m ...]

--all replays every row of finite_all_143_rows.log, beside this script, and compares each with the
recorded line: the exact maximum M, the maximising centre and every_r_proved must match (the timing
is not compared). Give budgets after --all to replay only those rows. The whole log, all 143 rows,
took 8,606 s of computation when it was made (the sum of its secs fields), about 2.4 hours on one
core; m = 6 alone takes under a minute. Exit 0 when every replayed row matches; 1 on any mismatch;
2 if it could not run: the log missing, not 143 rows or unreadable, no row with the budgets given,
or the arguments wrong.
"""
import os
import sys
import time
from fractions import Fraction

N_ROWS = 143
USAGE = "usage: release_boundary_rows.py m N [N ...]  or  release_boundary_rows.py --all [m ...]"

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))           # step3_dp.py lives one folder up
from step3_dp import scaled_max  # noqa: E402  released producer, unchanged


def row(m, N):
    t0 = time.time()
    best, L, word = scaled_max(N, m, N + 1)
    M = Fraction(best, L)
    r = N + 1
    holds_all_r = (M <= 1) or (r * (M - 1) < M)
    line = (f"m={m} N={N} n={N+1} r={r} M={M} ({float(M):.6f}) word={''.join('01X'[c] for c in word)} "
            f"every_r_proved={holds_all_r}")
    return line, time.time() - t0


def replay_all(budgets):
    path = os.path.join(HERE, "finite_all_143_rows.log")
    try:
        recorded = [ln.rstrip("\n") for ln in open(path) if ln.strip()]
    except OSError as exc:
        print(f"NOT RUN: {exc}")
        return 2
    if len(recorded) != N_ROWS:
        print(f"NOT RUN: {path} has {len(recorded)} rows, not {N_ROWS}")
        return 2
    wanted = []
    seen = []
    for k, ln in enumerate(recorded, 1):
        try:
            fields = dict(f.split("=", 1) for f in ln.split() if "=" in f)
            m, N = int(fields["m"]), int(fields["N"])
        except (KeyError, ValueError):
            print(f"NOT RUN: line {k} of {path} has no readable m and N")
            return 2
        seen.append((m, N))
        if not budgets or m in budgets:
            wanted.append((m, N, ln.rsplit(" secs=", 1)[0]))
    # 2026-10-24: the row COUNT was pinned but not WHICH cells, and the cells to replay are read
    # from the log itself. So replacing one row with a copy of another kept the count at 143 and
    # passed: the duplicate was replayed twice and the displaced cell was never computed at all.
    # A log of 143 rows describes 143 distinct cells, so requiring that is free and closes it.
    if len(set(seen)) != len(seen):
        dup = sorted({p for p in seen if seen.count(p) > 1})
        print(f"NOT RUN: {path} has {len(recorded)} rows but only {len(set(seen))} distinct "
              f"(m, N) cells; repeated: {dup}")
        print("         A row replaced by a copy of another keeps the count and silently drops "
              "a cell.")
        return 2
    if not wanted:
        print("NOT RUN: no recorded row has the budgets given")
        return 2
    mismatches = 0
    for m, N, expected in wanted:
        line, secs = row(m, N)
        same = line == expected
        mismatches += not same
        print(f"{'ok      ' if same else 'MISMATCH'} {line} secs={secs:.0f}", flush=True)
        if not same:
            print(f"         recorded {expected}", flush=True)
    print(f"replayed {len(wanted)} of {len(recorded)} recorded rows, mismatches {mismatches}")
    return 1 if mismatches else 0


def main(argv):
    if len(argv) >= 2 and argv[1] == "--all":
        try:
            budgets = {int(a) for a in argv[2:]}
        except ValueError:
            print(USAGE)
            return 2
        return replay_all(budgets)
    try:
        m, lengths = int(argv[1]), [int(a) for a in argv[2:]]
    except (IndexError, ValueError):
        print(USAGE)
        return 2
    if not lengths:
        print(USAGE)
        return 2
    for N in lengths:
        line, secs = row(m, N)
        print(f"{line} secs={secs:.0f}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
