#!/usr/bin/env python3
"""Tests for offtarget_budget.py.

The load-bearing test is reproduce_published_table: the tool must regenerate the paper's
Table 1 for the exa-cel worked example, byte for byte against table_main.tex. If that fails,
the tool and the paper disagree and one of them is wrong.

Run: python3 test_offtarget_budget.py
"""
import json
import os
import subprocess
import sys
import tempfile
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from offtarget_budget import budget_rows, render_tex  # noqa: E402

N_EXACEL = 28175
PUBLISHED = {  # k: (sites to validate, i.i.d. draws), from Table 1 of the paper
    1: (26767, 84404),
    2: (21875, 42201),
    5: (12699, 16880),
    10: (7293, 8439),
    25: (3181, 3375),
}
FAILURES = []


def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def reproduce_published_rows():
    rows = {r["k"]: r for r in budget_rows(N_EXACEL, sorted(PUBLISHED), "0.05")}
    for k, (sites, iid) in PUBLISHED.items():
        check(f"k={k} sites to validate == {sites:,}",
              rows[k]["sites_to_validate"] == sites, rows[k]["sites_to_validate"])
        check(f"k={k} i.i.d. draws == {iid:,}",
              rows[k]["iid_draws"] == iid, rows[k]["iid_draws"])
    check("k=1 is 95.00% of N", f"{100 * rows[1]['fraction_of_N']:.2f}" == "95.00")
    check("k=10 is 25.88% of N", f"{100 * rows[10]['fraction_of_N']:.2f}" == "25.88")


def reproduce_published_table():
    want = open(os.path.join(HERE, "table_main.tex")).read().strip()
    got = render_tex(budget_rows(N_EXACEL, sorted(PUBLISHED), "0.05"), N_EXACEL, "0.05").strip()
    check("render_tex matches table_main.tex exactly", got == want)
    if got != want:
        for a, b in zip(got.splitlines(), want.splitlines()):
            if a != b:
                print(f"        got  {a!r}\n        want {b!r}")


def monotonicity():
    rows = budget_rows(N_EXACEL, [1, 2, 5, 10, 25, 50], "0.05")
    sites = [r["sites_to_validate"] for r in rows]
    check("the budget falls strictly in k", all(a > b for a, b in zip(sites, sites[1:])), sites)
    tight = budget_rows(N_EXACEL, [10], "0.01")[0]["sites_to_validate"]
    loose = budget_rows(N_EXACEL, [10], "0.05")[0]["sites_to_validate"]
    check("a smaller delta1 costs more sites", tight > loose, f"{tight} vs {loose}")


def exactness():
    """m* must be the smallest m that works: m*-1 must fail."""
    from m0_search import m0_hypergeom
    for k in (1, 7, 25):
        m = m0_hypergeom(N_EXACEL, k, F(1, 20))
        num, den = 1, 1
        for i in range(k):
            num *= (N_EXACEL - (m - 1) - i)
            den *= (N_EXACEL - i)
        check(f"k={k}: m*-1 does not reach the level", F(num, den) > F(1, 20))


def rejects_bad_input():
    bad = [
        ("N below 1", dict(n_sites=0, k_values=[1], delta1="0.05")),
        ("k above N", dict(n_sites=10, k_values=[11], delta1="0.05")),
        ("k below 1", dict(n_sites=10, k_values=[0], delta1="0.05")),
        ("delta1 == 0", dict(n_sites=10, k_values=[1], delta1="0")),
        ("delta1 == 1", dict(n_sites=10, k_values=[1], delta1="1")),
    ]
    for name, kw in bad:
        try:
            budget_rows(**kw)
            check(f"rejects {name}", False, "no error raised")
        except (ValueError, ZeroDivisionError):
            check(f"rejects {name}", True)


def delta1_is_exact():
    """A float delta1 must be refused: 0.05 is not 1/20 in binary."""
    try:
        budget_rows(N_EXACEL, [1], 0.05)
        check("refuses a float delta1", False, "accepted a float")
    except Exception:
        check("refuses a float delta1", True)


def cli_runs_and_is_consistent():
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "b.json")
        r = subprocess.run(
            [sys.executable, os.path.join(HERE, "offtarget_budget.py"),
             "--n", str(N_EXACEL), "--k", "1", "10", "--delta1", "0.05", "--out-json", out],
            capture_output=True, text=True)
        check("CLI exits 0", r.returncode == 0, r.stderr.strip()[:200])
        if r.returncode == 0:
            payload = json.load(open(out))
            got = {row["k"]: row["sites_to_validate"] for row in payload["rows"]}
            check("CLI JSON agrees with the published table",
                  got == {1: PUBLISHED[1][0], 10: PUBLISHED[10][0]}, got)
            check("CLI stdout states the level is repeated-sampling, not posterior",
                  "none" in r.stdout and "(A_k)" in r.stdout)


if __name__ == "__main__":
    for fn in (reproduce_published_rows, reproduce_published_table, monotonicity, exactness,
               rejects_bad_input, delta1_is_exact, cli_runs_and_is_consistent):
        print(f"{fn.__name__}:")
        fn()
    print()
    if FAILURES:
        print(f"FAILED {len(FAILURES)}: " + ", ".join(FAILURES))
        sys.exit(1)
    print("all tests passed")
