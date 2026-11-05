#!/usr/bin/env python3
"""Re-derive the population law from the deposited executed records, and refuse a corrupted one.

The law: for single RNA-bulge expansion the 2016 route's word set is complete exactly when
g[0]==g[1] and g[-1]==g[-2], so it is incomplete on 15/16 of a uniform four-letter guide
population at every length n>=4.

This checker never trusts a stated fraction. For each deposited record it recomputes the predicate
count over the whole guide space and requires the executed `short` count to equal it, and requires
that the record reported no guide where execution and predicate disagreed.

  python3 -I verify_population_law.py             # check every deposited record
  python3 -I verify_population_law.py --control   # same, with one record corrupted; must FAIL
"""
import itertools
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def predicate_short_count(n):
    """Count guides the predicate calls incomplete, by enumeration over the whole space."""
    short = 0
    for letters in itertools.product("ACGT", repeat=n):
        complete = letters[0] == letters[1] and letters[-1] == letters[-2]
        short += not complete
    return short


def check(records, corrupt_first=False):
    failures = []
    for i, (path, rec) in enumerate(sorted(records)):
        rec = dict(rec)
        if corrupt_first and i == 0:
            rec["short"] = rec["short"] - 1          # remove one short guide
            rec["short_fraction"] = f"{rec['short']}/{rec['guides']}"
        n = rec["n"]
        if n < 4:
            failures.append(f"{path}: n={n} is outside the claimed range n>=4")
            continue
        expect = predicate_short_count(n)
        if rec["short"] != expect:
            failures.append(f"{path}: executed short={rec['short']} but the predicate gives {expect}")
        if rec["short"] + rec["complete"] != 4 ** n:
            failures.append(f"{path}: guides do not sum to 4^{n}")
        if rec["guides_where_execution_disagrees_with_predicate"] != 0:
            failures.append(f"{path}: record reports "
                            f"{rec['guides_where_execution_disagrees_with_predicate']} disagreements")
        if 16 * rec["short"] != 15 * rec["guides"]:
            failures.append(f"{path}: {rec['short_fraction']} is not 15/16")
    return failures


def main():
    control = "--control" in sys.argv[1:]
    records = []
    for path in sorted(HERE.glob("executed_population_n*.json")):
        records.append((path.name, json.loads(path.read_text())))
    if not records:
        print("no deposited records found")
        return 2
    failures = check(records, corrupt_first=control)
    for line in failures:
        print("FAIL ", line)
    if control:
        if failures:
            print(f"control fired: removing one short guide from {sorted(records)[0][0]} is caught")
            return 0
        print("CONTROL DEAD: a corrupted record passed")
        return 1
    if failures:
        return 1
    for name, rec in sorted(records):
        print(f"[ok  ] {name}: n={rec['n']} {rec['short_fraction']} = {rec['short_percent']}%, "
              f"0 disagreements, executed")
    print(f"[ok  ] {len(records)} executed record(s), each equal to the predicate count over all 4^n guides")
    return 0


if __name__ == "__main__":
    sys.exit(main())
