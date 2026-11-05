#!/usr/local/bin/python3
"""Execute deliberately false certificates/maps and require their rejection."""
import copy
import csv
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys

import verify_finite as v

ROOT = Path(__file__).resolve().parent
CONTROL = ROOT / "controls"


def subprocess_control(name, arguments, expected_text):
    p = subprocess.run([sys.executable, "-B", str(ROOT / "verify_finite.py"), *arguments],
                       cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (CONTROL / (name + ".log")).write_text(p.stdout + p.stderr)
    assert p.returncode == 1, (name, "control unexpectedly accepted", p.returncode)
    assert expected_text in p.stderr, (name, p.stderr)
    return {"name": name, "status": "REFUTED", "expected": "REJECTED",
            "actual": "REJECTED", "exit_code": p.returncode,
            "diagnostic": p.stderr.strip()}


def assertion_control(name, description, check):
    try:
        check()
    except AssertionError as e:
        return {"name": name, "status": "REFUTED", "expected": "REJECTED",
                "actual": "REJECTED", "mutation": description, "diagnostic": str(e)}
    raise AssertionError(f"Control unexpectedly accepted: {name}")


def main():
    CONTROL.mkdir(exist_ok=True)
    source = json.loads((ROOT / "certificates.json").read_text())
    outcomes = []
    bad = copy.deepcopy(source)
    c = next(c for c in bad["cases"] if (c["N"], c["r"], c["m"]) == (3, 2, 0))
    for entry in c["packing"]:
        entry["weight_per_word"] = str(2 * F(entry["weight_per_word"]))
    path = CONTROL / "bad_packing.json"
    path.write_text(json.dumps(bad, indent=2) + "\n")
    outcomes.append(subprocess_control("inflated_packing", ["--cert", str(path)],
                                       "packing ball load exceeds 1"))
    bad = copy.deepcopy(source)
    c = next(c for c in bad["cases"] if (c["N"], c["r"], c["m"]) == (3, 2, 0))
    c["covering"] = []
    path = CONTROL / "bad_covering.json"
    path.write_text(json.dumps(bad, indent=2) + "\n")
    outcomes.append(subprocess_control("empty_cover", ["--cert", str(path)], "cover short: 0"))
    with (ROOT / "tau_values.csv").open() as f:
        table = list(csv.DictReader(f))
    c = next(c for c in table if (c["N"], c["r"], c["m"]) == ("3", "2", "1"))
    c["tau_exact"] = "3"
    path = CONTROL / "bad_tau_values.csv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=table[0].keys())
        writer.writeheader()
        writer.writerows(table)
    outcomes.append(subprocess_control("injected_increase", ["--table", str(path)],
                                       "MONOTONICITY FAIL N=3, r=2, m=1: 2 < 3"))
    def wrong_direction():
        values = {(c["N"], c["r"], c["m"]): F(c["tau"]) for c in source["cases"]}
        for m in range(1, 4):
            v.require(values[3, 2, m-1] <= values[3, 2, m],
                      f"reversed inequality rejected at m={m}: "
                      f"{values[3,2,m-1]} > {values[3,2,m]}")
    outcomes.append(assertion_control("reversed_monotonicity",
                    "Assert tau(m-1) <= tau(m) on the actual exact N=3,r=2 row", wrong_direction))
    pats = v.patterns(3, 1)
    x, y = (1, 0, 0), (0, 0, 0)
    assert v.canonical_map(x, pats, 1) == y
    def identity_map():
        bad_image = x
        v.require(min(v.distance(bad_image, p) for p in pats) <= 0,
                  "identity mutant leaves shell word 100 outside U_0")
    outcomes.append(assertion_control("uncorrected_shell",
                    "Fix a shell word instead of correcting it", identity_map))
    def overwritten_collision():
        weights = {y: F(1, 2), x: F(1, 2)}
        # This source is feasible: its entire mass is exactly one.
        assert sum(weights.values(), F(0)) == 1
        pushed = {}
        for point, weight in weights.items():
            pushed[v.canonical_map(point, pats, 1)] = weight
        v.require(sum(pushed.values(), F(0)) == sum(weights.values(), F(0)),
                  "overwrite mutant retains 1/2 of source mass 1")
    outcomes.append(assertion_control("overwrite_pushforward",
                    "Overwrite colliding weights instead of summing each fibre", overwritten_collision))
    def same_radius_preimage():
        centre = (0, 1, 0)
        assert v.distance(y, centre) <= 1
        v.require(v.distance(x, centre) <= 1,
                  "same-radius mutant fails: d(f(100),010)=1 but d(100,010)=2")
    outcomes.append(assertion_control("same_radius_inclusion",
                    "Replace B_(m-1)(c) by B_m(c) in the preimage claim", same_radius_preimage))
    assert len(outcomes) == 7 and all(o["actual"] == "REJECTED" for o in outcomes)
    report = {"status": "VERIFIED-COMPUTATIONALLY", "controls": outcomes,
              "expected_rejections": 7, "actual_rejections": 7, "unexpected_acceptances": 0}
    (ROOT / "controls.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
