#!/usr/bin/env python3
"""Reproduce every figure in the appendix candidate-site table from one deposited file.

Input: producers_2026-10-20/PANELREF_GRCh38_NGG.json (deposited).
Output: the five-number summaries the table prints for the two count rows and for the
paired per-guide factor row. Pass the file path as the only argument.

The table's printed values use Python's round(), which breaks a half toward the even
digit. Two cells sit exactly on a half (12,122.5 and 17,821.5) and both match.
"""
import json
import statistics
import sys


def five(values, digits=0):
    ordered = sorted(values)
    q1, median, q3 = statistics.quantiles(ordered, n=4, method="inclusive")
    out = [ordered[0], q1, median, q3, ordered[-1]]
    return [round(v, digits) if digits else round(v) for v in out]


def main(path):
    rows = json.load(open(path))["rows"]
    sites = lambda k: [r["sites"][k] for r in rows]
    expected = {
        "N at <= 4 mismatches": [23, 107, 178, 299, 4542],
        "N at <= 6 mismatches": [3103, 8756, 12122, 17822, 103955],
        "paired factor N6/N4": [22.9, 54.6, 70.0, 86.9, 185.6],
    }
    got = {
        "N at <= 4 mismatches": five(sites("4")),
        "N at <= 6 mismatches": five(sites("6")),
        "paired factor N6/N4": five([r["sites"]["6"] / r["sites"]["4"] for r in rows], 1),
    }
    bad = 0
    print("guides read: %d" % len(rows))
    for name, want in expected.items():
        mark = "ok  " if got[name] == want else "FAIL"
        bad += got[name] != want
        print("  [%s] %-22s %s" % (mark, name, got[name]))
    if bad:
        print("GATE FAIL: %d row(s) differ from the published table." % bad)
        return 1
    print("All three table rows reproduce from the deposited panel file.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
