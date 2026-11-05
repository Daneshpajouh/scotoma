#!/usr/local/bin/python3
"""Independent exact verifier: standard library only, all 4**N words/centres.

No LP solver and no orbit incidence matrix is trusted. Certificate weights are
expanded to words. Integer Walsh convolution evaluates every Hamming-ball load.
All feasibility, mass, equality and monotonicity decisions are exact.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction as F
from itertools import product
import json
from math import comb, lcm
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def distance(x, y):
    return sum(a != b for a, b in zip(x, y))


def encode(word):
    n = 0
    for a in word:
        n = 4 * n + a
    return n


def patterns(N, r):
    L = r - 1
    return [tuple([0] * j + [1] * (L-j) + [0] * (N-L)) for j in range(r)]


def key_for_word(w, L):
    p = tuple(a if a in (0, 1) else 2 for a in w[:L])
    p2 = tuple({0: 1, 1: 0, 2: 2}[a] for a in p[::-1])
    return min(p, p2), sum(a != 0 for a in w[L:])


def expand(entries, N, r, words):
    table = {}
    for e in entries:
        p, h = tuple(e["prefix"]), e["suffix_nonzero"]
        require(len(p) == r-1 and all(a in (0, 1, 2) for a in p), "invalid prefix")
        require(isinstance(h, int) and 0 <= h <= N-r+1, "invalid suffix count")
        k = (p, h)
        require(k not in table, "duplicate orbit entry")
        canonical = key_for_word(p + (1,) * h + (0,) * (N-r+1-h), r-1)
        require(canonical == k, "noncanonical orbit entry")
        table[k] = F(e["weight_per_word"])
        require(table[k] >= 0, "negative weight")
    den = lcm(*(v.denominator for v in table.values()))
    nums = {k: v.numerator * (den // v.denominator) for k, v in table.items()}
    vector = [nums.get(key_for_word(w, r-1), 0) for w in words]
    return vector, den


def fwht(values):
    a = list(values)
    step = 1
    while step < len(a):
        for start in range(0, len(a), 2 * step):
            for j in range(start, start + step):
                x, y = a[j], a[j+step]
                a[j], a[j+step] = x+y, x-y
        step *= 2
    return a


def kernel_eigenvalues(N, m):
    # A nontrivial character sums to -1 over the three nonzero letters;
    # the trivial character sums to 3. Thus the generating polynomial for
    # frequency support t is (1+3z)**(N-t) * (1-z)**t.
    return [sum((-1)**j * comb(t, j) * comb(N-t, k-j) * 3**(k-j)
                for k in range(m+1)
                for j in range(max(0, k-(N-t)), min(t, k)+1))
            for t in range(N+1)]


def support_size(index):
    s = 0
    while index:
        s += (index & 3) != 0
        index >>= 2
    return s


def ball_sums(vector, N, m):
    eig = kernel_eigenvalues(N, m)
    transformed = fwht(vector)
    multiplied = [a * eig[support_size(i)] for i, a in enumerate(transformed)]
    inverse = fwht(multiplied)
    q = 4**N
    require(all(a % q == 0 for a in inverse), "nonintegral inverse convolution")
    return [a // q for a in inverse]


def check_kernel():
    # Direct Hamming-distance summation checks the full transform, not merely
    # transform/inverse consistency. Nonuniform integer weights prevent an
    # accidentally constant kernel from passing.
    N = 3
    words = list(product(range(4), repeat=N))
    v = [((17*i+5) % 31) - 11 for i in range(4**N)]
    for m in range(N+1):
        direct = [sum(v[j] for j, x in enumerate(words) if distance(c, x) <= m)
                  for c in words]
        require(ball_sums(v, N, m) == direct, f"convolution kernel failed at {m=}")
    return {"status": "VERIFIED-COMPUTATIONALLY", "radii": 4,
            "centre_checks_against_direct_summation": 256}


def canonical_map(word, pats, m):
    ds = [distance(word, p) for p in pats]
    nearest_distance = min(ds)
    require(nearest_distance <= m, "map source outside U_m")
    if nearest_distance <= m-1:
        return word
    j = ds.index(nearest_distance)
    k = next(k for k, (a, b) in enumerate(zip(word, pats[j])) if a != b)
    y = list(word)
    y[k] = pats[j][k]
    return tuple(y)


def check_case(case, words):
    N, r, m = case["N"], case["r"], case["m"]
    require(3 <= N <= 6 and 1 <= r <= N+1 and 0 <= m <= N, "invalid cell")
    pats = patterns(N, r)
    delta = [min(distance(w, p) for p in pats) for w in words]
    inside = [i for i, d in enumerate(delta) if d <= m]
    packing, pd = expand(case["packing"], N, r, words)
    covering, cd = expand(case["covering"], N, r, words)
    require(all(not v or delta[i] <= m for i, v in enumerate(packing)),
            f"{N,r,m}: packing outside union")
    loads = ball_sums(packing, N, m)
    coverage = ball_sums(covering, N, m)
    require(F(max(loads), pd) <= F(1),
            f"{N,r,m}: packing ball load exceeds 1: {F(max(loads), pd)}")
    require(F(min(coverage[i] for i in inside), cd) >= F(1),
            f"{N,r,m}: cover short: {F(min(coverage[i] for i in inside), cd)}")
    ptotal, ctotal, target = F(sum(packing), pd), F(sum(covering), cd), F(case["tau"])
    require(ptotal == ctotal == target,
            f"{N,r,m}: exact primal/dual objective mismatch: {ptotal}, {ctotal}, {target}")
    push_summary = None
    if m:
        pushed = [0] * len(words)
        for i in inside:
            y = canonical_map(words[i], pats, m)
            yi = encode(y)
            require(delta[yi] <= m-1, "map image outside smaller union")
            require(distance(words[i], y) <= 1, "map moves too far")
            require(delta[i] >= m or y == words[i], "map fails to fix inner union")
            pushed[yi] += packing[i]
        require(sum(pushed) == sum(packing), "pushforward loses mass")
        require(all(not v or delta[i] <= m-1 for i, v in enumerate(pushed)),
                "pushforward places mass outside target")
        push_loads = ball_sums(pushed, N, m-1)
        require(F(max(push_loads), pd) <= F(1), "pushforward infeasible")
        push_summary = {"mass": str(ptotal), "max_ball_load": str(F(max(push_loads), pd)),
                        "centres_checked": len(words)}
    return {"N": N, "r": r, "m": m, "tau": str(target),
            "status": "VERIFIED-COMPUTATIONALLY", "words_in_union": len(inside),
            "centres_checked": len(words), "max_packing_load": str(F(max(loads), pd)),
            "min_cover": str(F(min(coverage[i] for i in inside), cd)),
            "pushforward": push_summary}


def check_monotonicity(values):
    comparisons = []
    for N in range(3, 7):
        for r in range(1, N+2):
            row = [F(values[N, r, m]) for m in range(N+1)]
            require(row[0] == r, f"{N,r}: zero-radius endpoint wrong")
            require(row[-1] == 1, f"{N,r}: saturated endpoint wrong")
            for m in range(1, N+1):
                require(row[m-1] >= row[m],
                        f"MONOTONICITY FAIL {N=}, {r=}, {m=}: {row[m-1]} < {row[m]}")
                comparisons.append({"N": N, "r": r, "m": m,
                                    "decrease": str(row[m-1] - row[m])})
    require(len(comparisons) == 104, "wrong number of monotonicity comparisons")
    return comparisons


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cert", type=Path, default=ROOT / "certificates.json")
    parser.add_argument("--table", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    kernel = check_kernel()
    data = json.loads(args.cert.read_text())
    cases = data["cases"]
    keys = [(c["N"], c["r"], c["m"]) for c in cases]
    expected = {(N, r, m) for N in range(3, 7) for r in range(1, N+2) for m in range(N+1)}
    require(len(keys) == len(set(keys)) == 126 and set(keys) == expected,
            "certificate grid is incomplete or duplicated")
    if args.table:
        with args.table.open() as f:
            rows = list(csv.DictReader(f))
        vals = {(int(e["N"]), int(e["r"]), int(e["m"])): F(e["tau_exact"]) for e in rows}
        require(len(vals) == len(rows) == 126 and set(vals) == expected, "table grid wrong")
        # Test inequalities before checking correspondence to the certificates,
        # so an injected increase is diagnosed by the actual monotonicity check.
        comparisons = check_monotonicity(vals)
        require(all(vals[k] == F(c["tau"]) for k, c in zip(keys, cases)),
                "table differs from certified optima")
    summaries = []
    for N in range(3, 7):
        words = list(product(range(4), repeat=N))
        for c in cases:
            if c["N"] == N:
                summary = check_case(c, words)
                summaries.append(summary)
                print(f"EXACT PASS N={N} r={c['r']} m={c['m']} tau={c['tau']}", flush=True)
    values = {(c["N"], c["r"], c["m"]): F(c["tau"]) for c in cases}
    comparisons = check_monotonicity(values)
    report = {"status": "VERIFIED-COMPUTATIONALLY", "cells": len(cases),
              "rows": 22, "monotonicity_comparisons": len(comparisons),
              "strict_decreases": sum(F(c["decrease"]) > 0 for c in comparisons),
              "equalities": sum(F(c["decrease"]) == 0 for c in comparisons),
              "violations": 0, "kernel_check": kernel,
              "all_centre_packing_checks": sum(c["centres_checked"] for c in summaries),
              "pushforwards_checked": sum(c["pushforward"] is not None for c in summaries),
              "pushforward_all_centre_checks": sum(c["pushforward"]["centres_checked"]
                  for c in summaries if c["pushforward"] is not None),
              "cases": summaries, "comparisons": comparisons,
              "elapsed_seconds": time.monotonic()-start}
    if args.report:
        require(args.report.resolve().is_relative_to(ROOT), "report must remain inside J4")
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"EXACT VERIFICATION PASS: 126 optima, 104 comparisons, "
          f"104 mass-preserving feasible pushforwards; {time.monotonic()-start:.3f}s", flush=True)


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, ValueError, KeyError) as error:
        print(f"REJECTED: {error}", file=sys.stderr, flush=True)
        sys.exit(1)
