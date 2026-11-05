#!/usr/bin/env python3
"""Producer: private witnesses for every single omission of a deletion word.

A private witness for omitting deletion-word i at budget m is a DNA word at
Hamming distance <= m from that word and at distance > m from every other
deletion word of the same guide. Its existence proves that particular
(r-1)-subset of the classical library does not cover Um.

This is a necessary condition for K = r among library subsets. The paper's K
ranges over unrestricted centres; that is Theorems 5-7. Report negatives.

    python3 -I private_witnesses.py --budget 4 [--out PATH]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def runs(g: str) -> int:
    return 1 + sum(g[i] != g[i + 1] for i in range(len(g) - 1))


def deletion_words_ordered(g: str):
    words = []
    i = 0
    n = len(g)
    while i < n:
        j = i
        while j + 1 < n and g[j + 1] == g[i]:
            j += 1
        words.append(g[:i] + g[i + 1 :])
        i = j + 1
    return words


def ham(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def private_witness(words, omit_idx: int, m: int):
    T = words[omit_idx]
    others = [words[i] for i in range(len(words)) if i != omit_idx]
    if all(ham(T, o) > m for o in others):
        return T, "self"
    alphabet = "ACGT"
    L = len(T)
    uni = [p for p in range(L) if all(o[p] == T[p] for o in others)]
    dist0 = [ham(T, o) for o in others]
    if len(uni) >= m:
        w = list(T)
        for p in uni[:m]:
            for ch in alphabet:
                if ch != T[p]:
                    w[p] = ch
                    break
        wit = "".join(w)
        if ham(wit, T) <= m and all(ham(wit, o) > m for o in others):
            return wit, "unanimous"
    w = list(T)
    changed = set()
    dists = dist0[:]
    for _ in range(m):
        best = None
        for p in range(L):
            if p in changed:
                continue
            for ch in alphabet:
                if ch == w[p]:
                    continue
                nd = []
                for i, o in enumerate(others):
                    d = dists[i]
                    if o[p] == w[p]:
                        d += 1
                    if o[p] == ch:
                        d -= 1
                    nd.append(d)
                bad = sum(1 for d in nd if d <= m)
                key = (bad, -min(nd), p, ch, nd)
                if best is None or key < best:
                    best = key
        if best is None:
            break
        bad, _, p, ch, nd = best
        w[p] = ch
        changed.add(p)
        dists = nd
        if bad == 0:
            return "".join(w), "greedy"
    wit = "".join(w)
    if ham(wit, T) <= m and all(ham(wit, o) > m for o in others):
        return wit, "greedy"
    return None, "none"


def verify_witness(words, omit_idx, wit, m) -> bool:
    T = words[omit_idx]
    if ham(wit, T) > m:
        return False
    return all(ham(wit, words[i]) > m for i in range(len(words)) if i != omit_idx)


def load_guides():
    bench_path = ROOT / "exactness_dividend" / "published_benchmark_guides.txt"
    bench = [
        (f"b{i}", g.strip().upper(), "benchmark")
        for i, g in enumerate(bench_path.read_text().splitlines())
        if g.strip()
    ]
    panel = []
    seen = set()
    with (ROOT / "tool_conformance" / "per_guide_counts.tsv").open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row["tool"].startswith("Cas-OFFinder-bulge"):
                continue
            g = row["guide"].upper()
            if g in seen:
                continue
            seen.add(g)
            panel.append((row["guide_id"], g, "panel"))
    return bench + panel


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", type=int, default=4)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)
    m = args.budget
    if m < 0:
        raise SystemExit("REJECT: budget must be a nonnegative integer, got %r" % (m,))
    rows = []
    miss = 0
    for gid, g, src in load_guides():
        words = deletion_words_ordered(g)
        r = len(words)
        if r != runs(g):
            raise SystemExit("run/word mismatch on %s" % gid)
        ok = 0
        wits = []
        for i in range(r):
            wit, how = private_witness(words, i, m)
            if wit and verify_witness(words, i, wit, m):
                ok += 1
                wits.append({"omit": i, "witness": wit, "method": how})
            else:
                miss += 1
                wits.append({"omit": i, "witness": None, "method": how})
        rows.append(
            {
                "source": src,
                "guide_id": gid,
                "guide": g,
                "n": len(g),
                "r": r,
                "m": m,
                "omissions_with_private_witness": ok,
                "all_omissions_witnessed": ok == r,
                "witnesses": wits,
            }
        )
    summary = {
        "definition": (
            "Private witness: DNA word within distance m of an omitted deletion "
            "word and >m from every other deletion word. Shows that (r-1)-subset "
            "of the classical library fails to cover Um. Not a proof of K over "
            "unrestricted centres."
        ),
        "m": m,
        "guides": len(rows),
        "fully_witnessed": sum(1 for r in rows if r["all_omissions_witnessed"]),
        "miss_events": miss,
        "status": "PASS" if miss == 0 else "NEGATIVE_PARTIAL",
    }
    out = args.out or (HERE / "private_witnesses_m%d.json" % m)
    out.write_text(json.dumps({"summary": summary, "rows": rows}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("wrote", out)
    return 0 if miss == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
