#!/usr/bin/env python3
"""Checker for private_witnesses.py output. Rejects corrupted witnesses.

    python3 -I private_witnesses.py --budget 4 --out /tmp/pw.json
    python3 -I verify_private_witnesses.py /tmp/pw.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from private_witnesses import (  # noqa: E402
    deletion_words_ordered,
    ham,
    runs,
    verify_witness,
)

ALPHA = set("ACGT")


def fail(msg):
    print("REJECT:", msg)
    return 1


def main(argv):
    if len(argv) < 2:
        print("usage: verify_private_witnesses.py OUT.json", file=sys.stderr)
        return 2
    data = json.loads(Path(argv[1]).read_text())
    summary = data["summary"]
    rows = data["rows"]
    if summary.get("m") is None or summary.get("guides") != len(rows):
        return fail("summary/rows length mismatch")
    m = summary["m"]
    bad = 0
    seen_ids = set()
    seen_guides = set()
    fully = 0
    for row in rows:
        gid = row["guide_id"]
        g = row["guide"]
        if gid in seen_ids or g in seen_guides:
            print("REJECT: duplicate guide identity", gid)
            bad += 1
            continue
        seen_ids.add(gid)
        seen_guides.add(g)
        words = deletion_words_ordered(g)
        if len(words) != row["r"] or runs(g) != row["r"]:
            print("REJECT: r mismatch", gid)
            bad += 1
            continue
        if len(row["witnesses"]) != row["r"]:
            print("REJECT: witness list length", gid)
            bad += 1
            continue
        omits = [wrec["omit"] for wrec in row["witnesses"]]
        if sorted(omits) != list(range(row["r"])):
            print("REJECT: omit indices not exactly 0..r-1 once each", gid)
            bad += 1
            continue
        ok = 0
        for wrec in row["witnesses"]:
            wit = wrec.get("witness")
            if wit is None:
                continue
            if not isinstance(wit, str) or len(wit) != len(g) - 1 or set(wit) - ALPHA:
                print("REJECT: witness alphabet/length", gid, wrec["omit"])
                bad += 1
                continue
            if not verify_witness(words, wrec["omit"], wit, m):
                print("REJECT: forged witness", gid, wrec["omit"])
                bad += 1
            else:
                ok += 1
        if ok != row["omissions_with_private_witness"]:
            print("REJECT: count corruption", gid, ok, row["omissions_with_private_witness"])
            bad += 1
        if (ok == row["r"]) != bool(row["all_omissions_witnessed"]):
            print("REJECT: flag corruption", gid)
            bad += 1
        if ok == row["r"]:
            fully += 1
    if summary.get("fully_witnessed") != fully:
        return fail("summary fully_witnessed %s != recomputed %s" % (summary.get("fully_witnessed"), fully))
    # planted expectation on deposited panels at m=4
    if m == 4 and summary.get("status") == "PASS":
        if summary["fully_witnessed"] != summary["guides"]:
            return fail("PASS status but not all witnessed")
        if summary["guides"] not in (171,):
            return fail("unexpected guide count %s" % summary["guides"])
    print(
        "ACCEPT: %d guides, %d fully witnessed, %d reject events"
        % (summary["guides"], summary["fully_witnessed"], bad)
    )
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
