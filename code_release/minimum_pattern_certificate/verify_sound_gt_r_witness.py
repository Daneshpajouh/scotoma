#!/usr/bin/env python3
"""Checker for sound_gt_r_witness.json. Rejects corruption."""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from sound_gt_r_witness import analyze, ball, deletion_words, runs  # noqa: E402


def fail(msg):
    print("REJECT:", msg)
    return 1


def main(argv):
    path = Path(argv[1] if len(argv) > 1 else HERE / "sound_gt_r_witness.json")
    data = json.loads(path.read_text())
    p = data["primary_witness"]
    g, m = p["guide"], p["m"]
    fresh = analyze(g, m)
    for k in ("r", "|sound|", "Um_proper", "sound_gt_r", "outside_companion_soundness_range"):
        if fresh[k] != p[k]:
            return fail("primary field %s mismatch fresh=%s file=%s" % (k, fresh[k], p[k]))
    if not p["sound_gt_r"] or not p["outside_companion_soundness_range"] or not p["Um_proper"]:
        return fail("primary fails breakthrough predicates")
    # verify cover really covers
    words = deletion_words(g)
    Um = set()
    for w in words:
        Um |= ball(w, m)
    cover = p.get("explicit_rminus1_cover") or []
    if len(cover) != p["r"] - 1:
        return fail("cover size != r-1")
    u = set()
    for c in cover:
        u |= ball(c, m)
    if not (Um <= u):
        return fail("claimed r-1 cover does not cover Um")
    # control inside theorem
    ctrl = data.get("control_inside_theorem_AACCG_m2")
    if ctrl:
        if ctrl["|sound|"] != ctrl["r"] or ctrl["outside_companion_soundness_range"]:
            return fail("inside-theorem control corrupted")
    if data.get("status") != "PASS":
        return fail("status not PASS")
    print("ACCEPT: |sound|>r witness outside m<=n-3; Um proper; r-1 cover verified; control OK")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
