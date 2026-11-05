#!/usr/bin/env python3
"""Checker for family_tightness.json. Recomputes from raw guide strings; rejects corruption."""
from __future__ import annotations
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from family_tightness import analyze, alternating, find_rminus1_cover  # noqa

def fail(msg):
    print("REJECT:", msg)
    return 1

def main(argv):
    path = Path(argv[1] if len(argv) > 1 else HERE / "family_tightness.json")
    data = json.loads(path.read_text())
    s = data["settled"]
    # Recompute ACA
    fresh = analyze("ACA", 1)
    p = s["minimal_witness"]
    for k in ("r", "|sound|", "Um_proper", "sound_gt_r", "|Um|", "space"):
        if fresh[k] != p[k]:
            return fail("ACA %s mismatch" % k)
    if not (p["sound_gt_r"] and p["Um_proper"] and p["at_m_equals_n_minus_2"]):
        return fail("ACA predicates")
    ctrl = s["minimal_control"]
    fresh_c = analyze("ACA", 0)
    if fresh_c["|sound|"] != ctrl["|sound|"] or ctrl["|sound|"] != ctrl["r"]:
        return fail("ACA control")
    # AACCG secondary
    a = analyze("AACCG", 3)
    sec = s["secondary_AACCG"]
    if a["|sound|"] != sec["|sound|"] or not sec["sound_gt_r"]:
        return fail("AACCG")
    # Family flags
    if not s.get("family_alternating_AC_n_3_to_6_all_outside_sound_gt_r"):
        return fail("family outside flag")
    if not s.get("family_alternating_AC_n_3_to_6_all_inside_sound_eq_r"):
        return fail("family inside flag")
    # Recompute each family member
    for mem in data["family_rows"]:
        g = mem["outside"]["guide"]
        n = mem["outside"]["n"]
        fo = analyze(g, n - 2)
        if not (fo["sound_gt_r"] and fo["Um_proper"]):
            return fail("family member outside failed n=%d" % n)
        fi = analyze(g, n - 3)
        if fi["|sound|"] != fi["r"]:
            return fail("family member inside failed n=%d" % n)
        if g != alternating(n):
            return fail("guide not alternating")
    # Cover
    cover = s.get("minimal_rminus1_cover") or []
    if len(cover) != p["r"] - 1:
        return fail("cover size")
    # verify cover covers Um
    from family_tightness import ball, dels
    words = dels("ACA")
    Um = set()
    for w in words:
        Um |= ball(w, 1)
    u = set()
    for c in cover:
        u |= ball(c, 1)
    if not (Um <= u):
        return fail("cover does not cover Um")
    if data.get("status") != "PASS":
        return fail("status")
    # Negatives must record infinite-family open
    negs = data.get("negatives") or []
    if not any(n.get("result", "").startswith("NEGATIVE") for n in negs):
        return fail("missing negative for infinite family")
    print("ACCEPT: m<=n-3 tight; ACA+AACCG; family n=3..6 alternating; cover OK; infinite open recorded")
    return 0

if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
