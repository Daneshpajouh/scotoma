#!/usr/bin/env python3
"""Checker for tightness_hypothesis.json. Rejects corruption."""
from __future__ import annotations
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from tightness_hypothesis import analyze_idx  # noqa: E402

def fail(msg):
    print("REJECT:", msg)
    return 1

def main(argv):
    path = Path(argv[1] if len(argv) > 1 else HERE / "tightness_hypothesis.json")
    data = json.loads(path.read_text())
    p = data["primary_minimal_witness"]
    fresh = analyze_idx(p["guide"], p["m"])
    for k in ("r", "|sound|", "Um_proper", "sound_gt_r", "outside_companion_soundness_range", "|Um|", "space"):
        if fresh[k] != p[k]:
            return fail("primary %s mismatch fresh=%s file=%s" % (k, fresh[k], p[k]))
    if p["guide"] != "ACA" or p["n"] != 3 or p["m"] != 1:
        return fail("primary is not the minimal ACA witness")
    if not (p["sound_gt_r"] and p["Um_proper"] and p["outside_companion_soundness_range"]):
        return fail("primary fails tightness predicates")
    ctrl = data["control_at_theorem_boundary"]
    if ctrl["|sound|"] != ctrl["r"] or ctrl["sound_gt_r"] or ctrl["m"] != 0:
        return fail("boundary control corrupted")
    if data.get("inside_boundary_failures_count", 1) != 0:
        return fail("inside-boundary failures nonempty")
    cover = data.get("primary_explicit_rminus1_cover") or []
    if len(cover) != p["r"] - 1:
        return fail("cover size")
    # verify cover
    from tightness_hypothesis import ball_indices, encode, deletion_words
    words = deletion_words(p["guide"])
    L = p["n"] - 1
    Um = set()
    for w in words:
        Um |= ball_indices(encode(w), L, p["m"])
    u = set()
    for c in cover:
        u |= ball_indices(encode(c), L, p["m"])
    if not (Um <= u):
        return fail("claimed r-1 cover does not cover Um")
    if data.get("status") != "PASS":
        return fail("status")
    # census: every inside cell has 0 sound_gt_r
    for cell in data.get("census", []):
        if cell["m_mode"].startswith("inside") and cell["sound_gt_r_count"] != 0:
            return fail("census inside cell has sound_gt_r")
    print("ACCEPT: hypothesis m<=n-3 is tight; ACA witness; 0 inside failures; r-1 cover OK")
    return 0
if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
