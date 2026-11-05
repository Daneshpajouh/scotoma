#!/usr/bin/env python3
"""Plant corruptions of tightness_hypothesis.json; checker must reject."""
from __future__ import annotations
import json, subprocess, sys, tempfile
from pathlib import Path
HERE = Path(__file__).resolve().parent
PYEXE = sys.executable

def main():
    src = HERE / "tightness_hypothesis.json"
    data = json.loads(src.read_text())
    plants = []
    d1 = json.loads(json.dumps(data)); d1["primary_minimal_witness"]["|sound|"] = 3; plants.append(("drop_sound_gt", d1))
    d2 = json.loads(json.dumps(data)); d2["primary_minimal_witness"]["guide"] = "AACCG"; d2["primary_minimal_witness"]["n"] = 5; d2["primary_minimal_witness"]["m"] = 3; plants.append(("not_minimal_ACA", d2))
    d3 = json.loads(json.dumps(data)); d3["inside_boundary_failures_count"] = 2; plants.append(("fake_inside_fail", d3))
    d4 = json.loads(json.dumps(data)); d4["primary_explicit_rminus1_cover"] = ["AA"]; plants.append(("bad_cover", d4))
    bad = 0
    for name, payload in plants:
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump(payload, fh); path = fh.name
        proc = subprocess.run([PYEXE, "-I", "-B", str(HERE / "verify_tightness_hypothesis.py"), path], capture_output=True, text=True)
        ok = proc.returncode != 0
        print("%-20s exit=%s %s" % (name, proc.returncode, "REJECTED" if ok else "ACCEPTED_BUG"))
        bad += not ok
    print("%d plants, %d unexpected" % (len(plants), bad))
    return 1 if bad else 0
if __name__ == "__main__":
    raise SystemExit(main())
