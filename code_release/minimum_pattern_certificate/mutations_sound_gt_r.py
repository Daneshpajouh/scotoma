#!/usr/bin/env python3
"""Plant corruptions of sound_gt_r_witness.json; checker must reject."""
from __future__ import annotations
import json, subprocess, sys, tempfile
from pathlib import Path
HERE = Path(__file__).resolve().parent
PY = sys.executable

def main():
    src = HERE / "sound_gt_r_witness.json"
    if not src.exists():
        subprocess.check_call([PY, "-I", "-B", str(HERE / "sound_gt_r_witness.py")], cwd=HERE)
    data = json.loads(src.read_text())
    plants = []
    d1 = json.loads(json.dumps(data)); d1["primary_witness"]["|sound|"] = d1["primary_witness"]["r"]; plants.append(("drop_sound_gt", d1))
    d2 = json.loads(json.dumps(data)); d2["primary_witness"]["explicit_rminus1_cover"] = ["AAAA"]; plants.append(("bad_cover", d2))
    d3 = json.loads(json.dumps(data)); d3["primary_witness"]["m"] = 2; plants.append(("move_inside_theorem", d3))
    bad = 0
    for name, payload in plants:
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump(payload, fh); path = fh.name
        proc = subprocess.run([PY, "-I", "-B", str(HERE / "verify_sound_gt_r_witness.py"), path], capture_output=True, text=True)
        ok = proc.returncode != 0
        print("%-24s exit=%s %s" % (name, proc.returncode, "REJECTED" if ok else "ACCEPTED_BUG"))
        bad += not ok
    print("%d plants, %d unexpected" % (len(plants), bad))
    return 1 if bad else 0
if __name__ == "__main__":
    sys.exit(main())
