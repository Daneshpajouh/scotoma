#!/usr/bin/env python3
"""Plant corruptions; checker must reject each."""
from __future__ import annotations
import json, subprocess, sys, tempfile
from pathlib import Path
HERE = Path(__file__).resolve().parent
PY = sys.executable

def main():
    src = HERE / "family_tightness.json"
    data = json.loads(src.read_text())
    plants = []
    d = json.loads(json.dumps(data)); d["settled"]["minimal_witness"]["|sound|"] = 3; plants.append(("drop_sound", d))
    d = json.loads(json.dumps(data)); d["settled"]["family_alternating_AC_n_3_to_6_all_outside_sound_gt_r"] = False; plants.append(("family_flag", d))
    d = json.loads(json.dumps(data)); d["settled"]["minimal_rminus1_cover"] = ["AA"]; plants.append(("bad_cover", d))
    d = json.loads(json.dumps(data)); d["settled"]["secondary_AACCG"]["|sound|"] = 3; d["settled"]["secondary_AACCG"]["sound_gt_r"] = False; plants.append(("aaccg_drop", d))
    d = json.loads(json.dumps(data)); d["negatives"] = []; plants.append(("drop_negatives", d))
    d = json.loads(json.dumps(data)); d["status"] = "FAIL"; plants.append(("status_fail", d))
    bad = 0
    for name, payload in plants:
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump(payload, fh); path = fh.name
        proc = subprocess.run([PY, "-I", "-B", str(HERE / "verify_family_tightness.py"), path], capture_output=True, text=True)
        ok = proc.returncode != 0
        print("%-16s exit=%s %s" % (name, proc.returncode, "REJECTED" if ok else "ACCEPTED_BUG"))
        if not ok:
            print(proc.stdout)
            bad += 1
    print("%d plants, %d unexpected" % (len(plants), bad))
    return 1 if bad else 0

if __name__ == "__main__":
    raise SystemExit(main())
