#!/usr/bin/env python3
"""Plant corrupted private-witness outputs; checker must reject each."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable


def main():
    raw = HERE / "private_witnesses_m4.json"
    if not raw.exists():
        subprocess.check_call(
            [PY, "-I", "-B", str(HERE / "private_witnesses.py"), "--budget", "4", "--out", str(raw)],
            cwd=HERE,
        )
    data = json.loads(raw.read_text())
    plants = []
    # 1. flip a witness letter so it leaves the omitted ball
    d1 = json.loads(json.dumps(data))
    row = d1["rows"][0]
    w = row["witnesses"][0]["witness"]
    row["witnesses"][0]["witness"] = ("T" if w[0] != "T" else "A") + w[1:]
    plants.append(("corrupt_witness_ball", d1))
    # 2. claim all witnessed while zeroing one witness
    d2 = json.loads(json.dumps(data))
    d2["rows"][0]["witnesses"][0]["witness"] = None
    # leave counts claiming success
    plants.append(("drop_witness_keep_counts", d2))
    # 3. inflate fully_witnessed in summary
    d3 = json.loads(json.dumps(data))
    d3["summary"]["fully_witnessed"] = d3["summary"]["guides"] + 5
    plants.append(("inflate_summary", d3))

    # 4. duplicate first guide across the panel (identity collapse)
    d4 = json.loads(json.dumps(data))
    d4["rows"] = [json.loads(json.dumps(d4["rows"][0])) for _ in d4["rows"]]
    plants.append(("repeat_first_guide", d4))
    # 5. append a base to every witness (length corruption)
    d5 = json.loads(json.dumps(data))
    for row in d5["rows"]:
        for w in row["witnesses"]:
            if w.get("witness"):
                w["witness"] = w["witness"] + "A"
    plants.append(("append_base_to_every_witness", d5))
    # 6. repeat first omission index for every slot
    d6 = json.loads(json.dumps(data))
    for row in d6["rows"]:
        first = row["witnesses"][0]
        row["witnesses"] = [json.loads(json.dumps(first)) for _ in row["witnesses"]]
    plants.append(("repeat_first_omission", d6))


    bad = 0
    for name, payload in plants:
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump(payload, fh)
            path = fh.name
        proc = subprocess.run(
            [PY, "-I", "-B", str(HERE / "verify_private_witnesses.py"), path],
            cwd=HERE,
            capture_output=True,
            text=True,
        )
        ok = proc.returncode != 0
        print("%-28s exit=%s %s" % (name, proc.returncode, "REJECTED" if ok else "ACCEPTED_BUG"))
        if not ok:
            bad += 1
            print(proc.stdout, proc.stderr)
    print("%d plants, %d unexpected accepts" % (len(plants), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
