#!/usr/local/bin/python3
"""Exhaustive integer audit of the map and every centre/point inclusion pair.

This audit uses explicit Hamming distances and no LP certificates or transforms.
All generated files stay beside the script. Run with /usr/local/bin/python3 -B.
"""
from itertools import product
import json
from pathlib import Path
import time

import numpy as np

ROOT = Path(__file__).resolve().parent


def encode(w):
    i = 0
    for a in w:
        i = 4*i + int(a)
    return i


def main():
    start = time.monotonic()
    rows = []
    for N in range(3, 7):
        words = np.array(list(product(range(4), repeat=N)), dtype=np.uint8)
        D = np.zeros((len(words), len(words)), dtype=np.uint8)
        for k in range(N):
            D += words[:, k, None] != words[None, :, k]
        for r in range(1, N+2):
            L, C = r-1, N-r+1
            pats = np.array([encode([0]*j + [1]*(L-j) + [0]*C) for j in range(r)])
            ds = D[pats]
            delta = ds.min(axis=0)
            for m in range(1, N+1):
                inside = np.flatnonzero(delta <= m)
                image = inside.copy()
                shell_positions = np.flatnonzero(delta[inside] == m)
                tied_shell_sources = 0
                for a in shell_positions:
                    x = int(inside[a])
                    nearest = np.flatnonzero(ds[:, x] == m)
                    tied_shell_sources += len(nearest) >= 2
                    p = int(pats[nearest[0]])
                    k = int(np.flatnonzero(words[x] != words[p])[0])
                    image[a] = x + (int(words[p, k])-int(words[x, k])) * 4**(N-1-k)
                assert np.all(delta[image] <= m-1)
                assert np.all(D[inside, image] <= 1)
                inner = delta[inside] <= m-1
                assert np.all(image[inner] == inside[inner])
                # Test every raw centre, including every 2/3 letter assignment.
                violations = 0
                for begin in range(0, len(words), 128):
                    is_preimage = D[begin:begin+128, image] <= m-1
                    outside = D[begin:begin+128, inside] > m
                    violations += int(np.count_nonzero(is_preimage & outside))
                assert violations == 0
                row = {"N": N, "r": r, "m": m, "C": C,
                       "source_points": len(inside), "shell_points": len(shell_positions),
                       "tied_shell_sources": int(tied_shell_sources),
                       "all_raw_centres": len(words),
                       "centre_point_pairs": len(words) * len(inside),
                       "pairs_with_centres_containing_2_or_3": (4**N-2**N)*len(inside),
                       "inclusion_violations": violations,
                       "status": "VERIFIED-COMPUTATIONALLY"}
                rows.append(row)
                print(f"MAP PASS {N=} {r=} {m=} {C=} "
                      f"pairs={row['centre_point_pairs']} ties={tied_shell_sources}", flush=True)
    report = {"status": "VERIFIED-COMPUTATIONALLY", "map_instances": len(rows),
              "source_points_checked": sum(r["source_points"] for r in rows),
              "centre_point_pairs": sum(r["centre_point_pairs"] for r in rows),
              "pairs_with_centres_containing_2_or_3": sum(
                  r["pairs_with_centres_containing_2_or_3"] for r in rows),
              "tied_shell_sources": sum(r["tied_shell_sources"] for r in rows),
              "C_zero_map_instances": sum(r["C"] == 0 for r in rows),
              "C_zero_tied_shell_sources": sum(r["tied_shell_sources"] for r in rows if r["C"] == 0),
              "inclusion_violations": 0, "cases": rows,
              "elapsed_seconds": time.monotonic()-start}
    assert len(rows) == 104
    assert report["tied_shell_sources"] > 0 and report["C_zero_tied_shell_sources"] > 0
    (ROOT / "map_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "cases"}), flush=True)


if __name__ == "__main__":
    main()
