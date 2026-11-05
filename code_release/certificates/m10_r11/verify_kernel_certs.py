#!/usr/local/bin/python3
"""Integer verification of the kernel-cover certificates.

Every comparison uses Python integers. A nonzero exit status means the
certificate does not prove that every nonnegative integer solution has
sum at least 11.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    A = np.load(HERE / "A_reduced.npy")
    rhs = [int(v) for v in np.load(HERE / "rhs_reduced.npy").tolist()]
    pack = np.load(HERE / "pack_matrix.npy")
    DPROB = json.loads((HERE / "threshold.json").read_text())["dprob"]
    bad = 0
    for r, b in enumerate(rhs):
        g = DPROB // b
        if g * b != DPROB:
            bad += 1
            continue
        if not np.array_equal(A[r].astype(np.int64) * np.int64(g), pack[r].astype(np.int64)):
            bad += 1
    if bad:
        raise SystemExit(f"reduced matrix does not match pack_matrix on {bad} rows")

    index = json.loads((HERE / "cert_index.json").read_text())
    D = int(index["D"])
    rows, cols = A.shape
    for cert in index["certs"]:
        z = [int(v) for v in np.load(HERE / f"cert_{cert['name']}_z.npy").tolist()]
        if len(z) != rows or any(v < 0 for v in z):
            raise SystemExit(f"{cert['name']}: bad z")
        banned = set(cert["banned"])
        subtracted = cert["subtracted"]
        acc = [0] * cols
        for r, zr in enumerate(z):
            if zr == 0:
                continue
            row = A[r]
            for c in range(cols):
                a = int(row[c])
                if a:
                    acc[c] += zr * a
        worst = max(acc[c] for c in range(cols) if c not in banned)
        if worst > D:
            raise SystemExit(f"{cert['name']}: column product {worst} exceeds {D}")
        dot = 0
        for r, zr in enumerate(z):
            if zr == 0:
                continue
            gap = rhs[r]
            for j in subtracted:
                gap -= int(A[r, j])
            dot += zr * gap
        multiple = int(cert["multiple"])
        if dot != int(cert["signed_dot"]):
            raise SystemExit(f"{cert['name']}: signed dot {dot} != stored {cert['signed_dot']}")
        if dot <= multiple * D:
            raise SystemExit(f"{cert['name']}: dot {dot} does not exceed {multiple}*D")
        # Free variables must exceed `multiple`, so an integer free sum is
        # at least multiple+1. Adding the subtracted columns gives the total.
        total = (multiple + 1) + len(subtracted)
        if total < 11:
            raise SystemExit(f"{cert['name']}: bound {total} is below 11")
        print(f"{cert['name']}: worst {worst} <= {D}, dot {dot} > {multiple}*D, integer sum >= {total}")
    print("KERNEL_CERTS_OK")


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    main()
