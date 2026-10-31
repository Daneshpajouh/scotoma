"""Checker for resolved near-cut packings.

The checker is written from the shell definitions. It recomputes every shell
size, rejects a shell whose edit count is not m, rebuilds the integer
coefficients, and reruns the exhaustive load program. Acceptance of a closed
row requires the computed maximum numerator to equal the stored normalizer
and the mass at that normalizer to exceed r-1. A frontier row is accepted
when the same identities hold and the mass does not exceed r-1; that is a
measurement, not a closure.

The planted-overload control doubles every positive weight and keeps the
stored normalizer. The outside-support control adds the word of all twos,
which lies outside every radius-m ball of a binary pattern when N > m.
Both must be rejected.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
from fractions import Fraction
from math import lcm
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from support import layout, shell_size

ROOT = Path(__file__).resolve().parent
BIN = ROOT / "dp_resolved"
LOG = ROOT / "logs"
CERT = ROOT


def edit_count(row: dict) -> int:
    return (
        sum(s != 0 for s in row["guard"])
        + sum(s != 0 for s in row["resL"])
        + sum(s != 0 for s in row["resR"])
        + row["t"]
        + row["u"]
    )


def as_shell(row: dict):
    return (
        row["j"],
        tuple(row["guard"]),
        tuple(row["resL"]),
        tuple(row["resR"]),
        row["t"],
        row["u"],
        row["fl"],
        row["fr"],
    )


def hamming(a: str, b: list[int]) -> int:
    return sum(int(x) != y for x, y in zip(a, b))


def pattern(j: int, N: int) -> list[int]:
    return [1 if i < j else 0 for i in range(N)]


def check_shape(data: dict, rows: list[dict]) -> list[str]:
    errors = []
    N, m, depth = data["N"], data["m"], data["depth"]
    if data["free_guard"] not in (False, True):
        errors.append("free_guard flag missing")
    for row in rows:
        if not 0 <= row["j"] < data["r"] or row["a"] < 0:
            errors.append("invalid anchor or negative shell weight")
        lay = layout(row["j"], N, depth)
        if list(row["guard"]) != [0] * len(lay["guard"]) and not data["free_guard"]:
            errors.append(f"guard edited at j={row['j']} while the guard is fixed")
        if len(row["guard"]) != len(lay["guard"]) or len(row["resL"]) != len(lay["resL"]) or len(row["resR"]) != len(lay["resR"]):
            errors.append(f"layout mismatch at j={row['j']}")
        if any(s not in (0, 1, 2) for s in row["guard"] + row["resL"] + row["resR"]):
            errors.append(f"bad symbol at j={row['j']}")
        if edit_count(row) != m:
            errors.append(f"edit count {edit_count(row)} is not budget {m}")
        if row["fl"] > row["t"] or row["fr"] > row["u"] or min(row["t"], row["u"], row["fl"], row["fr"]) < 0:
            errors.append(f"bad far counts at j={row['j']}")
        if row["fl"] > data["flip_cap"] or row["fr"] > data["flip_cap"]:
            errors.append(f"flip cap exceeded at j={row['j']}")
        size = shell_size(N, m, depth, as_shell(row))
        if size <= 0:
            errors.append(f"empty shell at j={row['j']}")
        elif data["H"] and int(data["H"]) % size != 0:
            errors.append(f"denominator {data['H']} is not a multiple of size {size}")
    for item in data.get("outside_words", []):
        word = item["word"]
        if item.get("a", 0) <= 0:
            continue
        dist = min(hamming(word, pattern(j, data["N"])) for j in range(data["r"]))
        if dist > data["m"]:
            errors.append(f"outside support: min distance {dist} exceeds budget {data['m']}")
    return errors


def run_dp(data: dict, rows: list[dict]) -> dict:
    N = data["N"]
    lines = [f"{N} {data['m']} {data['depth']} {data['flip_cap']} {int(data['free_guard'])} {len(rows)}"]
    H = int(data["H"])
    coeffs = []
    for row in rows:
        size = shell_size(N, data["m"], data["depth"], as_shell(row))
        coeff = row["a"] * (H // size)
        coeffs.append(coeff)
        g, L, R = row["guard"], row["resL"], row["resR"]
        parts = [str(row["j"]), str(len(g)), str(len(L)), str(len(R))]
        parts += [str(s) for s in g + L + R]
        parts += [str(row["t"]), str(row["u"]), str(row["fl"]), str(row["fr"]), str(coeff)]
        lines.append(" ".join(parts))
    path = LOG / "verify.in"
    path.write_text("\n".join(lines) + "\n")
    env = os.environ.copy()
    env["DP_THREADS"] = "6"
    proc = subprocess.run([str(BIN), "max", str(path)], capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"dp exit {proc.returncode}")
    out = {}
    for line in proc.stdout.splitlines():
        parts = line.split()
        if parts and parts[0] in {"maximum", "witness", "classes", "seconds"}:
            out[parts[0]] = parts[1]
    out["coeffs"] = coeffs
    return out


def check_file(path: Path, require_closed: bool | None = None) -> dict:
    data = json.loads(path.read_text())
    rows = data["masses"]
    errors = check_shape(data, rows)
    if errors:
        return {"file": path.name, "passed": False, "errors": errors}
    H = 1
    for row in rows:
        H = lcm(H, shell_size(data["N"], data["m"], data["depth"], as_shell(row)))
    if H != int(data["H"]):
        errors.append(f"recomputed H {H} != stored {data['H']}")
    ans = run_dp(data, rows)
    maximum = int(ans["maximum"])
    normalizer = int(data["normalizer"])
    if maximum != normalizer:
        errors.append(f"computed maximum {maximum} != stored normalizer {normalizer}")
    if int(ans["classes"]) != 3 ** data["N"] or int(ans["classes"]) != data["classes"]:
        errors.append(f"class count {ans['classes']} != 3^{data['N']}")
    mass = Fraction(H * sum(row["a"] for row in rows), normalizer)
    stored = Fraction(data["total_mass"])
    if mass != stored:
        errors.append(f"recomputed mass {mass} != stored {stored}")
    closed = mass > data["r"] - 1
    if require_closed is True and not closed:
        errors.append(f"mass {mass} does not exceed {data['r'] - 1}")
    if require_closed is False and closed:
        errors.append("frontier row unexpectedly clears the target")
    if data["above_target"] != closed:
        errors.append("above_target flag disagrees with the mass")
    return {
        "file": path.name,
        "passed": not errors,
        "mass": str(mass),
        "max_load": f"{maximum}/{normalizer}",
        "closed": closed,
        "classes": int(ans["classes"]),
        "witness": ans["witness"],
        "errors": errors,
    }


def overload(data: dict) -> dict:
    forged = json.loads(json.dumps(data))
    for row in forged["masses"]:
        row["a"] *= 2
    path = LOG / "overload.json"
    path.write_text(json.dumps(forged) + "\n")
    return check_file(path)


def outside(data: dict) -> dict:
    forged = json.loads(json.dumps(data))
    forged["outside_words"] = [{"word": "2" * data["N"], "a": 1}]
    path = LOG / "outside.json"
    path.write_text(json.dumps(forged) + "\n")
    return check_file(path)


def main() -> None:
    try:
        os.nice(10)
    except OSError:
        pass
    files = sorted(CERT.glob("packing_N*_m10_*.json"))
    # Skip forged files if any were copied into certificates.
    files = [p for p in files if "overload" not in p.name and "outside" not in p.name]
    if {p.name for p in files} != {'packing_N22_r14_m10_d2_c2_all.json', 'packing_N24_r25_m10_d2_c1_all.json', 'packing_N22_r13_m10_d2_c2_all.json', 'packing_N23_r24_m10_d2_c1_all.json', 'packing_N22_r16_m10_d2_c2_all.json', 'packing_N22_r12_m10_d2_c2_all.json', 'packing_N26_r27_m10_d2_c0_all.json', 'packing_N22_r15_m10_d2_c2_all.json', 'packing_N25_r26_m10_d2_c0_all.json'}:
        raise RuntimeError("certificate set differs from the documented package")
    if not files:
        print(json.dumps({"passed": False, "errors": ["no certificates"]}))
        raise SystemExit(1)
    reports = []
    for p in files:
        print(json.dumps({"checking": p.name}), flush=True)
        reports.append(check_file(p, require_closed=True))
    base = json.loads(files[0].read_text())
    over = overload(base)
    out = outside(base)
    summary = {
        "certificates": len(reports),
        "certificates_passed": all(r["passed"] for r in reports),
        "overload_rejected": not over["passed"],
        "outside_rejected": not out["passed"],
        "overload_errors": over["errors"],
        "outside_errors": out["errors"],
    }
    summary["passed"] = (
        summary["certificates_passed"] and summary["overload_rejected"] and summary["outside_rejected"]
    )
    for row in reports:
        print(json.dumps(row), flush=True)
    print(json.dumps(summary), flush=True)
    print(json.dumps({"reports": reports}), flush=True)
    if not summary["passed"]:
        raise SystemExit(1)


def release_main():
    global LOG, BIN
    if sys.argv[1:] not in ([], ["verify"]):
        raise SystemExit("usage: python3 -I verify.py verify")
    with tempfile.TemporaryDirectory(prefix=".verify_", dir=ROOT) as work:
        LOG = Path(work)
        env = os.environ.copy()
        env["TMPDIR"] = work
        BIN = LOG / "dp_resolved"
        subprocess.run([os.environ.get("CXX", "c++"), "-O3", "-std=c++17", "-pthread",
                        str(ROOT / "dp_resolved.cpp"), "-o", str(BIN)], check=True, env=env)
        from calibrate_precompute import run as calibrate
        calibrate(BIN, LOG)
        return main()

if __name__ == "__main__":
    raise SystemExit(release_main())
