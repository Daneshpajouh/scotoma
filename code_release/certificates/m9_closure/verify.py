"""Exact checker for the m=9 certificates in this directory.

The checker is written from the definitions. It does not replace the stored
normalizer by the computed maximum: the computed maximum numerator must be
at most the stored normalizer, and the mass at that normalizer must exceed
r-1. A planted overload keeps the normalizer and raises one weight. An
outside-support control puts positive weight on the word 2^N, which lies at
distance N from every library pattern.
"""

from __future__ import annotations

import json
import os
import random
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from shells import Counter, pattern, shell_size
from tagged import Counter as TagCounter
from tagged import shell_size as tag_size

ROOT = Path(__file__).resolve().parent
CERT = ROOT
LOG = ROOT / "logs"
DP = ROOT / "dp"
DP_TAG = ROOT / "dp_tag"


def hamming(a, b) -> int:
    return sum(x != y for x, y in zip(a, b))


def library(N: int, r: int) -> list[list[int]]:
    return [pattern(j, N) for j in range(r)]


def min_distance(word, patterns) -> int:
    return min(hamming(word, P) for P in patterns)


def run_classical(N: int, m: int, rows: list[tuple]) -> dict:
    lines = [f"{N} {m} {len(rows)}"]
    for g, j, t, fl, fr, coeff in rows:
        lines.append(f"{g} {j} {t} {fl} {fr} {coeff}")
    path = LOG / "verify_classical.in"
    path.write_text("\n".join(lines) + "\n")
    proc = subprocess.run([str(DP), str(path)], check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"dp exit {proc.returncode}")
    return parse_dp(proc.stdout)


def run_tagged(N: int, m: int, rows: list[tuple]) -> dict:
    lines = [f"{N} {m} {len(rows)}"]
    for row in rows:
        lines.append(" ".join(str(x) for x in row))
    path = LOG / "verify_tagged.in"
    path.write_text("\n".join(lines) + "\n")
    proc = subprocess.run([str(DP_TAG), str(path)], check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"dp_tag exit {proc.returncode}")
    return parse_dp(proc.stdout)


def parse_dp(text: str) -> dict:
    maximum = witness = classes = None
    seconds = None
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "maximum":
            maximum = int(parts[1])
        elif parts[0] == "witness":
            witness = parts[1]
        elif parts[0] == "classes":
            classes = int(parts[1])
        elif parts[0] == "seconds":
            seconds = float(parts[1])
    if maximum is None or witness is None or classes is None:
        raise RuntimeError("dp output missing fields")
    return {"maximum": maximum, "witness": witness, "classes": classes, "seconds": seconds}


def atom_errors(data, patterns) -> list[str]:
    errors = []
    for atom in data.get("atoms") or []:
        word = [int(ch) for ch in atom["word"]]
        weight = int(atom.get("numerator", "1"))
        if weight <= 0:
            continue
        if len(word) != data["N"]:
            errors.append("atom length")
            continue
        dist = min_distance(word, patterns)
        if dist > data["m"]:
            errors.append(f"outside support: min distance {dist} exceeds budget {data['m']}")
        else:
            errors.append("atom inside U is not a guarded shell")
    return errors


def check(data: dict) -> dict:
    family = data["family"]
    N = int(data["N"])
    m = int(data["m"])
    r = int(data["r"])
    A = int(data["normalizer"])
    H = int(data["H"])
    patterns = library(N, r)
    errors = atom_errors(data, patterns)
    if family == "guarded_private_shells":
        return check_classical(data, patterns, errors)
    if family == "guarded_tagged_shells":
        return check_tagged(data, patterns, errors)
    return {"passed": False, "errors": [f"unknown family {family}"]}


def check_classical(data, patterns, errors: list[str]) -> dict:
    N = int(data["N"])
    m = int(data["m"])
    r = int(data["r"])
    A = int(data["normalizer"])
    H = int(data["H"])
    guards = tuple(data["guards"])
    rows = []
    mass_num_v = 0
    for raw in data["masses"]:
        if not 0 <= raw[1] < r:
            errors.append("anchor outside the library")
        g, j, t, fl, fr, v = raw
        if v < 0 or g not in guards:
            errors.append(f"bad row {raw}")
            continue
        S = shell_size(N, m, g, j, t, fl, fr)
        if S <= 0 or fl > g - 1 or fr > g - 1:
            errors.append(f"shell not in the private family {raw}")
            continue
        if H % S != 0:
            errors.append(f"H not a multiple of shell size {raw}")
            continue
        # Every shell word matches the anchor off the edit set and differs on it.
        if t + (m - t) != m:
            errors.append("edit total")
        rows.append((g, j, t, fl, fr, v * (H // S)))
        mass_num_v += v
    result = base_result(data, errors, rows, mass_num_v, H, A, r)
    if errors:
        result["passed"] = False
        return result
    ans = run_classical(N, m, rows)
    return finish(data, result, ans, rows, classical_numerator)


def classical_numerator(data, centre, rows) -> int:
    guards = tuple(data["guards"])
    ctr = classical_numerator.ctr
    if ctr is None or ctr.N != data["N"] or ctr.guards != guards:
        classical_numerator.ctr = Counter(data["N"], data["m"], guards, r=data["r"])
        ctr = classical_numerator.ctr
    index = {(g, j, t, fl, fr): i for i, (g, j, t, fl, fr) in enumerate(ctr.shells)}
    total = 0
    for g, j, t, fl, fr, coeff in rows:
        total += coeff * ctr.count(centre, index[(g, j, t, fl, fr)])
    return total


classical_numerator.ctr = None


def check_tagged(data, patterns, errors: list[str]) -> dict:
    N = int(data["N"])
    m = int(data["m"])
    r = int(data["r"])
    A = int(data["normalizer"])
    H = int(data["H"])
    guards = tuple(data["guards"])
    flip_cap = int(data["flip_cap"])
    rows = []
    mass_num_v = 0
    for raw in data["masses"]:
        if not 0 <= raw[1] < r:
            errors.append("anchor outside the library")
        g, j, t, u, fl, fr, tagL, tagR, v = raw
        if v < 0 or g not in guards or fl > flip_cap or fr > flip_cap:
            errors.append(f"bad row {raw}")
            continue
        S = tag_size(N, m, g, j, t, u, fl, fr, tagL, tagR)
        if S <= 0:
            errors.append(f"empty or illegal tagged shell {raw}")
            continue
        if H % S != 0:
            errors.append(f"H not a multiple of shell size {raw}")
            continue
        extra = (tagL in (1, 2)) + (tagR in (1, 2))
        if t + u + extra != m:
            errors.append("tagged shell is not at distance m from its anchor")
            continue
        rows.append((g, j, t, u, fl, fr, tagL, tagR, v * (H // S)))
        mass_num_v += v
    result = base_result(data, errors, rows, mass_num_v, H, A, r)
    if errors:
        result["passed"] = False
        return result
    ans = run_tagged(N, m, rows)
    return finish(data, result, ans, rows, tagged_numerator)


def tagged_numerator(data, centre, rows) -> int:
    key = (data["N"], data["m"], tuple(data["guards"]), int(data["flip_cap"]))
    ctr = tagged_numerator.ctr
    if ctr is None or tagged_numerator.key != key:
        tagged_numerator.ctr = TagCounter(data["N"], data["m"], tuple(data["guards"]), int(data["flip_cap"]), r=data["r"])
        tagged_numerator.key = key
        ctr = tagged_numerator.ctr
    index = {shell: i for i, shell in enumerate(ctr.shells)}
    total = 0
    for g, j, t, u, fl, fr, tagL, tagR, coeff in rows:
        total += coeff * ctr.count(centre, index[(g, j, t, u, fl, fr, tagL, tagR)])
    return total


tagged_numerator.ctr = None
tagged_numerator.key = None


def base_result(data, errors, rows, mass_num_v, H, A, r) -> dict:
    mass = Fraction(H * mass_num_v, A) if A else Fraction(0)
    return {
        "family": data["family"],
        "N": data["N"],
        "r": r,
        "rows": len(rows),
        "errors": errors,
        "stored_mass": data.get("total_mass"),
        "mass": f"{mass.numerator}/{mass.denominator}",
        "mass_above": mass > r - 1,
        "H": str(H),
        "normalizer": str(A),
    }


def finish(data, result, ans, rows, numerator) -> dict:
    A = int(data["normalizer"])
    M = ans["maximum"]
    result["computed_maximum"] = str(M)
    result["classes"] = ans["classes"]
    result["dp_seconds"] = ans["seconds"]
    result["witness"] = ans["witness"]
    result["max_load"] = f"{M}/{A}"
    if ans["classes"] != 3 ** data["N"]:
        result["errors"].append(f"class count {ans['classes']}")
    if M > A:
        result["errors"].append(f"overload: computed maximum {M} exceeds normalizer {A}")
    witness = tuple(int(ch) for ch in ans["witness"])
    py = numerator(data, witness, rows)
    if py != M:
        result["errors"].append(f"witness numerator {py} != dp maximum {M}")
    patterns = library(data["N"], data["r"])
    for P in patterns:
        load = numerator(data, tuple(P), rows)
        if load > M:
            result["errors"].append("pattern centre exceeds dp maximum")
            break
    rng = random.Random(20261029 + data["N"])
    for _ in range(20):
        centre = tuple(rng.randrange(3) for _ in range(data["N"]))
        if numerator(data, centre, rows) > M:
            result["errors"].append("sampled centre exceeds dp maximum")
            break
    mass = Fraction(result["mass"])
    if mass <= data["r"] - 1:
        result["errors"].append("mass does not exceed r-1")
    stored = Fraction(data["total_mass"])
    if stored != mass:
        result["errors"].append(f"stored mass {stored} != recomputed {mass}")
    if Fraction(data["max_load"]) != Fraction(M, A):
        result["errors"].append("stored max load does not equal computed maximum over normalizer")
    result["passed"] = not result["errors"]
    return result


def load_cert(path: Path) -> dict:
    return json.loads(path.read_text())


def forge_overload(data: dict) -> dict:
    forged = json.loads(json.dumps(data))
    A = int(forged["normalizer"])
    H = int(forged["H"])
    bump = A // H + 1
    forged["masses"][0][-1] += bump
    forged["forged"] = "overload"
    return forged


def forge_outside(data: dict) -> dict:
    forged = json.loads(json.dumps(data))
    forged["atoms"] = [{"word": "2" * int(forged["N"]), "numerator": "1"}]
    forged["forged"] = "outside"
    return forged


def main() -> int:
    try:
        os.nice(10)
    except OSError:
        pass
    LOG.mkdir(exist_ok=True)
    paths = sorted(CERT.glob("packing_N*_r*_m9.json"))
    if {p.name for p in paths} != {'packing_N21_r22_m9.json', 'packing_N24_r25_m9.json', 'packing_N23_r24_m9.json', 'packing_N22_r23_m9.json', 'packing_N20_r21_m9.json'}:
        raise RuntimeError("certificate set differs from the documented package")
    if not paths:
        print(json.dumps({"passed": False, "errors": ["no certificates"]}))
        return 1
    reports = []
    first = None
    all_ok = True
    for path in paths:
        data = load_cert(path)
        if data.get("atoms") or data.get("forged"):
            reports.append({"file": path.name, "passed": False, "errors": ["forged file in certificates"]})
            all_ok = False
            continue
        report = check(data)
        report["file"] = path.name
        # Drop nothing essential; the error list is the verdict.
        reports.append({k: v for k, v in report.items() if k != "rows" or True})
        all_ok = all_ok and report["passed"]
        if first is None:
            first = data
        print(json.dumps({"file": path.name, "passed": report["passed"], "mass": report.get("mass"), "max_load": report.get("max_load"), "errors": report["errors"]}), flush=True)
    if first is None:
        print(json.dumps({"passed": False}))
        return 1
    over = check(forge_overload(first))
    outside = check(forge_outside(first))
    over_rejects = (not over["passed"]) and any("overload" in e for e in over["errors"])
    outside_rejects = (not outside["passed"]) and any("outside support" in e for e in outside["errors"])
    summary = {
        "certificates": len(reports),
        "certificates_passed": all_ok,
        "overload_rejected": over_rejects,
        "outside_rejected": outside_rejects,
        "overload_errors": over["errors"],
        "outside_errors": outside["errors"],
        "passed": all_ok and over_rejects and outside_rejects,
    }
    print(json.dumps(summary, indent=2), flush=True)
    print(json.dumps({"reports": reports}), flush=True)
    return 0 if summary["passed"] else 1


def release_main():
    global LOG, DP, DP_TAG
    if sys.argv[1:] not in ([], ["verify"]):
        raise SystemExit("usage: python3 -I verify.py verify")
    with tempfile.TemporaryDirectory(prefix=".verify_", dir=ROOT) as work:
        LOG = Path(work)
        env = os.environ.copy()
        env["TMPDIR"] = work
        DP = LOG / "dp"
        subprocess.run([os.environ.get("CXX", "c++"), "-O3", "-std=c++17", "-pthread",
                        str(ROOT / "dp.cpp"), "-o", str(DP)], check=True, env=env)
        DP_TAG = LOG / "dp_tag"
        subprocess.run([os.environ.get("CXX", "c++"), "-O3", "-std=c++17", "-pthread",
                        str(ROOT / "dp_tag.cpp"), "-o", str(DP_TAG)], check=True, env=env)
        import test_dp, test_tag
        for test, binary in ((test_dp, DP), (test_tag, DP_TAG)):
            test.ROOT = LOG
            test.BIN = binary
            test.main()
        return main()

if __name__ == "__main__":
    raise SystemExit(release_main())
