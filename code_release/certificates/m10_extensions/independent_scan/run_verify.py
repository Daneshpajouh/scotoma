"""Fail-closed verification and Slurm handoff entry points. No producer imports."""
import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

from measure import (COMPLEMENT, Measure, find_unit_witness, intersection, load,
                     read_certificate, resource_estimate, wrong_input_controls)

BASE = Path(__file__).resolve().parent


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(data, indent=2)+"\n")
    tmp.replace(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def manifest():
    items = json.loads((BASE/"inputs/manifest.json").read_text())
    for item in items:
        if sha(BASE/item["certificate"]) != item["sha256"]:
            raise ValueError("staged input hash mismatch")
        if "producer_result" in item and sha(BASE/item["producer_result"]) != item["producer_result_sha256"]:
            raise ValueError("staged RESULT hash mismatch")
    return {x["name"]: x for x in items}


def forge_files(item, measure, witness, outdir):
    """Round-trip two actual wrong JSON inputs through the reader/checker."""
    control = wrong_input_controls(measure, witness)
    raw = json.loads((BASE/item["certificate"]).read_text())
    p, v = witness
    index = next(i for i, row in enumerate(raw["packing"])
                 if Fraction(row["mass"]) > 0 and intersection(measure.C, measure.m,
                    row["prefix"].translate(COMPLEMENT), row["suffix_weight"], p, v) > 0)
    control["raw_row_index"] = index
    scaled = json.loads(json.dumps(raw))
    scaled["packing"][index]["mass"] = str(Fraction(scaled["packing"][index]["mass"])*Fraction(1000000001, 1000000000))
    scaled_path = outdir/"scaled.certificate.json"
    write_json(scaled_path, scaled)
    parsed, _, _ = read_certificate(scaled_path)
    excess = load(parsed, p, v)-1
    if excess <= 0:
        raise AssertionError("wrong serialized scaled input accepted")
    moved = json.loads(json.dumps(raw))
    moved["packing"][index]["prefix"] = "X"*measure.L
    moved["packing"][index]["suffix_weight"] = measure.C
    moved_path = outdir/"outside.certificate.json"
    write_json(moved_path, moved)
    try:
        read_certificate(moved_path)
    except ValueError as e:
        if "outside U" not in str(e):
            raise
    else:
        raise AssertionError("wrong serialized outside input accepted")
    control["serialized_scaled_excess"] = str(excess)
    control["files"] = {p.name: sha(p) for p in (scaled_path, moved_path)}
    return control


def audit(item, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    measure, data, digest = read_certificate(BASE/item["certificate"])
    if (measure.N, measure.r, measure.m) != (item["N"], item["r"], item["m"]):
        raise ValueError("manifest cell mismatch")
    result = json.loads((BASE/item["producer_result"]).read_text()) if "producer_result" in item else None
    mass = measure.mass
    if item["closure_required"] and mass <= measure.r-1:
        raise ValueError("mass does not close the target cell")
    if result is not None and Fraction(result["tau_lower"]) != mass:
        raise ValueError("recomputed mass disagrees with claimed lower bound")
    witness = find_unit_witness(measure, result)
    info = {"name": item["name"], "certificate_sha256": digest,
            "N": measure.N, "r": measure.r, "m": measure.m,
            "total_mass": str(mass), "mass_decimal": format(float(mass), ".15g"),
            "mass_minus_r_minus_1": str(mass-(measure.r-1)),
            "mass_exceeds_r_minus_1": mass > measure.r-1,
            "support_classes": len(measure.terms), "support_inside_U": True,
            "support_min_distances": sorted({__import__("measure").min_distance(p, w) for p, w, _ in measure.terms}),
            "resource_estimate": resource_estimate(measure),
            "global_maximum_status": "NOT_EXECUTED"}
    if result is not None:
        info["producer_embedded_certificate_sha256"] = result.get("certificate_sha256")
        info["producer_embedded_hash_matches_decision_file"] = result.get("certificate_sha256") == digest
    if witness:
        info["independently_counted_unit_witness"] = {"prefix": witness[0], "suffix_weight": witness[1], "exact_load": "1"}
        info["negative_controls"] = forge_files(item, measure, witness, outdir/"controls")
    elif item["closure_required"]:
        raise ValueError("could not independently establish a unit-load seed for controls")
    write_json(outdir/"audit.json", info)
    print(f"AUDIT {item['name']}: mass={info['mass_decimal']} support={len(measure.terms)} seed={witness}", flush=True)
    return measure, info


def verify(item, outdir, binary, threads, maxbytes, seconds):
    started = time.monotonic()
    measure, info = audit(item, outdir)
    inp, native_path = outdir/"density.txt", outdir/"native.json"
    (outdir/"result.json").unlink(missing_ok=True)
    D = measure.write(inp)
    # Never allow a stale successful native result to survive an interrupted rerun.
    native_path.unlink(missing_ok=True)
    with (outdir/"native.log").open("w") as log:
        cmd = [str(binary), str(inp), str(native_path), str(threads), str(maxbytes), "auto"]
        print("NATIVE", " ".join(cmd), flush=True)
        subprocess.run(cmd, stdout=log, stderr=log, check=True, timeout=seconds)
    native = json.loads(native_path.read_text())
    if native["algorithm"] != "ZX2_exact_radial_tensor_v1" or native["classes_checked"] != 3**measure.L*(measure.C+1):
        raise ValueError("incomplete or wrong native scan")
    for key, value in (("N", measure.N), ("r", measure.r), ("m", measure.m)):
        if native[key] != value:
            raise ValueError("native dimensions disagree")
    if int(native["denominator"]) != D or Fraction(int(native["mass_numerator"]), D) != measure.mass:
        raise ValueError("native measure mismatch")
    maximum = Fraction(int(native["maximum_numerator"]), D)
    arg = native["argmax"]
    direct = load(measure, arg["prefix"], arg["suffix_weight"])
    if direct != maximum:
        raise ValueError("direct exact kernel disagrees at native argmax")
    if maximum != 1:
        raise ValueError(f"maximum is {maximum}, expected exactly one")
    witness = (arg["prefix"], arg["suffix_weight"])
    info["negative_controls_at_global_argmax"] = forge_files(item, measure, witness, outdir/"global_controls")
    info.update(global_maximum_status="EXHAUSTIVE_EXACT_PASS", exact_maximum=str(maximum),
                argmax=arg, argmax_in_certificate_coordinates={"prefix": arg["prefix"].translate(COMPLEMENT),
                    "suffix_weight": arg["suffix_weight"]},
                classes_checked=native["classes_checked"], backend=native["backend"],
                density_input_sha256=sha(inp), native_result_sha256=sha(native_path),
                binary_sha256=sha(binary), elapsed_seconds=time.monotonic()-started,
                verifier_identity=gate_identity(),
                slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                slurm_array_task_id=os.environ.get("SLURM_ARRAY_TASK_ID"),
                closure_verified=bool(item["closure_required"]),
                execution_state="H3" if item["closure_required"] else "CALIBRATION_PASS")
    write_json(outdir/"result.json", info)
    print(f"PASS {item['name']}: exact maximum 1 at {arg}; {native['classes_checked']} classes", flush=True)
    return info


def gate_identity():
    names = ("measure.py", "exact_tensor.cpp", "run_verify.py", "test_kernel.py", "inputs/manifest.json")
    return {name: sha(BASE/name) for name in names}


def check_gate(path):
    gate = json.loads(Path(path).read_text())
    if gate.get("status") != "SMOKE_PASS" or gate.get("identity") != gate_identity():
        raise ValueError("full runs require a passing smoke for these exact sources and inputs")
    for name in ("gap_19_3_11", "calibration_19_12_10", "cell_20_11_10"):
        row = gate["results"][name]
        result_path = Path(path).parent/name/"result.json"
        if sha(result_path) != row["result_sha256"]:
            raise ValueError("smoke result hash mismatch")
        d = json.loads(result_path.read_text())
        if d.get("global_maximum_status") != "EXHAUSTIVE_EXACT_PASS" or d.get("exact_maximum") != "1":
            raise ValueError("invalid smoke calibration result")
    return gate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("audit", "smoke", "cell", "collect"))
    parser.add_argument("--name")
    parser.add_argument("--binary", type=Path, default=BASE/"bin/exact_tensor")
    parser.add_argument("--output", type=Path, default=BASE/"runs")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--max-bytes", type=int, default=1073741824)
    parser.add_argument("--seconds", type=int, default=1100)
    parser.add_argument("--gate", type=Path, default=BASE/"runs/smoke/SMOKE_PASS.json")
    args = parser.parse_args()
    items = manifest()
    if args.phase == "audit":
        rows = [audit(item, args.output/item["name"])[1] for item in items.values()]
        write_json(args.output/"audit_summary.json", rows)
    elif args.phase == "smoke":
        gatepath = args.output/"SMOKE_PASS.json"
        gatepath.unlink(missing_ok=True)
        from test_kernel import test_all
        test_all(args.binary, args.output/"toy_tests")
        rows = {}
        for name in ("gap_19_3_11", "calibration_19_12_10", "cell_20_11_10"):
            result = verify(items[name], args.output/name, args.binary, args.threads, args.max_bytes, args.seconds)
            if name == "gap_19_3_11" and result["total_mass"] != "32/11":
                raise ValueError("gap calibration mass changed")
            if name == "calibration_19_12_10" and abs(Fraction(result["total_mass"])-Fraction("10.6014609356")) > Fraction(1, 10**10):
                raise ValueError("upper-row calibration mass disagrees with known answer")
            rows[name] = {"result_sha256": sha(args.output/name/"result.json")}
        write_json(gatepath, {"status": "SMOKE_PASS", "identity": gate_identity(), "results": rows,
                             "slurm_job_id": os.environ.get("SLURM_JOB_ID")})
    elif args.phase == "cell":
        check_gate(args.gate)
        if args.name not in items or not items[args.name]["closure_required"]:
            raise ValueError("unknown target cell")
        (args.output/args.name/"result.json").unlink(missing_ok=True)
        verify(items[args.name], args.output/args.name, args.binary, args.threads, args.max_bytes, args.seconds)
    elif args.phase == "collect":
        check_gate(args.gate)
        rows = []
        for item in items.values():
            if not item["closure_required"]:
                continue
            directory = args.gate.parent if item["name"] == "cell_20_11_10" else args.output
            path = directory/item["name"]/"result.json"
            result = json.loads(path.read_text())
            if result.get("certificate_sha256") != item["sha256"] or result.get("execution_state") != "H3":
                raise ValueError("missing or mismatched H3 result")
            if result.get("verifier_identity") != gate_identity() or sha(path.parent/"native.json") != result["native_result_sha256"]:
                raise ValueError("returned result/source hash mismatch")
            measure, _, _ = read_certificate(BASE/item["certificate"])
            native = json.loads((path.parent/"native.json").read_text())
            if Fraction(int(native["maximum_numerator"]), int(native["denominator"])) != 1:
                raise ValueError("returned maximum is not exactly one")
            if native["classes_checked"] != 3**measure.L*(measure.C+1) or native["ambient_words"] != 4**measure.N:
                raise ValueError("returned enumeration is incomplete")
            arg = native["argmax"]
            if load(measure, arg["prefix"], arg["suffix_weight"]) != 1:
                raise ValueError("returned argmax fails exact direct check")
            rows.append(result)
        write_json(args.output/"H3_summary.json", {"status": "SIX_H3_VERIFIED", "results": rows})


if __name__ == "__main__":
    main()
