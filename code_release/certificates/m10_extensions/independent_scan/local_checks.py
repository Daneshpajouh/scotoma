"""Small N=19 gap calibration and hostile checks of parser/arithmetic guards."""
from fractions import Fraction
from itertools import product
import json
import os
from pathlib import Path
import subprocess

from measure import (Measure, decode, load, read_certificate, resource_estimate)
from run_verify import check_gate, manifest, verify, write_json
from test_kernel import native

BASE = Path(__file__).resolve().parent


def main():
    output = BASE/"local/checks"
    output.mkdir(parents=True, exist_ok=True)
    binary = BASE/"bin/exact_tensor"
    items = manifest()
    item = items["gap_19_3_11"]
    result = verify(item, output/item["name"], binary, 1, 1073741824, 120)
    assert result["total_mass"] == "32/11"
    measure, data, _ = read_certificate(BASE/item["certificate"])
    gap_native = json.loads((output/item["name"]/"native.json").read_text())
    for w in range(measure.C+1):
        for index in range(3**measure.L):
            exact = load(measure, decode(index, measure.L), w)
            assert exact == Fraction(int(gap_native["loads"][w*3**measure.L+index]), int(gap_native["denominator"]))
    assert max(load(measure, decode(i, measure.L), w) for i in range(3**measure.L) for w in range(measure.C+1)) == 1
    print("GAP CALIBRATION PASS: mass 32/11, exact max 1, all 162 classes agree", flush=True)
    # Production parser rejection tests use actual serialized certificates.
    mutations = []
    for tag, key, value in (("float_mass", "mass", 0.25), ("negative_mass", "mass", "-1/7"),
                            ("bad_symbol", "prefix", "??"), ("bad_weight", "suffix_weight", 18),
                            ("float_weight", "suffix_weight", 1.0)):
        forged = json.loads(json.dumps(data))
        forged["packing"][0][key] = value
        path = output/f"{tag}.json"
        write_json(path, forged)
        try:
            read_certificate(path)
        except (ValueError, TypeError):
            mutations.append({"case": tag, "status": "REJECTED"})
        else:
            raise AssertionError(f"malformed input accepted: {tag}")
    # Native engine rejects local parallelism, oversized allocation, and 128-bit overflow.
    huge = Measure(2, 2, 1, (("0", 0, Fraction(2**400+1)),)).validate()
    inp = output/"huge_density.txt"
    huge.write(inp)
    for tag, args in (("wide_bound", ["1", "1073741824", "wide"]),
                      ("local_threads", ["2", "1073741824", "auto"]),
                      ("local_memory", ["1", "2147483648", "auto"])):
        with (output/f"{tag}.log").open("w") as log:
            ran = subprocess.run([str(binary), str(inp), str(output/f"{tag}.native.json"), *args],
                                 stdout=log, stderr=log, timeout=30)
        assert ran.returncode != 0, tag
        mutations.append({"case": tag, "status": "REJECTED"})
    _, arbitrary = native(binary, huge, output/"automatic_big")
    assert arbitrary["backend"] == "cpp_int"
    assert Fraction(int(arbitrary["maximum_numerator"]), int(arbitrary["denominator"])) == huge.mass
    # No uncalibrated or source-mismatched full run can pass the smoke gate.
    for tag, contents in (("gate_not_passed", {"status": "FAILED"}),
                          ("gate_wrong_identity", {"status": "SMOKE_PASS", "identity": {}})):
        gate = output/f"{tag}.json"
        write_json(gate, contents)
        try:
            check_gate(gate)
        except ValueError:
            mutations.append({"case": tag, "status": "REJECTED"})
        else:
            raise AssertionError("invalid smoke gate accepted")
    bad_native = output/"outside_density.txt"
    bad_native.write_text("ZX2_DENSITY_V1 2 2 1 1 1\nX 1 1\n")
    with (output/"native_outside_support.log").open("w") as log:
        ran = subprocess.run([str(binary), str(bad_native), str(output/"bad_native.json"),
                              "1", "1073741824", "auto"], stdout=log, stderr=log, timeout=30)
    assert ran.returncode != 0
    mutations.append({"case": "native_outside_support", "status": "REJECTED"})
    write_json(output/"summary.json", {"gap_calibration": "PASS", "all_gap_classes_crosschecked": 162,
                                       "automatic_big_backend": "PASS",
                                       "rejection_tests": mutations})


if __name__ == "__main__":
    main()
