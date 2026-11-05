#!/usr/local/bin/python3 -I
"""Replay the deposit; --manuscript also recompiles local paper copies."""
import argparse
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript", action="store_true")
    args = parser.parse_args()
    support = runpy.run_path(str(HERE / "run_support.py"))
    run, python = support["run"], support["PYTHON"]
    stages = [("scoreboard", "scoreboard.py"), ("validate-patch", "validate_patch.py"),
              ("selftest", "selftest.py"), ("unit-tests", "test_conformance.py")]
    if args.manuscript:
        stages.append(("check-sentence", "check_sentence.py"))
    stages.extend((("render-docs", "render_docs.py"), ("verify-deposit", "verify_deposit.py")))
    for label, filename in stages:
        result = run(label, [python, "-I", HERE / filename], timeout=300 if label == "check-sentence" else 60)
        print(result.stdout[-1000:], flush=True)
    run("checksum-check", ["/usr/bin/shasum", "-a", "256", "-c", HERE / "MANIFEST.sha256"], cwd=HERE.parent)
    print("Deposit replay passed. Commands and exit codes: evidence/runs.jsonl")


if __name__ == "__main__":
    main()
