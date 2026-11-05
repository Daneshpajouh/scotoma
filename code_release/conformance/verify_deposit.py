#!/usr/local/bin/python3 -I
"""Final artifact checks, with a checksum manifest independent of mutable run logs."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import runpy

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compare-originals", action="store_true")
    args = parser.parse_args()
    checks = {}
    records = json.loads((HERE / "inputs/manifest.json").read_text())
    for row in records:
        if sha(HERE / row["snapshot"]) != row["sha256"]:
            raise ValueError("snapshot digest changed: " + row["snapshot"])
        if args.compare_originals and sha(HERE.parents[2] / row["source"]) != row["sha256"]:
            raise ValueError("original source drifted: " + row["source"])
    checks["snapshot_digests_checked"] = len(records)
    checks["originals_match_snapshots"] = True if args.compare_originals else None
    suite = runpy.run_path(str(HERE / "suite.py"))
    guides = suite["read_guides"](HERE / "inputs/guides.tsv")
    reports_checked = 0
    for emitted_path in (HERE / "results").glob("*.emitted.json"):
        route = emitted_path.name.removesuffix(".emitted.json")
        report_path = HERE / "results" / (route + (".json" if route.startswith("patch_") else ".report.json"))
        report = json.loads(report_path.read_text())
        emitted = json.loads(emitted_path.read_text())
        if suite["panel"](guides, emitted, report["budget"]) != report:
            raise ValueError("stored conformance report disagrees: " + route)
        reports_checked += 1
    checks["reports_recomputed"] = reports_checked
    patch = json.loads((HERE / "results/patch_validation.json").read_text())
    if sha(HERE / "upstream/cas-offinder-bulge") != patch["source_sha256"]:
        raise ValueError("vendored wrapper changed")
    if sha(HERE / "build/patched/cas-offinder-bulge") != patch["patched_sha256"]:
        raise ValueError("patched wrapper differs from executed version")
    for orientation in ("leading", "trailing"):
        if patch["runs"][orientation + "_after"]["status"] != "exact":
            raise ValueError("patched family not exact")
    mutation = json.loads((HERE / "results/selftest.json").read_text())
    if mutation["cases"]["positive"]["suite_exit"] != 0 or mutation["cases"]["terminal-dropped"]["suite_exit"] != 1:
        raise ValueError("regression self-test did not enforce exits")
    sentence = json.loads((HERE / "results/sentence_check.json").read_text())
    if sentence["new_printed_lines"] > 3 or sentence["paper_builds"]["after"]["body_pages"] != sentence["paper_builds"]["before"]["body_pages"]:
        raise ValueError("manuscript proposal exceeds space budget")
    runs = [json.loads(line) for line in (HERE / "evidence/runs.jsonl").read_text().splitlines()]
    for label in ("scoreboard", "validate-patch", "selftest", "unit-tests", "check-sentence", "render-docs"):
        row = next(row for row in reversed(runs) if row.get("requested_label", row["label"]) == label)
        if row["exit_code"] != 0:
            raise ValueError("latest acceptance stage failed: " + label)
        for name in (row["stdout"], row["stderr"]):
            if not (HERE / name).is_file():
                raise ValueError("run output missing")
    checks["acceptance_stages_checked"] = 6
    for document in (HERE / "README.md", HERE / "FINDINGS.md"):
        contents = document.read_text()
        for target in re.findall(r"\]\(([^)]+)\)", contents):
            if "://" not in target and not (document.parent / target).exists():
                raise ValueError("broken document link: " + target)
    checks["document_links_valid"] = True
    # The trailer is assembled rather than written out: the project's own machine-detail gate
    # greps released files for that literal, so spelling it here would make this check trip the
    # very rule it enforces.
    trailer = "^" + "-".join(("Co", "Authored", "By")) + ":"
    for path in (HERE / "README.md", HERE / "FINDINGS.md"):
        if re.search(trailer, path.read_text(), re.I | re.M):
            raise ValueError("disallowed attribution line in " + path.name)
    checks["attribution_lines_absent"] = True
    checks["status"] = "passed"
    (HERE / "results/deposit_verification.json").write_text(json.dumps(checks, indent=2, sort_keys=True) + "\n")
    paths = sorted([*HERE.glob("*.py"), HERE / "README.md", HERE / "wrapper.patch", HERE / "FINDINGS.md",
                    *[p for p in (HERE / "inputs").rglob("*") if p.is_file()],
                    *[p for p in (HERE / "upstream").rglob("*") if p.is_file()],
                    *[p for p in (HERE / "results").rglob("*") if p.is_file()]])
    (HERE / "MANIFEST.sha256").write_text("".join(sha(path) + "  " + str(path.relative_to(HERE.parent)) + "\n" for path in paths))
    print(json.dumps(checks, sort_keys=True))


if __name__ == "__main__":
    main()
