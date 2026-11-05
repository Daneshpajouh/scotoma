#!/usr/local/bin/python3 -I
"""Build and apply the minimal diff locally; execute both wrappers in dry-run mode."""
import difflib
import hashlib
import json
from pathlib import Path
import runpy
import shutil

HERE = Path(__file__).resolve().parent
SUPPORT = runpy.run_path(str(HERE / "run_support.py"))
SUITE = runpy.run_path(str(HERE / "suite.py"))
run, PYTHON = SUPPORT["run"], SUPPORT["PYTHON"]


def main():
    original = HERE / "upstream/cas-offinder-bulge"
    source = original.read_text()
    old = "for i in range(1, len(target)-bulge_size):"
    new = "for i in range(len(target)-bulge_size+1):"
    if source.count(old) != 2:
        raise ValueError("unexpected vendored RNA loop count")
    revised = source.replace(old, new)
    patch = "".join(difflib.unified_diff(source.splitlines(keepends=True), revised.splitlines(keepends=True),
                                       fromfile="a/cas-offinder-bulge", tofile="b/cas-offinder-bulge", n=3))
    patch_path = HERE / "wrapper.patch"
    patch_path.write_text(patch)
    build = HERE / "build/patched"
    build.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(original, build / original.name)
    run("patch-dry-run", ["/usr/bin/patch", "--dry-run", "-p1", "-i", patch_path], cwd=build)
    run("patch-apply", ["/usr/bin/patch", "-p1", "-i", patch_path], cwd=build)
    if (build / original.name).read_text() != revised:
        raise ValueError("applied patch does not equal proposed source")
    guides = SUITE["read_guides"](HERE / "inputs/guides.tsv")
    expected_union = set().union(*(SUITE["REF"]["required_words"](g) for g in guides.values()))
    reports = {}

    def execute(label, script, selected, leading=False):
        directory = HERE / "build/dry_runs" / label
        directory.mkdir(parents=True, exist_ok=True)
        input_path = directory / "input.txt"
        pattern = "CCN" + "N" * 20 if leading else "N" * 20 + "NGG"
        queries = ["NNN" + g if leading else g + "NNN" for g in selected.values()]
        input_path.write_text(str(HERE / "inputs/fixture_rna.fa") + "\n" + pattern + " 0 1\n" +
                              "".join(query + " 4\n" for query in queries))
        run(label, [PYTHON, "-I", script, input_path, "D", directory / "unused-output.txt"])
        lines = (directory / "input_bulge.txt").read_text().splitlines()
        words = []
        for line in lines[2:]:
            p, budget = line.split()
            if budget != "4":
                raise ValueError("wrapper changed mismatch budget")
            if len(p) == 23 and ((leading and p.endswith("N")) or (not leading and p.startswith("N"))):
                word = p[3:-1] if leading else p[1:-3]
                if len(word) != 19 or set(word) - set("ACGT"):
                    raise ValueError("malformed projected wrapper deletion word")
                words.append(word)
        return words

    for leading, orientation in ((False, "trailing"), (True, "leading")):
        for version, script in (("before", original), ("after", build / original.name)):
            emitted = {gid: execute(f"wrapper-{orientation}-{version}-{gid}", script, {gid: g}, leading)
                       for gid, g in guides.items()}
            report = SUITE["panel"](guides, emitted, budget=4)
            key = f"{orientation}_{version}"
            reports[key] = {k: v for k, v in report.items() if k != "per_guide"}
            (HERE / f"results/patch_{key}.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            emit_path = HERE / f"results/patch_{key}.emitted.json"
            emit_path.write_text(json.dumps(emitted, indent=2, sort_keys=True) + "\n")
            run(f"suite-patch-{key}", [PYTHON, "-I", HERE / "suite.py", "--guides", HERE / "inputs/guides.tsv",
                                      "--emitted", emit_path], expected=0 if version == "after" else 1)
            if version == "before":
                deposited = json.loads((HERE / "results/wrapper.emitted.json").read_text())
                if any(set(emitted[gid]) != set(deposited[gid]) for gid in guides):
                    raise ValueError("fresh unmodified wrapper differs from deposited per-guide outputs")
            pooled = execute(f"wrapper-fixture-{orientation}-{version}", script, guides, leading)
            reports[key]["fixture_pooled_emitted"] = len(set(pooled))
            reports[key]["fixture_pooled_required"] = len(expected_union)
            reports[key]["fixture_pooled_missing"] = len(expected_union - set(pooled))
            reports[key]["fixture_pooled_invalid"] = len(set(pooled) - expected_union)
            if set(pooled) != set().union(*(set(words) for words in emitted.values())):
                raise ValueError("pooled emission differs from per-guide union")
            if version == "after" and (report["status"] != "exact" or set(pooled) != expected_union):
                raise ValueError("patch not exact")
    lines = patch.splitlines()
    added = sum(line.startswith("+") and not line.startswith("+++") for line in lines)
    removed = sum(line.startswith("-") and not line.startswith("---") for line in lines)
    result = {"added_lines": added, "removed_lines": removed, "modified_source_lines": removed,
              "unified_diff_physical_lines": len(lines), "source_sha256": hashlib.sha256(original.read_bytes()).hexdigest(),
              "patched_sha256": hashlib.sha256((build / original.name).read_bytes()).hexdigest(),
              "scope": "single-deletion word generation only; dry-run, no native search", "runs": reports}
    (HERE / "results/patch_validation.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
