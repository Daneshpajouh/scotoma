#!/usr/local/bin/python3 -I
"""Mutate the actual producer, then require the unmodified CLI suite to reject it."""
import ast
import difflib
import json
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent
SUPPORT = runpy.run_path(str(HERE / "run_support.py"))
run, PYTHON = SUPPORT["run"], SUPPORT["PYTHON"]
GUIDE = "GAGTAGCGCGAGCACAGCTA"


def main():
    source = (HERE / "reference.py").read_text()
    tree = ast.parse(source)
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "required_words")
    statement = next(node for node in function.body if isinstance(node, ast.Return))
    statement.value = ast.BinOp(left=statement.value, op=ast.Sub(),
                               right=ast.parse("{guide[:-1]}", mode="eval").body)
    mutant = ast.unparse(ast.fix_missing_locations(tree)) + "\n"
    directory = HERE / "build/terminal_regression"
    directory.mkdir(parents=True, exist_ok=True)
    mutated_path = directory / "reference.py"
    mutated_path.write_text(mutant)
    (directory / "mutation.patch").write_text("".join(difflib.unified_diff(
        source.splitlines(keepends=True), mutant.splitlines(keepends=True),
        fromfile="reference.py", tofile="mutated/reference.py")))
    reports = {}
    for label, producer, expected in (("positive", HERE / "reference.py", 0), ("terminal-dropped", mutated_path, 1)):
        produced = run(f"selftest-producer-{label}", [PYTHON, "-I", producer, GUIDE])
        certificate = json.loads(produced.stdout)
        emitted_path = directory / f"{label}.json"
        emitted_path.write_text(json.dumps(certificate["words"]) + "\n")
        checked = run(f"selftest-suite-{label}", [PYTHON, "-I", HERE / "suite.py", "--guide", GUIDE,
                                               "--emitted", emitted_path], expected=expected)
        report = json.loads(checked.stdout)
        reports[label] = {"producer_exit": produced.returncode, "suite_exit": checked.returncode,
                          "r": certificate["r"], "emitted_words": len(certificate["words"]),
                          "suite_report": report}
    bad = reports["terminal-dropped"]["suite_report"]["per_guide"]["guide"]
    if bad["missing_words"] != [GUIDE[:-1]] or bad["status"] != "incomplete":
        raise ValueError("suite did not identify precisely the planted terminal word")
    result = {"status": "passed", "mutation": "required_words(g) minus {g[:-1]}", "cases": reports}
    (HERE / "results/selftest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "passed", "positive_suite_exit": reports["positive"]["suite_exit"],
                      "mutant_suite_exit": reports["terminal-dropped"]["suite_exit"],
                      "required": bad["r"], "mutant_emitted": bad["valid_distinct"],
                      "missing_word": bad["missing_words"][0]}, sort_keys=True))


if __name__ == "__main__":
    main()
