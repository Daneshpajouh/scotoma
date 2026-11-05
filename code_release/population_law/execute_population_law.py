#!/usr/bin/env python3
"""Execute the vendored 2016 cas-offinder-bulge source, as shipped, once per guide, over every
guide of a given length, and record how many it leaves short of the required single-RNA-deletion
word set.

Nothing here models the generator. The generator is run in its own dry-run mode ('D'), which
needs no genome and no OpenCL, and the expanded pattern file it writes is read back.

Usage:  python3 -I execute_population_law.py <n> [workdir]
"""
import itertools
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WRAPPER = HERE.parent / "conformance/upstream/cas-offinder-bulge"
FIXTURE = HERE.parent / "conformance/inputs/fixture_rna.fa"
BUDGET = "4"


def required_words(guide):
    """Every distinct word reachable by deleting exactly one base. The definition, nothing else."""
    return {guide[:i] + guide[i + 1:] for i in range(len(guide))}


def wrapper_words(guide, n, work):
    """Run the shipped wrapper in dry-run mode and read back what it expanded to."""
    inp = work / "input.txt"
    inp.write_text(f"{FIXTURE}\n{'N' * n}NGG 0 1\n{guide}NNN {BUDGET}\n")
    done = subprocess.run([sys.executable, "-I", str(WRAPPER), str(inp), "D",
                           str(work / "unused.txt")],
                          capture_output=True, text=True, timeout=120)
    if done.returncode != 0:
        raise SystemExit(f"wrapper exit {done.returncode} on {guide}: {done.stderr[:400]}")
    words = set()
    for line in (work / "input_bulge.txt").read_text().splitlines()[2:]:
        pattern, budget = line.split()
        if budget != BUDGET:
            raise SystemExit("wrapper changed the mismatch budget")
        if len(pattern) == n + 3 and pattern.startswith("N"):
            word = pattern[1:-3]
            if len(word) == n - 1 and not set(word) - set("ACGT"):
                words.add(word)
    return words


def main():
    n = int(sys.argv[1])
    work = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / f"run_n{n}"
    work.mkdir(parents=True, exist_ok=True)
    if not WRAPPER.is_file():
        raise SystemExit(f"vendored wrapper not found at {WRAPPER}")

    short = complete = 0
    disagreements = []
    for letters in itertools.product("ACGT", repeat=n):
        guide = "".join(letters)
        need = required_words(guide)
        got = wrapper_words(guide, n, work)
        if not got <= need:
            raise SystemExit(f"wrapper emitted a word outside the required set on {guide}")
        is_short = got != need
        # the predicate, from two character comparisons and no search at all
        predicted_complete = guide[0] == guide[1] and guide[-1] == guide[-2]
        if is_short == predicted_complete:
            disagreements.append({"guide": guide, "executed_short": is_short,
                                  "predicate_says_complete": predicted_complete,
                                  "required": sorted(need), "emitted": sorted(got)})
        short += is_short
        complete += not is_short

    total = short + complete
    result = {
        "n": n,
        "guides": total,
        "short": short,
        "complete": complete,
        "short_fraction": f"{short}/{total}",
        "short_percent": round(100.0 * short / total, 4),
        "predicate_expectation": total * 15 // 16,
        "guides_where_execution_disagrees_with_predicate": len(disagreements),
        "disagreements": disagreements[:5],
        "evidence_class": "the shipped 2016 wrapper source executed per guide in its own dry-run mode",
        "route": "cas-offinder-bulge v1.2 beta (2016-07-06); CRISPRseek backend 2.4.1",
    }
    out = HERE / f"executed_population_n{n}.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 1 if disagreements else 0


if __name__ == "__main__":
    sys.exit(main())
