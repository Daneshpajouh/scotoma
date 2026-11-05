#!/usr/bin/env python3
"""Build deliberately corrupted copies of this directory's inputs and run the checker on each.

    python3 -I mutations.py

Every corruption below is invalid for a reason visible in the deposited files alone. A
corruption the checker accepts is a defect in the checker, so this script exits non-zero if
any corruption passes.
"""
import csv, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, os.pardir, "exactness_dividend", "published_benchmark_guides.txt")
PANEL = os.path.join(HERE, os.pardir, "tool_conformance", "per_guide_counts.tsv")
VERIFY = os.path.join(HERE, "verify_minimum_pattern_certificate.py")


def run_verify(env_bench, env_panel):
    env = os.environ.copy()
    # checker reads fixed relative paths; we swap files in a temp mirror of code_release layout
    return subprocess.run(
        [sys.executable, "-I", "-B", VERIFY],
        cwd=HERE, env=env, capture_output=True, text=True,
    )


def with_swapped(bench_text=None, panel_rows=None, panel_fields=None):
    """Copy the two inputs aside, write corruptions, run, restore."""
    bench_bak = BENCH + ".bak_mut"
    panel_bak = PANEL + ".bak_mut"
    shutil.copy2(BENCH, bench_bak)
    shutil.copy2(PANEL, panel_bak)
    try:
        if bench_text is not None:
            open(BENCH, "w").write(bench_text)
        if panel_rows is not None:
            with open(PANEL, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=panel_fields, delimiter="\t",
                                   lineterminator="\n")
                w.writeheader()
                w.writerows(panel_rows)
        proc = subprocess.run([sys.executable, "-I", "-B", VERIFY],
                              cwd=HERE, capture_output=True, text=True)
        return proc.returncode == 0, proc.stdout + proc.stderr
    finally:
        shutil.move(bench_bak, BENCH)
        shutil.move(panel_bak, PANEL)


def load_panel():
    with open(PANEL, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows, list(rows[0].keys())


def drop_one_benchmark_guide():
    lines = [l for l in open(BENCH) if l.strip()]
    return with_swapped(bench_text="".join(lines[:-1]) + ("\n" if lines else ""))


def duplicate_all_benchmark_to_one():
    lines = [l for l in open(BENCH) if l.strip()]
    one = lines[0]
    return with_swapped(bench_text=one * len(lines))


def alter_panel_r_column():
    rows, fields = load_panel()
    for r in rows:
        r["r"] = str(int(r["r"]) + 1)
    return with_swapped(panel_rows=rows, panel_fields=fields)


def swap_two_panel_r_values():
    rows, fields = load_panel()
    # exchange r on EVERY row of the first two distinct guide_ids (a single-row
    # swap is invisible: later tool rows for the same guide overwrite the panel map)
    ids = []
    for r in rows:
        if r["guide_id"] not in ids:
            ids.append(r["guide_id"])
        if len(ids) == 2:
            break
    a_id, b_id = ids
    a_r = next(r["r"] for r in rows if r["guide_id"] == a_id)
    b_r = next(r["r"] for r in rows if r["guide_id"] == b_id)
    for r in rows:
        if r["guide_id"] == a_id:
            r["r"] = b_r
        elif r["guide_id"] == b_id:
            r["r"] = a_r
    return with_swapped(panel_rows=rows, panel_fields=fields)


def zero_patterns_emitted_on_incomplete_tool():
    rows, fields = load_panel()
    for r in rows:
        if r["tool"].startswith("Cas-OFFinder-bulge") and r["patterns_emitted"].strip().isdigit():
            r["patterns_emitted"] = "0"
    return with_swapped(panel_rows=rows, panel_fields=fields)


def inflate_benchmark_with_non_dna():
    lines = [l for l in open(BENCH) if l.strip()]
    lines[0] = "ACGTACGTACGTACGTACG N\n"  # illegal letter / length
    return with_swapped(bench_text="".join(lines))


CASES = [
    ("drop_one_benchmark_guide", drop_one_benchmark_guide, False),
    ("duplicate_all_benchmark_to_one", duplicate_all_benchmark_to_one, False),
    ("alter_panel_r_column", alter_panel_r_column, False),
    ("swap_two_panel_r_values", swap_two_panel_r_values, False),
    ("zero_patterns_emitted_on_incomplete_tool", zero_patterns_emitted_on_incomplete_tool, False),
    ("inflate_benchmark_with_non_dna", inflate_benchmark_with_non_dna, False),
]


def main():
    # clean run must pass first
    clean = subprocess.run([sys.executable, "-I", "-B", VERIFY],
                           cwd=HERE, capture_output=True, text=True)
    if clean.returncode != 0:
        print("clean verifier failed before mutations")
        print(clean.stdout)
        print(clean.stderr)
        return 1
    print("%-42s %-8s %s" % ("corruption", "verdict", "expected"))
    rejected = accepted = wrong = 0
    for name, fn, expect_pass in CASES:
        passed, out = fn()
        rejected += not passed
        accepted += passed
        ok = passed == expect_pass
        wrong += not ok
        print("%-42s %-8s %s%s" % (
            name, "accepted" if passed else "rejected",
            "accepted" if expect_pass else "rejected",
            "" if ok else "  UNEXPECTED"))
        if not ok:
            print(out[-1500:])
    print("%d corruptions built, %d rejected, %d accepted, %d unexpected"
          % (len(CASES), rejected, accepted, wrong))
    return 1 if wrong else 0


if __name__ == "__main__":
    sys.exit(main())
