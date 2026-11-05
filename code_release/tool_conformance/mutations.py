#!/usr/bin/env python3
"""Plant corrupted tool_conformance inputs; exit non-zero if an unmarked one is accepted."""
import csv, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def copy_tree(dst):
    for name in os.listdir(HERE):
        src = os.path.join(HERE, name)
        if name.endswith(".tsv") or name.endswith(".py") and name.startswith("verify"):
            shutil.copy(src, os.path.join(dst, name))
        elif name == "verify_tool_conformance.py":
            shutil.copy(src, os.path.join(dst, name))


def run(dst):
    return subprocess.run([sys.executable, "-I", "-B", "verify_tool_conformance.py"],
                          cwd=dst, capture_output=True, text=True).returncode


def two_gaps(dst):
    for fn in ("crispritz_2.9.0_rna1.targets.tsv", "crispritz_2.7.1_rna1.targets.tsv"):
        p = os.path.join(dst, fn)
        with open(p, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
            fields = list(rows[0].keys())
        for r in rows:
            r["DNA"] = r["DNA"] + "-"
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)


CASES = [("conformance_two_gaps", two_gaps, False)]
print("%-28s %-8s %s" % ("corruption", "verdict", "expected"))
unexpected = 0
for name, fn, expect_accept in CASES:
    d = tempfile.mkdtemp()
    try:
        copy_tree(d)
        fn(d)
        code = run(d)
        accepted = code == 0
        verdict = "accepted" if accepted else "rejected"
        expected = "accepted" if expect_accept else "rejected"
        bad = accepted != expect_accept
        unexpected += bad
        print("%-28s %-8s %s%s" % (name, verdict, expected, "   <-- UNEXPECTED" if bad else ""))
    finally:
        shutil.rmtree(d)
print("%d corruptions, %d unexpected" % (len(CASES), unexpected))
sys.exit(1 if unexpected else 0)
