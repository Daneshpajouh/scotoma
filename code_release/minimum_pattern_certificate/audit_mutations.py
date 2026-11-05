#!/usr/bin/env python3
"""Adversarial emitted sets for min_patterns.py --audit.

    python3 -I audit_mutations.py

Each fixture must produce the named verdict on at least one guide, and the exact
fixture must exit 0. A wrong verdict is a defect.
"""
import os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "audit_fixtures")
PROG = os.path.join(HERE, "min_patterns.py")
GUIDES = os.path.join(FIX, "guides.txt")

CASES = [
    ("emitted_exact.tsv", "EXACT", 0),
    ("emitted_incomplete.tsv", "INCOMPLETE", 1),
    ("emitted_redundant.tsv", "REDUNDANT", 1),
    ("emitted_invalid.tsv", "INVALID", 1),
]

def run(emitted):
    return subprocess.run(
        [sys.executable, "-I", "-B", PROG, GUIDES, "--budget", "4", "--audit", emitted],
        cwd=HERE, capture_output=True, text=True)

bad = 0
print("%-28s %-10s %-6s %s" % ("fixture", "want", "exit", "verdict_seen"))
for name, want, want_exit in CASES:
    path = os.path.join(FIX, name)
    proc = run(path)
    out = proc.stdout + proc.stderr
    seen = want in out
    ok = seen and (proc.returncode == want_exit)
    bad += not ok
    print("%-28s %-10s %-6s %s%s" % (
        name, want, proc.returncode, "yes" if seen else "NO",
        "" if ok else "  UNEXPECTED"))
    if not ok:
        print(out)
print("%d fixtures, %d unexpected" % (len(CASES), bad))
sys.exit(1 if bad else 0)
