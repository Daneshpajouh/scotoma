#!/usr/bin/env python3
"""Fixtures for the whole-token figure pattern in figure_lookup.py.
usage: test_boundaries.py [DIR]  (DIR holds the script; default: beside this file). Exit 0 if every case passes, 1 if not."""
import importlib.util, os, sys

sys.dont_write_bytecode = True
HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
CASES = [  # figure, text, should match
    ("46854", "fae46854a39f sha", False), ("978103", "chr978103", False),
    ("978103", "978103bp", False), ("978103", "x 978,103 y", True),
    ("978103", "x 978103 y", True), ("978103", "ends 978,103.", True),
    ("978103", "51978103", False), ("978103", "978,1030", False),
    ("978103", "978,103.5", False), ("202040", "of 202 040 total", True),
    ("202040", "202 040", True), ("734238", "4 734 238", False),
    ("4734", "4 734 238", False), ("1207622", "1207622 1207623", True),
    ("1207623", "1207622 1207623", True), ("1207622", "1207622\t1207623", True),
    ("5516.51", "5,516.51 s", True), ("5,516.51", "5516.51", True),
    ("5516.51", "5,516.519", False), ("5516.51", "15,516.51", False),
    ("20000", "id_20000", False), ("341716", "(341,716)", True),
    ("341716", "341,716,", True), ("100947", "|100,947|", True),
]


def load(name):
    spec = importlib.util.spec_from_file_location(name[:-3], os.path.join(HERE, name))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    fl = load("figure_lookup.py")
    tools = [("figure_lookup", fl.pattern)]
    bad = 0
    for name, pat in tools:
        for fig, text, want in CASES:
            if bool(pat(fig).search(text)) != want:
                bad += 1
                print("FAIL", name, repr(fig), repr(text), "expected", want)
    print("cases", len(CASES) * len(tools), "failures", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
