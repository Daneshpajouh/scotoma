#!/usr/bin/env python3
"""Recompute three bulge-census figures of the RECOMB paper from the deposited producers alone.

  24,385  loci reachable by both a DNA-bulge and an RNA-bulge pattern, complete position set
          (appendix_moved.tex:2387, sections/04_class.tex:103): |DNA| + |RNA| - |DNA or RNA|.
  3,834   single-bulge patterns, complete set, per-guide sum (appendix_moved.tex:465).
  3,581   single-bulge patterns, wrapper set, per-guide sum (appendix_moved.tex:465).

Usage: python3 -I census_figures.py [PRODUCERS_DIR]
Default PRODUCERS_DIR is data/producers_2026-10-20 under the current folder (run from code_release).
No paper value is typed here. Exit 0 on a clean run, 2 if an input check fails.
"""
import json
import os
import sys

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data", "producers_2026-10-20")


def load(name):
    with open(os.path.join(D, name)) as fh:
        return json.load(fh)


def per_guide_sum(idx):
    # Each pattern maps to the list of (guide index, guide, slot) entries that emit it,
    # so the per-guide sum counts a pattern once for every guide that emits it.
    return sum(len({entry[0] for entry in entries}) for entries in idx.values())


full = load("BULGE_CENSUS_size1_FULL_mm4.json")["totals"]
wrap = load("BULGE_CENSUS_size1_mm4.json")["totals"]
if full["ungapped"] != wrap["ungapped"]:
    print("INPUT CHECK FAILED: the two censuses do not share an ungapped universe")
    sys.exit(2)
for t in (full, wrap):
    if not t["add_any"] <= t["add_dna"] + t["add_rna"]:
        print("INPUT CHECK FAILED: union larger than the sum of its arms")
        sys.exit(2)

both_full = full["add_dna"] + full["add_rna"] - full["add_any"]
both_wrap = wrap["add_dna"] + wrap["add_rna"] - wrap["add_any"]
print(f"loci reachable by both arms, complete set: {both_full}")
print(f"loci reachable by both arms, wrapper set:  {both_wrap}")

f_dna = per_guide_sum(load("bulge_census_full_idx_dna.json"))
f_rna = per_guide_sum(load("bulge_census_full_idx_rna.json"))
w_dna = per_guide_sum(load("idx_dna.json"))
w_rna = per_guide_sum(load("bulge_census_idx_rna.json"))
print(f"single-bulge patterns, complete set, per-guide sum: {f_dna} DNA + {f_rna} RNA = {f_dna + f_rna}")
print(f"single-bulge patterns, wrapper set, per-guide sum:  {w_dna} DNA + {w_rna} RNA = {w_dna + w_rna}")
print(f"pattern deficit of the wrapper: {f_dna + f_rna - (w_dna + w_rna)}")
