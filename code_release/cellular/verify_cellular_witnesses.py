#!/usr/bin/env python3
"""The cellular claim, recomputed from sequences: only the terminal deletion reaches these loci.

    python3 -I verify_cellular_witnesses.py

Tracking-seq2 reports five loci of guide GACCCCCTCCACCCCGCCTC active in 293T cells. For each,
this recomputes the distance from every single-deletion word of the guide and from the ungapped
register in the supplied genomic window, and asserts what the paper says: the complete library
reaches the locus at four substitutions, the emitted libraries need five, the deleted position
that achieves four is unique and is position 20, and no ungapped NGG register in the window
comes close.

The amplicon validation is the source's own measurement and is not recomputed here; it is cited
as the authors' result, not as ours.
"""
import importlib.util
import json
import os
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SPACER = 20
GUIDE = "GACCCCCTCCACCCCGCCTC"
SITES = {
    "chr1:121395975-121395997": "TATTCCCCCCACCCCGCCT",
    "chr11:32073136-32073159": "ACCCCCCCCAACCCCGCCT",
    "chr18:57434112-57434134": "CTACCCCTCCACCCCGACT",
    "chr19:5536276-5536298": "AACCCACCCCACCCCACCT",
    "chrX:68433937-68433960": "TCTCCCCACCACCCCGCCT",
}
# The paper's five-row record table prints these sheet row indices, so they are pinned and
# checked here rather than only quoted: 2026-11-02, with the interior column beside them.
SHEET_ROWS = {
    "chr1:121395975-121395997": 181,
    "chr11:32073136-32073159": 422,
    "chr18:57434112-57434134": 1143,
    "chr19:5536276-5536298": 1194,
    "chrX:68433937-68433960": 2665,
}
VERDICTS = []


def claim(ok, label, got, want):
    VERDICTS.append(ok)
    print("  %-5s %-52s got %-18s want %s" % ("PASS" if ok else "FAIL", label, got, want))


spec = importlib.util.spec_from_file_location("mps", ROOT / "code_release/make_pattern_sets.py")
mps = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mps)

rc = lambda s: s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
hd = lambda a, b: sum(1 for x, y in zip(a, b) if x != y) if len(a) == len(b) else 10 ** 6

recs = json.load(open(HERE / "Tracking2_verified_witnesses.json"))
claim(len(recs) == 5, "deposited cellular witnesses", len(recs), 5)
claim({r["guide"] for r in recs} == {GUIDE}, "all five are one guide",
      len({r["guide"] for r in recs}), 1)
claim({r["location"] for r in recs} == set(SITES),
      "  and the pinned loci", "identical" if {r["location"] for r in recs} == set(SITES)
      else "differ", "identical")

complete = mps.rna_words(GUIDE, 1)
emitted = mps.sets_for(GUIDE, 1)["rna_emitted"]
claim(len(complete) > len(emitted), "the complete library exceeds the emitted one",
      "%d vs %d" % (len(complete), len(emitted)), "complete larger")

worst_pos1, best_nob = 0, 10 ** 6
for r in sorted(recs, key=lambda x: x["location"]):
    loc, core = r["location"], r["canonical_core"].upper()
    claim(core == SITES[loc], "  %.26s: sequence matches the pin" % loc,
          "identical" if core == SITES[loc] else core, "identical")
    # every single-deletion word, by position, recomputed here
    per = {i + 1: hd(GUIDE[:i] + GUIDE[i + 1:], core) for i in range(SPACER)}
    full = min(per.values())
    best = sorted(i for i in per if per[i] == full)
    wrap = min(hd(core, p) for p in emitted)
    nob = min((hd(GUIDE, s[i:i + 23][:SPACER]) for s in (r["window"], rc(r["window"]))
               for i in range(len(s) - 22) if s[i:i + 23][21:23] == "GG"), default=10 ** 6)
    claim(full == 4, "    complete library reaches it at four", full, 4)
    claim(best == [SPACER], "    and only by deleting position 20", best, [SPACER])
    interior = min(per[i] for i in range(2, SPACER))
    claim(interior == 5, "    interior deletions 2..19 reach five", interior, 5)
    claim(r["row"] == SHEET_ROWS[loc], "    deposited sheet row", r["row"], SHEET_ROWS[loc])
    claim(wrap == 5, "    every emitted word needs five", wrap, 5)
    claim(nob >= 7, "    no ungapped NGG register in the window is near", nob, ">= 7")
    claim(r["locus_min"]["full"] == full and r["locus_min"]["wrapper"] == wrap,
          "    agrees with the deposited record", "yes", "yes")
    worst_pos1 = max(worst_pos1, per[1])
    best_nob = min(best_nob, nob)

claim(worst_pos1 >= 6, "the position-1 deletion never reaches four", worst_pos1, ">= 6")

# ------------------------------------------------- the PAM scope, in full
# The paper calls this an NGG claim, and that qualification has to be checkable. Under NGG no
# locus is ungapped-reachable at four; under NGA two are; under any second-base-G PAM four are.
# An earlier version of the paper said "ungapped-reachable only under a non-NGG PAM", which
# reads as one locus and understated it.
PAM_EXPECT = {
    "chr1:121395975-121395997": (12, 15, 4),
    "chr11:32073136-32073159": (7, 14, 4),
    "chr18:57434112-57434134": (11, 4, 4),
    "chr19:5536276-5536298": (10, 4, 4),
    "chrX:68433937-68433960": (10, 5, 5),
}
ngg = nga = ngstar = 0
for r in sorted(recs, key=lambda x: x["location"]):
    w = r["window"].upper()
    best = {"GG": 99, "GA": 99, "G*": 99}
    for strand in (w, rc(w)):
        for i in range(len(strand) - 22):
            t23 = strand[i:i + 23]
            d = hd(GUIDE, t23[:SPACER])
            if t23[21:23] == "GG":
                best["GG"] = min(best["GG"], d)
            if t23[21:23] == "GA":
                best["GA"] = min(best["GA"], d)
            if t23[21] == "G":
                best["G*"] = min(best["G*"], d)
    got = (best["GG"], best["GA"], best["G*"])
    claim(got == PAM_EXPECT[r["location"]],
          "  %.26s: ungapped NGG / NGA / NG*" % r["location"], got,
          PAM_EXPECT[r["location"]])
    ngg += best["GG"] <= 4
    nga += best["GA"] <= 4
    ngstar += best["G*"] <= 4
claim(ngg == 0, "under NGG none of the five is ungapped at four", ngg, 0)
claim(nga == 2, "  under NGA two are, which the paper now states", nga, 2)
claim(ngstar == 4, "  and under any second-base-G PAM four are", ngstar, 4)
claim(best_nob >= 7, "the closest ungapped register over all five", best_nob, ">= 7")

# control: the PAM matrix above is a comparison against pinned values, so corrupt a window and
# confirm the comparison can fail. Without this the PAM block could pass on any input that
# happened to agree, and a reader could not tell a check from a restatement.
def pam_triple(window):
    best = {"GG": 99, "GA": 99, "G*": 99}
    for strand in (window, rc(window)):
        for i in range(len(strand) - 22):
            t23 = strand[i:i + 23]
            d = hd(GUIDE, t23[:SPACER])
            if t23[21:23] == "GG":
                best["GG"] = min(best["GG"], d)
            if t23[21:23] == "GA":
                best["GA"] = min(best["GA"], d)
            if t23[21] == "G":
                best["G*"] = min(best["G*"], d)
    return (best["GG"], best["GA"], best["G*"])


planted = 0
for r in sorted(recs, key=lambda x: x["location"]):
    w = r["window"].upper()
    # plant a perfect NGG match for the guide in the middle of the window
    cut = len(w) // 2
    forged = (w[:cut] + GUIDE + "AGG" + w[cut + 23:])[:len(w)]
    if pam_triple(forged) != PAM_EXPECT[r["location"]]:
        planted += 1
claim(planted == 5, "  control: a forged window stops reproducing the matrix", planted, 5)

# control: move one base of a core and the exclusivity must break
bad = 0
for loc, core in SITES.items():
    alt = core[:SPACER - 2] + ("A" if core[SPACER - 2] != "A" else "C") + core[SPACER - 1:]
    per = {i + 1: hd(GUIDE[:i] + GUIDE[i + 1:], alt) for i in range(SPACER)}
    if min(per.values()) == 4 and sorted(i for i in per if per[i] == 4) == [SPACER]:
        bad += 1
claim(bad == 0, "  control: altering a core base breaks the exclusivity", bad, 0)

n = VERDICTS.count(False)
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - n, n))
if n:
    sys.exit(1)
print("On all five cellular loci the complete library reaches four substitutions and the only")
print("deleted position that does so is 20, the position the audited emitters never visit.")
