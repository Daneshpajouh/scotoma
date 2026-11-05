#!/usr/bin/env python3
"""Replay the clinical-guide result from the four raw tool outputs in this directory.

    python3 -I verify_clinical_guide.py

Recomputes every number the paper states for guide sg1617 and prints a claim-by-claim
verdict. Exit status 0 only when every claim is reproduced. No genome, no network.

Locus convention, identical to the panel analysis: a locus is
(guide, chromosome, strand, PAM-side forward boundary), where the boundary is
Position + length of the ungapped DNA on the + strand and Position on the - strand,
zero-based, with no 3-bp merging.
"""
import csv
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SPACER = "CTAACAGTTGCTTTTATCAC"
VERDICTS = []


def claim(ok, label, got, want):
    VERDICTS.append(ok)
    print("  %-5s %-56s got %-16s want %s" % ("PASS" if ok else "FAIL", label, got, want))


def read(label):
    with open(os.path.join(HERE, label + ".tsv"), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def locus(r):
    pos, strand = int(r["Position"]), r["Direction"]
    return (r["#Id"], r["Chromosome"], strand,
            pos + len(r["DNA"].replace("-", "")) if strand == "+" else pos)


def del_pos(r):
    i = r["DNA"].find("-")
    return None if i < 0 else i + 1


rec = json.load(open(os.path.join(HERE, "run_record.json")))
claim(rec["spacer"] == SPACER, "guide spacer", rec["spacer"], SPACER)
claim(rec["staged_records"] == 25, "records scanned", rec["staged_records"], 25)
claim(rec["staged_bases"] == 3088286401, "assembled bases scanned",
      rec["staged_bases"], 3088286401)
claim(rec["reference_sha256"] == "9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7",
      "reference hash matches the build gate", rec["reference_sha256"][:16] + "...", "9cce8b9264...")
claim(rec["binary_sha256"]["released"].startswith("b8991a89"), "released binary hash",
      rec["binary_sha256"]["released"][:8], "b8991a89")
claim(rec["binary_sha256"]["corrected"].startswith("cee92311"), "corrected binary hash",
      rec["binary_sha256"]["corrected"][:8], "cee92311")

data = {l: read(l) for l in ("released_rna1", "corrected_rna1", "released_rna0", "corrected_rna0")}
rel, cor = {locus(r) for r in data["released_rna1"]}, {locus(r) for r in data["corrected_rna1"]}
claim(len(rel) == 4414, "loci the released build reports", len(rel), 4414)
claim(len(cor) == 4778, "loci the corrected build reports", len(cor), 4778)
claim(len(cor - rel) == 364, "loci only the corrected build reports", len(cor - rel), 364)
claim(len(rel - cor) == 0, "loci only the released build reports", len(rel - cor), 0)


def multiset(rows):
    return Counter((r["#Id"], r["Chromosome"], r["Position"], r["Direction"], r["DNA"],
                    r["Mismatches"]) for r in rows)


claim(multiset(data["released_rna0"]) == multiset(data["corrected_rna0"]),
      "bulge-free control identical as a multiset", len(data["released_rna0"]), "equal")

by = defaultdict(list)
for r in data["corrected_rna1"]:
    by[locus(r)].append(r)
cls, pos, mm, nongg = Counter(), Counter(), Counter(), 0
for k in cor - rel:
    for r in by[k]:
        cls[r["Bulge type"]] += 1
        pos[del_pos(r)] += 1
        mm[r["Mismatches"]] += 1
    if by[k][0]["DNA"][-2:] != "GG":
        nongg += 1
claim(list(cls) == ["RNA"], "every corrected-only alignment is an RNA bulge", dict(cls), "{'RNA': ...}")
claim(list(pos) == [20], "every corrected-only deletion is at position 20", dict(pos), "{20: ...}")
claim(list(mm) == ["4"], "every corrected-only alignment is at four mismatches", dict(mm), "{'4': ...}")
claim(nongg == 0, "corrected-only loci with a PAM that is not literally GG", nongg, 0)

# The mismatch shells. An earlier version of this file computed BOTH shells from the
# corrected build's alignments, so the released shell never read released mismatch
# values and the check could not fail; it also asserted an equality that is false on
# this data. Each build is now read from its own rows, and the real asymmetry is the
# claim being tested.
def best_own(rows):
    d = {}
    for r in rows:
        k, v = locus(r), int(r["Mismatches"])
        if v < d.get(k, 99):
            d[k] = v
    return d


rel_best, cor_best = best_own(data["released_rna1"]), best_own(data["corrected_rna1"])
rel_shell = {i: sum(1 for v in rel_best.values() if v == i) for i in range(5)}
cor_shell = {i: sum(1 for v in cor_best.values() if v == i) for i in range(5)}
claim(rel_shell == {0: 1, 1: 4, 2: 29, 3: 381, 4: 3999},
      "released mismatch shells, from released values", rel_shell,
      "{0: 1, 1: 4, 2: 29, 3: 381, 4: 3999}")
claim(cor_shell == {0: 1, 1: 4, 2: 30, 3: 429, 4: 4314},
      "corrected mismatch shells, from corrected values", cor_shell,
      "{0: 1, 1: 4, 2: 30, 3: 429, 4: 4314}")

# The ranking error, visible at the deployed budget alone: among loci BOTH builds
# report, some carry a released mismatch count higher than the truth.
shared = set(rel_best) & set(cor_best)
inflated = [k for k in shared if rel_best[k] > cor_best[k]]
claim(len(inflated) == 50, "shared loci whose released count is too high", len(inflated), 50)
claim({rel_best[k] - cor_best[k] for k in inflated} == {1},
      "every such inflation is exactly one mismatch",
      sorted({rel_best[k] - cor_best[k] for k in inflated}), "[1]")
true_counts = {i: sum(1 for k in inflated if cor_best[k] == i) for i in (2, 3)}
claim(true_counts == {2: 1, 3: 49},
      "their true counts match the budget-2 and budget-3 misses", true_counts, "{2: 1, 3: 49}")
claim(not [k for k in shared if rel_best[k] < cor_best[k]],
      "no locus where the released count is too low",
      len([k for k in shared if rel_best[k] < cor_best[k]]), 0)

chroms = Counter(k[1] for k in cor - rel)
claim(len(chroms) == 24, "records carrying at least one missed locus", len(chroms), 24)

ovh = len(cor - rel) / len(rel)
claim(abs(ovh - 0.082465) < 1e-5, "candidate overhead for this guide",
      "364/4414 = %.6f" % ovh, "0.082465")

# ---------------------------------------------------------------- budget sweep
def sweep(label):
    with open(os.path.join(HERE, "sweep_" + label + ".tsv"), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def best(rows):
    d = {}
    for r in rows:
        k, v = locus(r), int(r["Mismatches"])
        if v < d.get(k, 99):
            d[k] = v
    return d


print()
print("  budget sweep, loci under the same convention")
EXPECT = {0: (1, 1, 0), 1: (5, 5, 0), 2: (34, 35, 1), 3: (415, 464, 49)}
sweep_best = {}
for mm in (0, 1, 2, 3):
    rl, cr = best(sweep("mm%d_released_rna1" % mm)), best(sweep("mm%d_corrected_rna1" % mm))
    sweep_best[mm] = (rl, cr)
    co, ro = len(set(cr) - set(rl)), len(set(rl) - set(cr))
    er, ec, eco = EXPECT[mm]
    claim((len(rl), len(cr), co) == (er, ec, eco),
          "budget %d: released / corrected / corrected-only" % mm,
          (len(rl), len(cr), co), (er, ec, eco))
    claim(ro == 0, "budget %d: released-only loci" % mm, ro, 0)
    claim(multiset(sweep("mm%d_released_rna0" % mm)) == multiset(sweep("mm%d_corrected_rna0" % mm)),
          "budget %d: bulge-free control identical" % mm, "equal", "equal")

# The decisive mechanism claim: a locus the released build cannot report at budget b is
# reported at a higher budget with a mismatch count inflated by exactly one, never fewer.
rl4 = best(data["released_rna1"])
cr3 = sweep_best[3][1]
missed3 = set(cr3) - set(sweep_best[3][0])
claim(len(missed3) == 49, "loci invisible to the released build at budget 3", len(missed3), 49)
infl = [rl4[k] - cr3[k] for k in missed3 if k in rl4]
claim(len(infl) == 49 and set(infl) == {1},
      "each is reported at budget 4, inflated by exactly one mismatch",
      "%d loci, inflation set %s" % (len(infl), sorted(set(infl))), "49 loci, [1]")

cr2 = sweep_best[2][1]
missed2 = set(cr2) - set(sweep_best[2][0])
claim(len(missed2) == 1, "loci invisible to the released build at budget 2", len(missed2), 1)
k2 = next(iter(missed2))
claim(cr2[k2] == 2 and rl4[k2] == 3,
      "that locus: true count 2, best the released build can report 3",
      (cr2[k2], rl4[k2]), (2, 3))

# ------------------------------------------------- the alignments themselves, not their labels
# Every count above came from the Mismatches column. A reviewer rewrote one four-mismatch
# alignment into a perfect match over A, C, G, T while keeping its gap, PAM, coordinate and the
# label 4, and nothing here noticed, because no check compared the two strings. Each recorded
# alignment is now recomputed from its own crRNA and DNA.
ALPHA = set("ACGT")
alpha_bad = relabel = guide_bad = gapcount = pam_bad = len_bad = 0
# Main budget-4 outputs and every budget-0..3 sweep file. A reviewer rewrote a budget-2
# sweep alignment to a perfect match while keeping Mismatches=2, and separately set shared
# PAM bases to AA on budget-4 rows, both while leaving other files alone. String checks
# must cover every row that feeds a claim.

def check_rows(rows):
    global alpha_bad, relabel, guide_bad, gapcount, pam_bad, len_bad
    for r in rows:
        cr, dna = r["crRNA"].upper(), r["DNA"].upper()
        if set(dna) - ALPHA - {"-"} or set(cr) - ALPHA - {"N", "-"}:
            alpha_bad += 1
            continue
        if len(cr) != len(dna):
            len_bad += 1
            continue
        want_gaps = 1 if r["Bulge type"] == "RNA" else 0
        if dna.count("-") != want_gaps:
            gapcount += 1
            continue
        if cr[:len(cr) - 3].replace("-", "") != SPACER:
            guide_bad += 1
        a, b = cr[:len(cr) - 3], dna[:len(dna) - 3]
        got = sum(1 for x, y in zip(a, b) if x != y and x != "-" and y != "-" and x != "N")
        if got != int(r["Mismatches"]):
            relabel += 1
        pam = dna[len(dna) - 3:]
        if "-" in pam or pam[1:] != "GG":
            pam_bad += 1

for name in ("released_rna1", "corrected_rna1", "released_rna0", "corrected_rna0"):
    check_rows(read(name))
for mm in range(4):
    for lab in ("released_rna0", "released_rna1", "corrected_rna0", "corrected_rna1"):
        check_rows(sweep("mm%d_%s" % (mm, lab)))
claim(alpha_bad == 0, "every main and sweep alignment is over the four letters", alpha_bad, 0)
claim(len_bad == 0, "  with equal crRNA and DNA column lengths", len_bad, 0)
claim(gapcount == 0, "  with one gap on an RNA row and none on an ungapped row", gapcount, 0)
claim(guide_bad == 0, "  and carries this guide's own spacer", guide_bad, 0)
claim(relabel == 0, "  and its mismatch count recomputed from the two strings", relabel, 0)
claim(pam_bad == 0, "  and every row ends in an NGG PAM on DNA", pam_bad, 0)
# Coordinates must be genomic positions, not free integers. A reviewer shifted every Position
# by -10^12; locus keys and mismatch shells survived because the shift was uniform.
pos_bad = 0
for name in ("released_rna1", "corrected_rna1", "released_rna0", "corrected_rna0"):
    for r in read(name):
        if int(r["Position"]) < 0:
            pos_bad += 1
for mm in range(4):
    for lab in ("released_rna0", "released_rna1", "corrected_rna0", "corrected_rna1"):
        for r in sweep("mm%d_%s" % (mm, lab)):
            if int(r["Position"]) < 0:
                pos_bad += 1
claim(pos_bad == 0, "every main and sweep Position is a nonnegative genomic coordinate", pos_bad, 0)

bad = VERDICTS.count(False)
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - bad, bad))
if bad:
    sys.exit(1)
print("Every clinical-guide claim in the paper is reproduced from these outputs.")
