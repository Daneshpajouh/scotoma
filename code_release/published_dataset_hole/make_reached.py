#!/usr/bin/env python3
"""Reduce the raw whole-assembly scan to the two statistics the hole tables need.

Input   a directory holding queries.tsv, guides_meta.json and out/chr*.tsv, the per-record
        outputs of the scan for every distinct deletion word of the sixteen guides.
Output  reached_positions.tsv   one row per (guide, off-target 19-mer) the scan reached, with
                                the spacer deletion positions whose word reaches it
        reached_loci.tsv        per guide, distinct loci reached and loci reached only through
                                a word that sits at spacer position 1 or 20

A locus is (guide, record, position, strand). Everything the paper says about this scan is a
function of these two files; the 63 MB of raw tool output is needed only to rebuild them.
"""
import csv, glob, json, os, sys
from collections import defaultdict

W, OUT = sys.argv[1], sys.argv[2]
SPACER = 20
EMITTED = set(range(2, SPACER))

qpos, qguide = {}, {}
with open(os.path.join(W, "queries.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        qpos[r["query_id"]] = tuple(int(x) for x in r["positions"].split(","))
        qguide[r["query_id"]] = r["query_id"].split("_")[0]

seqpos, locpos = defaultdict(set), defaultdict(set)
rows = 0
files = sorted(glob.glob(os.path.join(W, "out", "chr*.tsv")))
for p in files:
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows += 1
            qid = r["#Id"]
            g, ps = qguide[qid], qpos[qid]
            seqpos[(g, r["DNA"][:19].upper())].update(ps)
            locpos[(g, r["Chromosome"], r["Position"], r["Direction"])].update(ps)

with open(os.path.join(OUT, "reached_positions.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["guide_id", "off_target_19mer", "positions"])
    for k in sorted(seqpos):
        w.writerow([k[0], k[1], ",".join(str(i) for i in sorted(seqpos[k]))])

# loci, endpoint-only loci, and which endpoint: the terminal one is the class that cannot be
# restored without moving the PAM boundary, so the split is what the paper's claim turns on
per = defaultdict(lambda: [0, 0, 0, 0])
for (g, _, _, _), ps in locpos.items():
    per[g][0] += 1
    if not (ps & EMITTED):
        per[g][1] += 1
        if ps == {1}:
            per[g][2] += 1
        elif ps == {SPACER}:
            per[g][3] += 1
with open(os.path.join(OUT, "reached_loci.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["guide_id", "loci", "unreachable_loci",
                "loci_position_1_only", "loci_position_20_only"])
    for g in sorted(per):
        w.writerow([g] + per[g])

print("record files            : %d" % len(files))
print("alignment rows read     : %d" % rows)
print("distinct (guide, 19mer) : %d" % len(seqpos))
print("distinct loci           : %d" % sum(v[0] for v in per.values()))
print("loci only at 1 or 20    : %d" % sum(v[1] for v in per.values()))
print("  of those, only at 1   : %d" % sum(v[2] for v in per.values()))
print("  of those, only at 20  : %d" % sum(v[3] for v in per.values()))
