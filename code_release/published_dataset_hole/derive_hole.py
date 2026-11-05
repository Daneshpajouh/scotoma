#!/usr/bin/env python3
"""Build the hole tables from the scan statistics and the two public CSVs.

    python3 -I derive_hole.py <scan_dir> <csv_dir> <out_dir>

<scan_dir> holds reached_positions.tsv, reached_loci.tsv and guides_meta.json; <csv_dir> holds
the two deposited CSVs of the published model. Two definitions do the work, and both are
stated here rather than assumed:

UNREACHABLE  a (guide, off-target 19-mer) pair the scan reached, no word of which sits at a
             spacer position the generator's loop visits; equivalently its position set misses
             2..19 entirely, so only a deletion at position 1 or 20 produces it.

HELD         a 19-mer occurs in the published files for a guide when it is a contiguous
             19-column window of any protospacer those files report for that guide, gap and
             pad characters removed. This is the conservative test. It counts a sequence as
             held even when the row holding it is ungapped or carries a DNA bulge, that is
             even when the files reach the site by an alignment other than the deletion one.
             The narrow alternative, counting only the 19-mers that the files' own gapped rows
             spell out, is also computed and reported, and it is the weaker test: it would
             call a site absent that the files plainly contain.
"""
import csv, hashlib, json, os, sys
from collections import Counter, defaultdict

SPACER, EMITTED = 20, set(range(2, 20))
SCAN, DATA, OUTDIR = sys.argv[1], sys.argv[2], sys.argv[3]
CSVS = {"CIRCLE_seq_10gRNA_wholeDataset.csv": ("sgRNA_seq", "off_seq"),
        "elevation_6gRNA_wholeDataset.csv": ("crRNA", "DNA")}

meta = json.load(open(os.path.join(SCAN, "guides_meta.json")))
spacer_of = {g: meta[g]["spacer"] for g in meta}
by_spacer = {meta[g]["spacer"]: g for g in meta}

held, strict, kinds = defaultdict(set), defaultdict(set), defaultdict(set)
for fn, (gc, oc) in CSVS.items():
    p = os.path.join(DATA, fn)
    print("%-40s sha256 %s" % (fn, hashlib.sha256(open(p, "rb").read()).hexdigest()))
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh):
            g, o = r[gc], r[oc]
            sp = g.replace("_", "").replace("-", "")[:SPACER]
            if sp not in by_spacer:
                continue
            gid = by_spacer[sp]
            proto = o.replace("-", "").replace("_", "").upper()[:-3]
            kind = {19: "gapped19", 20: "ungapped20", 21: "dnabulge21"}.get(len(proto))
            if kind is None:
                continue
            for i in range(len(proto) - 19 + 1):
                w = proto[i:i + 19]
                held[gid].add(w)
                kinds[(gid, w)].add(kind)
            if kind == "gapped19":
                strict[gid].add(proto)

reached = defaultdict(dict)
with open(os.path.join(SCAN, "reached_positions.tsv"), newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        reached[r["guide_id"]][r["off_target_19mer"]] = \
            {int(x) for x in r["positions"].split(",")}
loci = {r["guide_id"]: r for r in
        csv.DictReader(open(os.path.join(SCAN, "reached_loci.tsv")), delimiter="\t")}

rows, per, how = [], [], Counter()
for g in sorted(reached):
    comp = reached[g]
    unreach = sorted(s for s, ps in comp.items() if not (ps & EMITTED))
    for s in unreach:
        inf = s in held[g]
        rows.append([g, spacer_of[g], s, "yes" if inf else "no",
                     "+".join(sorted(kinds[(g, s)])) if inf else "-"])
        how["+".join(sorted(kinds[(g, s)])) if inf else "absent"] += 1
    words = meta[g]["words"]
    per.append(dict(
        guide_id=g, spacer=spacer_of[g], r=len(words),
        unsearched_words=sum(1 for _, ps in words.items() if not (set(ps) & EMITTED)),
        dataset_sequences=len(strict[g]), complete_sequences=len(comp),
        unreachable_sequences=len(unreach),
        unreachable_in_files=sum(1 for s in unreach if s in held[g]),
        unreachable_absent=sum(1 for s in unreach if s not in held[g]),
        absent_gapped_only=sum(1 for s in unreach if s not in strict[g]),
        loci=int(loci[g]["loci"]), unreachable_loci=int(loci[g]["unreachable_loci"])))

with open(os.path.join(OUTDIR, "unreachable_sequences.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["guide_id", "spacer", "off_target_19mer", "in_published_files", "held_by"])
    w.writerows(rows)
with open(os.path.join(OUTDIR, "per_guide_hole.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(per[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(per)

T = lambda k: sum(p[k] for p in per)
print()
for lab, v in (("guides", len(per)), ("distinct deletion words required", T("r")),
               ("words only at position 1 or 20", T("unsearched_words")),
               ("loci reached", T("loci")), ("loci only through those positions", T("unreachable_loci")),
               ("sequences, complete library", T("complete_sequences")),
               ("unreachable sequences", T("unreachable_sequences")),
               ("of those, in the published files", T("unreachable_in_files")),
               ("of those, absent outright", T("unreachable_absent")),
               ("absent under the narrow gapped-row test", T("absent_gapped_only"))):
    print("%-42s %d" % (lab, v))
print()
for k, n in how.most_common():
    print("  held by %-32s %d" % (k, n))
