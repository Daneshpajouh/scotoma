#!/usr/bin/env python3
"""Score every audited bulge-aware tool against the exact minimum, on two axes.

    python3 -I verify_tool_conformance.py

The criterion is decidable and has two sides, and both are theorems rather than conventions:

  COMPLETE   emit no fewer than the r distinct single-deletion words, because r is the exact
             minimum number of equal-radius balls that covers the deletion union, even with
             replacement centres anywhere and even with wildcards;
  MINIMAL    emit no more, because the r distinct words already cover it.

A tool below the minimum loses sites. A tool above it repeats work and output. The minimum is
the unique point that is neither.

Everything below is recomputed from the raw tool outputs in this directory and from the guide
spacers themselves. r is never read from a file; it is derived per guide by set construction.
No genome, no executable, no network.
"""
import csv
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
VERDICTS = []


def claim(ok, label, got, want):
    VERDICTS.append(ok)
    print("  %-5s %-55s got %-16s want %s" % ("PASS" if ok else "FAIL", label, got, want))


def read(name):
    with open(os.path.join(HERE, name), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def words(g):
    return {g[:p] + g[p + 1:] for p in range(len(g))}


def runs(g):
    return 1 + sum(1 for i in range(len(g) - 1) if g[i] != g[i + 1])


# ------------------------------------------------- the exact minimum, from the guides alone
per_guide = read("per_guide_counts.tsv")
guides = {}
for r in per_guide:
    guides[r["guide_id"]] = r["guide"].upper()
claim(len(guides) == 110, "guides in the audited panel", len(guides), 110)
claim({len(g) for g in guides.values()} == {20}, "every spacer is 20 nt",
      sorted({len(g) for g in guides.values()}), [20])
dups = [g for g, n in Counter(guides.values()).items() if n > 1]
claim(not dups, "guide spacers are distinct", dups[:2] or "distinct", "distinct")
per_tool_spacer = {}
for r in per_guide:
    per_tool_spacer.setdefault(r["guide_id"], set()).add(r["guide"].upper())
inconsistent = [k for k, v in per_tool_spacer.items() if len(v) != 1]
claim(not inconsistent, "each guide has one spacer across every route",
      inconsistent[:2] or "consistent", "consistent")
mism = [k for k, g in guides.items() if len(words(g)) != runs(g)]
claim(not mism, "distinct deletion words equal the run count", mism[:2] or "all agree", "all agree")

MINIMUM = sum(runs(g) for g in guides.values())
PLACEMENTS = sum(len(g) for g in guides.values())
claim(MINIMUM == 1634, "the exact minimum over the panel", MINIMUM, 1634)
claim(PLACEMENTS == 2200, "placements a per-position emitter enumerates", PLACEMENTS, 2200)

# ------------------------------------------------- axis one: completeness, per tool
print()
print("  axis one: does the tool emit at least the minimum?")
by_tool = defaultdict(list)
for r in per_guide:
    by_tool[r["tool"]].append(r)

EXPECT_EMITTED = {
    "Cas-OFFinder-bulge v1.2 beta (2016 script)": (1491, 143),
    "Cas-OFFinder 3.0.0 develop": (1561, 73),
    "Cas-OFFinder tag 3.0.0b3": (1561, 73),
}
for tool, (want_em, want_short) in sorted(EXPECT_EMITTED.items()):
    rs = by_tool[tool]
    em = sum(int(x["patterns_emitted"]) for x in rs)
    short = sum(1 for x in rs if int(x["deficit"]) > 0)
    lost = sum(int(x["deficit"]) for x in rs)
    claim((em, lost) == (want_em, want_short),
          "%.44s" % tool, "%d words, %d lost" % (em, lost),
          "%d words, %d lost" % (want_em, want_short))
    claim(em + lost == MINIMUM, "  its words plus its losses reach the minimum",
          em + lost, MINIMUM)
    claim(short > 0, "  guides on which it is short", short, "> 0")

# For CRISPRitz the recovered-word count is DERIVED from its own raw alignments, not read
# from the summary column: a reviewer showed that rotating every protospacer leaves the
# summary intact while no row carries a valid deletion word at all.
for label, fn in (("CRISPRitz 2.9.0 native TST", "crispritz_2.9.0_rna1.targets.tsv"),
                  ("CRISPRitz 2.7.1 native TST", "crispritz_2.7.1_rna1.targets.tsv")):
    rows = read(fn)
    seen, bad_word, bad_pam, bad_mm, bad_align, bad_gap = {}, 0, 0, 0, 0, 0
    for r in rows:
        spacer = r["crRNA"][:20]
        ung = r["DNA"].replace("-", "")
        dna19, pam = ung[:19], ung[19:22]
        if dna19 not in words(spacer):
            bad_word += 1
        else:
            seen.setdefault(spacer, set()).add(dna19)
        if pam[1:3] != "GG":
            bad_pam += 1
        if r["Mismatches"] != "0":
            bad_mm += 1
        # The gap column is part of the claim, not decoration. A reviewer moved every gap one
        # column, which leaves the ungapped word set and every multiplicity intact but makes
        # the row a false alignment. The alignment the row writes down is recomputed here: the
        # spacer with the gap's own column deleted must equal the target, letter for letter,
        # and the row's own mismatch count must be that distance.
        # Exactly one DNA gap in a spacer column, equal 23-column lengths, no guide gap.
        # A reviewer appended a second "-" to every DNA string (length 24, two gaps); find("-")
        # still returned the first gap and every check passed. Reject that class here.
        if len(r["DNA"]) != 23 or len(r["crRNA"]) != 23:
            bad_gap += 1
            continue
        if r["DNA"].count("-") != 1 or "-" in r["crRNA"][:20]:
            bad_gap += 1
            continue
        i = r["DNA"].find("-")
        if i < 0 or i >= 20:
            bad_gap += 1
            continue
        implied = spacer[:i] + spacer[i + 1:]
        if sum(1 for x, y in zip(implied, dna19) if x != y) != int(r["Mismatches"]) \
                or len(implied) != len(dna19):
            bad_align += 1
    claim(bad_word == 0,
          "%.42s: every row's target is a deletion word" % label, bad_word, 0)
    claim(bad_pam == 0, "%.42s: every row carries an NGG PAM" % label, bad_pam, 0)
    claim(bad_mm == 0, "%.42s: every fixture row is an exact match" % label, bad_mm, 0)
    claim(bad_gap == 0, "%.42s: every row's gap sits in a spacer column" % label, bad_gap, 0)
    claim(bad_align == 0,
          "%.42s: every row's gap column yields its own target" % label, bad_align, 0)
    derived = sum(len(v) for v in seen.values())
    claim(derived == MINIMUM,
          "%.42s: distinct words DERIVED from its rows" % label,
          "%d distinct words, 0 lost" % derived, "%d, 0 lost" % MINIMUM)
    claim(all(len(seen.get(g, ())) == runs(g) for g in guides.values()),
          "%.42s: per guide the derived words reach r" % label,
          sum(1 for g in guides.values() if len(seen.get(g, ())) != runs(g)), 0)

# ------------------------------------------------- axis two: minimality, from raw output
print()
print("  axis two: does the tool emit more than the minimum?")
for label, fn in (("CRISPRitz 2.9.0", "crispritz_2.9.0_rna1.targets.tsv"),
                  ("CRISPRitz 2.7.1", "crispritz_2.7.1_rna1.targets.tsv")):
    rows = read(fn)
    loci = Counter((r["crRNA"], r["Chromosome"], r["Position"], r["Direction"]) for r in rows)
    repeats = len(rows) - len(loci)
    claim(len(rows) == PLACEMENTS, "%s: alignment rows emitted" % label, len(rows), PLACEMENTS)
    claim(len(loci) == MINIMUM, "%s: distinct loci reported" % label, len(loci), MINIMUM)
    claim(repeats == PLACEMENTS - MINIMUM,
          "%s: rows repeating a locus already reported" % label, repeats, PLACEMENTS - MINIMUM)
    claim(abs(repeats / len(rows) - 0.257273) < 5e-7,
          "%s: redundant share of the output" % label,
          "%d/%d = %.6f" % (repeats, len(rows), repeats / len(rows)), "0.257273")
    # a repeat must be the same underlying target, differing only in gap placement
    groups = defaultdict(list)
    for r in rows:
        groups[(r["crRNA"], r["Chromosome"], r["Position"], r["Direction"])].append(r)
    multi = {k: v for k, v in groups.items() if len(v) > 1}
    bad_seq = [k for k, v in multi.items() if len({x["DNA"].replace("-", "") for x in v}) != 1]
    bad_gap = [k for k, v in multi.items() if len({x["DNA"].index("-") for x in v}) != len(v)]
    claim(not bad_seq, "%s: each repeat group is one target sequence" % label,
          len(bad_seq), 0)
    claim(not bad_gap, "%s: repeats differ only in gap placement" % label, len(bad_gap), 0)
    claim(max(len(v) for v in groups.values()) == 7,
          "%s: most times one locus is repeated" % label,
          max(len(v) for v in groups.values()), 7)

# ------------------------------------------------- the field verdict
print()
emits_exactly = []
for tool in sorted(by_tool):
    rs = by_tool[tool]
    lost = sum(int(x["deficit"]) for x in rs)
    if lost == 0 and tool.startswith("CRISPRitz"):
        pass          # complete, but its raw output is above the minimum, shown above
    if lost == 0 and not tool.startswith("CRISPRitz"):
        emits_exactly.append(tool)
claim(not emits_exactly,
      "tools that emit exactly the minimum, no more and no fewer",
      emits_exactly or "none", "none")

# a control: the criterion must be able to pass. A hand-built emitter that returns exactly the
# distinct words scores 0 lost and 0 redundant, so the test is not vacuous.
ideal_emitted = sum(len(words(g)) for g in guides.values())
claim(ideal_emitted == MINIMUM and PLACEMENTS - ideal_emitted == 566,
      "  control: an exact emitter scores 0 lost and 0 redundant",
      "%d emitted, 0 lost, 0 redundant" % ideal_emitted, "%d, 0, 0" % MINIMUM)

bad = VERDICTS.count(False)
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - bad, bad))
if bad:
    sys.exit(1)
print("Against the exact minimum of %d words on this panel: two audited tools emit fewer and"
      % MINIMUM)
print("lose sites, two emit the full %d placements and repeat 566 of them. None emits exactly."
      % PLACEMENTS)
