#!/usr/bin/env python3
"""What CRISPRseek's own bulge route emits, from the patterns it actually wrote.

    python3 -I verify_crisprseek_executed.py

The companion previously described this route by reading its source. This directory holds the
route executed: `run_crisprseek.R` loads CRISPRseek, takes the bulge generator the package
itself returns for its default backend, and calls it on `input.txt`, the 110 audited spacers at
one RNA bulge. The generator's 'D' mode writes its expanded pattern file and returns before any
binary runs, so the whole expansion is observable with no genome and no executable.

`input_bulge.txt` is that file, as written. Everything below is recomputed from it.
"""
import csv, hashlib, os, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.join(HERE, os.pardir, "tool_conformance", "per_guide_counts.tsv")
# the file the Bioconductor release ships, as installed and sourced by the R script here
WRAPPER_SHA = "15ba0d6470f630be7291c3d5e87f597df2b978cc1f1419d28352cf8ddbaf5f56"
VERDICTS = []


def claim(ok, label, got, want):
    VERDICTS.append(ok)
    print("  %-5s %-54s got %-16s want %s" % ("PASS" if ok else "FAIL", label, got, want))


def runs(g):
    return 1 + sum(1 for i in range(len(g) - 1) if g[i] != g[i + 1])


def words(g):
    return {g[:i] + g[i + 1:] for i in range(len(g))}


guides, deposited = {}, {}
with open(PANEL, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        guides.setdefault(r["guide_id"], r["guide"].upper())
        if r["tool"].startswith("Cas-OFFinder-bulge v1.2"):
            deposited[r["guide_id"]] = int(r["patterns_emitted"])
claim(len(guides) == 110, "audited spacers", len(guides), 110)

# ---- the fixture really is those spacers, in the generator's own input format
fx = [l.strip() for l in open(os.path.join(HERE, "input.txt")) if l.strip()]
claim(fx[1] == "NNNNNNNNNNNNNNNNNNNNNRG 0 1",
      "fixture asks for one RNA bulge and no DNA bulge", fx[1], "NNNN...RG 0 1")
asked = {l.split()[0][:-3] for l in fx[2:]}
claim(asked == set(guides.values()), "fixture carries exactly the audited spacers",
      "identical" if asked == set(guides.values()) else "differs", "identical")
mm = {l.split()[1] for l in fx[2:]}
claim(len(mm) == 1, "every line shares one mismatch value, as the format requires",
      sorted(mm), "one value")

# ---- what the generator wrote
lines = [l.strip() for l in open(os.path.join(HERE, "input_bulge.txt")) if l.strip()]
emitted_raw = [l.split()[0] for l in lines[2:]]
claim(len(emitted_raw) == len(set(emitted_raw)),
      "the generator wrote no pattern twice", len(emitted_raw) - len(set(emitted_raw)), 0)
ungapped = [p for p in emitted_raw if not p[:-3].startswith("N")]
deletion = [p[:-3].lstrip("N") for p in emitted_raw if p[:-3].startswith("N")]
claim(len(ungapped) == 110, "one ungapped pattern per spacer", len(ungapped), 110)
claim(all(len(w) == 19 for w in deletion),
      "every deletion pattern spells a 19-mer", Counter(map(len, deletion)).most_common(1), "19")

# ---- per guide, from the loop the shipped file implements
# A shared word must NOT be credited to every guide that could own it: five of these 19-mers
# are deletion words of more than one spacer, and crediting each to all of its owners
# overstates four guides and the total. The attribution is therefore by position, per guide.
owner = {}
for g, s in guides.items():
    for w in words(s):
        owner.setdefault(w, set()).add(g)
shared = {w for w, gs in owner.items() if len(gs) > 1}
claim(len(shared) == 5, "19-mers that are deletion words of several spacers", len(shared), 5)

WROTE = set(deletion)
per, missing_pos = {}, Counter()
for g, s in guides.items():
    visited = {s[:i] + s[i + 1:] for i in range(1, len(s) - 1)}   # range(1, len-bulge)
    per[g] = visited & words(s)
    for i in range(len(s)):
        w = s[:i] + s[i + 1:]
        if w not in per[g]:
            missing_pos[i + 1] += 1
    claim_ok = per[g] <= WROTE
    if not claim_ok:
        claim(False, "  %s: its words are all in the written file" % g, "missing", "present")
out = sum(1 for g in guides for w in per[g] if w not in WROTE)
claim(out == 0, "every word the loop reaches is in the written file", out, 0)
# Subset is not enough: a reviewer added one valid terminal word that the loop never reaches
# and the file still passed. The two sets must be equal, so an extra word is caught too.
reached_union = set().union(*per.values())
claim(WROTE == reached_union,
      "  and the written file contains nothing else",
      "%d written, %d reached" % (len(WROTE), len(reached_union)), "equal sets")
# The thresholds are part of what was run. Rewriting them all to zero changed the search and
# passed, because nothing read them.
budgets = {l.split()[1] for l in lines[2:] if len(l.split()) > 1}
claim(budgets == {"4"}, "  every written line carries the fixture's own budget",
      sorted(budgets), ["4"])

tot = sum(len(v) for v in per.values())
need = sum(runs(s) for s in guides.values())
short = [g for g in guides if len(per[g]) < runs(guides[g])]
claim(need == 1634, "distinct deletion words the panel requires", need, 1634)
claim(tot == 1491, "distinct deletion words this route emits", tot, 1491)
claim(len(short) == 98, "spacers served with fewer words than they need", len(short), 98)
dis = [g for g in guides if len(per[g]) != deposited[g]]
claim(not dis, "  agrees with the deposited standalone-script row on every spacer",
      len(dis) or "all 110 agree", "all 110 agree")

claim(set(missing_pos) == {1, 20},
      "the only spacer positions ever missed", sorted(missing_pos), [1, 20])
claim(missing_pos[1] == 70 and missing_pos[20] == 73,
      "  how many spacers lose each endpoint",
      (missing_pos[1], missing_pos[20]), (70, 73))
# The cross-check that catches a mis-attributed word: the two endpoint counts must add up to
# the shortfall exactly. An earlier pass credited each shared 19-mer to every spacer that could
# own it, which gave 69 at position 20 and a sum of 139 against a shortfall of 143.
claim(missing_pos[1] + missing_pos[20] == need - tot,
      "  and they sum to the shortfall",
      "%d + %d = %d" % (missing_pos[1], missing_pos[20], missing_pos[1] + missing_pos[20]),
      need - tot)

# ---- controls
mid = sum(1 for g in guides for i in range(1, 19)
          if guides[g][:i] + guides[g][i + 1:] not in per[g])
claim(mid == 0, "  control: no interior position is ever missed", mid, 0)
rescued = [g for g in guides if len(per[g]) == runs(guides[g])]
claim(len(rescued) == 12,
      "  control: spacers that lose nothing, by the boundary law", len(rescued), 12)

# ---- the version gate the package applies before it will run this route at all
# The gate splits the Python version and compares the minor part as a STRING, so "13" >= "5"
# is false and every Python from 3.10 on is refused as below 3.5.
def gate(v):
    p = v.split(".")
    return len(p) >= 2 and p[0] == "3" and p[1] >= "5"


claim(gate("3.9") and gate("3.5"), "the version gate admits Python 3.5 and 3.9",
      (gate("3.5"), gate("3.9")), (True, True))
claim(not gate("3.10") and not gate("3.13"),
      "  and refuses 3.10 and 3.13, comparing the minor part as text",
      (gate("3.10"), gate("3.13")), (False, False))

print()
w = os.environ.get("CRISPRSEEK_WRAPPER")
if w and os.path.exists(w):
    got = hashlib.sha256(open(w, "rb").read()).hexdigest()
    claim(got == WRAPPER_SHA, "the installed bundled wrapper matches the pinned digest",
          got[:12] + "...", WRAPPER_SHA[:12] + "...")
else:
    print("  NOTE  set CRISPRSEEK_WRAPPER to the installed")
    print("        extdata/Cas-OFFinder/cas-offinder-bulge to bind the digest %s" % WRAPPER_SHA[:12])

bad = VERDICTS.count(False)
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - bad, bad))
if bad:
    sys.exit(1)
print("Executed through R, CRISPRseek's default bulge route emits 1,491 of the 1,634 distinct")
print("deletion words this panel requires, falling short on 98 of 110 spacers, and every word")
print("it misses sits at spacer position 1 or 20.")
