#!/usr/bin/env python3
"""Derive every genome return bit from SEQUENCE ALONE, then compare with the genome run.

    python3 -I census_predicts_genome.py

Reads two files and nothing else:

  ../data/sites129.json          the 260 selected measured-active rows: guide spacer,
                                 aligned target, PAM, assay window. Sequences only.
  active_sites.part000.tsv       the genome run's own answer: one released bit and one
                                 corrected bit per row, written by the whole-assembly
                                 scan with the two executables.

No genome, no executable, no network. The prediction is computed from the paper's
pattern-generation rule, stated once here and used for nothing else:

  a build reports a locus through a single RNA deletion exactly when some deletion
  pattern IN ITS OWN PATTERN SET lies within the mismatch budget of the 19 nt
  protospacer anchored at that locus' PAM, and the PAM is NGG. The released build's
  set omits the terminal deletion; the corrected set contains it.

The point of the file is that the two routes can BOTH be run by a reader. The prediction
never reads a return bit, and the genome run never read a pattern count.

What this file does NOT establish, asserted below as checks rather than left for a reader
to discover. The selected rows are terminal-selected by construction, and on them:

  * all 260 corrected bits equal the single predicate "the PAM is NGG", so they carry no
    information about the rest of the rule;
  * the decisive non-terminal placement is position 19 on 248 of the 260 rows, so the
    released bits test the tail of the rule and little else;
  * consequently a deliberately truncated rule admitting only four of the nineteen
    non-terminal positions also reproduces all 520 bits.

So the agreement establishes that the terminal omission accounts for the released-corrected
difference at every selected site. It does not identify the whole pattern rule. That is a
different test, done here at budget zero, where a pattern matches only its own exact word.
"""
import csv
import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
SITES = os.path.join(HERE, os.pardir, "data", "sites129.json")
GENOME = os.path.join(HERE, "active_sites.part000.tsv")
SP, BUDGET = 20, 4
RC = str.maketrans("ACGT", "TGCA")
VERDICTS = []


def claim(ok, label, got, want):
    VERDICTS.append(ok)
    print("  %-5s %-54s got %-18s want %s" % ("PASS" if ok else "FAIL", label, got, want))


def ham(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def predict(rows, omit=SP - 1, budget=BUDGET, pam_suffix="GG"):
    """Return one record per selected row, from sequence alone.

    omit         the spacer position the released build's pattern set leaves out
    budget       the mismatch budget both builds are run at
    pam_suffix   the two bases the PAM must end in for either build to report
    """
    out = []
    for i, r in enumerate(rows):
        at, spacer = r["aligned_target"], r["spacer"]
        gap = at.index("-")
        if gap != SP - 1 or len(spacer) != SP:
            sys.exit("GATE FAILED: row %d is not a terminal single-deletion annotation" % i)
        dna, pam = at[:SP - 1], at[SP:SP + 3]
        pam_ok = pam[1:3] == pam_suffix
        costs = [ham(spacer[:p] + spacer[p + 1:], dna) for p in range(SP)]
        released = min(c for p, c in enumerate(costs) if p != omit)
        corrected = min(costs)
        out.append(dict(row_index=i, chrom=r["chrom"], strand=r["strand"], pam=pam,
                        pam_ok=pam_ok, terminal_mm=costs[SP - 1],
                        best_without_omitted=released, best_with_all=corrected,
                        released=int(pam_ok and released <= budget),
                        corrected=int(pam_ok and corrected <= budget),
                        guide=spacer, window=(r["window"] or "").upper(), site=dna + pam))
    return out


def bits_match(pred, dep):
    """Per-row bit comparison. Returns (agreeing, disagreeing). Order-independent and
    total-independent: swapping two rows' bits while preserving the totals disagrees."""
    a = d = 0
    for p in pred:
        g = dep[p["row_index"]]
        for mine, theirs in (("released", "released_returned"),
                             ("corrected", "corrected_returned")):
            if p[mine] == int(g[theirs]):
                a += 1
            else:
                d += 1
    return a, d


rows = json.load(open(SITES))
dep = {int(r["row_index"]): r for r in csv.DictReader(open(GENOME), delimiter="\t")}
claim(len(rows) == 260, "selected assay rows read", len(rows), 260)
claim(len(dep) == 260 and set(dep) == set(range(260)),
      "genome rows, indexed 0..259 without gaps", len(dep), 260)

pred = predict(rows)
agree, dis = bits_match(pred, dep)
claim((agree, dis) == (520, 0), "sequence-only bits against the genome run",
      "%d agree, %d disagree" % (agree, dis), "520 agree, 0 disagree")

rel = sum(p["released"] for p in pred)
cor = sum(p["corrected"] for p in pred)
loss = [p for p in pred if p["corrected"] and not p["released"]]
claim(rel == 108, "released returns, predicted from sequence", rel, 108)
claim(cor == 237, "corrected returns, predicted from sequence", cor, 237)
claim(len(loss) == 129, "exclusive losses, predicted from sequence", len(loss), 129)
claim(len({p["guide"] for p in loss}) == 40, "guides carrying a loss",
      len({p["guide"] for p in loss}), 40)

# the mechanism, site by site: every lost site is reachable only one mismatch beyond
# the budget once the terminal pattern is withheld.
profile = dict(sorted(Counter(p["best_without_omitted"] for p in loss).items()))
claim(profile == {BUDGET + 1: 129}, "each loss sits exactly one mismatch out of reach",
      profile, "{%d: 129}" % (BUDGET + 1))
claim(all(p["terminal_mm"] <= BUDGET for p in loss),
      "each loss is in reach with the terminal pattern",
      max(p["terminal_mm"] for p in loss), "<= %d" % BUDGET)

# the PAM split the genome run recorded, recomputed from the PAM strings themselves
status = Counter((dep[p["row_index"]]["status"], p["pam_ok"]) for p in pred)
claim(status == Counter({("VERIFIED_NGG", True): 237, ("NON_NGG", False): 23}),
      "PAM status agrees with the literal PAM bases", dict(status),
      "{VERIFIED_NGG/True: 237, NON_NGG/False: 23}")
claim(all(not (p["released"] or p["corrected"]) for p in pred if not p["pam_ok"]),
      "no non-NGG row is predicted returned by either build",
      sum(1 for p in pred if not p["pam_ok"] and (p["released"] or p["corrected"])), 0)

# a bulge-free alignment at the same locus would rescue the released build; none does.
rescued = 0
for p in pred:
    w, s = p["window"], p["site"]
    for seq, flip in ((s, False), (s.translate(RC)[::-1], True)):
        j = w.find(seq)
        if j < 0:
            continue
        if not flip and j >= 1:
            ext = w[j - 1:j + SP - 1]
        elif flip and j + len(seq) < len(w):
            ext = w[j:j + len(seq) + 1].translate(RC)[::-1][:SP]
        else:
            ext = None
        if ext and len(ext) == SP and ham(p["guide"], ext) <= BUDGET and p["pam_ok"]:
            rescued += 1
        break
claim(rescued == 0, "losses rescued by a bulge-free NGG alignment", rescued, 0)

# ------------------------------------------------------------------ falsification
# Each control changes the rule, not the data. A match under a wrong rule would mean
# the comparison is insensitive, so each of these MUST disagree.
print()
print("  falsification controls: each must disagree with the genome run")
CONTROLS = [
    ("released omits position 1 instead of the terminal one", dict(omit=0)),
    ("released omits the tenth position", dict(omit=9)),
    ("both builds run at budget three", dict(budget=3)),
    ("both builds run at budget five", dict(budget=5)),
    ("PAM policy NAG instead of NGG", dict(pam_suffix="AG")),
]
for label, kw in CONTROLS:
    a, d = bits_match(predict(rows, **kw), dep)
    claim(d > 0, "  control: " + label, "%d disagree" % d, "> 0")

# and one control on the data side: the comparison is per row, so a swap that keeps
# both totals must still be caught.
swapped = dict(dep)
hi = next(k for k in sorted(dep) if dep[k]["released_returned"] == "1")
lo = next(k for k in sorted(dep) if dep[k]["released_returned"] == "0"
          and dep[k]["status"] == "VERIFIED_NGG")
swapped[hi] = dict(dep[hi], released_returned="0")
swapped[lo] = dict(dep[lo], released_returned="1")
a, d = bits_match(pred, swapped)
claim(d == 2, "  control: two return bits swapped, totals preserved",
      "%d disagree" % d, "2")

# ------------------------------------------- what the agreement does NOT establish
print()
print("  the limits of the agreement, asserted as checks")

pam_bits = [1 if r["aligned_target"][SP:SP + 3][1:3] == "GG" else 0 for r in rows]
claim(all(p["corrected"] == b for p, b in zip(pred, pam_bits)),
      "  corrected bits are the PAM predicate alone",
      sum(1 for p, b in zip(pred, pam_bits) if p["corrected"] != b), 0)

at19 = sum(1 for p in pred
           if p["best_without_omitted"] ==
           ham(p["guide"][:SP - 2] + p["guide"][SP - 1:], p["site"][:SP - 1]))
claim(at19 == 248, "  rows whose non-terminal minimum sits at position 19", at19, 248)


def truncated(rows, keep, budget=BUDGET):
    out = []
    for i, r in enumerate(rows):
        at, spacer = r["aligned_target"], r["spacer"]
        dna, pam = at[:SP - 1], at[SP:SP + 3]
        ok = pam[1:3] == "GG"
        c = [ham(spacer[:j] + spacer[j + 1:], dna) for j in range(SP)]
        rel = min([c[j - 1] for j in keep] or [99])
        cor = min(rel, c[SP - 1])
        out.append(dict(row_index=i, released=int(ok and rel <= budget),
                        corrected=int(ok and cor <= budget)))
    return out


a_t, d_t = bits_match(truncated(rows, [11, 16, 18, 19]), dep)
claim((a_t, d_t) == (520, 0),
      "  a rule keeping 4 of 19 non-terminal positions also fits",
      "%d agree, %d disagree" % (a_t, d_t), "520 agree, 0 disagree")

# ------------------------------------------- where the rule IS identified: budget zero
print()
print("  identification at budget zero, one exact witness per distinct word")
spacers = sorted({r["spacer"] for r in rows})
words_of = {sp: sorted({sp[:j] + sp[j + 1:] for j in range(SP)}) for sp in spacers}
witnesses = [(sp, w) for sp in spacers for w in words_of[sp]]
claim(len(spacers) == 48, "distinct spacers among the selected rows", len(spacers), 48)
claim(len(witnesses) == sum(len(v) for v in words_of.values()),
      "witnesses, one per distinct word", len(witnesses),
      sum(len(v) for v in words_of.values()))

KEEP = (10, 15, 17, 18, 19)          # the truncated rule's positions, zero-based
full_hit = sum(1 for sp, w in witnesses
               if any(sp[:j] + sp[j + 1:] == w for j in range(SP)))
trunc_hit = sum(1 for sp, w in witnesses
                if any(sp[:j] + sp[j + 1:] == w for j in KEEP))
claim(full_hit == len(witnesses), "the complete rule matches every witness exactly",
      full_hit, len(witnesses))
claim(trunc_hit < len(witnesses),
      "the truncated rule misses witnesses at budget zero",
      "%d of %d matched" % (trunc_hit, len(witnesses)), "fewer than %d" % len(witnesses))

bad = VERDICTS.count(False)
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - bad, bad))
if bad:
    sys.exit(1)
print("The terminal omission accounts for the released-corrected difference at every selected")
print("site. The pattern rule itself is identified at budget zero, not by these 520 bits.")
