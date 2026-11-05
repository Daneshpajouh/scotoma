#!/usr/bin/env python3
"""The higher-budget panel claim, checked from the deposited records.

    python3 -I verify_budget_extend.py

Paper A reports that raising the substitution budget does not close the released build's
deficit: over the 110-guide panel and the whole GRCh38 primary assembly the corrected-only
counts at five, six and seven substitutions are 253,519, 1,546,525 and 6,275,358 loci, with
released-only zero throughout and every locus reached by the terminal deletion word.

The two 77 GB alignment tables those counts come from are not redistributed. What is deposited
is the output of two programs that read them, and they share no code:

  shard_<i>_record.json   written by panel_extend.py, the producer, which runs both published
                          binaries on the staged assembly and tallies loci per guide
  shard_<i>_VERIFY.json   written by verify_panel.py, an independent reader that re-derives the
                          locus convention from its definition and recomputes the terminal
                          attribution from each alignment's own genomic sequence, with a
                          forgery control and a corruption control

This checker re-derives every published total from the per-guide records, asserts that the two
independent programs agree guide by guide, and asserts that the controls in the verifier's own
output behaved as they must. A disagreement anywhere is a failure, not a warning.
"""
import copy
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
BUDGETS = ["5", "6", "7"]
PUBLISHED = {"5": 253519, "6": 1546525, "7": 6275358}
PUBLISHED_GUIDES_WITH_DEFICIT = 73
VERDICTS = []


def claim(ok, label, got, want):
    VERDICTS.append(bool(ok))
    print("  %-5s %-56s got %-16s want %s" % ("PASS" if ok else "FAIL", label, got, want))


def load():
    prod, ver = {}, {}
    for i in range(8):
        r = json.load(open(HERE / ("shard_%d_record.json" % i)))
        v = json.load(open(HERE / ("shard_%d_VERIFY.json" % i)))
        for g, d in r["guides"].items():
            prod[g] = {b: d["mm" + b] for b in BUDGETS}
        for b in BUDGETS:
            for g, d in v["budgets"][b]["per_guide"].items():
                ver.setdefault(g, {})[b] = d
        yield_v = v
        ver.setdefault("__controls__", []).append(
            (i, yield_v["verdict"], yield_v["controls"]))
    return prod, ver


prod, ver = load()
controls = ver.pop("__controls__")

claim(len(prod) == 110, "guides in the producer's records", len(prod), 110)
claim(len(ver) == 110, "guides in the verifier's records", len(ver), 110)
claim(set(prod) == set(ver), "the two programs cover the same guides",
      "identical" if set(prod) == set(ver) else "differ", "identical")

# guide-by-guide agreement between the two independent programs
fields = [("released", "released"), ("corrected", "corrected"),
          ("corrected_only", "corrected_only"), ("released_only", "released_only"),
          ("terminal_attributable", "terminal_attributable")]
disagree = []
for g in sorted(prod):
    for b in BUDGETS:
        for pk, vk in fields:
            if prod[g][b][pk] != ver[g][b][vk]:
                disagree.append((g, b, pk, prod[g][b][pk], ver[g][b][vk]))
claim(not disagree, "producer and independent verifier agree on every guide, every budget",
      "%d disagreement(s)" % len(disagree), 0)
if disagree:
    for d in disagree[:10]:
        print("      ", d)

# the published totals, re-derived here rather than read from the aggregate file
for b in BUDGETS:
    tot = sum(prod[g][b]["corrected_only"] for g in prod)
    claim(tot == PUBLISHED[b], "corrected-only total at budget %s" % b, tot, PUBLISHED[b])
    ro = sum(prod[g][b]["released_only"] for g in prod)
    claim(ro == 0, "  released-only total at budget %s" % b, ro, 0)
    ta = sum(prod[g][b]["terminal_attributable"] for g in prod)
    claim(ta == tot, "  and every one is reached by the terminal word", ta, tot)
    nd = sum(1 for g in prod if prod[g][b]["corrected_only"] > 0)
    claim(nd == PUBLISHED_GUIDES_WITH_DEFICIT, "  guides carrying a deficit at budget %s" % b,
          nd, PUBLISHED_GUIDES_WITH_DEFICIT)

# the deficit set is a property of the guide, not of the budget
sets = {b: {g for g in prod if prod[g][b]["corrected_only"] > 0} for b in BUDGETS}
same = sets["5"] == sets["6"] == sets["7"]
claim(same, "the same 73 guides carry the deficit at all three budgets",
      "identical" if same else "differ", "identical")
claim(len(prod) - len(sets["5"]) == 37, "  and 37 guides carry none",
      len(prod) - len(sets["5"]), 37)

# the aggregate file must agree with what was just re-derived, or it is stale
agg = json.load(open(HERE / "PANEL_AGGREGATE.json"))
rows = {str(r["budget"]).replace("mm", ""): r for r in agg["rows"]}
ok = all(rows[b]["corrected_only"] == PUBLISHED[b] for b in BUDGETS)
claim(ok, "the deposited aggregate agrees with the re-derivation",
      "agrees" if ok else "stale", "agrees")

# the independent verifier's own controls
for i, verdict, c in controls:
    claim(verdict == "VERIFY_OK", "shard %d verifier verdict" % i, verdict, "VERIFY_OK")
    claim(not c["C1_C2_violations"], "  shard %d: no malformed alignment rows" % i,
          len(c["C1_C2_violations"]), 0)
    claim(not c["C3_released_only"], "  shard %d: released set is a subset of corrected" % i,
          len(c["C3_released_only"]), 0)
    claim(not c["C4_first_word_explains_all"],
          "  shard %d: forgery control, first word explains nothing" % i,
          len(c["C4_first_word_explains_all"]), 0)
    c5 = c["C5_corruption"]
    claim(c5.get("rejected") is True and c5.get("corrupted_terminal_attributed_extra") == 0
          and c5.get("real_terminal_attributed", 0) > 0,
          "  shard %d: corruption control rejected a corrupted guide" % i,
          "%s attributed %s, corrupted %s" % (c5.get("guide"),
                                              c5.get("real_terminal_attributed"),
                                              c5.get("corrupted_terminal_attributed_extra")),
          "real > 0, corrupted 0")

# control: the agreement test above must be able to fail
forged = copy.deepcopy(prod)
g0 = sorted(forged)[0]
forged[g0]["7"]["corrected_only"] += 1
broke = any(forged[g][b][pk] != ver[g][b][vk]
            for g in forged for b in BUDGETS for pk, vk in fields)
claim(broke, "control: a single altered count breaks the agreement test",
      "broken" if broke else "held", "broken")

n = VERDICTS.count(False)
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - n, n))
if n:
    sys.exit(1)
print("At five, six and seven substitutions the released build misses 253,519, 1,546,525 and")
print("6,275,358 loci that the corrected build returns, it gains none, and every missed locus")
print("is reached by the terminal deletion word the released build never emits.")
