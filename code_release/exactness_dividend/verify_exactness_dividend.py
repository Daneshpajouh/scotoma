#!/usr/bin/env python3
"""How much of a per-position emitter's deletion work is provably redundant.

    python3 -I verify_exactness_dividend.py

Reads one file, `published_benchmark_guides.txt`, the 61 spacers of an independent published
tool-comparison benchmark, and nothing else. No genome, no executable, no network.

A per-position emitter enumerates one deletion placement per spacer position, n of them. The
number of DISTINCT deletion words is the run count r, and r is also the minimum number of
radius-m balls that cover the deletion union, so no family of any kind is smaller. The gap
n - r is therefore work no correct emitter needs to do, and it is not a heuristic estimate.

Provenance, because it matters and is not clean: the benchmark repository the preprint names is
not public under its authors' account. This guide list is the copy carried by a public
derivative repository that retains the original project file. The preprint states 61 sgRNAs and
this file holds exactly 61, all 20 nt, which is the only check available on the list itself.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VERDICTS = []


def claim(ok, label, got, want):
    VERDICTS.append(ok)
    print("  %-5s %-52s got %-16s want %s" % ("PASS" if ok else "FAIL", label, got, want))


guides = [l.strip().upper() for l in open(os.path.join(HERE, "published_benchmark_guides.txt"))
          if l.strip()]
claim(len(guides) == 61, "guides in the published benchmark", len(guides), 61)
claim({len(g) for g in guides} == {20}, "every spacer is 20 nt",
      sorted({len(g) for g in guides}), [20])
claim(set("".join(guides)) <= set("ACGT"), "spacers over the four letters only",
      sorted(set("".join(guides)) - set("ACGT")) or "clean", "clean")


def runs(g):
    return 1 + sum(1 for i in range(len(g) - 1) if g[i] != g[i + 1])


def words(g):
    return {g[:p] + g[p + 1:] for p in range(len(g))}


# the identity is re-derived per guide by set construction, not assumed
mismatched = [g for g in guides if len(words(g)) != runs(g)]
claim(not mismatched, "distinct deletion words equal the run count",
      mismatched[:2] or "all 61 agree", "all 61 agree")

placements = sum(len(g) for g in guides)
minimum = sum(runs(g) for g in guides)
redundant = placements - minimum
claim(placements == 1220, "placements a per-position emitter enumerates", placements, 1220)
claim(minimum == 891, "distinct words actually required", minimum, 891)
claim(redundant == 329, "provably redundant placements", redundant, 329)
claim(abs(redundant / placements - 0.269672) < 5e-7, "redundant share",
      "%d/%d = %.6f" % (redundant, placements, redundant / placements), "0.269672")

per = sorted(((len(g) - runs(g)), g) for g in guides)
claim((per[0][0], per[-1][0]) == (0, 12), "per-guide redundancy range",
      (per[0][0], per[-1][0]), (0, 12))
claim(per[-1][1] == "TTTTTTTGTTCCTCTTTTGG", "the worst guide",
      per[-1][1], "TTTTTTTGTTCCTCTTTTGG")

# exposure to the terminal omission: a terminal deletion is unique exactly when the last two
# spacer bases differ, so equal last bases make the terminal word reachable without it
at_risk = [g for g in guides if g[-1] != g[-2]]
claim(len(at_risk) == 47, "guides at risk under the terminal law", len(at_risk), 47)
claim(len(guides) - len(at_risk) == 14, "guides the law rescues",
      len(guides) - len(at_risk), 14)

# a control: the redundancy is a property of the guides, so a run-free guide must show none,
# and an all-one-letter guide must show the maximum
claim(len("ACGTACGTACGTACGTACGT") - runs("ACGTACGTACGTACGTACGT") == 0,
      "  control: an alternating spacer has no redundancy",
      len("ACGTACGTACGTACGTACGT") - runs("ACGTACGTACGTACGTACGT"), 0)
claim(len("A" * 20) - runs("A" * 20) == 19,
      "  control: a homopolymer spacer is 19 of 20 redundant",
      len("A" * 20) - runs("A" * 20), 19)

bad = VERDICTS.count(False)
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - bad, bad))
if bad:
    sys.exit(1)
print("On this published benchmark %d of %d enumerated placements, %.1f%%, are provably"
      % (redundant, placements, 100 * redundant / placements))
print("redundant: the exact minimum is %d distinct words." % minimum)
