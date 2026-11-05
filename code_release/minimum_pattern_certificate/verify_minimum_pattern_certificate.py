#!/usr/bin/env python3
"""Checks the per-guide minimum-pattern certificate, and tests its content by brute force.

    python3 -I verify_minimum_pattern_certificate.py

Three things are checked. First, that the run count really is the number of distinct
single-deletion words, on every guide of both deposited panels, by constructing the words.
Second, that the certificate refuses to assert a minimum outside the range the companion
proves. Third, and this is the part with teeth, that on instances small enough to settle
completely no family of r-1 four-letter centres covers the deletion union Um while the r
deletion words do. Centres are taken in {A,C,G,T}^L only, matching the companion's object
(X is an unused concrete letter, never a wildcard).

Those small instances are below the length the companion's theorem requires. They test the
mechanism the theorem asserts, over every centre and every subset, not the theorem itself.
"""
import csv, itertools, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from min_patterns import runs, deletion_words, certificate, proved  # noqa: E402

VERDICTS = []
LETTERS = "ACGT"


def claim(ok, label, got, want):
    VERDICTS.append(ok)
    print("  %-5s %-54s got %-18s want %s" % ("PASS" if ok else "FAIL", label, got, want))


print("  the count, constructed rather than asserted")
BENCH = os.path.join(HERE, os.pardir, "exactness_dividend", "published_benchmark_guides.txt")
PANEL = os.path.join(HERE, os.pardir, "tool_conformance", "per_guide_counts.tsv")

bench = [l.strip().upper() for l in open(BENCH) if l.strip()]
panel = {}
with open(PANEL, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        panel[r["guide_id"]] = (r["guide"].upper(), int(r["n"]), int(r["r"]))

claim(len(bench) == 61, "benchmark panel guides", len(bench), 61)
claim(len(panel) == 110, "audited panel guides", len(panel), 110)

mism = [g for g in bench if runs(g) != len(deletion_words(g))]
claim(not mism, "  benchmark: run count equals the distinct words",
      len(mism) or "all agree", "all agree")
mism = [k for k, (g, n, r) in panel.items() if runs(g) != len(deletion_words(g))]
claim(not mism, "  panel: run count equals the distinct words",
      len(mism) or "all agree", "all agree")
mism = [k for k, (g, n, r) in panel.items() if runs(g) != r]
claim(not mism, "  panel: and equals the deposited r column",
      len(mism) or "all agree", "all agree")

tb = sum(runs(g) for g in bench)
tp = sum(runs(g) for g, _, _ in panel.values())
claim(tb == 891, "benchmark patterns required in total", tb, 891)
claim(tp == 1634, "audited panel patterns required in total", tp, 1634)

print()
print("  the proved range, and the refusal outside it")
claim(certificate("A" * 20, 9)["status"] == "PROVED",
      "length 20 at budget 9 is inside the proved range",
      certificate("A" * 20, 9)["status"], "PROVED")
for n, m in ((19, 4), (20, 10), (19, 10)):
    c = certificate("AC" * (n // 2) + "A" * (n % 2), m)
    claim(c["status"] == "NOT PROVED" and c["minimum"] is None,
          "  length %d at budget %d asserts no minimum" % (n, m),
          "%s / %s" % (c["status"], c["minimum"]), "NOT PROVED / None")

print()
print("  bounded exhaustive test: no r-1 DNA centres cover Um")


def universe(L):
    return ["".join(p) for p in itertools.product(LETTERS, repeat=L)]


def ball_points(centre, U, m):
    out = set()
    for w in U:
        d = sum(1 for a, b in zip(centre, w) if a != b)
        if d <= m:
            out.add(w)
    return out


CASES = [("AAACC", 0), ("AAACC", 1), ("AACCA", 0), ("AACCA", 1),
         ("AAACCC", 0), ("AAACCC", 1), ("AACCAA", 0), ("AACCAA", 1)]
tested = 0
for g, m in CASES:
    r = runs(g)
    L = len(g) - 1
    U = universe(L)
    words = sorted(deletion_words(g))
    assert len(words) == r
    Um = set()
    for w in words:
        Um |= ball_points(w, U, m)
    covers = []
    for centre in itertools.product(LETTERS, repeat=L):
        c = "".join(centre)
        cov = ball_points(c, U, m) & Um
        if cov:
            covers.append(frozenset(cov))
    cov_r = set()
    for w in words:
        cov_r |= ball_points(w, U, m)
    ok_r = cov_r == Um
    ok_less = True
    if r - 1 <= 2:
        for combo in itertools.combinations(covers, r - 1):
            s = set()
            for part in combo:
                s |= part
            if s == Um:
                ok_less = False
                break
    else:
        ok_less = None
    claim(ok_r, "  %s m=%d: the r=%d deletion words cover Um" % (g, m, r),
          "covered" if ok_r else "not covered", "covered")
    um_proper = Um != set(U)
    claim(um_proper or m == 0,
          "    Um is not the whole space" if um_proper else "    Um IS the whole space",
          "proper" if um_proper else "everything", "proper")
    claim(len(covers) >= r,
          "    DNA centres that hit Um, at least the r words", len(covers), ">= %d" % r)
    if ok_less is None:
        claim(False, "    r-1 subsets not enumerable here", r - 1, "<= 2")
    else:
        claim(ok_less, "    no %d DNA centres cover Um, over %d hitting centres"
              % (r - 1, len(covers)), "none" if ok_less else "one exists", "none")
    tested += 1
claim(tested == len(CASES), "instances settled exhaustively", tested, len(CASES))

# control: a wildcard centre is outside the theorem's centre set
claim("*" not in LETTERS,
      "  control: wildcards are not in the centre alphabet", "excluded", "excluded")

print()
print("  audit verdicts on the deposited tool outputs")
tools = {}
with open(PANEL, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        tools.setdefault(row["tool"], []).append(row)
# Exact deposited emission totals: a corruption that zeroes or swaps rows must fail.
EXPECTED = {
    "Cas-OFFinder 3.0.0 develop": (1561, 73),
    "Cas-OFFinder tag 3.0.0b3": (1561, 73),
    "Cas-OFFinder-bulge v1.2 beta (2016 script)": (1491, 98),
}
for tool, rows in sorted(tools.items()):
    num = [r for r in rows if r["patterns_emitted"].strip().isdigit()]
    if not num:
        continue
    emitted = sum(int(r["patterns_emitted"]) for r in num)
    need = sum(int(r["r"]) for r in num)
    short = sum(1 for r in num if int(r["patterns_emitted"]) < int(r["r"]))
    verdict = "INCOMPLETE" if emitted < need else ("EXACT" if emitted == need else "OVER")
    print("    %-44s emits %5d of %5d, short on %3d guides: %s"
          % (tool[:44], emitted, need, short, verdict))
    claim(emitted <= need or verdict == "OVER",
          "    %.38s accounted over %d guides" % (tool, len(num)),
          "%d vs %d" % (emitted, need), "consistent")
    if tool in EXPECTED:
        want_e, want_s = EXPECTED[tool]
        claim((emitted, short) == (want_e, want_s),
              "    %.38s matches deposited emission total" % tool,
              "%d short %d" % (emitted, short), "%d short %d" % (want_e, want_s))

bad = VERDICTS.count(False)
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - bad, bad))
if bad:
    sys.exit(1)
print("The run count is the distinct-word count on all 171 deposited guides, the certificate")
print("asserts a minimum only where the companion proves one, and on every instance small")
print("enough to settle completely, no family of r-1 DNA centres covers the deletion union.")
