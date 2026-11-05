#!/usr/bin/env python3
"""A printed private witness at budget six on the clinical guide sg1617.

Why this exists. Paper B proves K = r for every guide of length at least twenty and every
budget at most nine, with the competing centres unrestricted. An external reviewer asked for
the claim to be made tangible in the body: one concrete set of r-1 centres, at a budget in
{6,7,8,9}, with a PRINTED word of U_m that those r-1 radius-m balls do not cover, on a guide
whose off-target list a reader already cares about. This file is that artifact and its checks.

The guide is sg1617, the BCL11A +58 erythroid enhancer guide used in the exa-cel trials, and the
same guide Cancellieri et al. searched at up to six mismatches with two DNA/RNA bulges.

Everything below is recomputed from the guide string alone. Nothing is read from a stored list.
"""
from itertools import combinations
import sys

GUIDE = "CTAACAGTTGCTTTTATCAC"
BUDGET = 6
WITNESS = "ACCGGCGTTGCTTTTATCA"
SIG = "ACGT"
claims = [0, 0]

def claim(ok, label, got, want):
    claims[0] += 1
    if ok: claims[1] += 1
    print("  %-5s %-62s got %-26s want %s" % ("PASS" if ok else "FAIL", label, got, want))
    if not ok: claim.failed = True
claim.failed = False

def hd(a, b):
    assert len(a) == len(b)
    return sum(x != y for x, y in zip(a, b))

n = len(GUIDE); N = n - 1
# the library, rebuilt from the guide: distinct words from deleting one base. Deletions inside a
# run of equal bases coincide, so the count is the run count.
library = sorted({GUIDE[:i] + GUIDE[i + 1:] for i in range(n)})
runs = []
i = 0
while i < n:
    j = i
    while j + 1 < n and GUIDE[j + 1] == GUIDE[i]: j += 1
    runs.append(GUIDE[i:j + 1]); i = j + 1

print("guide %s, length %d, budget %d" % (GUIDE, n, BUDGET))
claim(len(GUIDE) == 20, "the guide has length twenty, so the floor theorem applies", len(GUIDE), 20)
claim(6 <= BUDGET <= 9, "the budget lies in the range only this paper's floor certifies", BUDGET, "6..9")
claim(len(library) == len(runs), "distinct deletion words equal the run count", len(library), len(runs))
r = len(library)
print("  r = %d runs: %s" % (r, " ".join(runs)))

# the omitted centre is the terminal deletion, the one the released build never emits
terminal = GUIDE[:N]
claim(terminal in library, "the terminal deletion word is in the library", terminal, "present")
others = [w for w in library if w != terminal]
claim(len(others) == r - 1, "the competing set is the library minus that one word", len(others), r - 1)

# THE CLAIM
d_own = hd(WITNESS, terminal)
d_other = min(hd(WITNESS, w) for w in others)
claim(len(WITNESS) == N, "the witness has the deletion-word length", len(WITNESS), N)
claim(set(WITNESS) <= set(SIG), "the witness is over the DNA alphabet", "".join(sorted(set(WITNESS))), "subset of ACGT")
claim(d_own <= BUDGET, "the witness lies in U_6, within the budget of a library word", d_own, "<= %d" % BUDGET)
claim(d_other > BUDGET, "no other library word covers it at that budget", d_other, "> %d" % BUDGET)
claim(d_other == BUDGET + 1, "and it attains the triangle bound exactly", d_other, BUDGET + 1)

# C1. control: the triangle bound. The nearest other word is at distance one, so BUDGET+1 is the
# most any witness within BUDGET of the terminal word can achieve. The witness is extremal.
nearest = min(hd(terminal, w) for w in others)
claim(nearest == 1, "the nearest competing word is adjacent to the terminal word", nearest, 1)
claim(d_other <= BUDGET + nearest, "the achieved separation respects the triangle bound", d_other, "<= %d" % (BUDGET + nearest))

# C2. control: a forged witness must lose the property. Change one coordinate at a time; the
# property must fail for at least one such change, so the check is not satisfied by any string.
broke = 0
for c in range(N):
    for L in SIG:
        if L == WITNESS[c]: continue
        f = WITNESS[:c] + L + WITNESS[c + 1:]
        if not (hd(f, terminal) <= BUDGET and min(hd(f, w) for w in others) > BUDGET):
            broke += 1
claim(broke > 0, "control: single-coordinate forgeries of the witness lose the property", broke, "> 0")

# C3. control: the property must fail if the omitted word is put back
with_term = others + [terminal]
claim(min(hd(WITNESS, w) for w in with_term) <= BUDGET,
      "control: restoring the omitted word covers the witness", min(hd(WITNESS, w) for w in with_term), "<= %d" % BUDGET)

# C4. control: the search was exhaustive over this construction, so the printed witness is not
# merely one lucky find. Recompute the best achievable separation over all six-coordinate flips.
best = -1
for S in combinations(range(N), BUDGET):
    w = list(terminal)
    for c in S:
        w[c] = max((L for L in SIG if L != terminal[c]),
                   key=lambda L: sum(1 for x in others if x[c] != L))
    w = "".join(w)
    if hd(w, terminal) <= BUDGET:
        best = max(best, min(hd(w, x) for x in others))
claim(best == d_other, "exhaustive over all six-flip constructions, the printed witness is optimal", best, d_other)

print()
print("%d claims checked, %d reproduced, %d failed" % (claims[0], claims[1], claims[0] - claims[1]))
if claim.failed:
    sys.exit(1)
print("At budget six on sg1617 the %d library words other than the terminal deletion leave" % (r - 1))
print("%s uncovered, at distance %d from the omitted word and %d from every other." % (WITNESS, d_own, d_other))
