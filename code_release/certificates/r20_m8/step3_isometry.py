"""
Step 3 check, independent of any reviewer's arithmetic.

For a length-n guide w with r runs, its distinct single-deletion descendants (length n-1) are one
per run. Claim: an isometry of Hamming space (a coordinate permutation followed by a per-coordinate
alphabet permutation) maps them exactly onto the staircase prefix {p_0, ..., p_(r-1)},
p_i = 1^i 0^(n-1-i), with p_i the image of the descendant from deleting inside run i.

Construction tested: coordinate k of a descendant is "switching" iff w[k] != w[k+1]; there are r-1
of them. Put switching coordinates first, in order, then the constant ones. On a switching
coordinate send w[k] -> 1 and w[k+1] -> 0; on a constant coordinate send w[k] (= w[k+1]) -> 0.
Extend each to a permutation of {A,C,G,T}, so the map is an isometry of all of Z_4^(n-1).

Checks per guide: (1) exactly r distinct descendants; (2) image equals the staircase prefix exactly,
in run order; (3) every pairwise Hamming distance preserved (d(desc_a, desc_b) = |a - b|).
CONTROL: the same map without the coordinate permutation must FAIL on guides whose switching
coordinates are not already the first r-1.
"""
import itertools, random, sys

def runs(w):
    return 1 + sum(w[k] != w[k + 1] for k in range(len(w) - 1))

def descendants_by_run(w):
    out, seen = [], set()
    for j in range(len(w)):
        d = w[:j] + w[j + 1:]
        if d not in seen:
            seen.add(d); out.append(d)
    return out  # in deletion order, which is run order

def iso_map(w, permute=True):
    n1 = len(w) - 1
    sw = [k for k in range(n1) if w[k] != w[k + 1]]
    co = [k for k in range(n1) if w[k] == w[k + 1]]
    order = sw + co if permute else list(range(n1))
    def f(d):
        out = []
        for k in order:
            if w[k] != w[k + 1]:
                out.append("1" if d[k] == w[k] else ("0" if d[k] == w[k + 1] else "X"))
            else:
                out.append("0" if d[k] == w[k] else "X")
        return "".join(out)
    return f

def ham(a, b):
    return sum(x != y for x, y in zip(a, b))

def check(w, permute=True):
    r = runs(w); n1 = len(w) - 1
    ds = descendants_by_run(w)
    if len(ds) != r:
        return False, "descendant count %d != runs %d" % (len(ds), r)
    f = iso_map(w, permute)
    img = [f(d) for d in ds]
    stair = ["1" * i + "0" * (n1 - i) for i in range(r)]
    if img != stair:
        return False, "image is not the staircase prefix"
    for a in range(r):
        for b in range(r):
            if ham(ds[a], ds[b]) != abs(a - b):
                return False, "distance d(%d,%d)=%d not %d" % (a, b, ham(ds[a], ds[b]), abs(a - b))
    return True, "ok"

guides = [l.split()[-1].strip().upper() for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
guides = [g for g in guides if set(g) <= set("ACGT")]
random.seed(5)
rand = ["".join(random.choice("ACGT") for _ in range(20)) for _ in range(20000)]
adversarial = ["A" * 20, "AC" * 10, "ACGT" * 5, "A" * 19 + "C", "C" + "A" * 19, "AAAACCCCGGGGTTTTAAAA",
               "ACACACACACGGGGGGGGGG", "ATATATATATATATATATAT", "AACCAACCAACCAACCAACC"]
for name, S in (("panel", guides), ("random 20-mers", rand), ("adversarial", adversarial)):
    bad = [(g, why) for g in S for ok, why in [check(g)] if not ok]
    rs = sorted({runs(g) for g in S})
    print(f"{name:15s}: {len(S):>6} guides, run counts {rs[0]}..{rs[-1]}, failures {len(bad)}")
    for g, why in bad[:3]:
        print("   FAIL", g, why)
# CONTROL: without the coordinate permutation it must fail on most guides
ctrl_fail = sum(1 for g in guides if not check(g, permute=False)[0])
print(f"CONTROL (no coordinate permutation) on the panel: {ctrl_fail} of {len(guides)} fail, as they must")
assert ctrl_fail > 0
# A claim supplied from an independent route, recomputed here rather than taken on report
w = "01" * 10
d18, d19 = w[:18] + w[19:], w[:19]
print("the supplied case (01)^10: delete pos 19 ->", w[:18] + w[19:], " delete pos 20 ->", w[:19], " distance", ham(d18, d19))
