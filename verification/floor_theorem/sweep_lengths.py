# The witness family of the scattered-transversal theorem at EVERY budget m = 1..5 and every length
# n = 20..61 (N = n - 1 = 19..60), with the largest library a guide of that length can have, N + 1
# centres P_k = 1^k 0^(N-k). verify_construction.py checks m = 5 at n = 20 only, and verify_m4.py
# searches for a support at each budget rather than running the paper's own letters, so before this
# file the letters at m = 1 and m = 3 ran nowhere. Here they are read from the paper's table and used
# as printed: offsets 4, 6, 9, 12, 16 write letters 2, 3, 2, 3, 2, and budget m uses the first m.
# The cycle tau(k, s) = 1 + ((k + s - 1) mod 19) is fixed at 19 for every length, as in the paper.
OFFSETS = [4, 6, 9, 12, 16]
LETTERS = [2, 3, 2, 3, 2]
LENGTHS = range(19, 61)  # N

def tau(k, s):
    return 1 + ((k + s - 1) % 19)

def P(k, N):
    return [1 if i <= k else 0 for i in range(1, N + 1)]

def W(k, N, support):
    w = P(k, N)
    for s, a in support:
        w[tau(k, s) - 1] = a
    return w

def d(x, y):
    return sum(a != b for a, b in zip(x, y))

def measure(N, support):
    R = N + 1
    Ps = [P(k, N) for k in range(R)]
    Ws = [W(k, N, support) for k in range(R)]
    own = sorted({d(Ws[k], Ps[k]) for k in range(R)})
    pair, priv = {}, {}
    for i in range(R):
        for j in range(R):
            if j == i:
                continue
            t = abs(j - i)
            priv[t] = min(priv.get(t, N + 1), d(Ws[i], Ps[j]))
            if j > i:
                pair[t] = min(pair.get(t, N + 1), d(Ws[i], Ws[j]))
    return Ws, own, pair, priv

def windows_exact(Ws, m):
    # Every four consecutive witnesses pairwise at distance exactly 2m + 1, checked at m = 5 only.
    # Adjacent witnesses are capped at 2m + 1 by the triangle inequality, so those three pairs are
    # forced; the other three are a measured property of this family. A counting argument sharper
    # than the paper's N >= 4m - 4 says the six pairs of a window sum to at most 66 at N = 19, which
    # would force all six to 11 there, but that argument is not in the paper and nothing here relies
    # on it. This file measures the property; it does not prove it.
    bad = 0
    for k in range(len(Ws) - 3):
        if any(d(Ws[k + a], Ws[k + b]) != 2 * m + 1 for a in range(4) for b in range(a + 1, 4)):
            bad += 1
    return len(Ws) - 3, bad

def check(m, support, label, quiet=False):
    # The cyclic table: per-shift minima at N = 19, where the chain fills the cycle. The paper says
    # this is a lower bound at every length; shifts of 19 or more are bounded by t - m (pairwise)
    # and t (privacy) instead.
    _, _, cyc_pair, cyc_priv = measure(19, support)
    fails = []
    worst_pair = worst_priv = None
    windows = [0, 0]
    for N in LENGTHS:
        Ws, own, pair, priv = measure(N, support)
        if own != [m]:
            fails.append("N=%d own distance %s, not %d" % (N, own, m))
        for t in pair:
            lo_pair = cyc_pair[t] if t < 19 else t - m
            lo_priv = cyc_priv[t] if t < 19 else t
            if pair[t] < lo_pair:
                fails.append("N=%d t=%d pairwise %d below %d" % (N, t, pair[t], lo_pair))
            if priv[t] < lo_priv:
                fails.append("N=%d t=%d privacy %d below %d" % (N, t, priv[t], lo_priv))
        mp, mv = min(pair.values()), min(priv.values())
        if mp <= 2 * m:
            fails.append("N=%d minimum pairwise %d, needs > %d" % (N, mp, 2 * m))
        if mv <= m:
            fails.append("N=%d minimum privacy %d, needs > %d" % (N, mv, m))
        worst_pair = mp if worst_pair is None else min(worst_pair, mp)
        worst_priv = mv if worst_priv is None else min(worst_priv, mv)
        if m == 5:
            nwin, bad = windows_exact(Ws, m)
            windows[0] += nwin
            windows[1] += bad
            if bad:
                fails.append("N=%d %d of %d four-witness windows not all at 11" % (N, bad, nwin))
    if not quiet:
        print("%-28s pairwise min %2d (needs > %2d)  privacy min %2d (needs > %d)  %s"
              % (label, worst_pair, 2 * m, worst_priv, m, "PASS" if not fails else "FAIL"))
        if m == 5:
            print("   four-witness windows with all six pairs at exactly 11: %d of %d"
                  % (windows[0] - windows[1], windows[0]))
    return cyc_pair, cyc_priv, fails

print("Lengths n = %d..%d, N + 1 centres at each; offsets %s, letters %s"
      % (LENGTHS[0] + 1, LENGTHS[-1] + 1, OFFSETS, LETTERS))
ok = True
tables = {}
for m in range(1, 6):
    support = list(zip(OFFSETS[:m], LETTERS[:m]))
    cp, cv, fails = check(m, support, "m = %d  %s" % (m, " ".join("%d->%d" % x for x in support)))
    tables[m] = (cp, cv)
    for f in fails[:5]:
        print("   ", f)
    ok &= not fails

# The paper prints the nineteen cyclic cases at m = 4 and m = 5 (appendix, "The nineteen cases in
# full"). Compared value by value, both rows at both budgets.
PRINTED = {
    4: ([9, 9, 9, 11, 10, 10, 13, 12, 13, 14, 14, 14, 14, 14, 16, 16, 17, 18, 15],
        [5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 12, 12, 13, 14, 15, 16, 17, 18, 19]),
    5: ([11, 11, 11, 11, 11, 11, 11, 13, 14, 14, 15, 12, 15, 15, 17, 16, 17, 18, 14],
        [6, 7, 8, 8, 9, 9, 10, 11, 11, 12, 13, 13, 14, 15, 16, 16, 17, 18, 19]),
}
for m, (pp, pv) in PRINTED.items():
    cp, cv = tables[m]
    got_p = [cp[t] for t in range(1, 20)]
    got_v = [cv[t] for t in range(1, 20)]
    same = got_p == pp and got_v == pv
    print("printed nineteen-case table at m = %d reproduced: %s" % (m, same))
    if not same:
        print("    pairwise measured", got_p)
        print("    privacy  measured", got_v)
    ok &= same

# Controls, run every time: each must make the check fail, or this file could not fail either.
# 1. Four positions at budget five (offset 16 dropped): own distance 4, minimum pair 9.
# 2. Letter 3 moved from offset 12 to offset 7, an adversarial move used on verify_construction.py.
# 3. The window check alone, on the paper's m = 5 family with one witness overwritten at one more
#    coordinate, so that some window holds a pair at 12.
_, _, f1 = check(5, list(zip([4, 6, 9, 12], [2, 3, 2, 3])), "", quiet=True)
_, _, f2 = check(5, [(4, 2), (6, 3), (9, 2), (7, 3), (16, 2)], "", quiet=True)
Ws5 = measure(19, list(zip(OFFSETS, LETTERS)))[0]
Ws5[2] = list(Ws5[2])
Ws5[2][[p for p in range(19) if Ws5[2][p] in (0, 1)][0]] = 3
_, bad3 = windows_exact(Ws5, 5)
controls = [("support of four at m = 5", bool(f1)), ("letter 3 at offset 7", bool(f2)),
            ("one extra overwrite, window check", bad3 > 0)]
for name, failed in controls:
    print("control %-34s %s" % (name + ":", "fails as it must" if failed else "PASSED, so the check is blind"))
    ok &= failed

print("\nVERDICT:", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
