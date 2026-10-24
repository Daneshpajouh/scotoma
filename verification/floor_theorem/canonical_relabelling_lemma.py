"""Is the canonical relabelling a structural fact or a property of our panel?

Claim to test: for ANY guide, each coordinate of the deletion-pattern matrix carries
at most two distinct letters and switches between them at most once as the index runs
over the patterns. If true for every guide, the appendix states a lemma, not a data check.
"""
import itertools, random
A = "ACGT"

def runs(g):
    out = []; i = 0
    while i < len(g):
        j = i
        while j + 1 < len(g) and g[j + 1] == g[i]: j += 1
        out.append((i, j)); i = j + 1
    return out

def patterns(g):
    return [g[:s] + g[s + 1:] for s, _e in runs(g)]

def bad(g):
    P = patterns(g)
    if len(P) < 2: return None
    N = len(P[0])
    for c in range(N):
        col = [p[c] for p in P]
        if len(set(col)) > 2: return ("3+ letters", c, "".join(col))
        sw = sum(1 for t in range(1, len(col)) if col[t] != col[t - 1])
        if sw > 1: return ("2+ switches", c, "".join(col))
    return None

total_failures = 0
print("exhaustive over all guides, lengths 2..8")
tot = 0; fails = 0
for n in range(2, 9):
    for t in itertools.product(A, repeat=n):
        g = "".join(t); tot += 1
        b = bad(g)
        if b: fails += 1; print("  FAIL", g, b)
print("  checked %d guides, failures %d" % (tot, fails))
total_failures += fails

print("random guides at the deployed length 20")
random.seed(11)
tot = 0; fails = 0
for _ in range(20000):
    g = "".join(random.choice(A) for _ in range(20)); tot += 1
    b = bad(g)
    if b: fails += 1; print("  FAIL", g, b)
print("  checked %d guides, failures %d" % (tot, fails))
total_failures += fails

print("adversarial: homopolymer, alternating, single run")
for g in ("A"*20, "AC"*10, "A"*19+"C", "ACGT"*5, "AAC"+"G"*17):
    print("  %-22s runs=%2d  %s" % (g, len(runs(g)), bad(g) or "ok"))
    if bad(g): total_failures += 1

# A verdict and an exit code. Until 2026-10-23 every script in this folder printed its checks
# and exited 0 whatever they said, so a reader running one in CI saw a pass even when the check
# had failed. A verifier that cannot fail verifies nothing.
print("\nVERDICT:", "PASS" if total_failures == 0 else "FAIL",
      "(%d guide(s) failed the lemma)" % total_failures)
raise SystemExit(0 if total_failures == 0 else 1)
