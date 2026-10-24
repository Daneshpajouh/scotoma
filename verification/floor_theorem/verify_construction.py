# Verify the cyclic-shift construction in the canonical (threshold) coordinates: witnesses
# obtained from the thresholds by writing letters 2 and 3 on a support shifted cyclically
# with the centre index, support S = {4, 9, 16} for letter 2 and {6, 12} for letter 3.
# P_k = 1^k 0^(N-k), 1-indexed coordinates 1..N, alphabet {0,1,2,3}.
N = 19
A = {4, 9, 16}
B = {6, 12}
S = sorted(A | B)

def tau(k, s):
    return 1 + ((k + s - 1) % 19)

def P(k):
    return [1 if i <= k else 0 for i in range(1, N + 1)]

def W(k):
    w = P(k)
    for s in A: w[tau(k, s) - 1] = 2
    for s in B: w[tau(k, s) - 1] = 3
    return w

def d(x, y):
    return sum(a != b for a, b in zip(x, y))

R = 19  # centres k = 0..R-1
print("N =", N, " centres =", R, " support S =", S)

# 1. distance from each witness to its own centre
dd = {d(W(k), P(k)) for k in range(R)}
print("d(W_k, P_k) values:", sorted(dd), " (claimed 5)")

# 2. pairwise witness distances, by index difference
by_h = {}
for i in range(R):
    for j in range(i + 1, R):
        by_h.setdefault(j - i, []).append(d(W(i), W(j)))
claimed = [11,11,11,11,11,11,11,13,14,14,15,12,15,15,17,16,17,18]
print("\n h   min  max   claimed  match")
ok = True
for h in range(1, 19):
    v = by_h.get(h, [])
    if not v: continue
    c = claimed[h - 1]
    m_ = min(v); M_ = max(v)
    good = (m_ == c)
    ok &= good
    print(" %2d  %3d  %3d   %4d     %s" % (h, m_, M_, c, "yes" if good else "NO"))
print("\nall claimed minima reproduced:", ok)
print("overall minimum pairwise witness distance:", min(min(v) for v in by_h.values()))

# 3. privacy: is W_k farther than m from every other centre?
# This family is the BUDGET-FIVE one, so only m = 5 is a pass condition. The m = 4 line is printed
# for contrast and its "False" is expected: a budget-five witness sits at distance 5 from its own
# centre, which is more than 4. verify_m4.py carries the budget-four family. Reported separately so
# that a reader, and the exit code below, do not read an expected contrast as a failure.
own = max(d(W(k), P(k)) for k in range(R))
mn = min(min(v) for v in by_h.values())
for m in (4, 5):
    bad = [(k, j, d(W(k), P(j))) for k in range(R) for j in range(R)
           if j != k and d(W(k), P(j)) <= m]
    label = "pass condition" if m == 5 else "contrast only, see verify_m4.py"
    print("\nm = %d (%s):  d(W_k,P_k) max = %d (needs <= m: %s)" % (m, label, own, own <= m))
    print("   violations of d(W_k,P_j) > m for j != k:", len(bad))
    print("   min pairwise witness distance %d > 2m = %d : %s" % (mn, 2 * m, mn > 2 * m))
    if m == 5:
        privacy_ok = not bad
        own_ok = own <= m
        scattered_ok = mn > 2 * m

# The verdict, and an exit code that carries it. Until 2026-10-23 every script in this folder
# printed its booleans and exited 0 whatever they said: an adversarial run moved letter 3 from
# {6, 12} to {6, 7}, which makes two printed checks False, and this script still exited 0. Anyone
# running it in CI or a loop saw a pass. A verifier that cannot fail verifies nothing.
good = ok and own_ok and privacy_ok and scattered_ok
print("\nVERDICT:", "PASS" if good else "FAIL")
if not good:
    print("  minima reproduced:", ok, " own distance <= 5:", own_ok,
          " privacy:", privacy_ok, " 2m-scattered:", scattered_ok)
raise SystemExit(0 if good else 1)
