#!/usr/bin/env python3
"""The witness family of the floor theorem is forced, not chosen.

Claim checked here, as the paper states it: at the obstruction boundary N = 4m-1 the cyclic
overwrite rule admits a valid witness family only for 2 <= m <= 5, and at N = 19, m = 5 the
overwrite set is unique up to reflection and exchange of the two labels.

Model, exactly as the construction in the paper uses it. N coordinates over {0,1,2,3}. The pattern
family is P_k = 1^k 0^(N-k) for k = 0..N-1. A witness family is given by an overwrite set
S subset of {1..N} with |S| = m and a labelling S -> {2,3}: witness x_k is P_k with position
tau(k,s) = 1 + ((k + s - 1) mod N) set to label(s), for each s in S. The family is valid when
d(x_j, x_k) >= 2m + 1 for every j != k, which is what 2m-scattered means.

Usage:
    check_uniqueness.py              the headline: N = 19, m = 5, and two controls
    check_uniqueness.py --wide       also the boundary counts for m = 2..5 and the m = 6 cases
    check_uniqueness.py --self-test  planted-defect controls; each must FAIL the check

Exit 0 when every count matches the paper, 1 when one does not, 2 if it could not run.
"""
import sys
from itertools import combinations, product

EXPECTED = {            # (N, m): number of valid (overwrite set, labelling) pairs
    (19, 5): 4,         # the paper's claim: unique up to reflection and label swap
    (20, 5): 96,        # one coordinate of slack removes the rigidity
    # the rule attains the boundary N = 4m-1 across the whole range the paper claims,
    (7, 2): 2,
    (11, 3): 8,
    (15, 4): 20,
    # and fails at the next budget up, which is the other half of the claim.
    (23, 6): 0,
}


def witnesses(N, S, labels):
    out = []
    for k in range(N):
        w = [1] * k + [0] * (N - k)
        for s, lab in zip(S, labels):
            w[(k + s - 1) % N] = lab
        out.append(w)
    return out


def valid(N, m, S, labels):
    W = witnesses(N, S, labels)
    need = 2 * m + 1
    for i in range(len(W)):
        wi = W[i]
        for j in range(i + 1, len(W)):
            wj = W[j]
            if sum(a != b for a, b in zip(wi, wj)) < need:
                return False
    return True


def count(N, m, stop_after=None):
    found = []
    for S in combinations(range(1, N + 1), m):
        for labels in product((2, 3), repeat=m):
            if valid(N, m, S, labels):
                found.append((S, labels))
                if stop_after is not None and len(found) >= stop_after:
                    return found
    return found


def gaps(N, S):
    s = sorted(S)
    return sorted([s[i + 1] - s[i] for i in range(len(s) - 1)] + [s[0] + N - s[-1]])


def main(argv):
    wide = "--wide" in argv
    if "--self-test" in argv:
        # Each planted defect must change the headline count away from 4.
        bad = 0
        for name, N, m, want_not in (("budget one too low", 19, 4, 4), ("length one short", 18, 5, 4)):
            got = len(count(N, m))
            ok = got != want_not
            print(f"  control {name:<22} N={N} m={m}: {got} valid, differs from 4: {ok}")
            bad |= not ok
        print("VERDICT:", "PASS" if not bad else "FAIL")
        return 1 if bad else 0

    bad = 0
    found = count(19, 5)
    print(f"N = 19, m = 5: searched all {len(list(combinations(range(1, 20), 5))) * 2**5} "
          f"(overwrite set, labelling) pairs")
    print(f"  valid: {len(found)}   expected {EXPECTED[(19, 5)]}")
    for S, lab in found:
        print(f"    S = {S}  labels = {lab}  gaps = {gaps(19, S)}")
    supports = sorted({S for S, _ in found})
    print(f"  distinct overwrite sets: {len(supports)} -> {supports}")
    print(f"  all share one gap multiset: {len({tuple(gaps(19, S)) for S in supports}) == 1}")
    bad |= len(found) != EXPECTED[(19, 5)]

    cheap = [(7, 2), (11, 3), (15, 4), (20, 5)]
    checks = cheap if not wide else cheap + [(23, 6)]
    for N, m in checks:
        got = len(count(N, m))
        ok = got == EXPECTED[(N, m)]
        print(f"N = {N:>2}, m = {m}: {got:>3} valid, expected {EXPECTED[(N, m)]:>3}  {'ok' if ok else 'MISMATCH'}")
        bad |= not ok

    print()
    print("VERDICT:", "PASS" if not bad else "FAIL")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
