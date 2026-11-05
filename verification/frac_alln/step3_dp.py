#!/usr/bin/env python3
"""Step 3, independently: does ANY centre collect more than one unit of the private-sphere weight?

This is the step we re-ran independently, with separately written code, because it
is an exact computation rather than a reading proof. Steps 1 and 2 I have already reproduced in
privsphere.py: privacy is a one-line argument, the total is 16,633,344, every pattern carries
exactly 1, and the reported best non-pattern centre gives exactly 29/34.

WHY A MAXIMISATION IS POSSIBLE AT ALL. The ball count at a centre c depends on c only through, for
each pattern j, the four summaries (a, b_diff, b_same, e) derived in privsphere.py. Those depend on
c only through whether each coordinate is 0, 1, or one of the two unused letters -- which of 2 or 3
never matters, since the count C(t_a,k) is the same either way. So the search space is words over
{0,1,X} of length 19, which is 3^19 = 1,162,261,467: too many to enumerate, few enough to DP.

THE DECOMPOSITION. Write Z_{<=j} and O_{<=j} for the number of 0s and 1s in the first j
coordinates, and Z, O for the totals. Then, before removing the two excluded coordinates,
  b_diff_j = Z_{<=j} + (O - O_{<=j}),   b_same_j = O_{<=j} + (Z - Z_{<=j}),
  a_j = X - (X among coordinates j, j+1),  e_j = [c_j != 1] + [c_{j+1} != 0].
So pattern j's contribution is a function of (Z_{<=j}, O_{<=j}, c_j, c_{j+1}) once the TOTALS are
fixed. Fix the totals in an outer loop and the objective becomes a sum of per-position local terms,
which a max-plus DP over (position, Z_{<=p}, O_{<=p}, c_p) maximises exactly.

EXACT ARITHMETIC. Spheres have size C(17,6)*2^6 = 792,064 or C(18,6)*2^6 = 1,188,096, whose lcm is
2,376,192. Scaling by 3 and 2 respectively makes every term an integer and feasibility becomes
"total <= 2,376,192". No floating point anywhere.

THE CONTROL THAT MAKES THE ANSWER WORTH ANYTHING. The DP is checked against BRUTE FORCE over the
whole space at sizes where brute force is possible. If the DP's maximum differs from the true
maximum at (N,m) = (7,2), (8,2) or (8,3), its answer at (19,6) is discarded. A DP whose answer
cannot be validated is not evidence.
"""
import itertools
import os
import sys
from fractions import Fraction
from math import comb, gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # privsphere.py beside this file, any working directory
from privsphere import ball_weight, pattern, sphere_size, allowed  # noqa: E402


def term(a, b_diff, b_same, e, m):
    """#{x in S_j : d(c,x) <= m} from the four summaries."""
    tot = 0
    for t_a in range(0, min(a, m) + 1):
        for t_d in range(0, min(b_diff, m - t_a) + 1):
            t_s = m - t_a - t_d
            if t_s < 0 or t_s > b_same:
                continue
            base = a + b_diff + e - t_a + t_s
            for k in range(0, t_a + 1):
                if base + k <= m:
                    tot += (comb(a, t_a) * comb(b_diff, t_d) * comb(b_same, t_s)
                            * comb(t_a, k) * 2 ** (t_d + t_s))
    return tot


def scaled_max(N, m, r, verbose=False):
    """Exact maximum of the scaled ball weight over all centres, by DP. Returns (max, lcm)."""
    sizes = [sphere_size(j, N, m) for j in range(r)]
    L = sizes[0]
    for s in sizes[1:]:
        L = L * s // gcd(L, s)
    scale = [L // s for s in sizes]

    # j-term as a function of (Z_le_j, O_le_j, cj, cj1) given totals. Coordinates are 1-indexed;
    # c is indexed 0..N-1. Letters: 0, 1, 2 meaning "one of the unused letters".
    def jterm(j, Zle, Ole, cj, cj1, Z, O, X):
        A = set(allowed(j, N))
        excl = [p for p in (j, j + 1) if 1 <= p <= N]
        # letters at the excluded coordinates
        lets = {}
        if j in excl:
            lets[j] = cj
        if (j + 1) in excl:
            lets[j + 1] = cj1
        aX = sum(1 for p in excl if lets[p] == 2)
        a = X - aX
        bd = Zle + (O - Ole)
        bs = Ole + (Z - Zle)
        e = 0
        for p in excl:
            Pp = 1 if p <= j else 0
            lp = lets[p]
            lv = 2 if lp == 2 else lp          # unused letter is never 0 or 1
            if lv != Pp:
                e += 1
            # remove this coordinate from b_diff / b_same since it is not allowed
            if lp == 0:
                if Pp == 0:
                    bs -= 1
                else:
                    bd -= 1
            elif lp == 1:
                if Pp == 1:
                    bs -= 1
                else:
                    bd -= 1
        if a < 0 or bd < 0 or bs < 0:
            return None
        if a + bd + bs != len(A):
            return None
        return term(a, bd, bs, e, m) * scale[j]

    best = None
    best_word = None
    for Z in range(N + 1):
        for O in range(N + 1 - Z):
            X = N - Z - O
            # DP over positions 1..N. State: (Zle, Ole, letter at this position).
            # Pattern j's term needs c_j and c_{j+1}, and prefix counts at j, so it is added once
            # we have placed coordinate j+1. j = 0 needs only c_1; j = N needs only c_N.
            # states[(Zle_at_p, Ole_at_p, c_p)] = best accumulated
            states = {}
            for c1 in (0, 1, 2):
                z = 1 if c1 == 0 else 0
                o = 1 if c1 == 1 else 0
                if z > Z or o > O or (c1 == 2 and X == 0):
                    continue
                # j = 0: prefix counts are 0, c_j does not exist, c_{j+1} = c_1
                t = jterm(0, 0, 0, None, c1, Z, O, X)
                if t is None:
                    continue
                key = (z, o, c1)
                if key not in states or states[key][0] < t:
                    states[key] = (t, (c1,))
            for p in range(2, N + 1):
                nxt = {}
                for (z, o, cp), (acc, word) in states.items():
                    for cn in (0, 1, 2):
                        z2 = z + (1 if cn == 0 else 0)
                        o2 = o + (1 if cn == 1 else 0)
                        x2 = p - z2 - o2
                        if z2 > Z or o2 > O or x2 > X or x2 < 0:
                            continue
                        # pattern j = p-1 is completed now: prefix counts are at j = p-1,
                        # which is (z, o) BEFORE adding coordinate p... careful: (z,o) are counts
                        # over coordinates 1..p-1, which is exactly prefix j = p-1.
                        t = jterm(p - 1, z, o, cp, cn, Z, O, X)
                        if t is None:
                            continue
                        key = (z2, o2, cn)
                        val = acc + t
                        if key not in nxt or nxt[key][0] < val:
                            nxt[key] = (val, word + (cn,))
                states = nxt
                if not states:
                    break
            # j = N: prefix counts over all N coordinates, c_N is the last letter, no c_{N+1}
            for (z, o, cp), (acc, word) in states.items():
                if z != Z or o != O:
                    continue
                t = jterm(N, z, o, cp, None, Z, O, X)
                if t is None:
                    continue
                val = acc + t
                if best is None or val > best:
                    best, best_word = val, word
    return best, L, best_word


def brute_max(N, m, r):
    best, bw = None, None
    for cc in itertools.product(range(4), repeat=N):
        w = ball_weight(list(cc), N, m, r)
        if best is None or w > best:
            best, bw = w, cc
    return best, bw


def main():
    # CONTROL: the DP must equal brute force wherever brute force is possible.
    for (N, m) in ((7, 2), (8, 2), (8, 3)):
        r = N + 1
        dp, L, word = scaled_max(N, m, r)
        bf, bw = brute_max(N, m, r)
        dpf = Fraction(dp, L)
        status = 'ok' if dpf == bf else 'MISMATCH'
        print(f'control (N,m)=({N},{m}) r={r}: DP {dpf} vs brute force {bf}  -> {status}')
        if dpf != bf:
            print(f'  DP word {word}, brute-force centre {bw}')
            print('CONTROL FAILED. The DP is wrong and its answer at (19,6) is discarded.')
            return 1
    print('controls: the DP reproduces the exact maximum at every size where brute force is '
          'possible.\n')

    N, m, r = 19, 6, 20
    dp, L, word = scaled_max(N, m, r)
    w = Fraction(dp, L)
    print(f'N={N}, m={m}, r={r}')
    print(f'  lcm of sphere sizes          : {L:,}')
    print(f'  exact maximum ball weight    : {w} = {float(w):.9f}')
    print(f'  attained at (0/1/X word)     : {"".join("01X"[c] for c in word)}')
    print(f'  feasible (no ball above 1)   : {w <= 1}')
    if w <= 1:
        print(f'\nSo every covering code of the union has at least r = {r} centres, and the r '
              f'patterns attain it: the covering number is exactly {r} at m = {m}.')
    else:
        print(f'\nINFEASIBLE: some ball carries {w} > 1, so this weighting does NOT prove the '
              f'floor at m = {m}.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
