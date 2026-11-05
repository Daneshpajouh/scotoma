#!/usr/bin/env python3
"""The step-3 maximum AT A GIVEN RUN COUNT, which is what the rescaling table needs.

WHY THIS SCRIPT EXISTS. scaled_max(N, m, r) in step3_dp.py is only ever correct at r = N+1: it
builds `sizes` with r entries and then indexes scale[j] for j = 0..N, so any r < N+1 raises
IndexError rather than answering a smaller-run-count question. The rescaling table in the appendix
prints a largest-proved r per row beside a maximum M, and the M it prints is the r = N+1 maximum.
Those are different numbers above m = 7, where the maximum depends on the run count.

WHAT THE RIGHT OBJECT IS. A spacer with r runs contributes the patterns D_k for k in its own
threshold set, an arbitrary r-subset of {0,...,N}. The quantity that proves the floor for every
such spacer is

    M*_r = max over r-subsets S, max over centres c, of  sum_{j in S} w_j(c),
           where w_j(c) = |B_m(c) ∩ S_j| / |S_j|.

For a FIXED centre the best subset is the r largest terms, so

    M*_r = max_c ( sum of the r largest of w_0(c),...,w_N(c) ).

A top-r selection is not a sum of per-position terms, so the existing DP cannot maximise it
directly. The standard dual does:

    sum of top-r of {v_j}  =  min over lambda of [ r*lambda + sum_j max(0, v_j - lambda) ],

so for EVERY lambda,  M*_r  <=  r*lambda + max_c sum_j max(0, w_j(c) - lambda),
and the right-hand maximum IS a sum of per-position terms, so the same DP computes it exactly.
Each lambda therefore gives a rigorous upper bound on M*_r, and the smallest one found is reported.
A row claiming r is proved needs an upper bound strictly below r/(r-1).

The other direction is cheap: any centre exhibited by the DP gives a LOWER bound on M*_r, namely
the top-r sum of its own weights. A lower bound at or above r/(r-1) shows r is NOT proved, which is
what justifies calling the row's r the largest.

All arithmetic is exact integer arithmetic over the common denominator.
"""
import os
import sys
from fractions import Fraction
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # sibling files, any working directory
from privsphere import sphere_size, allowed, ball_count, ball_weight   # noqa: E402
from step3_dp import term                                   # noqa: E402


def _lcm_sizes(N, m):
    sizes = [sphere_size(j, N, m) for j in range(N + 1)]
    L = sizes[0]
    for s in sizes[1:]:
        L = L * s // gcd(L, s)
    return sizes, L, [L // s for s in sizes]


def dp_shifted(N, m, lam, want_word=False):
    """max_c sum_j max(0, w_j(c) - lam), exact, as a Fraction. lam a Fraction >= 0."""
    sizes, L, scale = _lcm_sizes(N, m)
    a, b = lam.numerator, lam.denominator
    aL = a * L                                   # lam in units of 1/(b*L) is a*L

    def f(j, Zle, Ole, cj, cj1, Z, O, X):
        """b*scale_j*term_j - a*L, clipped at 0; None when the state is impossible."""
        A = allowed(j, N)
        excl = [p for p in (j, j + 1) if 1 <= p <= N]
        lets = {}
        if j in excl:
            lets[j] = cj
        if (j + 1) in excl:
            lets[j + 1] = cj1
        aX = sum(1 for p in excl if lets[p] == 2)
        aa = X - aX
        bd = Zle + (O - Ole)
        bs = Ole + (Z - Zle)
        e = 0
        for p in excl:
            Pp = 1 if p <= j else 0
            lp = lets[p]
            lv = 2 if lp == 2 else lp
            if lv != Pp:
                e += 1
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
        if aa < 0 or bd < 0 or bs < 0:
            return None
        if aa + bd + bs != len(A):
            return None
        v = b * term(aa, bd, bs, e, m) * scale[j] - aL
        return v if v > 0 else 0

    best, best_word = None, None
    for Z in range(N + 1):
        for O in range(N + 1 - Z):
            X = N - Z - O
            states = {}
            for c1 in (0, 1, 2):
                z = 1 if c1 == 0 else 0
                o = 1 if c1 == 1 else 0
                if z > Z or o > O or (c1 == 2 and X == 0):
                    continue
                t = f(0, 0, 0, None, c1, Z, O, X)
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
                        t = f(p - 1, z, o, cp, cn, Z, O, X)
                        if t is None:
                            continue
                        key = (z2, o2, cn)
                        val = acc + t
                        if key not in nxt or nxt[key][0] < val:
                            nxt[key] = (val, word + (cn,))
                states = nxt
                if not states:
                    break
            for (z, o, cp), (acc, word) in states.items():
                if z != Z or o != O:
                    continue
                t = f(N, z, o, cp, None, Z, O, X)
                if t is None:
                    continue
                val = acc + t
                if best is None or val > best:
                    best, best_word = val, word
    return Fraction(best, b * L), best_word


def exhibit_weights(word, N, m):
    """Exact w_j for the centre the DP exhibited, j = 0..N. Letters 0,1,2 -> 0,1,2."""
    c = [w if w != 2 else 2 for w in word]
    ws = [Fraction(ball_count(c, j, N, m), sphere_size(j, N, m)) for j in range(N + 1)]
    # Two invariants that hold by construction, asserted because a violated one means the
    # instrument is wrong rather than the answer surprising. w_j is a fraction of a sphere, so it
    # cannot exceed 1, and the total over r patterns cannot exceed r.
    assert all(0 <= w <= 1 for w in ws), 'w_j outside [0,1]: %s' % [str(w) for w in ws if not 0 <= w <= 1]
    assert sum(ws) <= len(ws), 'sum of w_j exceeds the pattern count'
    return ws


def check_dp(N, m):
    """The DP at lambda = 0 must equal the independent per-pattern sum, exactly.

    Not a formality: it is the only check that does not share the max-plus decomposition. 2026-10-22
    I read an apparent violation of both invariants above at m=11 and went looking for a bug in
    term(). There was none. My own diagnostic printed str(x)[:14], which turned 6694579/1584128
    into 6694579/158412 and a maximum of 606003/1810432 into 606003/18104, so a 4.23 read as 42.26
    and a 0.33 as 33.5. TRUNCATING A FRACTION IN DIAGNOSTIC OUTPUT CAN MANUFACTURE AN IMPOSSIBLE
    VALUE. Print exact fractions, and when an invariant appears to break, re-measure before
    diagnosing.
    """
    v, word = dp_shifted(N, m, Fraction(0))
    indep = ball_weight([w for w in word], N, m, N + 1)
    return v == indep, v, indep


def topr_lower(word, N, m, r):
    ws = sorted(exhibit_weights(word, N, m), reverse=True)
    return sum(ws[:r])
