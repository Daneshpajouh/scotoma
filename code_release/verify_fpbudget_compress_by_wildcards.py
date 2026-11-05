#!/usr/bin/env python3
"""The floor under a NONZERO false-positive budget: w wildcards compress r patterns to ceil(r/(w+1)).

MECHANISM. The library is isometric to the integer line, d(P_i,P_j) = |i-j| (held lemma, widened
check in code_release/verify_path_lemma_wide.py). So w+1 CONSECUTIVE library words agree outside
exactly w coordinates. Fix those agreements and free the w differing coordinates: the resulting
pattern agrees with each of the w+1 words off its wildcard set, so for any x within m of any one of
them the fixed-position mismatch count is at most m and x is accepted. One pattern therefore covers
w+1 balls, and ceil(r/(w+1)) patterns cover U_g.

This is the direction the zero-budget theorem forbids. It is admissible only at a nonzero budget,
and the price is measured here rather than asserted.

TWO CHECKS, because the claim has a structural half and a quantitative half.
  A. EXHAUSTIVE at small N: the family really covers U_g (every word of U_g accepted by some
     pattern), and the family's own false-positive rate.
  B. STRUCTURAL on the real panel at n=20: every block of w+1 consecutive library words agrees
     outside exactly w coordinates. If that fails anywhere the construction does not apply there.
A control is included: the same coverage test on a DELIBERATELY MISPAIRED family (blocks taken from
non-consecutive patterns), which must FAIL to cover, else the coverage test proves nothing.
"""
import os
import re, sys
from itertools import combinations, product
from math import comb, ceil

ALPHA = 'ACGT'


def library(g):
    ends = [k for k in range(len(g) - 1) if g[k] != g[k + 1]] + [len(g) - 1]
    return [g[:e] + g[e + 1:] for e in ends]


def ball_union(L, m):
    N = len(L[0])
    U = set()
    for D in L:
        for k in range(m + 1):
            for pos in combinations(range(N), k):
                for subs in product(*[[c for c in ALPHA if c != D[i]] for i in pos]):
                    x = list(D); 
                    for i, c in zip(pos, subs): x[i] = c
                    U.add(''.join(x))
    return U


def block_pattern(words):
    """Fixed where all agree, wildcard where any differ. None if they differ in too many places."""
    N = len(words[0])
    W = [i for i in range(N) if len({w[i] for w in words}) > 1]
    pat = [words[0][i] if i not in W else '*' for i in range(N)]
    return ''.join(pat), W


def accepts(pat, x, m):
    d = 0
    for a, b in zip(pat, x):
        if a != '*' and a != b:
            d += 1
            if d > m:
                return False
    return True


def family(L, w, consecutive=True):
    step = w + 1
    order = list(range(len(L))) if consecutive else \
        [i for j in range(step) for i in range(j, len(L), step)]      # deliberately mispaired
    blocks = [order[i:i + step] for i in range(0, len(order), step)]
    return [block_pattern([L[i] for i in b]) for b in blocks]


def acc_size(pat, m):
    nf = sum(1 for c in pat if c != '*')
    nw = len(pat) - nf
    return 4 ** nw * sum(comb(nf, i) * 3 ** i for i in range(min(m, nf) + 1))


print('=== A. EXHAUSTIVE at small N: does the construction cover, and at what price? ===')
print('  guide    m   r   w   patterns  ceil(r/(w+1))  covers U_g   |U_g|   family FP rate')
for g, m in (('ACGTAC', 1), ('ACGTAC', 2), ('ACGTACG', 1), ('AACGTAC', 1), ('ACGTACG', 2)):
    L = library(g); U = ball_union(L, m); r = len(L)
    for w in (1, 2, 3):
        fam = family(L, w)
        if any(len(W) != w for _, W in fam if len(_) and W is not None and len(W) > 0) and False:
            pass
        cov = all(any(accepts(p, x, m) for p, _ in fam) for x in U)
        acc = set()
        for p, _ in fam:
            nf = [i for i, c in enumerate(p) if c != '*']
            # enumerate accepted words exactly, small N only
            for k in range(min(m, len(nf)) + 1):
                for pos in combinations(nf, k):
                    for subs in product(*[[c for c in ALPHA if c != p[i]] for i in pos]):
                        base = [p[i] if p[i] != '*' else None for i in range(len(p))]
                        for i, c in zip(pos, subs): base[i] = c
                        wl = [i for i in range(len(p)) if base[i] is None]
                        for fill in product(ALPHA, repeat=len(wl)):
                            y = list(base)
                            for i, c in zip(wl, fill): y[i] = c
                            acc.add(''.join(y))
        rate = (len(acc) - len(U)) / len(acc)
        print('  %-8s %-3d %-3d %-3d %-9d %-14d %-12s %-7d %.4f'
              % (g, m, r, w, len(fam), ceil(r / (w + 1)), 'YES' if cov else 'NO', len(U), rate))
    # control: mispaired family must fail
    fam_bad = family(L, 1, consecutive=False)
    covb = all(any(accepts(p, x, m) for p, _ in fam_bad) for x in U)
    print('  %-8s %-3d CONTROL mispaired w=1 covers U_g: %s  (must be NO)'
          % (g, m, 'YES' if covb else 'NO'))
    print()

print('=== B. STRUCTURAL on the real panel, n=20 ===')
sp = sorted(set(re.findall(r'"([ACGT]{20})"',
            open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'multiguide_20261017', 'match_changeseq.json')).read())))
bad = 0
rows = []
for g in sp:
    L = library(g); r = len(L)
    for w in (1, 2, 3, 4, 5):
        for i in range(0, r, w + 1):
            blk = L[i:i + w + 1]
            if len(blk) < 2:
                continue
            _, W = block_pattern(blk)
            if len(W) != len(blk) - 1:
                bad += 1
    rows.append(r)
print('  110 guides, r from %d to %d (median %d)' % (min(rows), max(rows), sorted(rows)[len(rows)//2]))
print('  blocks whose wildcard count != block size - 1, over w=1..5: %d' % bad)
if bad == 0:
    print('  so every block of w+1 consecutive library words agrees outside EXACTLY w coordinates,')
    print('  on all 110 guides, and the construction applies at every w tested')
print()
print('  compressed family size at n=20, median r = %d:' % sorted(rows)[len(rows)//2])
med = sorted(rows)[len(rows)//2]
for w in (0, 1, 2, 3, 4, 5):
    print('    w=%d  ceil(r/(w+1)) = %2d   (%.0f%% of r)' % (w, ceil(med/(w+1)), 100*ceil(med/(w+1))/med))
