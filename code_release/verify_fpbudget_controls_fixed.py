#!/usr/bin/env python3
"""The coverage control, redone, after the first one could not fail.

WHAT WENT WRONG. I tested coverage against a family built from NON-CONSECUTIVE blocks, expecting it
to fail. It covered. The reason is structural, not incidental: a block of non-consecutive library
words differs in MORE coordinates, so the block pattern frees MORE wildcards and accepts MORE. A
mispaired family is strictly more permissive, so it can only cover better. The control was incapable
of failing, which makes the coverage result it accompanied unevidenced.

TWO CONTROLS THAT CAN FAIL.
  1. DROP-ONE. Remove one pattern from the family. The words near the orphaned library ball must
     then be uncovered, so coverage must be NO. If it is still YES the family was redundant and the
     count ceil(r/(w+1)) is not tight.
  2. MISPAIRING, scored on the right quantity. Mispairing cannot hurt coverage, but it must hurt the
     FALSE-POSITIVE RATE, which is what the budget actually spends. Consecutive blocks must be
     strictly better than mispaired ones.
"""
from itertools import combinations, product
from math import comb, ceil
from verify_fpbudget_compress_by_wildcards import library, ball_union, block_pattern, accepts, family, ALPHA


def accepted_set(fam, m, N):
    acc = set()
    for p, _ in fam:
        nf = [i for i, c in enumerate(p) if c != '*']
        wl = [i for i, c in enumerate(p) if c == '*']
        for k in range(min(m, len(nf)) + 1):
            for pos in combinations(nf, k):
                for subs in product(*[[c for c in ALPHA if c != p[i]] for i in pos]):
                    base = list(p)
                    for i, c in zip(pos, subs):
                        base[i] = c
                    for fill in product(ALPHA, repeat=len(wl)):
                        y = list(base)
                        for i, c in zip(wl, fill):
                            y[i] = c
                        acc.add(''.join(y))
    return acc


print('=== CONTROL 1: drop one pattern, coverage must become NO ===')
print('  guide    m   w   full covers   dropped covers   uncovered words')
ok1 = True
for g, m in (('ACGTAC', 1), ('ACGTAC', 2), ('ACGTACG', 1), ('AACGTAC', 1)):
    L = library(g); U = ball_union(L, m); N = len(L[0])
    for w in (1, 2):
        fam = family(L, w)
        full = all(any(accepts(p, x, m) for p, _ in fam) for x in U)
        short = fam[:-1]
        miss = [x for x in U if not any(accepts(p, x, m) for p, _ in short)]
        print('  %-8s %-3d %-3d %-13s %-16s %d'
              % (g, m, w, 'YES' if full else 'NO', 'NO' if miss else 'YES', len(miss)))
        if not full or not miss:
            ok1 = False
print('  verdict: %s' % ('drop-one always breaks coverage, so the family is not redundant'
                         if ok1 else 'FAILED, a dropped family still covered: the count is not tight'))
print()

print('=== CONTROL 2: mispairing must cost false-positive RATE ===')
print('  guide    m   w   consecutive rate   mispaired rate   consecutive better?')
ok2 = True
for g, m in (('ACGTAC', 1), ('ACGTAC', 2), ('ACGTACG', 1), ('AACGTAC', 1), ('ACGTACG', 2)):
    L = library(g); U = ball_union(L, m); N = len(L[0])
    for w in (1, 2):
        a = accepted_set(family(L, w, True), m, N)
        b = accepted_set(family(L, w, False), m, N)
        ra = (len(a) - len(U)) / len(a)
        rb = (len(b) - len(U)) / len(b)
        better = ra < rb
        print('  %-8s %-3d %-3d %-18.4f %-16.4f %s'
              % (g, m, w, ra, rb, 'YES' if better else 'NO'))
        if not better:
            ok2 = False
print('  verdict: %s' % ('consecutive blocks are strictly cheaper, as the mechanism predicts'
                         if ok2 else 'FAILED, mispairing was not more expensive'))
