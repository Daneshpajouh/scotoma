#!/usr/bin/env python3
"""Re-earn every figure in the appendix's one-mismatch paragraph, by exhaustion.

Published claims checked here:
  1. For every spacer of length n = 3..10 over {A,C,G,T}, 1,398,080 in total, the omitted i=0 word
     g_2..g_n and the emitted i=1 word g_1 g_3..g_n are at Hamming distance 1[g_1 != g_2] <= 1.
  2. The distance is 0 in exactly 349,520 of those cases, which is one quarter, and those are
     exactly the spacers with g_1 == g_2.
  3. The symmetric statement at the PAM-proximal end holds with indicator 1[g_n != g_{n-1}].
  4. The bound is attained, with an exhibited (spacer, target) pair at every length.
  5. Among targets within budget D of the OMITTED word, the fraction outside budget D of the
     EMITTED one is 75.0 / 50.9 / 28.7 per cent at D = 1, 2, 3 for n = 6, and 69.2 / 40.3 / 15.4
     for n = 5. This is the paragraph's honest limit and the reason it is not a second rescue.

Nothing here is asserted: every number is recomputed and compared to the published figure, and the
script exits non-zero on any disagreement. Runtime is a few seconds; no data files and no network.
"""
import itertools, sys

A = 'ACGT'
PUBLISHED = {
    'spacers_n3_to_n10': 1398080,
    'coincidences': 349520,
    # CONDITIONAL on g_1 != g_2, i.e. restricted to the spacers the omission can hurt at all.
    'pct_exposed': {(6, 1): 75.0, (6, 2): 50.9, (6, 3): 28.7,
                    (5, 1): 69.2, (5, 2): 40.3, (5, 3): 15.4},
    # UNCONDITIONAL, degenerate quarter included. Both are published, because the gate as first
    # written compared my computation against my published value and BOTH carried my convention,
    # so a convention drift was invisible to it. A peer recomputing unconditionally got 56.2%
    # against the published 75.0% and would have had a defensible error report. Checking both,
    # and checking that they stand in the ratio the degenerate quarter forces, closes that.
    'pct_all': {(6, 1): 56.2, (6, 2): 38.2, (6, 3): 21.5,
                (5, 1): 51.9, (5, 2): 30.2, (5, 3): 11.6},
}


def dH(a, b):
    if len(a) != len(b):
        sys.exit('GATE FAILED: comparing words of different lengths, %d and %d; the two deletion '
                 'words must both have length n-1 or the proposition is about the wrong objects'
                 % (len(a), len(b)))
    return sum(x != y for x, y in zip(a, b))


def main():
    total = coincide = 0
    worst = 0
    tight_at = {}
    for n in range(3, 11):
        for g in itertools.product(A, repeat=n):
            g = ''.join(g)
            total += 1
            # distal: omitted deletes g_1, emitted (i=1) deletes g_2
            d0 = dH(g[1:], g[0] + g[2:])
            # proximal: omitted deletes g_n, emitted (i=n-2) deletes g_{n-1}
            d1 = dH(g[:-1], g[:-2] + g[-1])
            if d0 != (0 if g[0] == g[1] else 1):
                sys.exit('GATE FAILED: distal distance %d on %s does not equal 1[g_1 != g_2]' % (d0, g))
            if d1 != (0 if g[-1] == g[-2] else 1):
                sys.exit('GATE FAILED: proximal distance %d on %s does not equal 1[g_n != g_(n-1)]'
                         % (d1, g))
            worst = max(worst, d0, d1)
            if d0 == 0:
                coincide += 1
            if n not in tight_at and d0 == 1:
                # the omitted word itself is at distance 0 from omitted and 1 from emitted
                tight_at[n] = (g, g[1:])

    print('exhaustive over n = 3..10: %d spacers, max distance %d, %d coincidences'
          % (total, worst, coincide))
    if total != PUBLISHED['spacers_n3_to_n10']:
        sys.exit('GATE FAILED: %d spacers, appendix says %d'
                 % (total, PUBLISHED['spacers_n3_to_n10']))
    if coincide != PUBLISHED['coincidences']:
        sys.exit('GATE FAILED: %d coincidences, appendix says %d'
                 % (coincide, PUBLISHED['coincidences']))
    if coincide * 4 != total:
        sys.exit('GATE FAILED: coincidences are %d of %d, not one quarter' % (coincide, total))
    if worst != 1:
        sys.exit('GATE FAILED: max distance %d; the +1 bound is not what this proves' % worst)
    print('  [ok  ] 1,398,080 spacers, bound 1, 349,520 coincidences = exactly one quarter')
    g, t = tight_at[6]
    print('  [ok  ] tightness exhibited, e.g. spacer %s with target %s at distance 0 from the '
          'omitted word and 1 from the emitted' % (g, t))

    fails = []
    for n in (5, 6):
        for D in (1, 2, 3):
            need = tot_all = tot_exp = 0
            for g in itertools.product(A, repeat=n):
                g = ''.join(g)
                p_om, p_em = g[1:], g[0] + g[2:]
                degenerate = (p_om == p_em)       # g_1 == g_2: one word, not two
                for t in itertools.product(A, repeat=n - 1):
                    t = ''.join(t)
                    if dH(t, p_om) <= D:
                        tot_all += 1
                        if not degenerate:
                            tot_exp += 1
                        if dH(t, p_em) > D:
                            if degenerate:
                                sys.exit('GATE FAILED: a spacer with g_1 == g_2 contributed to the '
                                         'numerator, which is impossible: its two words coincide')
                            need += 1
            for label, tot, table in (('exposed', tot_exp, 'pct_exposed'),
                                      ('all    ', tot_all, 'pct_all')):
                pct = 100.0 * need / tot
                want = PUBLISHED[table][(n, D)]
                ok = abs(pct - want) < 0.05
                print('  [%s] n=%d budget %d, %s: %d of %d need the extra mismatch, %.1f%% '
                      '(appendix says %.1f%%)'
                      % ('ok  ' if ok else 'FAIL', n, D, label, need, tot, pct, want))
                if not ok:
                    fails.append('n=%d D=%d %s: measured %.2f%%, appendix says %.1f%%'
                                 % (n, D, label.strip(), pct, want))
            # the two denominators must differ by exactly 4/3, since the degenerate quarter can
            # never contribute to the numerator. This is the invariant that makes a silent
            # convention swap impossible rather than merely unlikely.
            if tot_all * 3 != tot_exp * 4:
                fails.append('n=%d D=%d: denominators %d and %d are not in the ratio 4/3 that the '
                             'degenerate quarter forces' % (n, D, tot_all, tot_exp))
    if fails:
        print('GATE FAILED: ' + '; '.join(fails))
        return 2
    print('\nCLEAN: the appendix paragraph\'s bound, its degenerate cases, its tightness and its '
          'honest limit all reproduce by exhaustion.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
