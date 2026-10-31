#!/usr/bin/env python3
"""The 2.46% median and 2.88% tolerance bound: the method, checked, and reproduced when given its data.

WHY THIS EXISTS. The paper states, twice, "a simultaneous order-statistic tolerance bound across the
nested family of thresholds from zero to four mismatches, at family-wise alpha = 0.05, places that
share at or below 2.88%". On 2026-10-24 a search established that **no code anywhere, published or
local, computed it**: `git grep` over the published tree and a grep over every local .py for
`simultaneous`, `family-wise` and `order-statistic` returned nothing, the figure was absent from the
110-row figure manifest, and no gate bound it.

The provenance was found on 2026-10-25 in
`PORTFOLIO/06_markdown_manuscripts/crispr_pangenome/CORRECTED_RESULTS_2026-08-18.md`, which records
the method exactly: order statistic, n = 459 haplotypes, simultaneous over the pre-specified nested
family mm <= 0..4 at family-wise alpha = 0.05, **per-member alpha = 0.01, k = 456**, with this table:

    mm<=0  median 0.0000  bound 0.0027
    mm<=1  median 0.0139  bound 0.0279
    mm<=2  median 0.0206  bound 0.0254
    mm<=3  median 0.0238  bound 0.0279
    mm<=4  median 0.0246  bound 0.0288     <- the paper's 2.46% and 2.88%

WHAT IS CHECKED WITHOUT THE DATA. The method itself, which is arithmetic and needs no inputs:
  * the k-th order statistic of n observations is a one-sided distribution-free prediction bound with
    coverage k/(n+1), so k must satisfy k/(n+1) >= 1 - alpha;
  * k = 456 gives 456/460 = 0.991304 >= 0.99, and k = 455 gives 0.989130, which fails. So the stated
    k is the TIGHTEST valid choice and not an arbitrary one;
  * Bonferroni over the five nested thresholds gives 0.05/5 = 0.01, which is the stated per-member
    alpha.

WHAT NEEDS THE DATA. Reproducing the five bounds requires the 459 per-haplotype shares, one per
threshold. Those are not in this repository and were not located on fir or narval. Give this script
a JSON holding them and it recomputes every row; without it the script exits 2 NOT RUN rather than
pretending.
"""
import json, os, sys

N, K, FAMILY_ALPHA, THRESHOLDS = 459, 456, 0.05, 5
PUBLISHED = {0: (0.0000, 0.0027), 1: (0.0139, 0.0279), 2: (0.0206, 0.0254),
             3: (0.0238, 0.0279), 4: (0.0246, 0.0288)}


def check_method():
    per_member = FAMILY_ALPHA / THRESHOLDS
    cov = K / (N + 1)
    ok = cov >= 1 - per_member
    tight = (K - 1) / (N + 1) < 1 - per_member
    print('method, which needs no data:')
    print('  family-wise alpha %.2f over %d nested thresholds -> per-member %.3f'
          % (FAMILY_ALPHA, THRESHOLDS, per_member))
    print('  coverage of the k-th order statistic = k/(n+1) = %d/%d = %.6f' % (K, N + 1, cov))
    print('  needs >= %.2f : %s' % (1 - per_member, 'OK' if ok else 'FAILS'))
    print('  k-1 = %d would give %.6f : %s, so k = %d is the TIGHTEST valid index'
          % (K - 1, (K - 1) / (N + 1), 'fails as it must' if tight else 'ALSO VALID, so k is not tight', K))
    return ok and tight


def check_data(path):
    d = json.load(open(path))
    shares = d.get('shares') if isinstance(d, dict) else d
    if not isinstance(shares, dict):
        print('NOT RUN: expected {"shares": {"<mm>": [459 floats], ...}}')
        return 2
    bad = []
    print()
    print('  mm  n     median      published   k-th (k=%d)  published   agree' % K)
    for mm in sorted(PUBLISHED):
        v = sorted(shares.get(str(mm), []))
        if len(v) != N:
            print('  %-3d NOT RUN: %d values, expected %d' % (mm, len(v), N)); return 2
        med = v[len(v) // 2] if len(v) % 2 else (v[len(v)//2 - 1] + v[len(v)//2]) / 2
        bound = v[K - 1]
        pm, pb = PUBLISHED[mm]
        agree = abs(med - pm) < 5e-5 and abs(bound - pb) < 5e-5
        if not agree:
            bad.append(mm)
        print('  %-3d %-5d %-11.4f %-11.4f %-12.4f %-11.4f %s'
              % (mm, len(v), med, pm, bound, pb, 'yes' if agree else 'NO'))
    if bad:
        print('\n[FAIL] %d threshold(s) do not reproduce: %s' % (len(bad), bad)); return 1
    print('\nCLEAN: every published median and bound recomputes from the shares')
    return 0


if __name__ == '__main__':
    good = check_method()
    if not good:
        print('[FAIL] the stated method does not check out'); sys.exit(1)
    print('  method OK')
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        sys.exit(check_data(sys.argv[1]))
    print()
    print('NOT RUN against data: the 459 per-haplotype shares are not in this repository and were')
    print('not located on fir or narval. Pass a JSON of them as the first argument to reproduce the')
    print('five published rows. The method above is checked regardless.')
    sys.exit(2)
