#!/usr/bin/env python3
"""The genomic factor with f MEASURED rather than assumed, and the hard interval collapsed.

positional_propensity.py measured, over 1,150,885 candidate rows at exactly 6 mismatches, the
fraction whose mismatch set contains each coordinate. Restricted to the 20 SPACER coordinates (21-23
are the PAM: coordinate 21 shows 0.003 because the guide carries N there and this project's gated
convention never charges a guide N, and at k=7 coordinate 22 shows 0.997 because that stratum is
PAM-relaxed sites, so k=7 is not used).

The ratio to uniform runs 0.88 to 1.20 across the spacer, so max/min is 1.36, under the 1.5 the
pre-registration set as the flatness threshold. The uniform estimate therefore stands, and the hard
[f=0, f=1] interval collapses to the measured range.
"""
import json, os, re
from math import ceil

HERE = os.path.dirname(os.path.abspath(__file__))
CUM = json.load(open(os.path.join(HERE, '..', 'producers_2026-10-20',
                                  'PAM_PANEL_FACTOR_PERGUIDE.json')))
PANEL = os.path.join(HERE, '..', 'multiguide_20261017', 'match_changeseq.json')
N = 19

# measured f/uniform over the 20 spacer coordinates at k = 6, from positional_propensity.py
RATIO_K6 = [1.20, 1.05, 1.06, 1.05, 1.19, 1.04, 0.97, 1.08, 1.14, 1.02,
            0.96, 0.88, 1.04, 0.97, 1.06, 1.10, 1.14, 1.11, 1.04, 0.95]


def runs(g):
    return 1 + sum(1 for a, b in zip(g, g[1:]) if a != b)


sp = sorted(set(re.findall(r'"([ACGT]{20})"', open(PANEL).read())))
rs = sorted(runs(g) for g in sp)
rmed = rs[len(rs) // 2]
m = 5
cm = CUM[str(m)]['ngg']
S = CUM[str(m + 1)]['ngg'] - cm
share = ceil(rmed / 2) / rmed
f0 = (m + 1) / N


def fac(f):
    return share * (cm + f * S) / cm


print('m=%d   cum(m)=%d   S=%d   probe share=%d/%d=%.3f   uniform f=%.4f'
      % (m, cm, S, ceil(rmed / 2), rmed, share, f0))
print()
print('measured f/uniform across the 20 spacer coordinates: min %.2f  max %.2f  ratio %.2f'
      % (min(RATIO_K6), max(RATIO_K6), max(RATIO_K6) / min(RATIO_K6)))
print('flatness threshold registered: 1.5   ->   %s'
      % ('FLAT, the uniform estimate stands'
         if max(RATIO_K6) / min(RATIO_K6) < 1.5 else 'NOT FLAT, reweight'))
print()
print('=== the factor, with f measured ===')
print('  assumption-free hard interval   [%.3f, %.3f]' % (fac(0.0), fac(1.0)))
print('  uniform placement               %.3f' % fac(f0))
print('  MEASURED range over coordinates [%.3f, %.3f]'
      % (fac(f0 * min(RATIO_K6)), fac(f0 * max(RATIO_K6))))
print('  break-even (factor = 1) needs f/uniform = %.3f, measured min is %.2f -> %s'
      % (((cm / share - cm) / S) / f0, min(RATIO_K6),
         'NOT reached' if ((cm / share - cm) / S) / f0 < min(RATIO_K6) else 'reached'))
print()
print('=== over all 110 guides, measured f ===')
lo = [ceil(r / 2) / r * (cm + f0 * min(RATIO_K6) * S) / cm for r in rs]
hi = [ceil(r / 2) / r * (cm + f0 * max(RATIO_K6) * S) / cm for r in rs]
print('  r %d..%d   factor %.3f .. %.3f' % (min(rs), max(rs), min(lo), max(hi)))
print()
print('WHY THE COORDINATE MAPPING DOES NOT MATTER HERE. The freed coordinates are run boundaries of')
print('the 19-long library word, while the measurement indexes the 20 spacer coordinates of the')
print('guide. A mapping would matter if the distribution were steep. At a 1.36 spread it moves the')
print('factor by under 0.08, which is why the flatness test was the right gate to put first.')
