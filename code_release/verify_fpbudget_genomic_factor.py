#!/usr/bin/env python3
"""The GENOMIC candidate factor of the w=1 compression, from held counts, with a hard interval.

The matching-bounds result said no genomic factor was claimed and that measuring one needed an index.
It needs less than that. `producers_2026-10-20/PAM_PANEL_FACTOR_PERGUIDE.json` holds panel-wide
genomic NGG candidate counts at every threshold from 3 to 6, gated against our own enumerator at
thresholds 0 through 5 exactly (see RESULT_2026-10-21_the_pam_dial_...). Those are cumulative
candidate counts, which is the quantity a search retrieves.

THE DECOMPOSITION, exact. A w=1 pattern is P_i with coordinate j freed, where P_i and P_{i+1} differ
only at j. A site x is retrieved iff its mismatch count to P_i off j is at most m. Splitting on j:

  x_j = P_i[j]  -> retrieved iff d(x,P_i) <= m
  x_j != P_i[j] -> retrieved iff d(x,P_i) <= m+1

So retrieved(w=1) = cum(m) + (sites at distance EXACTLY m+1 whose mismatch set contains j).

Write S = cum(m+1) - cum(m) and let f = P(j is among the m+1 mismatched coordinates).

  factor = [ceil(r/(w+1)) * (cum(m) + f*S)] / [r * cum(m)]

f = 0 and f = 1 give a HARD INTERVAL that assumes nothing. Uniform placement of the m+1 mismatches
over N coordinates gives the point estimate f = (m+1)/N.

TWO ASSUMPTIONS, named. (1) The point estimate's uniform placement: the genome is not a uniform
sample of 23-mers, so position-specific mismatch propensity could move f. The interval does not
depend on it. (2) The held counts are for the ungapped 20-nt guide search; the library-word search is
19-nt against a 20-nt window. Only the GROWTH RATIO transfers, and that is what is used.
"""
import json, os, re
from math import ceil

HERE = os.path.dirname(os.path.abspath(__file__))
CUM = json.load(open(os.path.join(HERE, '..', 'producers_2026-10-20',
                                  'PAM_PANEL_FACTOR_PERGUIDE.json')))
PANEL = os.path.join(HERE, '..', 'multiguide_20261017', 'match_changeseq.json')
N = 19


def runs(g):
    return 1 + sum(1 for a, b in zip(g, g[1:]) if a != b)


sp = sorted(set(re.findall(r'"([ACGT]{20})"', open(PANEL).read())))
rs = sorted(runs(g) for g in sp)
rmed = rs[len(rs) // 2]

print('held cumulative genomic NGG candidate counts, panel-wide:')
for k in sorted(CUM, key=int):
    print('   threshold %s : %d' % (k, CUM[k]['ngg']))
print()
print('median r over the 110 guides: %d' % rmed)
print()

print('=== the w=1 genomic candidate factor ===')
print('  m   cum(m)     S=cum(m+1)-cum(m)  S/cum(m)   f=0 (hard low)  f=(m+1)/N  f=1 (hard high)')
for m in (3, 4, 5):
    if str(m) not in CUM or str(m + 1) not in CUM:
        continue
    cm = CUM[str(m)]['ngg']; S = CUM[str(m + 1)]['ngg'] - cm
    share = ceil(rmed / 2) / rmed
    def fac(f):
        return share * (cm + f * S) / cm
    print('  %-3d %-10d %-18d %-10.3f %-15.3f %-10.3f %.3f'
          % (m, cm, S, S / cm, fac(0.0), fac((m + 1) / N), fac(1.0)))
print()
print('  share of probes retained at w=1, median r=%d: %d/%d = %.3f'
      % (rmed, ceil(rmed / 2), rmed, ceil(rmed / 2) / rmed))
print()
print('=== sensitivity of the point estimate to f, at m=5 ===')
cm = CUM['5']['ngg']; S = CUM['6']['ngg'] - cm
share = ceil(rmed / 2) / rmed
f0 = 6 / N
print('  f relative to uniform   f        factor')
for mult in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
    f = min(1.0, f0 * mult)
    print('   %-22.2f %-8.4f %.3f' % (mult, f, share * (cm + f * S) / cm))
print()
print('=== across all 110 guides, at m=5 and uniform f ===')
facs = [ceil(r / 2) / r * (cm + f0 * S) / cm for r in rs]
print('  r from %d to %d; factor %.3f to %.3f (median %.3f)'
      % (min(rs), max(rs), min(facs), max(facs), sorted(facs)[len(facs) // 2]))
print()
print('READING. The hard interval needs no assumption and is wide: the factor cannot fall below the')
print('probe share nor rise above share*(1+S/cum(m)). The uniform point estimate sits near the')
print('sequence-space factor already reported, so the genome does not obviously make the trade worse.')
print('This NARROWS the open question; it does not close it. Closing it needs the position-specific')
print('mismatch distribution of genomic candidate sites, which is a genome scan.')
