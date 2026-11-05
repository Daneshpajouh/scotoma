#!/usr/bin/env python3
"""Positional mismatch propensity of candidate sites, from the held deposit's alignments.

Pre-registered in PREREG_2026-10-24_the_positional_mismatch_propensity.md.

Decides f, the chance that a given coordinate is among the m+1 mismatches of a distance-(m+1)
candidate. The genomic factor's point estimate assumed f = (m+1)/N; break-even sits at half that.

Mismatch convention is this project's own, already gated on 12,463,279 rows: a guide 'N' is a
wildcard and is never charged; a target 'N' is charged. Rows with a gap on either side are skipped,
because this is about the mismatch-only geometry.
"""
import csv, subprocess, sys, os
from collections import Counter

# The dataset archive is NOT in this repository. Give its path as the first argument, or set
# CFD_DATASETS_ZIP. A hardcoded session path made the released copy unrunnable for a reader.
ZIP = (sys.argv[1] if len(sys.argv) > 1
       else os.environ.get('CFD_DATASETS_ZIP', 'datasets.zip'))
if not os.path.exists(ZIP):
    print('NOT RUN: datasets.zip not found. Pass its path as the first argument,')
    print('         or set CFD_DATASETS_ZIP. Looked for: %s' % os.path.abspath(ZIP))
    raise SystemExit(2)          # 0 passed, 1 FAILED, 2 NOT RUN
MEMBER = 'datasets/CrisprNet.csv'
csv.field_size_limit(10 ** 7)

p = subprocess.Popen(['/usr/bin/unzip', '-p', ZIP, MEMBER],
                     stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
r = csv.DictReader(p.stdout)
if 'Align.sgRNA' not in (r.fieldnames or []):
    print('NOT RUN: Align.sgRNA absent; columns are %s' % (r.fieldnames or [])[:8])
    sys.exit(2)

rows = kept = 0
lens = Counter()
# pos_hits[k][i] = rows at exactly k mismatches whose mismatch set contains coordinate i
pos_hits = {}
kcount = Counter()

for row in r:
    rows += 1
    gs = row.get('Align.sgRNA') or ''
    ts = row.get('Align.off-target') or ''
    if not gs or not ts or len(gs) != len(ts):
        continue
    if '-' in gs or '-' in ts:
        continue
    lens[len(gs)] += 1
    mm = [i for i, (a, b) in enumerate(zip(gs, ts))
          if a != 'N' and a != b]
    k = len(mm)
    if k == 0 or k > 8:
        continue
    kept += 1
    kcount[k] += 1
    d = pos_hits.setdefault(k, Counter())
    for i in mm:
        d[i] += 1

print('rows read %d   mismatch-only rows kept %d' % (rows, kept))
print('alignment lengths seen: %s' % dict(lens.most_common(4)))
print()
N = lens.most_common(1)[0][0] if lens else 0
print('=== fraction of rows at exactly k mismatches whose mismatch set contains coordinate i ===')
print('   (uniform placement would give k/N at every i)')
for k in sorted(kcount):
    if kcount[k] < 2000:
        continue
    d = pos_hits[k]
    fr = [d.get(i, 0) / kcount[k] for i in range(N)]
    unif = k / N
    lo, hi = min(fr), max(fr)
    print()
    print(' k=%d  n=%d  uniform=%.4f   observed min %.4f  max %.4f  ratio %.2f'
          % (k, kcount[k], unif, lo, hi, (hi / lo if lo > 0 else float('inf'))))
    print('   i :  ' + ' '.join('%5d' % (i + 1) for i in range(N)))
    print('   f :  ' + ' '.join('%5.3f' % v for v in fr))
    print('   /u:  ' + ' '.join('%5.2f' % (v / unif) for v in fr))
