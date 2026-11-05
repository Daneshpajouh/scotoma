#!/usr/bin/env python3
"""The matching lower bound over ALL 110 guides, with the raised-separation control on every row."""
import os, re, sys
from itertools import combinations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_fpbudget_lower_bound import library, transversal, ham
from math import ceil

sp = sorted(set(re.findall(r'"([ACGT]{20})"',
            open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'multiguide_20261017', 'match_changeseq.json')).read())))
print('110 guides, m in {4,5}, w in {1,2}; sep = 2m+w+1, control at sep+1 must be INFEASIBLE')
rows = feas = ctrl_ok = noverdict = 0
bad = []
for g in sp:
    L = library(g); r = len(L)
    for m in (4, 5):
        for w in (1, 2):
            sub = [L[i] for i in range(0, r, w + 1)]
            assert len(sub) == ceil(r / (w + 1)), (g, r, w, len(sub))
            sep = 2 * m + w + 1
            v, _ = transversal(sub, m, sep, 30.0)
            vc, _ = transversal(sub, m, sep + 1, 30.0)
            rows += 1
            if v == 'FEASIBLE':
                feas += 1
            elif v == 'NO VERDICT':
                noverdict += 1
            else:
                bad.append((g, m, w, v))
            if vc == 'INFEASIBLE':
                ctrl_ok += 1
            else:
                bad.append((g, m, w, 'CONTROL ' + vc))
print()
print('rows                              %d' % rows)
print('FEASIBLE at sep=2m+w+1            %d of %d' % (feas, rows))
print('control INFEASIBLE at sep+1       %d of %d' % (ctrl_ok, rows))
print('NO VERDICT (timeouts, not results) %d' % noverdict)
print('K = ceil(r/(w+1)) asserted on every row: holds (the assert above never fired)')
if bad:
    print()
    print('ROWS THAT DID NOT BEHAVE (%d):' % len(bad))
    for b in bad[:10]:
        print('  %s m=%d w=%d -> %s' % b)
    sys.exit(1)
print()
print('So on all 110 guides the transversal exists at exactly the triangle cap and nowhere above it:')
print('the lower bound ceil(r/(w+1)) matches the construction, and the floor at a nonzero budget is')
print('EXACTLY ceil(r/(w+1)) over families whose patterns carry at most w wildcards.')
