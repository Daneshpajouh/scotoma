"""Gate: the run count the index subset reference exact duals reach must equal what the certificates prove.

The appendix prints a boundary table and the abstract quotes a panel-guide count derived from it.
Both are downstream of data/floor_certificates/INDEX.json, so both can drift silently when a cell is
added or withdrawn. This recomputes them.

A cell counts as reached only if its stored mass actually CLEARS the criterion. The criterion is
K >= ceil(tau*), so a dual proves K = r only when its mass EXCEEDS r-1; a mass that merely rounds up
to r is not enough, and testing >= r instead of > r-1 is an error this project has made and paid
for. The test here is on exact rationals, never on a float.
"""
import json, os, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
IDX  = os.environ.get('CEILING_INDEX',
                      os.path.join(HERE, 'data', 'floor_certificates', 'INDEX.json'))

# the panel's homopolymer-run distribution, from bulge_census_full_idx_dna.json
RUNS = {9: 1, 10: 2, 11: 2, 12: 10, 13: 11, 14: 19, 15: 16, 16: 28, 17: 15, 18: 4, 19: 1, 20: 1}
PANEL = sum(RUNS.values())

QUOTED_REACH_M8   = 16      # appendix table, unrestricted column; abstract "up to $12$"
# The ABSTRACT no longer quotes this seat's own reach: Theorem thm:eight settles all 110 at m=8
# from the r=20 dual, which verify_floor_new_results.py gates. What this gate still guards is the
# corroborating run of cells this seat proved independently, which the appendix quotes, so the two
# numbers are about different claims and must not be conflated again.
QUOTED_PANEL_SETTLED = 89   # appendix: this seat's own cells corroborate 89 of the 110

# 2026-10-25: m=9 and m=10 were printed in the appendix and gated NOWHERE. Nine verified cells were
# added to the index and this gate's output did not change by one character, which is the silent
# drift its own docstring warns about, left open in the two rows that moved today.
QUOTED_REACH_M9        = 16   # historical index subset; the full budget-nine proof has separate replays
QUOTED_REACH_M10       = 10   # appendix unrestricted-dual column, was 7
QUOTED_PANEL_SETTLED_M9 = 89  # historical index subset, not the current universal budget-nine claim

try:
    idx = json.load(open(IDX))
except FileNotFoundError:
    # See verify_ceiling_certificate.py: a traceback reads as broken software, a sentence reads as
    # a missing file. Overridable so the index can sit anywhere a reader puts it.
    print("GATE NOT RUN: the certificate index is not at %s. Set CEILING_INDEX to its path, or\n"
          "run from a tree that contains data/floor_certificates/INDEX.json.\n"
          "This is an incomplete release, not a finding about the claim." % IDX, file=sys.stderr)
    sys.exit(2)   # 0 passed, 1 FAILED, 2 NOT RUN. sys.exit(str) would have exited 1, i.e. FAILED.
bad = []
proved = {}
for c in idx['cells']:
    if c.get('status') != 'PROVED': continue
    # The panel is 20-mers, so N=19. Cells at other spacer lengths are a robustness check and must
    # not be counted into the panel reach; without this filter an N=20 cell would inflate the m=8 row
    # and the guide count derived from it.
    if c.get('N', 19) != 19: continue
    mass = Fr(c['exact_total_mass'])
    if not mass > c['r'] - 1:
        bad.append('m=%d r=%d is marked PROVED with mass %s, which does not exceed r-1=%d'
                   % (c['m'], c['r'], c['exact_total_mass'], c['r']-1))
        continue
    if -(-mass.numerator // mass.denominator) != c['r']:
        bad.append('m=%d r=%d: ceil of %s is not r' % (c['m'], c['r'], c['exact_total_mass']))
        continue
    proved.setdefault(c['m'], set()).add(c['r'])

# the reach at a budget is the top of the INITIAL RUN of proved cells, not the largest one proved.
# A gap would mean some run count below the maximum is unproved, and "reaches r" would be false.
reach = {}
for mm, rs in proved.items():
    k = min(rs)
    while k + 1 in rs: k += 1
    reach[mm] = k
    if set(range(min(rs), k+1)) != rs:
        bad.append('m=%d: proved cells %s are not contiguous, so a single reach figure is wrong'
                   % (mm, sorted(rs)))

got8 = reach.get(8)
print('  m=8 proved run counts %s -> reach %s   index subset reference %d   %s'
      % (sorted(proved.get(8, ())), got8, QUOTED_REACH_M8,
         'ok' if got8 == QUOTED_REACH_M8 else 'MISMATCH'))
if got8 != QUOTED_REACH_M8:
    bad.append('m=8 reach is %s, index subset reference %d' % (got8, QUOTED_REACH_M8))

settled = sum(v for k, v in RUNS.items() if got8 is not None and k <= got8)
print('  panel guides settled at m=8: %d of %d   index subset reference %d   %s'
      % (settled, PANEL, QUOTED_PANEL_SETTLED,
         'ok' if settled == QUOTED_PANEL_SETTLED else 'MISMATCH'))
if settled != QUOTED_PANEL_SETTLED:
    bad.append('panel settled is %d, index subset reference %d' % (settled, QUOTED_PANEL_SETTLED))
if PANEL != 110:
    bad.append('run distribution sums to %d, not the 110 guides' % PANEL)

# 2026-10-24: this gate took its SCOPE FROM ITS INPUT. The reach is the top of the run starting at
# min(rs), the lowest cell that happens to be in the index, not at the panel's lowest run count,
# while `settled` above counts every guide whose run count is at or below that reach. Deleting the
# m=8 cells r=7 to 10 left the run 11 to 16, reach 16, settled 89: both figures matched the paper
# and the gate exited 0, with the three guides at run counts 9 and 10 counted as settled and no
# certificate proving them. Removing evidence made the gate agree with the paper.
# The two must be consistent: every run count the panel actually contains, up to the reach, needs
# a proved cell, or the settled total is counting guides nothing settles.
if got8 is not None:
    needed = {k for k in RUNS if k <= got8}
    missing = sorted(needed - proved.get(8, set()))
    print('  panel run counts at or below the reach, all proved: %s'
          % ('yes' if not missing else 'NO, missing %s' % missing))
    if missing:
        bad.append('run counts %s are at or below the m=8 reach of %d and the panel contains '
                   'them, but no certificate proves them, so the %d settled guides are not all '
                   'settled' % (missing, got8, settled))

# The appendix also prints the INCREMENT at each newly reached run count, and those three numbers
# have to add up to the jump in the total. A gate on the total alone passes while the breakdown is
# wrong, which is how "11 at 11, a further 10 at 12 and 11 more at 13" got written for a jump of 23.
QUOTED_FROM = 3                      # proved region before this round, run counts up to 10
QUOTED_INCREMENTS = {11: 2, 12: 10, 13: 11, 14: 19, 15: 16, 16: 28}
inc_ok = all(RUNS.get(k) == v for k, v in QUOTED_INCREMENTS.items())
jump = sum(QUOTED_INCREMENTS.values())
print('  increments %s sum to %d; %d + %d = %d   index subset reference %d   %s'
      % (QUOTED_INCREMENTS, jump, QUOTED_FROM, jump, QUOTED_FROM + jump, QUOTED_PANEL_SETTLED,
         'ok' if inc_ok and QUOTED_FROM + jump == QUOTED_PANEL_SETTLED else 'MISMATCH'))
for k, v in QUOTED_INCREMENTS.items():
    if RUNS.get(k) != v:
        bad.append('increment at run count %d quoted as %d, distribution says %s'
                   % (k, v, RUNS.get(k)))
if QUOTED_FROM + jump != QUOTED_PANEL_SETTLED:
    bad.append('increments sum to %d from %d, giving %d, but the total quoted is %d'
               % (jump, QUOTED_FROM, QUOTED_FROM + jump, QUOTED_PANEL_SETTLED))
if QUOTED_FROM != sum(v for k, v in RUNS.items() if k <= 10):
    bad.append('the prior region is quoted as %d guides, distribution gives %d'
               % (QUOTED_FROM, sum(v for k, v in RUNS.items() if k <= 10)))

# The appendix quotes five cells at three spacer lengths as a robustness check. They are excluded
# from the panel reach above, so nothing else would notice if one were withdrawn or its mass moved.
LENGTH_CELLS = {(N, 8, r) for N in (19, 20, 21) for r in (11, 12)}
have = {(c.get('N', 19), c['m'], c['r']) for c in idx['cells'] if c.get('status') == 'PROVED'}
missing = LENGTH_CELLS - have
print('  length-robustness cells present: %d of %d   %s'
      % (len(LENGTH_CELLS) - len(missing), len(LENGTH_CELLS), 'ok' if not missing else 'MISSING %s' % sorted(missing)))
if missing:
    bad.append('appendix quotes length-robustness cells that are not in the index: %s' % sorted(missing))

# The same three checks the m=8 row gets, applied to the two rows the appendix also prints. Written
# as a loop so a future row cannot be added to the paper and left ungated the way these two were.
for mm, quoted, quoted_settled in ((9, QUOTED_REACH_M9, QUOTED_PANEL_SETTLED_M9),
                                   (10, QUOTED_REACH_M10, None)):
    got = reach.get(mm)
    print('  m=%d proved run counts %s -> reach %s   index subset reference %d   %s'
          % (mm, sorted(proved.get(mm, ())), got, quoted,
             'ok' if got == quoted else 'MISMATCH'))
    if got != quoted:
        bad.append('m=%d reach is %s, index subset reference %d' % (mm, got, quoted))
    if got is not None:
        need = {k for k in RUNS if k <= got}
        miss = sorted(need - proved.get(mm, set()))
        if miss:
            bad.append('run counts %s are at or below the m=%d reach of %d and the panel contains '
                       'them, but no certificate proves them' % (miss, mm, got))
        st = sum(v for k, v in RUNS.items() if k <= got)
        if quoted_settled is None:
            print('    panel guides settled at m=%d: %d of %d   (not quoted in the paper)'
                  % (mm, st, PANEL))
        else:
            print('    panel guides settled at m=%d: %d of %d   index subset reference %d   %s'
                  % (mm, st, PANEL, quoted_settled,
                     'ok' if st == quoted_settled else 'MISMATCH'))
            if st != quoted_settled:
                bad.append('panel settled at m=%d is %d, index subset reference %d'
                           % (mm, st, quoted_settled))

if bad:
    print('\nCEILING REACH GATE FAILURES:')
    for b in bad: print('  [FAIL] %s' % b)
    sys.exit(1)
print('the reach quoted in the paper is what the exact certificates prove')
