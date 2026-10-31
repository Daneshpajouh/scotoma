#!/usr/bin/env python3
"""Gate: the paper's structural claim about CFD's gap entries, recomputed from the shipped table.

The paper states that CFD "does score bulges", that its mismatch table "carries a gap penalty at
every guide position from 2 to 20 in both orientations", and that position 1 is the only position
without one. That claim is load-bearing twice: it corrects a standing asset that said the opposite,
and it is the premise of the measured asymmetry between the two bulge orientations.

WHY THE TABLE IS NOT IN THIS RELEASE. It is CRISPRme's file, not ours, so it is not redistributed.
Point this gate at a local copy instead:

    python3 code_release/verify_cfd_gap_table.py path/to/mismatch_score.pkl

or set CFD_MM_PKL. With no path it reports that it did not run, which is an invocation state and exits 2, never 0, so no runner can read a skipped check as a pass and
not a finding, the same contract verify_missed_activity.py uses for the CHANGE-seq deposit.

A PICKLE IS EXECUTABLE. This reads it with an unpickler that refuses every global, so a tampered
file raises instead of running. Do not replace that with a bare pickle.load.
"""
import io, os, pickle, re, sys

CLAIM_GAP_ENTRIES = 152
CLAIM_POSITIONS = list(range(2, 21))


class _Restricted(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError('refused global %s.%s' % (module, name))


def audit(path):
    with open(path, 'rb') as fh:
        mm = _Restricted(io.BytesIO(fh.read())).load()
    if not isinstance(mm, dict):
        return ['table is %s, expected dict' % type(mm).__name__], {}
    gap = [k for k in mm if '-' in k]
    guide_side = [k for k in gap if re.match(r'^r-:d', k)]
    target_side = [k for k in gap if re.match(r'^r.:d-', k)]

    def pos(keys):
        out = set()
        for k in keys:
            try:
                out.add(int(k.split(',')[1]))
            except (IndexError, ValueError):
                pass
        return sorted(out)

    import statistics as _st
    gv = sorted(mm[k] for k in guide_side)
    tv = sorted(mm[k] for k in target_side)
    got = {'entries': len(mm), 'gap': len(gap), 'guide_side': len(guide_side),
           'target_side': len(target_side), 'pos_guide': pos(guide_side),
           'pos_target': pos(target_side),
           'med_guide': _st.median(gv) if gv else None,
           'med_target': _st.median(tv) if tv else None,
           'zero_guide': sum(1 for x in gv if x == 0),
           'zero_target': sum(1 for x in tv if x == 0)}
    bad = []
    if got['gap'] != CLAIM_GAP_ENTRIES:
        bad.append('%d gap entries, paper states %d' % (got['gap'], CLAIM_GAP_ENTRIES))
    for side in ('guide', 'target'):
        if got['pos_' + side] != CLAIM_POSITIONS:
            bad.append('%s-side gap positions %s, paper states 2..20'
                       % (side, got['pos_' + side]))
    if not got['guide_side'] or not got['target_side']:
        bad.append('one orientation is empty, so "both orientations" is false')
    # The MECHANISM behind the measured asymmetry, from the table alone and with no indexing
    # convention involved: the DNA family is harsher and far more often exactly zero, and a
    # multiplicative zero annihilates a score whatever every other position says.
    if got['med_guide'] is not None and got['med_target'] is not None:
        if not got['med_guide'] < got['med_target']:
            bad.append('the guide-side (DNA) median %.4f is not harsher than the target-side '
                       '(RNA) median %.4f, so the stated mechanism does not hold in this table'
                       % (got['med_guide'], got['med_target']))
        if got['zero_guide'] <= got['zero_target']:
            bad.append('the guide-side (DNA) family has %d zero entries against the target-side '
                       '%d; the paper says the DNA family is the one more often annihilating'
                       % (got['zero_guide'], got['zero_target']))
    return bad, got


def main():
    path = (sys.argv[1] if len(sys.argv) > 1 else os.environ.get('CFD_MM_PKL'))
    if not path or not os.path.exists(path):
        print('NOT CHECKED HERE: CRISPRme\'s mismatch_score.pkl was not supplied, so the gap-table')
        print('  claim was not recomputed. This is an invocation state, not a finding. Pass the')
        print('  path as an argument or set CFD_MM_PKL. The file is not redistributed with this')
        print('  release because it is not ours.')
        return 2
    bad, got = audit(path)
    print('  table entries            %d' % got.get('entries', -1))
    print('  gap entries              %d   paper states %d' % (got.get('gap', -1), CLAIM_GAP_ENTRIES))
    print('  guide-side / target-side %d / %d' % (got.get('guide_side', -1), got.get('target_side', -1)))
    print('  positions, guide side    %s' % (got.get('pos_guide'),))
    print('  positions, target side   %s' % (got.get('pos_target'),))
    print('  position 1 present       %s' % (1 in (got.get('pos_guide') or []) or
                                             1 in (got.get('pos_target') or [])))
    print('  median factor, DNA / RNA %.4f / %.4f   (smaller is harsher)'
          % (got.get('med_guide') or -1, got.get('med_target') or -1))
    print('  entries exactly zero     DNA %d, RNA %d' % (got.get('zero_guide', -1),
                                                         got.get('zero_target', -1)))
    if bad:
        print('\nCFD GAP TABLE GATE FAILURES:')
        for b in bad:
            print('  [FAIL] %s' % b)
        return 1
    print('\nCFD does score bulges: the shipped table carries a gap penalty at every guide position')
    print('from 2 to 20 in both orientations, and position 1 is the only one without.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
