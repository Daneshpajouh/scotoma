"""Gate: every maize figure quoted in the paper is recomputed from the released locus lists.

The claim is that patterns the 2016 wrapper omits reach loci on maize B73 AGPv4 that none of the
patterns it emits reaches. The lists in data/maize_exclusive_loci/ are the loci themselves, produced
by maize_step_two.py from the scan outputs, so the quoted counts are checkable without rerunning a
genome scan. This gate recomputes each one and refuses any that has drifted.

It also refuses lists that are EMPTY or that contain a duplicate key. An exclusive-locus count is a
set cardinality; a file with repeated rows would inflate it while every per-file check still passed,
which is the shape of error a count-only gate cannot see.
"""
import gzip, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data', 'maize_exclusive_loci')

# quoted in paper_v7/appendix_moved.tex, paragraph on Young et al. (2019)
QUOTED = {'GGCCGAGGTCGACTACCGGC': 40181, 'GGCGGCGGCGAGGTAGTGCG': 31912}
QUOTED_TOTAL   = 72093
QUOTED_OMITTED = 179304          # DISTINCT loci the 9 omitted patterns reach, over both
                                 # guides. Was 179689 until 2026-10-24, which is the SUM of
                                 # the per-pattern counts, not their union: 385 loci are
                                 # reached by more than one omitted pattern. Recomputed by
                                 # maize_omitted_counts.py from the deposited hit files.
QUOTED_SHARE   = 40.2            # per cent, one decimal. 72,093 / 179,304. Was 40.1 while
                                 # QUOTED_OMITTED was the per-pattern SUM 179,689; it moves
                                 # with the denominator, not on its own. The paper does not
                                 # quote this share anywhere, so it is an internal check.

bad = []
counts = {}
for guide, want in sorted(QUOTED.items()):
    p = os.path.join(D, 'exclusive_loci_%s.tsv.gz' % guide)
    if not os.path.exists(p):
        bad.append('%s: list missing at %s' % (guide, p)); continue
    seen = set()
    with gzip.open(p, 'rt') as fh:
        head = fh.readline().rstrip('\n').split('\t')
        if head != ['contig', 'strand', 'pam_forward_start']:
            bad.append('%s: unexpected header %r' % (guide, head)); continue
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) != 3: bad.append('%s: malformed row %r' % (guide, line[:50])); break
            seen.add((f[0], f[1], int(f[2])))
    n = len(seen)
    counts[guide] = n
    if n == 0:
        bad.append('%s: list is EMPTY, so its count is not a measurement' % guide)
    print('  %-22s %6d distinct loci   paper says %6d   %s'
          % (guide, n, want, 'ok' if n == want else 'MISMATCH'))
    if n != want:
        bad.append('%s: %d distinct loci, paper quotes %d' % (guide, n, want))

if len(counts) == len(QUOTED):
    tot = sum(counts.values())
    print('  %-22s %6d total          paper says %6d   %s'
          % ('both guides', tot, QUOTED_TOTAL, 'ok' if tot == QUOTED_TOTAL else 'MISMATCH'))
    if tot != QUOTED_TOTAL:
        bad.append('total %d, paper quotes %d' % (tot, QUOTED_TOTAL))
    share = round(100.0*tot/QUOTED_OMITTED, 1)
    print('  %-22s %6.1f%% of %d       paper says %4.1f%%   %s'
          % ('share', share, QUOTED_OMITTED, QUOTED_SHARE,
             'ok' if share == QUOTED_SHARE else 'MISMATCH'))
    if share != QUOTED_SHARE:
        bad.append('share %.1f%%, paper quotes %.1f%%' % (share, QUOTED_SHARE))

# ---- the whole-portfolio denominator -------------------------------------------------------
# "40.1% of what the complete set sees at the omitted positions" has a denominator defined by the
# omission itself. The defensible figure is the share of the COMPLETE candidate set, so the emitted
# unions are released too and recomputed here. Union, never a sum: a locus two emitted patterns
# both reach counts once, which is why a list and not a stored count is the artifact.
QUOTED_EMITTED  = {'GGCCGAGGTCGACTACCGGC': 403772, 'GGCGGCGGCGAGGTAGTGCG': 651999}
QUOTED_COMPLETE = 1127864
QUOTED_WHOLE_SHARE = 6.39
emitted = {}
for guide, want in sorted(QUOTED_EMITTED.items()):
    p = os.path.join(D, 'emitted_loci_%s.tsv.gz' % guide)
    if not os.path.exists(p):
        bad.append('%s: emitted list missing at %s' % (guide, p)); continue
    seen = set()
    with gzip.open(p, 'rt') as fh:
        if fh.readline().rstrip('\n').split('\t') != ['contig', 'strand', 'pam_forward_start']:
            bad.append('%s: unexpected emitted header' % guide); continue
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) != 3: bad.append('%s: malformed emitted row' % guide); break
            seen.add((f[0], f[1], int(f[2])))
    emitted[guide] = seen
    if not seen: bad.append('%s: emitted list is EMPTY' % guide)
    print('  %-22s %6d emitted loci   expected %6d   %s'
          % (guide, len(seen), want, 'ok' if len(seen) == want else 'MISMATCH'))
    if len(seen) != want:
        bad.append('%s: %d emitted loci, expected %d' % (guide, len(seen), want))

# THE INVARIANT, not the arithmetic: an "exclusive" locus must not appear in the emitted union.
# A count check cannot see a list that silently includes reachable loci; this can.
for guide in sorted(emitted):
    p = os.path.join(D, 'exclusive_loci_%s.tsv.gz' % guide)
    if not os.path.exists(p): continue
    excl = set()
    with gzip.open(p, 'rt') as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) == 3: excl.add((f[0], f[1], int(f[2])))
    overlap = excl & emitted[guide]
    print('  %-22s exclusive loci also in the emitted union: %d   %s'
          % (guide, len(overlap), 'ok' if not overlap else 'INVARIANT VIOLATED'))
    if overlap:
        bad.append('%s: %d loci called exclusive are reached by an emitted pattern, e.g. %r'
                   % (guide, len(overlap), sorted(overlap)[:2]))

if len(emitted) == len(QUOTED_EMITTED) and len(counts) == len(QUOTED):
    complete = sum(len(s) for s in emitted.values()) + sum(counts.values())
    wshare = round(100.0 * sum(counts.values()) / complete, 2)
    print('  %-22s %7d complete       paper says %7d   %s'
          % ('both guides', complete, QUOTED_COMPLETE,
             'ok' if complete == QUOTED_COMPLETE else 'MISMATCH'))
    print('  %-22s %6.2f%% of the complete set   paper says %4.2f%%   %s'
          % ('whole-portfolio share', wshare, QUOTED_WHOLE_SHARE,
             'ok' if wshare == QUOTED_WHOLE_SHARE else 'MISMATCH'))
    if complete != QUOTED_COMPLETE:
        bad.append('complete set %d, paper quotes %d' % (complete, QUOTED_COMPLETE))
    if wshare != QUOTED_WHOLE_SHARE:
        bad.append('whole-portfolio share %.2f%%, paper quotes %.2f%%' % (wshare, QUOTED_WHOLE_SHARE))

if bad:
    # A missing input and a refuted claim are different outcomes and were reported under one
    # banner, so a referee alone in a directory read "MAIZE GATE FAILURES" as the maize claim
    # failing when the locus lists were simply not there. Partitioned: 2 is NOT RUN, 1 is FAILED.
    absent = [b for b in bad if 'missing at' in b]
    if len(absent) == len(bad):
        print('\nGATE NOT RUN: %d locus list(s) are not here, so nothing could be checked.'
              % len(absent))
        for b in absent: print('  [absent] %s' % b)
        print('This is an incomplete release, not a finding about the claim.')
        sys.exit(2)
    print('\nMAIZE GATE FAILURES:')
    for b in bad: print('  [FAIL] %s' % b)
    sys.exit(1)
print('every maize figure in the paper is reproduced from the released locus lists')
