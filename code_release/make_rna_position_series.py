#!/usr/bin/env python3
"""Produce the twenty per-position RNA-bulge activity counts the appendix prints.

The appendix sentence at appendix_moved.tex:1764 lists measured active sites on the RNA side by
deleted spacer position, over all twenty positions. Nothing in the release produced that series:
make_sites129.py measures the same rows but keeps only the terminal position, because that is the
one it needs, and the other nineteen counts were never captured.

This script is make_sites129.py's selection with the terminal filter removed and a histogram in
its place. Every other rule is taken from that script unchanged, so the terminal cell of this
series is the same quantity that script's requirement 1 reports:

  * a row counts as active when its read count is above zero;
  * it carries a single RNA bulge when it has exactly one gap, in the target and not the guide;
  * mismatches are capped at four;
  * the deleted guide base's position is counted in real, ungapped guide bases.

Two checks run on the result, both against figures the release already holds, so the series cannot
drift from them silently:

  * the twenty counts must sum to rna_bulge_actives, 17,267, from verify_missed_actives.py;
  * position 1 must be 0 and the terminal position must be 583, that script's rna_at_pos_0 and
    rna_at_terminal.

Usage, the same stdin convention as make_sites129.py:

    cat <the CHANGE-seq deposit CSV> | python3 make_rna_position_series.py

Exits non-zero if any check fails.
"""
import csv
import sys

N = 20  # spacer length; the deposit's guides are 20 nt plus a 3 nt PAM

# From code_release/verify_missed_actives.py, DEPOSIT_DERIVED.
EXPECTED_TOTAL = 17267
EXPECTED_AT_POS_1 = 0
EXPECTED_AT_TERMINAL = 583


def main():
    counts = [0] * N
    rows = active = bulge = single_rna = 0
    csv.field_size_limit(1 << 24)
    for r in csv.DictReader(sys.stdin):
        rows += 1
        try:
            reads = float((r.get('reads') or '0').strip())
            nb = int(float((r.get('Align.#Bulges') or '0').strip()))
            nm = int(float((r.get('Align.#Mismatches') or '0').strip()))
        except (TypeError, ValueError):
            continue
        if reads <= 0:
            continue
        active += 1
        if nb < 1:
            continue
        bulge += 1

        g = (r.get('Align.sgRNA') or '').strip().upper()
        t = (r.get('Align.off-target') or '').strip().upper()
        if nb != 1 or t.count('-') != 1 or g.count('-') != 0 or nm > 4:
            continue
        single_rna += 1

        pos = len(g[:t.index('-')].replace('-', ''))
        if 0 <= pos < N:
            counts[pos] += 1

    if rows == 0:
        sys.exit('GATE FAILED: no rows read from stdin. This is an invocation error, not a '
                 'finding. Pipe the CHANGE-seq deposit in.')

    print('deposit rows %d -> measured active %d -> bulge-bearing %d -> single RNA bulge '
          'within 4 mm %d' % (rows, active, bulge, single_rna))
    print('by deleted spacer position, g1 first:')
    print('  ' + ', '.join('%d' % c for c in counts))

    bad = []
    total = sum(counts)
    if total != EXPECTED_TOTAL:
        bad.append('the series sums to %d, not the released rna_bulge_actives %d'
                   % (total, EXPECTED_TOTAL))
    if counts[0] != EXPECTED_AT_POS_1:
        bad.append('position 1 is %d, not the released rna_at_pos_0 %d'
                   % (counts[0], EXPECTED_AT_POS_1))
    if counts[N - 1] != EXPECTED_AT_TERMINAL:
        bad.append('the terminal position is %d, not the released rna_at_terminal %d'
                   % (counts[N - 1], EXPECTED_AT_TERMINAL))
    if bad:
        print('\nRNA-POSITION-SERIES GATE FAILURES')
        for b in bad:
            print('  [FAIL] %s' % b)
        return 1
    print('  [ok  ] sums to the released total %d, with %d at g1 and %d at the terminal position'
          % (total, counts[0], counts[N - 1]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
