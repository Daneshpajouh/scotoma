#!/usr/bin/env python3
"""The paper's "up to two bases" census row, recomputed from the deposited cluster measurement.

PROVENANCE, stated first because it is the point. The counts come from job 212667 on Fir, run by
`taxonomy_pair.py` over the panel's bulge census outputs `bulge_census_full` and `bulge2_census`,
and the ungapped set `cof_h2h/out_ngg_mm6.txt`, all under that job's scratch directory on Fir.
Those raw outputs run to roughly 2 GB and are **not redistributed**; the summary they produced,
`data/TAXONOMY_PAIR_mm4.json`, is. That run carried three gates and printed that all three passed:
the ungapped denominator reproduced 34,652, and the size-1 and size-2 strict additions reproduced
1,212,681 and 2,712,081, both described in the original producer as independently measured.

WHAT THIS SCRIPT RECOMPUTES, as opposed to reads. The distinction matters because the manifest's own
rule is that a value a script merely hardcodes does not earn a row.

  RECOMPUTED from the deposited counts:
    * the strict total,    1,212,681 + 2,712,081                  = 3,924,762
    * the strict union,    that total - the strict intersection   = 3,780,370
    * all four panel factors, (base + added) / base
    * the merged-specific intersection implied by the merged union
  READ from the deposited file, and NOT recomputable from anything smaller:
    * the merged union itself, 3,713,396

That last point is the finding this deposit exists to support. The merged union is **not**
inclusion-exclusion on the strict intersection: that arithmetic gives 3,713,381, fifteen short. The
merged sets have their own intersection, 144,377, and recovering it requires the locus sets rather
than any summary, which is why the figure needed a direct measurement and why it had no producer
until now.

Exit 0 = every paper value checks. 1 = a disagreement. 2 = NOT RUN.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'data', 'TAXONOMY_PAIR_mm4.json')

PAPER = {
    'ungapped base':                      34652,
    'size-1 added, strict':               1212681,
    'size-1 added, merged':               1175323,
    'size-2 added, strict':               2712081,
    'strict intersection':                144392,
    'strict total of the two classes':    3924762,
    'up to two bases, union, strict':     3780370,
    'up to two bases, union, merged':     3713396,
}
PAPER_FACTORS = {'size-1 merged': 34.9, 'size-2 merged': 78.4,
                 'union strict': 110.1, 'union merged': 108.2}


def main():
    if not os.path.exists(SRC):
        print('NOT RUN: %s absent.' % SRC); return 2
    d = json.load(open(SRC, encoding='utf-8'))
    base = d['ungapped']
    s1, s2 = d['classes']['size1'], d['classes']['size2']
    inter = d['overlaps']['size1_size2']['intersection']
    uni = d['union_size1_size2']

    got = {
        'ungapped base':                   base,
        'size-1 added, strict':            s1['add_strict'],
        'size-1 added, merged':            s1['add_merged'],
        'size-2 added, strict':            s2['add_strict'],
        'strict intersection':             inter,
        # recomputed, not read
        'strict total of the two classes': s1['add_strict'] + s2['add_strict'],
        'up to two bases, union, strict':  s1['add_strict'] + s2['add_strict'] - inter,
        # read: see the docstring
        'up to two bases, union, merged':  uni['add_merged'],
    }

    bad = 0
    print('%-34s %12s %18s' % ('quantity', 'measured', 'paper'))
    for k in PAPER:
        v, want = got[k], PAPER[k]
        ok = (v == want)
        if not ok:
            bad += 1
        print('%-34s %12d   paper %12d   %s' % (k, v, want, 'ok' if ok else 'MISMATCH'))

    # the deposited strict union must equal the inclusion-exclusion of the deposited class counts,
    # or the summary is internally inconsistent and nothing from it may be quoted
    if uni['add_strict'] != got['up to two bases, union, strict']:
        print('GATE FAILED: the deposited strict union %d disagrees with inclusion-exclusion on '
              'the deposited class counts, %d. The summary is internally inconsistent.'
              % (uni['add_strict'], got['up to two bases, union, strict']))
        return 1

    print()
    for name, added, want in (('size-1 merged', s1['add_merged'], PAPER_FACTORS['size-1 merged']),
                              ('size-2 merged', s2['add_merged'], PAPER_FACTORS['size-2 merged']),
                              ('union strict', uni['add_strict'], PAPER_FACTORS['union strict']),
                              ('union merged', uni['add_merged'], PAPER_FACTORS['union merged'])):
        f = (base + added) / base
        ok = (round(f, 1) == want)
        if not ok:
            bad += 1
        print('panel factor, %-14s %9.4fx   paper %5.1fx   %s'
              % (name, f, want, 'ok' if ok else 'MISMATCH'))

    ie = s1['add_merged'] + s2['add_merged'] - inter
    implied = s1['add_merged'] + s2['add_merged'] - uni['add_merged']
    print()
    print('inclusion-exclusion on the STRICT intersection would give %d, which is %+d from the '
          'measured merged union' % (ie, uni['add_merged'] - ie))
    print('merged-specific intersection implied by the measurement: %d' % implied)

    print()
    if bad:
        print('FAILED: %d quantity/quantities disagree with the paper.' % bad); return 1
    print('every figure the paper prints from this census agrees with the deposited measurement')
    return 0


if __name__ == '__main__':
    sys.exit(main())
