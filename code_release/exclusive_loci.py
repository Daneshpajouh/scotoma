#!/usr/bin/env python3
"""Step two of the maize question: which candidate loci NO emitted pattern reaches.

Step one scans the genome for the patterns the standard position set OMITS (scan_pattern.py over
make_pattern_sets.py's omitted lists). Those hits are candidates, not the claim. A candidate is the
claim only if NO pattern the standard set DOES emit reaches the same locus within the same budget,
because a locus both sets reach is not lost by the omission.

THE DEFINITION, matching what the paper already does for the 129: a pattern reaches a locus if it
matches some window overlapping that locus within the budget. Patterns here have lengths 18 to 22,
so the windows differ by pattern, and every offset inside the flank is tried rather than assuming
the omitted pattern's own alignment.

THE CONTROL IS A NEGATIVE ONE and --control runs it: the guide's own on-target site must NOT be
reported exclusive, because the emitted set plainly reaches it. A step-two that reports everything
exclusive would look like a large result and be worthless, which is the failure this guards.
"""
import argparse, json, os, sys

COMP = str.maketrans('ACGTN', 'TGCAN')


def revcomp(s):
    return s.translate(COMP)[::-1]


def matches(pattern, window, budget):
    """Does `pattern` match some offset of `window` within `budget`? N in the pattern is free."""
    L = len(pattern)
    if len(window) < L:
        return False
    for off in range(len(window) - L + 1):
        seg = window[off:off + L]
        mm = 0
        for p, w in zip(pattern, seg):
            if p == 'N':
                continue
            if p != w:
                mm += 1
                if mm > budget:
                    break
        else:
            if mm <= budget:
                return True
    return False


def reached_by_any(patterns, window, budget):
    """The first emitted pattern that reaches this window, or None."""
    rc = revcomp(window)
    for p in patterns:
        if matches(p, window, budget) or matches(p, rc, budget):
            return p
    return None


def classify(candidates, emitted, budget):
    """candidates: list of dicts with 'window' (and whatever else). Returns (exclusive, reached)."""
    exclusive, reached = [], []
    for c in candidates:
        hit = reached_by_any(emitted, c['window'], budget)
        if hit is None:
            exclusive.append(c)
        else:
            reached.append(dict(c, reached_by=hit))
    return exclusive, reached


def control(patterns_json, budget=5, flank=26):
    """NEGATIVE control: the guide's own site must not come out exclusive."""
    d = json.load(open(patterns_json))
    fails = []
    for g, v in d.items():
        emitted = []
        for b in ('b1', 'b2'):
            if b in v:
                emitted += v[b]['rna_emitted'] + v[b]['dna_emitted']
        if not emitted:
            sys.exit('GATE FAILED: no emitted patterns for %s, so this control would compare '
                     'against an empty set and pass trivially.' % g)
        # the on-target window: the spacer itself with a real PAM and random-free flanks
        window = 'A' * flank + g + 'AGG' + 'A' * flank
        hit = reached_by_any(emitted, window, budget)
        ok = hit is not None
        print('  [%s] %s own site reached by an emitted pattern: %s'
              % ('ok  ' if ok else 'FAIL', g, hit if hit else 'NONE'))
        if not ok:
            fails.append(g)
        # and a POSITIVE direction: a window of pure filler must NOT be reached
        junk = 'A' * (flank * 2 + 23)
        hit2 = reached_by_any(emitted, junk, budget)
        ok2 = hit2 is None
        print('  [%s] %s pure-filler window NOT reached (budget %d)'
              % ('ok  ' if ok2 else 'FAIL', g, budget))
        if not ok2:
            fails.append(g + ':filler')
    if fails:
        print('GATE FAILED: %s. A step-two that cannot tell a reached locus from an unreached one '
              'would report a large and worthless result.' % ', '.join(fails))
        return 2
    print('CLEAN: own sites are reached, filler is not, so the classifier discriminates.')
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--patterns', required=True, help='make_pattern_sets.py --out-json')
    ap.add_argument('--candidates', help='JSON list of {contig,start,strand,window} from step one')
    ap.add_argument('--guide', help='which spacer the candidates belong to')
    ap.add_argument('--budget', type=int, default=5)
    ap.add_argument('--out-json')
    ap.add_argument('--control', action='store_true')
    a = ap.parse_args(argv)
    if a.control:
        return control(a.patterns, a.budget)
    if not a.candidates or not a.guide:
        ap.error('--candidates and --guide are required unless --control')
    d = json.load(open(a.patterns))
    if a.guide not in d:
        sys.exit('GATE FAILED: %s is not in %s' % (a.guide, a.patterns))
    emitted = []
    for b in ('b1', 'b2'):
        if b in d[a.guide]:
            emitted += d[a.guide][b]['rna_emitted'] + d[a.guide][b]['dna_emitted']
    cands = json.load(open(a.candidates))
    # LOCI, not windows: both strands can report one genomic position
    seen, uniq = set(), []
    for c in cands:
        k = (c.get('contig'), c.get('start'))
        if k in seen:
            continue
        seen.add(k)
        uniq.append(c)
    exclusive, reached = classify(uniq, emitted, a.budget)
    print('guide %s  emitted patterns %d  candidate rows %d  distinct loci %d'
          % (a.guide, len(emitted), len(cands), len(uniq)))
    print('  reached by an emitted pattern : %d' % len(reached))
    print('  EXCLUSIVE to the complete set : %d' % len(exclusive))
    if a.out_json:
        json.dump({'guide': a.guide, 'budget': a.budget, 'emitted_patterns': len(emitted),
                   'candidate_rows': len(cands), 'distinct_loci': len(uniq),
                   'exclusive': exclusive, 'reached': reached},
                  open(a.out_json, 'w'), indent=1)
        print('wrote %s' % a.out_json)
    return 0


if __name__ == '__main__':
    sys.exit(main())
