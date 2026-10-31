#!/usr/bin/env python3
"""The omitted maize patterns' own locus counts, recomputed from the deposited hit files.

WHY THIS EXISTS. maize_step_two.py prints the per-pattern counts and the EXCLUSIVE union, and
verify_maize_exclusive.py checks the exclusive figures against the deposited lists. Neither
recomputes two figures the appendix quotes: how many distinct loci the nine omitted patterns reach
at all, and the sum of the per-pattern exclusive counts for the second guide. Both were quoted
without a producer until 2026-10-24, and the first of them was wrong: 179,689 is the SUM over
patterns, while the distinct union is 179,304, the two differing by the 385 loci that more than one
omitted pattern reaches. The appendix paragraph makes exactly that union-versus-sum distinction for
the exclusive counts, and had not applied it to its own denominator.

WHAT IT READS. data/maize_omitted_hits/ holds the nine omitted patterns' raw scan hits, gzipped,
as scan_pattern.py wrote them on the cluster. data/maize_exclusive_loci/ holds the emitted-loci
union per guide, which is all the exclusive difference needs, so the 133 emitted per-pattern scans
do not have to be redistributed.

THE LOCUS KEY IS NOT RE-DERIVED HERE. `key_of` is spliced out of maize_step_two.py's source text
and its bytes are checked against the definition this script was written against before it is used.
The anchor is the part of this measurement most likely to be wrong, and a second hand-written copy
would be a second chance for the two to differ. Pattern length comes from the file name, as there.

CONTROL. Every deposited exclusive locus must appear in the omitted set it was differenced out of.
That compares this script's key against an artifact produced independently of it; if the PAM offset
or the pattern length were wrong, containment would break. A failure aborts before any count is
printed, because a count under a wrong key is not a measurement.

Exit 0 = computed and controlled. 1 = the control failed. 2 = NOT RUN, which is not a pass.
"""
import glob
import gzip
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HITS = os.path.join(HERE, 'data', 'maize_omitted_hits')
DEPOSIT = os.path.join(HERE, 'data', 'maize_exclusive_loci')
STEP2 = os.path.join(HERE, 'maize_step_two.py')

EXPECTED_KEY_OF = (
    "def key_of(contig, strand, start, L):\n"
    "    return (contig, strand, start + L if strand == '+' else start)\n")


def splice_key_of(src_path):
    src = open(src_path, encoding='utf-8').read()
    m = re.search(r"^def key_of\(.*?\n(?:[ \t].*\n|\n)*", src, re.M)
    if not m:
        return None, 'key_of not found in %s' % src_path
    text = m.group(0).rstrip() + '\n'
    if text != EXPECTED_KEY_OF:
        return None, ('key_of in %s is not the definition this script was written against.\n'
                      'found:\n%r\nexpected:\n%r' % (src_path, text, EXPECTED_KEY_OF))
    ns = {}
    exec(compile(text, src_path, 'exec'), ns)
    return ns['key_of'], None


def _open(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p, encoding='utf-8')


def read_keyed(path, key_of, L):
    out = set()
    with _open(path) as fh:
        head = fh.readline()
        if not head.startswith('contig\t'):
            raise ValueError('bad header in %s: %r' % (os.path.basename(path), head[:50]))
        for line in fh:
            f = line.rstrip('\n').split('\t')
            out.add(key_of(f[0], f[1], int(f[2]), L))
    return out


def read_deposit(path):
    keys = set()
    with _open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        if head != ['contig', 'strand', 'pam_forward_start']:
            raise ValueError('unexpected deposit header %r in %s' % (head, os.path.basename(path)))
        for line in fh:
            f = line.rstrip('\n').split('\t')
            keys.add((f[0], f[1], int(f[2])))
    return keys


def main():
    if not os.path.exists(STEP2):
        print('NOT RUN: %s absent, so the locus key cannot be spliced.' % STEP2); return 2
    key_of, err = splice_key_of(STEP2)
    if err:
        print('NOT RUN: ' + err); return 2
    print('key_of spliced from maize_step_two.py and byte-identical to the definition relied on here.')

    files = sorted(glob.glob(os.path.join(HITS, 'hits_*omitted*.tsv*')))
    if not files:
        print('NOT RUN: no omitted hit files in %s' % HITS); return 2

    per_guide = {}
    for p in files:
        n = os.path.basename(p)
        stem = n[len('hits_'):]
        for suf in ('.tsv.gz', '.tsv'):
            if stem.endswith(suf):
                stem = stem[:-len(suf)]
                break
        guide, pat = stem.split('_')[0], stem.split('_')[-1]
        try:
            s = read_keyed(p, key_of, len(pat))
        except ValueError as e:
            print('NOT RUN: %s' % e); return 2
        if not s:
            print('NOT RUN: %s is empty, so a count from it is not a measurement.' % n); return 2
        per_guide.setdefault(guide, {})[pat] = s
        print('  %-22s %-24s len %2d  %7d distinct loci' % (guide, pat, len(pat), len(s)))

    npat = sum(len(v) for v in per_guide.values())
    print('\npatterns loaded: %d  (the appendix says 9)' % npat)

    # CONTROL before any headline count
    for g in sorted(per_guide):
        dep = os.path.join(DEPOSIT, 'exclusive_loci_%s.tsv.gz' % g)
        if not os.path.exists(dep):
            print('NOT RUN: deposited %s absent, so the key cannot be controlled.'
                  % os.path.basename(dep)); return 2
        try:
            keys = read_deposit(dep)
        except ValueError as e:
            print('NOT RUN: %s' % e); return 2
        union = set().union(*per_guide[g].values())
        missing = keys - union
        print('  control %-22s deposited exclusive %6d, missing from the omitted set %d'
              % (g, len(keys), len(missing)))
        if missing:
            print('CONTROL FAILED: a deposited exclusive locus is absent from the omitted set it '
                  'was differenced out of, so the locus key is wrong and no count below may be '
                  'quoted.')
            return 1

    grand = set()
    for g in sorted(per_guide):
        u = set().union(*per_guide[g].values())
        grand |= u
        ssum = sum(len(s) for s in per_guide[g].values())
        print('  %-22s %d patterns  union %7d  sum %7d  (%d reached by more than one)'
              % (g, len(per_guide[g]), len(u), ssum, ssum - len(u)))
    gsum = sum(len(s) for v in per_guide.values() for s in v.values())
    print('\nOMITTED PATTERNS REACH %d DISTINCT LOCI' % len(grand))
    print('per-pattern counts sum to %d, so %d loci are reached by more than one omitted pattern'
          % (gsum, gsum - len(grand)))

    print('\nper-pattern EXCLUSIVE counts, against the deposited emitted-loci union:')
    for g in sorted(per_guide):
        ep = os.path.join(DEPOSIT, 'emitted_loci_%s.tsv.gz' % g)
        if not os.path.exists(ep):
            print('  %s: emitted list absent, skipped' % g); continue
        try:
            E = read_deposit(ep)
        except ValueError as e:
            print('NOT RUN: %s' % e); return 2
        only_u, only_sum = set(), 0
        for pat in sorted(per_guide[g]):
            o = per_guide[g][pat] - E
            only_sum += len(o)
            only_u |= o
            print('    %-24s exclusive %6d' % (pat, len(o)))
        print('  %-22s emitted %7d   exclusive union %6d   sum %6d   (%d shared)'
              % (g, len(E), len(only_u), only_sum, only_sum - len(only_u)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
