#!/usr/bin/env python3
"""Emit the EMITTED and OMITTED pattern sets for a spacer, per the paper's Proposition 1.

This is the input to the maize question, so getting the convention right is the whole job. The
2016 wrapper's loops are the authority:

  DNA, size b: the scheme has n+1 slots and exactly one is degenerate (it reproduces the no-bulge
    pattern), leaving n distinct patterns. The wrapper's `range(1, len(target))` covers slots
    1..n-1, so it emits n-1 of the n and OMITS one. For a trailing PAM slot 0 is the degenerate
    one, so the omitted distinct pattern is slot n.
  RNA, size b: the full set is the b-burst deletions at starts i = 1..n-b+1 (1-indexed), and the
    wrapper's `range(1, len(target)-b)` covers the interior. The i=1 word is omitted unless
    g_1 = g_{b+1}, and the i=n-b+1 word is omitted unless g_n = g_{n-b}, which is exactly the
    deficit law 1 + 1[g_1 != g_{b+1}] + 1[g_n != g_{n-b}].

THE CONTROL IS A PUBLISHED FIGURE, not a self-check: on the 110-guide panel the terminal indicator
must split 73 exposed against 37 protected, which is the paper's own number, and the per-guide
deficits must equal the deficit law for every guide and every b. --control runs both.
"""
import argparse, json, os, sys

N_SPACER = 20


def rna_words(g, b):
    """The distinct b-burst-deletion words of g, with the start index that produced each."""
    seen = {}
    for i in range(len(g) - b + 1):            # 0-indexed start
        w = g[:i] + g[i + b:]
        seen.setdefault(w, i)                  # first start wins; collisions are the law's point
    return seen


def sets_for(g, b):
    """(emitted, omitted) pattern lists for one spacer at bulge size b, as strings over ACGT/N."""
    n = len(g)
    # ---- RNA arm: deletions from the guide
    words = rna_words(g, b)
    interior = {g[:i] + g[i + b:] for i in range(1, n - b)}     # the wrapper's own range
    rna_emitted = sorted(interior)
    rna_omitted = sorted(w for w in words if w not in interior)
    # ---- DNA arm: an inserted contiguous block, one wildcard per inserted base
    # SLOT 0 IS DEGENERATE and is not a distinct pattern of this class. In the wrapper's padded
    # encoding N^(bmax-b) N^b g N^L collapses to N^bmax g N^L, the no-bulge pattern, so of the n+1
    # slots exactly n are distinct. The control caught this: counting slot 0 as an omission gave a
    # deficit of 4 against the law's 3 on every guide, because it charged both slot 0 and slot n.
    dna = {}
    for slot in range(1, n + 1):               # slot 0 excluded: degenerate under padding
        dna.setdefault(g[:slot] + 'N' * b + g[slot:], slot)
    dna_emitted = sorted(w for w, s in dna.items() if s <= n - 1)
    dna_omitted = sorted(w for w, s in dna.items() if s == n)
    return dict(rna_emitted=rna_emitted, rna_omitted=rna_omitted,
                dna_emitted=dna_emitted, dna_omitted=dna_omitted,
                rna_distinct=len(words), dna_distinct=len(dna))


def deficit_law(g, b):
    """1 + 1[g_1 != g_{b+1}] + 1[g_n != g_{n-b}], the paper's closed form."""
    return 1 + (1 if g[0] != g[b] else 0) + (1 if g[-1] != g[-1 - b] else 0)


def control(panel_json):
    """Two published checks. Exits non-zero on either."""
    rows = json.load(open(panel_json))['rows']
    spacers = [r['spacer'] for r in rows]
    if len(spacers) != 110:
        sys.exit('GATE FAILED: expected the 110-guide panel, got %d spacers' % len(spacers))
    # 1. the deficit law must equal the generated omission count, every guide, every b
    bad = []
    for g in spacers:
        for b in (1, 2, 3):
            s = sets_for(g, b)
            got = len(s['rna_omitted']) + len(s['dna_omitted'])
            want = deficit_law(g, b)
            if got != want:
                bad.append((g, b, got, want))
    print('  [%s] deficit law vs generated omissions: %d disagreement(s) over %d guides x 3 sizes'
          % ('ok  ' if not bad else 'FAIL', len(bad), len(spacers)))
    if bad:
        for t in bad[:4]:
            print('        %s b=%d generated %d, law %d' % t)
    # 2. the terminal indicator must split the panel 73 / 37, which is the paper's own figure
    exposed = sum(1 for g in spacers if g[-1] != g[-2])
    protected = len(spacers) - exposed
    ok2 = (exposed, protected) == (73, 37)
    print('  [%s] terminal indicator splits the panel %d exposed / %d protected (paper says 73/37)'
          % ('ok  ' if ok2 else 'FAIL', exposed, protected))
    if bad or not ok2:
        print('GATE FAILED: the generator does not reproduce a published figure, so its pattern '
              'sets are not trustworthy. Fix this before any scan uses them.')
        return 2
    print('CLEAN: the generator reproduces the deficit law and the published 73/37 split.')
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--spacer', action='append', default=[],
                    help='repeatable; 20 nt over ACGT')
    ap.add_argument('--bulge-sizes', default='1,2')
    ap.add_argument('--out-json')
    ap.add_argument('--control', metavar='PANEL_JSON',
                    help='run the two published checks against the 110-guide panel')
    a = ap.parse_args(argv)
    if a.control:
        return control(a.control)
    if not a.spacer:
        ap.error('give --spacer at least once, or --control')
    sizes = [int(x) for x in a.bulge_sizes.split(',')]
    out = {}
    for g in a.spacer:
        g = g.strip().upper()
        if len(g) != N_SPACER or set(g) - set('ACGT'):
            sys.exit('GATE FAILED: %r is not 20 nt over ACGT' % g)
        out[g] = {}
        tot_om = 0
        for b in sizes:
            s = sets_for(g, b)
            out[g]['b%d' % b] = s
            nom = len(s['rna_omitted']) + len(s['dna_omitted'])
            tot_om += nom
            print('%s  b=%d  RNA %2d distinct, %2d emitted, %d omitted | DNA %2d distinct, %2d '
                  'emitted, %d omitted | deficit law %d'
                  % (g, b, s['rna_distinct'], len(s['rna_emitted']), len(s['rna_omitted']),
                     s['dna_distinct'], len(s['dna_emitted']), len(s['dna_omitted']),
                     deficit_law(g, b)))
        out[g]['omitted_total_across_sizes'] = tot_om
        print('%s  patterns never emitted across sizes %s: %d' % (g, sizes, tot_om))
    if a.out_json:
        json.dump(out, open(a.out_json, 'w'), indent=1)
        print('wrote %s' % a.out_json)
    return 0


if __name__ == '__main__':
    sys.exit(main())
