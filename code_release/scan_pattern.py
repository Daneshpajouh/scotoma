#!/usr/bin/env python3
"""Scan a reference for ONE bulge pattern: any length, and a contiguous wildcard block.

WHY THIS EXISTS. enumerate_sites.py refuses anything that is not a 20-mer over ACGT (its line
"the spacer must be 20 nt over ACGT"), and a bulge pattern is not one: an RNA bulge of size b
gives an (n-b)-mer and a DNA bulge of size b gives an (n+b)-mer carrying one wildcard. So the
claim "the standard position set omits patterns, and here are the loci only the complete set
reaches" cannot be checked with the tools that ship. This scans one pattern.

IT REUSES enumerate_sites.py's OWN encoding, reverse complement, contig reader and PAM test by
importing them, rather than reimplementing any of it, so the two agree by construction on
everything except the length and wildcard handling this file adds.

THE CONTROL IS THE POINT, and --self-test runs it: on a 20-mer with no wildcard this file and
enumerate_sites.py must return the identical window set. A bulge scanner that has not been held
against the shipped scanner on a case the shipped scanner can do is how a confident wrong number
gets published.

Usage:
    python3 scan_pattern.py --self-test --ref-dir DIR            # the control, run this first
    python3 scan_pattern.py --ref-dir DIR --pattern ACGT...N --max-mismatches 5 --out-tsv F
"""
import argparse, os, sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import enumerate_sites as ES          # noqa: E402  the shipped scanner, reused not reimplemented

WILD = 4                              # encode() maps anything outside ACGT to >= 4


def encode_pattern(p):
    """Encode a pattern over ACGT plus N. N becomes WILD and is free at that position."""
    p = p.strip().upper()
    out = np.empty(len(p), dtype=np.int8)
    for i, ch in enumerate(p):
        if ch == 'N':
            out[i] = WILD
        elif ch in 'ACGT':
            out[i] = 'ACGT'.index(ch)
        else:
            sys.exit('GATE FAILED: pattern character %r is neither ACGT nor N' % ch)
    # 2026-10-22, corrected the same day it was written: the original guard allowed exactly ONE
    # wildcard, "because one bulge contributes exactly one". That is true only at b=1. A DNA bulge
    # of size b inserts b CONTIGUOUS bases, so its pattern carries b contiguous wildcards, and
    # Young et al. searched with up to 2 bulges. The guard would have silently blocked half the
    # patterns the claim is about. What actually characterises a bulge pattern is CONTIGUITY, not
    # a count of one, so that is what is enforced.
    w = np.nonzero(out == WILD)[0]
    if w.size and int(w[-1] - w[0]) != int(w.size) - 1:
        sys.exit('GATE FAILED: the wildcards are not contiguous (positions %s). A bulge inserts a '
                 'contiguous block, so a scattered wildcard set is not a bulge pattern and this '
                 'scan would not mean what the claim says it means.' % list(map(int, w)))
    return out


def scan_contig(enc, pat, max_mm, CH=8_000_000):
    """(strand, 1-based forward start, window, mismatches) for one encoded contig.

    Mirrors ES.scan_encoded's geometry exactly: window is the protospacer plus three PAM bases,
    the PAM test is NGG or NAG on the last two, and the start is the leftmost forward base of the
    window on both strands. The two differences from ES are deliberate and are the reason this
    file exists: L is the pattern's own length rather than 20, and wildcard positions are free.
    """
    L = int(pat.size)
    free = (pat == WILD)
    n_free = int(free.sum())
    out = []
    if enc.size < L + 3:
        return out
    # 2026-10-22: CHUNKED. The unchunked version materialised (proto != pat) over the whole contig,
    # which is 5.9 GB of booleans on the 312 MB maize chromosome and made the scan memory-bound
    # rather than compute-bound. Chunks of CH bases cap that at a few hundred MB. Consecutive
    # chunks OVERLAP by L+2 so no window straddling a boundary is lost, which is exactly where an
    # off-by-one would hide; --self-test compares against the unchunked enumerate_sites.py and
    # would catch one.
    for strand, seq in (("+", enc), ("-", ES.revcomp(enc))):
        n = int(seq.size)
        base = 0
        while base + L + 3 <= n:
            stop = min(n, base + CH + L + 2)
            sub = seq[base:stop]
            win = np.lib.stride_tricks.sliding_window_view(sub, L + 3)
            proto = win[:, :L]
            acgt = (win < 4).all(axis=1)
            p2, p3 = win[:, L + 1], win[:, L + 2]
            pam_ok = ((p2 == 2) | (p2 == 0)) & (p3 == 2)      # NGG or NAG, as ES does
            diff = (proto != pat)
            if n_free:
                diff = diff & ~free                           # a wildcard never mismatches
            mm = diff.sum(axis=1)
            hit = acgt & pam_ok & (mm <= max_mm)
            for i in np.nonzero(hit)[0]:
                i = int(i)
                g = base + i                                  # index within this strand
                w = ''.join(ES._DEC[win[i]])
                start = g + 1 if strand == "+" else int(enc.size) - (g + L + 3) + 1
                out.append((strand, start, w, int(mm[i])))
            if stop >= n:
                break
            base += CH
    return out


def self_test(ref_dir, spacer=None, max_mm=3, chunk=8_000_000, max_bases=40_000_000):
    """On a 20-mer with no wildcard this file must equal enumerate_sites.py, window for window.

    IT RUNS ON A BOUNDED SUBSET OF THE REFERENCE, SMALLEST CONTIGS FIRST. The scanner it
    compares against, enumerate_sites.py, is deliberately NOT chunked, so on a 312 Mb chromosome it
    materialises about 6 GB of booleans and gets OOM-killed. The control's job is to establish that
    the two agree, which needs real hits and not the whole genome; max_bases caps it. The existing
    zero-window guard still refuses a comparison with nothing in it, so bounding the input cannot
    make this vacuous.

    THE CHUNK SIZE IS SETTABLE SO THIS CAN EXERCISE THE BOUNDARY. With the default chunk a small
    test reference fits in ONE chunk, so the overlap logic is never run and the control passes
    without testing the thing most likely to be wrong. Pass --chunk small enough to force several
    chunks over whatever reference is used."""
    spacer = spacer or ES.SPACER_EXACEL
    pat = encode_pattern(spacer)
    mine, theirs = [], []
    items = sorted(ES.iter_contigs(ref_dir), key=lambda t: len(t[1]))
    used = 0
    kept = []
    for name, raw in items:
        if used + len(raw) > max_bases and kept:
            break
        kept.append((name, raw)); used += len(raw)
    print('         control input: %d of %d contigs, %d bases (cap %d)'
          % (len(kept), len(items), used, max_bases))
    for name, raw in kept:
        enc = ES.encode(raw)
        mine += [(name,) + r for r in scan_contig(enc, pat, max_mm, CH=chunk)]
        theirs += [(name,) + r for r in ES.scan_encoded(enc, ES.encode(spacer.encode()), max_mm)]
    same = sorted(mine) == sorted(theirs)
    print('CONTROL: %d windows from this file, %d from enumerate_sites.py at <=%d mismatches'
          % (len(mine), len(theirs), max_mm))
    print('         chunk size %d, so the overlap logic %s exercised'
          % (chunk, 'WAS' if chunk < sum(1 for _ in []) + max(
              (len(r) for _, r in ES.iter_contigs(ref_dir)), default=0) else 'was NOT'))
    # 2026-10-22: the first run of this control returned 0 windows on BOTH sides and printed
    # IDENTICAL, which compares nothing. That is the same vacuous-pass defect two other scripts in
    # this folder were repaired for the same day, and it appeared here within minutes of my writing
    # the guard for them. A control that can pass without comparing anything is not a control.
    if not mine and not theirs:
        sys.exit('GATE FAILED: both scanners returned ZERO windows, so this comparison establishes\n'
                 '  nothing. Give --pattern a spacer that actually occurs in --ref-dir within the\n'
                 '  mismatch budget, or point --ref-dir at a reference that contains one.')
    print('IDENTICAL WINDOW SETS: %s' % same)
    if not same:
        a, b = set(map(tuple, mine)), set(map(tuple, theirs))
        print('  only here  (%d): %s' % (len(a - b), sorted(a - b)[:3]))
        print('  only there (%d): %s' % (len(b - a), sorted(b - a)[:3]))
        sys.exit('GATE FAILED: this scanner disagrees with the shipped one on a case the shipped '
                 'one can do, so its bulge results are not trustworthy. Fix this before using it.')
    print('The two agree, so the only behaviour this file adds is length and wildcard handling.')
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--ref-dir', required=True)
    ap.add_argument('--pattern', help='over ACGT plus a contiguous block of N')
    ap.add_argument('--max-mismatches', type=int, default=5)
    ap.add_argument('--out-tsv')
    ap.add_argument('--chunk', type=int, default=8_000_000,
                    help='bases per scan chunk; set small in --self-test to exercise the overlap')
    ap.add_argument('--control-bases', type=int, default=40_000_000,
                    help='cap on bases used by --self-test; the reference scanner is unchunked')
    ap.add_argument('--self-test', action='store_true',
                    help='compare against enumerate_sites.py on a 20-mer; run this first')
    a = ap.parse_args(argv)
    if a.self_test:
        # The 3-mismatch cap here was arbitrary and it interacts badly with a bounded input: on
        # 40 Mb of maize a 20-mer has an expected well under one hit at <=3, which would trip the
        # zero-window guard and abort the job for the right reason at the wrong time. The guard
        # stays; the budget is now the caller's.
        return self_test(a.ref_dir, a.pattern, a.max_mismatches, a.chunk, a.control_bases)
    if not a.pattern:
        ap.error('--pattern is required unless --self-test')
    pat = encode_pattern(a.pattern)
    rows, contigs = [], 0
    for name, raw in ES.iter_contigs(a.ref_dir):
        contigs += 1
        for r in scan_contig(ES.encode(raw), pat, a.max_mismatches, CH=a.chunk):
            rows.append((os.path.splitext(name)[0],) + r)
        print('[scan] %-14s loci so far: %8d' % (name, len(rows)), flush=True)
    print('pattern %s  length %d  wildcards %d  contigs %d  loci <=%d mm: %d'
          % (a.pattern, len(a.pattern), a.pattern.upper().count('N'), contigs,
             a.max_mismatches, len(rows)))
    if a.out_tsv:
        with open(a.out_tsv, 'w') as fh:
            fh.write('contig\tstrand\tstart_1based_forward\twindow\tmismatches\n')
            for r in rows:
                fh.write('\t'.join(map(str, r)) + '\n')
        print('wrote %s' % a.out_tsv)
    return 0


if __name__ == '__main__':
    sys.exit(main())
