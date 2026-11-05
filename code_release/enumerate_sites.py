#!/usr/bin/env python3
"""Enumerate the candidate off-target windows of a spacer in a reference genome.

A candidate window is a stretch of 20 nt (the protospacer) followed by a 3-nt PAM, read on
either strand, such that

  * the PAM is NGG or NAG,
  * the protospacer differs from the spacer in at most MAX_MISMATCHES positions
    (substitutions only, no bulges),
  * the window contains only A, C, G, T (soft-masked lowercase counts as the same base).

The unit is the window, that is a (strand, position) pair: two windows at different loci
that carry the same 23-mer are two sites, and one 23-mer may occupy several. The window with
zero mismatches is the on-target site, which is cleaved by design and is not a candidate
off-target; callers that count off-target sites must exclude it (see summarise()).

Usage:
  python3 enumerate_sites.py --ref-dir DIR --out-json summary.json [--out-sequences seqs.txt]
                             [--out-windows windows.tsv]
where DIR holds one file per contig (.fa, .fasta or .seq; an optional FASTA header line is
skipped and line breaks are ignored). The module has no dependency beyond numpy.
"""
import argparse
import json
import os
import sys

import numpy as np

SPACER_EXACEL = "CTAACAGTTGCTTTTATCAC"
PROTOSPACER_LEN = 20
MAX_MISMATCHES = 6

_ENC = np.full(256, 4, dtype=np.uint8)           # 0..3 = A, C, G, T; 4 = anything else
for _i, _b in enumerate(b"ACGTacgt"):
    _ENC[_b] = _i % 4
_COMP = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
_DEC = np.array(list("ACGTN"))


def encode(seq):
    """Encode a bytes sequence as 0..4 (A, C, G, T, other)."""
    return _ENC[np.frombuffer(seq, dtype=np.uint8)]


def revcomp(enc):
    return _COMP[enc[::-1]]


def scan_encoded(enc, spacer_enc, max_mismatches=MAX_MISMATCHES):
    """Return a list of (strand, start, window23, mismatches) for one encoded contig.

    strand is '+' or '-'; start is the 1-based position on the forward strand of the first
    (leftmost) base covered by the 23-base window, so the window covers start .. start + 22
    on the forward strand for both strands. The window23 string is read 5' to 3' on its own
    strand (the reverse complement of the forward bases when strand is '-')."""
    L = len(spacer_enc)
    n_contig = int(enc.size)
    out = []
    if enc.size < L + 3:
        return out
    for strand, seq in (("+", enc), ("-", revcomp(enc))):
        win = np.lib.stride_tricks.sliding_window_view(seq, L + 3)
        proto = win[:, :L]
        acgt = (win < 4).all(axis=1)
        p2, p3 = win[:, L + 1], win[:, L + 2]
        pam_ok = ((p2 == 2) | (p2 == 0)) & (p3 == 2)        # NGG or NAG
        mm = (proto != spacer_enc).sum(axis=1)
        for i in np.where(acgt & pam_ok & (mm <= max_mismatches))[0]:
            fwd_start = int(i) + 1 if strand == "+" else n_contig - int(i) - (L + 3) + 1
            out.append((strand, fwd_start, "".join(_DEC[win[i]]), int(mm[i])))
    return out


def iter_contigs(ref_dir):
    for fn in sorted(os.listdir(ref_dir)):
        if fn.startswith(".") or not fn.endswith((".seq", ".fa", ".fasta")):
            continue
        with open(os.path.join(ref_dir, fn), "rb") as fh:
            raw = fh.read()
        if raw.startswith(b">"):
            nl = raw.find(b"\n")
            raw = raw[nl + 1:] if nl >= 0 else b""
        yield fn, raw.replace(b"\n", b"").replace(b"\r", b"")


def summarise(records, spacer, max_mismatches=MAX_MISMATCHES):
    """Counts from scan records. off_target_windows excludes the zero-mismatch window(s)."""
    hist = {}
    seqs, seqs_ngg = set(), set()
    for _strand, _start, w23, mm in records:
        hist[mm] = hist.get(mm, 0) + 1
        seqs.add(w23)
        if w23[PROTOSPACER_LEN + 1:] == "GG":
            seqs_ngg.add(w23)
    total = sum(hist.values())
    return {
        "spacer": spacer,
        "max_mismatches": max_mismatches,
        "mismatch_histogram_windows": {str(k): hist[k] for k in sorted(hist)},
        "windows_total": total,
        "on_target_windows": hist.get(0, 0),
        "off_target_windows": total - hist.get(0, 0),
        "distinct_23mers": len(seqs),
        "distinct_23mers_ngg_only": len(seqs_ngg),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref-dir", required=True)
    ap.add_argument("--spacer", default=SPACER_EXACEL)
    ap.add_argument("--max-mismatches", type=int, default=MAX_MISMATCHES)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-sequences", help="write the distinct 23-mers, one per line")
    ap.add_argument("--out-windows", help="write every window: contig, strand, 1-based forward start, "
                                          "23-mer 5' to 3' on its strand, mismatches")
    a = ap.parse_args(argv)

    spacer_enc = encode(a.spacer.encode())
    if spacer_enc.size != PROTOSPACER_LEN or (spacer_enc > 3).any():
        sys.exit("the spacer must be 20 nt over ACGT")
    records, contigs, rows = [], 0, []
    for name, raw in iter_contigs(a.ref_dir):
        contigs += 1
        recs = scan_encoded(encode(raw), spacer_enc, a.max_mismatches)
        records.extend(recs)
        rows.extend((os.path.splitext(name)[0],) + r for r in recs)
        print(f"[scan] {name:<12} windows so far: {len(records):>8,}", flush=True)
    summary = summarise(records, a.spacer, a.max_mismatches)
    summary["contigs_scanned"] = contigs
    with open(a.out_json, "w") as fh:
        json.dump(summary, fh, indent=1)
    if a.out_windows:
        with open(a.out_windows, "w") as fh:
            fh.write("contig\tstrand\tstart_1based_forward\twindow_5to3\tmismatches\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
    if a.out_sequences:
        with open(a.out_sequences, "w") as fh:
            fh.write("\n".join(sorted({r[2] for r in records})) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
