# The candidate set for the worked example

The paper's worked example needs one input that is not derived from first principles: the
number N of candidate off-target sites of the exa-cel guide. This note says exactly how it is
obtained.

## Definition

The spacer is `CTAACAGTTGCTTTTATCAC` (sg1617). A *window* is a stretch of 20 nt (the
protospacer) followed by a 3-nt PAM, read on either strand of GRCh38, such that

* the PAM is NGG or NAG,
* the protospacer differs from the spacer in at most six positions (substitutions only; no
  bulges),
* the window contains only A, C, G or T (lowercase soft-masked bases count as the same base).

The unit is the window, that is a (strand, position) pair, not the 23-mer sequence: two
windows at different loci with the same 23-mer are two sites. The reference is the 25 primary
contigs chr1 to chr22, X, Y and M.

## Counts

| Quantity | Value |
|---|---|
| Windows returned by the scan | 28,176 |
| of which with zero mismatches (the on-target site) | 1 |
| **Candidate off-target sites N** | **28,175** |
| Distinct 23-mers (NGG and NAG) | 26,335 |
| Distinct 23-mers (NGG only) | 10,442 |

Windows by mismatch count: 0: 1, 2: 1, 3: 18, 4: 362, 5: 3,185, 6: 24,609. Every bin has at
least as many windows as distinct sequences because a sequence can occur at several loci.

The on-target window is cleaved by design and is not a candidate off-target, so
N = 28,175. The site count, not the sequence count, is the denominator: gamma = k/N is smaller
with the larger N, so the requirement is larger, which is the conservative direction.

## Reproduction

`enumerate_sites.py` implements the scan (about four minutes and several GB of memory for the
whole genome). The output reported here was produced by two independent implementations of
the scan and agrees exactly: the same window histogram, and a byte-identical sorted list of the
26,335 distinct 23-mers, SHA-256
`24ee0b98a956aab31c3a66bce8ff631e25d9a85a1693064e15babedc8f7a3d33`.
`test_enumerate_sites.py` checks the scan against a window-by-window reference implementation
on sequences with planted windows.

## What the count does and does not cover

The enumeration is reference-based. It does not contain sites created by genetic variants
(variant-aware searches add them, and one has been reported for this guide; see the paper), it
has no bulges, and it stops at six mismatches. The certified statement in the paper is about
the enumerated sites only, and it is conditional on every cleavable site being among them.

## The on-target window, the reference and the site list

* The single zero-mismatch window is `CTAACAGTTGCTTTTATCACAGG` on chromosome 2, minus strand,
  forward-strand positions 60,495,261 to 60,495,283 (GRCh38), in the BCL11A locus. It is excluded
  by coordinates; the histogram has no other perfect match.
* Four windows overlap a neighbour on the same strand (distance below 23 positions); windows
  are not merged, and merging would lower N by two.
* `data/exacel_windows.tsv.gz` lists every window (contig, strand, 1-based forward start of
  the 23-base window, the window read 5' to 3' on its strand, mismatches), 28,176 rows;
  `data/exacel_enumeration_summary.json` records its SHA-256, the size and SHA-256 of each of
  the 25 reference contig files (raw sequence, one file per contig; sizes equal the GRCh38
  primary-assembly contig lengths), and the counts.
