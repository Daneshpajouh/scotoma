# Comparator locus agreement, 110-guide panel on hg38: the record behind the locus counts

Scientific Director seat, 2026-10-23. Re-run at 12:32 PT (outputs written 12:32:48) from the saved tool outputs; every output file
is byte-identical to the first run of 06:17 PT.

## What this folder records

CRISPRitz v2.9.0 against SWOffinder v1.0.1 on the hg38 primary chromosomes (hg38.fa SHA-256
5be01555d98347fdb3714dc84c6f77c9d8bc774adcf32c6f7a8fa06f5baf5e51) with the 110-guide panel
(sgrna_panel110.txt, SHA-256 4be44a5da3ac83942f715f2eddb08be12e976235ed8f20b220290e03d840ff88).

Scope A: at most 4 mismatches and at most one bulge of either kind. SWOffinder --maxE 5 --maxM 4
--maxMB 4 --maxB 1. CRISPRitz -mm 4 on its bulge-1 index, as the union of -bDNA 1 -bRNA 0 and
-bDNA 0 -bRNA 1.

| scope A | count |
|---|---:|
| SWOffinder loci | 1,207,623 |
| CRISPRitz loci | 1,207,622 |
| found by both | 1,207,622 |
| shared loci with the same best alignment | 1,207,622 of 1,207,622 |
| SWOffinder only | 1 (chr7, + strand, PAM coordinate 61,063,449; the site contains an N) |
| CRISPRitz only | 0 |
| guides whose locus counts differ | 1 of 110, by one locus |

Scope B (one RNA bulge, no DNA bulge): CRISPRitz 978,104 loci; SWOffinder rows filtered to no DNA
bulge 973,505, a lower bound because SWOffinder keeps one alignment per locus; 4,599 CRISPRitz only.

## Files

- `idx1_panel_compare_out.txt`: the printed result above.
- `idx1_panel_per_guide.tsv`: loci per guide for each tool. The columns sum to 1,207,622 and 1,207,623.
- `idx1_panel_either_swo_only.tsv`, `idx1_panel_either_crz_only.tsv`: scope A differences (1 row, 0 rows).
- `idx1_panel_rna_only_crz_only.tsv`, `idx1_panel_rna_only_swo_only.tsv`: scope B differences (4,599 rows, 0 rows).
- `compare_panel_idx1.py`: the script that wrote them. It reads the CRISPRitz outputs and the SWOffinder
  locus keys at the relative paths in its source.
- `inputs.sha256`: SHA-256 of all 51 inputs the script read. The 50 CRISPRitz files match, hash for
  hash, `expected_crispritz.sha256` in benchmark rerun kit v5. The SWOffinder key file (panel_keys.tsv,
  112 MB) is not in this folder; its hash is listed.

## Not in this folder

The raw tool outputs (about 209 MB of CRISPRitz files and the SWOffinder CSVs) and the SWOffinder key
file. The rerun kit regenerates and checks the raw outputs.
