# Bulge shares, in vitro and in cells

`share_strata.py` reproduces the bulge shares in the appendix section "Bulge
policy, priced on both axes", and every count behind them, from three
third-party tables. It needs only the standard library and any Python 3.9 or
later. Run it once per table:

    python3 verification/bulge_shares/share_strata.py CHANGEseq_CR_Lazzarotto_2020_dataset.csv
    python3 verification/bulge_shares/share_strata.py GUIDEseq_CR_Lazzarotto_2020_dataset.csv
    python3 verification/bulge_shares/share_strata.py NewGUIDEseq_CR_Lazzarotto_2020_dataset.csv

The script knows each table by its SHA-256 and checks that table's figures. A
table can also be piped straight out of the archive, without unpacking it:

    unzip -p datasets.zip datasets/GUIDEseq/include_on_targets/GUIDEseq_CR_Lazzarotto_2020_dataset.csv \
        | python3 verification/bulge_shares/share_strata.py /dev/stdin

## What it checks

| Figure | Counts | Basis |
|---|---|---|
| 55.7% | 23,692 of 42,572 | bulged active sites at 4 or fewer mismatches, over all active sites at 4 or fewer mismatches |
| 45.4% | 15,698 of 34,578 | matched on total edits: bulged sites at 3 or fewer mismatches plus bulge-free sites at 4 or fewer |
| 11.7% | 23,692 of 202,040 | the published analysis of the same deposit, which counts every active site, on-targets included |
| 25.8% | 312 of 1,207 | GUIDE-seq, in cells: bulged active sites at 4 or fewer mismatches, over all active sites at 4 or fewer mismatches |
| 27.6% | 147 of 533 | NewGUIDE-seq, in cells: the same rule |

It also checks the bulged active rows at 0 to 4 mismatches (9, 246, 2,715,
12,728 and 7,994 in vitro; 2, 47, 139, 95 and 29 in GUIDE-seq; 3, 14, 57, 57
and 16 in NewGUIDE-seq) and the single-bulge rows in each table (1,749,342,
739,570 and 438,375).

An active row is a row with reads above zero. Every share excludes the
zero-edit active rows, which are the on-targets, one per guide with a measured
on-target: 109 in vitro, 58 in GUIDE-seq, 20 in NewGUIDE-seq. Kept in, the
shares are 55.5% and 45.3% in vitro, 24.7% and 26.6% in cells, so the exclusion
is a condition of reproducing the printed figures. The script checks both
versions.

Every share counts the active sites that carry a bulge and no mismatch: 9 in
vitro, 2 in GUIDE-seq, 3 in NewGUIDE-seq. Dropped, the shares are 55.6%, 25.7%
and 27.2%. The 23,692 that the published analysis counts includes them, so
they stay in, on the same rule for all three tables.

## The input

The tables are not released here. They are the CHANGE-seq, GUIDE-seq and
NewGUIDE-seq tables of Lazzarotto et al. (2020) as redistributed with
CRISPR-Bulge:

- repository: `OrensteinLab/CRISPR-Bulge` on GitHub
- archive: `files/datasets.zip`, stored in Git LFS, 524,400,344 bytes,
  SHA-256 `f892f70ba4ac3b05b03b2171b4ad38746630de08ad630e650f355dd61203eab0`
- tables inside it, all under `datasets/<assay>/include_on_targets/`:
  - `CHANGEseq_CR_Lazzarotto_2020_dataset.csv`, 878,760,836 bytes,
    SHA-256 `68f3e1d7eb9e9e18e3ddc859752be8d156378a3e681b0fff329fe679366adc11`
  - `GUIDEseq_CR_Lazzarotto_2020_dataset.csv`, 381,626,484 bytes,
    SHA-256 `a3398cba2ef452ec37ebebda38a44aa62615c66165a43699691ca210f467bcd2`
  - `NewGUIDEseq_CR_Lazzarotto_2020_dataset.csv`, 207,094,458 bytes,
    SHA-256 `ac9ca55a384a9af5ba37d667766e6723a7c211b5fa1d54875e450db23a13a23b`

The Zenodo record 10.5281/zenodo.10902445 holds only the 134-byte LFS pointer
for the archive, not the archive itself, so fetch it with Git LFS. The archive
holds other tables with near-identical names; the three in
`include_on_targets` above are the ones. The script prints the SHA-256 of what
it read and fails on any other file: exit 1 if it parses as such a table, exit
2 if not. `--table changeseq`, `guideseq` or `newguideseq` names the claims to
check against a file that is none of the three, and that file still fails the
digest check.

## Reading the output

One line per check, `ok` or `NO`, then `VERDICT: PASS` or `VERDICT: FAIL`. Exit
0 only on PASS; 1 on FAIL; 2 if the file could not be read or parsed as such a
table, in which case no digest is printed. Controls run every time. Keeping the
on-targets must move every share, to the values above. Dropping the bulged
sites at zero mismatches must move every share, to 55.6%, 25.7% and 27.2%. And
in vitro, each neighbouring rule for the total-edit share (bulged sites to 4
mismatches, or bulge-free sites to 3 or to 5) must give a count other than
15,698 of 34,578, so that the match shows which rule was used.
