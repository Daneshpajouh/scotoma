# In vitro bulge shares

`share_strata.py` reproduces the two in vitro bulge shares in the appendix
section "Bulge policy, priced on both axes", and every count behind them, from
one third-party table. It needs only the standard library and any Python 3.9
or later.

    python3 verification/bulge_shares/share_strata.py CHANGEseq_CR_Lazzarotto_2020_dataset.csv

## What it checks

| Figure | Counts | Basis |
|---|---|---|
| 55.7% | 23,692 of 42,572 | bulged active sites at 4 or fewer mismatches, over all active sites at 4 or fewer mismatches |
| 45.4% | 15,698 of 34,578 | matched on total edits: bulged sites at 3 or fewer mismatches plus bulge-free sites at 4 or fewer |
| 11.7% | 23,692 of 202,040 | the published analysis of the same deposit, which counts every active site, on-targets included |

It also checks the bulged active rows at 0 to 4 mismatches (9, 246, 2,715,
12,728 and 7,994) and the 1,749,342 single-bulge rows in the table.

An active row is a row with reads above zero. Both shares exclude the 109
zero-edit active rows, which are the on-targets, one per guide with a measured
on-target. Kept in, the shares are 55.5% and 45.3%, so the exclusion is a
condition of reproducing the printed figures. The script checks both versions.

## The input

The table is not released here. It is the CHANGE-seq table of Lazzarotto et al.
(2020) as redistributed with CRISPR-Bulge:

- repository: `OrensteinLab/CRISPR-Bulge` on GitHub
- archive: `files/datasets.zip`, stored in Git LFS, 524,400,344 bytes,
  SHA-256 `f892f70ba4ac3b05b03b2171b4ad38746630de08ad630e650f355dd61203eab0`
- table inside it: `datasets/CHANGEseq/include_on_targets/CHANGEseq_CR_Lazzarotto_2020_dataset.csv`,
  878,760,836 bytes, SHA-256 `68f3e1d7eb9e9e18e3ddc859752be8d156378a3e681b0fff329fe679366adc11`

The Zenodo record 10.5281/zenodo.10902445 holds only the 134-byte LFS pointer
for the archive, not the archive itself, so fetch it with Git LFS. The archive
also holds GUIDE-seq tables with near-identical names; the CHANGE-seq table in
`include_on_targets` is the one. The script prints the SHA-256 of what it read
and fails on any other file: exit 1 if it parses as this table, exit 2 if not.

## Reading the output

One line per check, `ok` or `NO`, then `VERDICT: PASS` or `VERDICT: FAIL`. Exit
0 only on PASS; 1 on FAIL; 2 if the file could not be read or parsed as this
table, in which case no digest is printed. Two controls run every time. Keeping
the on-targets must move both shares, to 55.5% and 45.3%. And each neighbouring
rule for the total-edit share (bulged sites to 4 mismatches, or bulge-free
sites to 3 or to 5) must give a count other than 15,698 of 34,578, so that the
match shows which rule was used.
