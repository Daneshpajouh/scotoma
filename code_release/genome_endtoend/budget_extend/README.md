# The deficit at the budgets the literature reports

Paper A states the released build's genome-scale deficit at four substitutions. This directory
carries the same comparison at five, six and seven, over the same 110-guide panel and the same
GRCh38 primary assembly, at the same PAM policy and bulge setting.

| budget | released | corrected | corrected only | released only | terminal attributable | guides with a deficit | share of corrected |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 5,941,903 | 6,195,422 | 253,519 | 0 | 253,519 (100%) | 73 of 110 | 4.092% |
| 6 | 32,769,159 | 34,315,684 | 1,546,525 | 0 | 1,546,525 (100%) | 73 of 110 | 4.507% |
| 7 | 144,444,173 | 150,719,531 | 6,275,358 | 0 | 6,275,358 (100%) | 73 of 110 | 4.164% |

The same 73 guides carry the deficit at all three budgets and the same 37 carry none, which is a
property of the guide's run structure and not of the budget.

## What this does and does not say

It says the deficit does not close as the budget rises, and that every missed locus is one the
terminal deletion word reaches. It does **not** say that a locus missed at four is still missed
at five. Those are different sets: adjacent deletion words differ in one coordinate, so the
radius-four union of the complete family is contained in the radius-five union of the emitted
family, and a locus lost at one budget can be recovered at the next by a different placement.
The counts above are each measured at their own budget for that reason.

## Two programs, no shared code

`panel_extend.py` is the producer. It runs both published binaries, gated on their SHA-256 and
on the staged reference's record count and byte total, shards the panel eight ways and tallies
loci per guide. Its output is `shard_<i>_record.json`.

`verify_panel.py` is an independent reader of the same alignment tables. It re-derives the locus
convention from its definition, guide, chromosome, strand and PAM-side forward boundary, zero
based and unmerged, and recomputes the terminal attribution from each alignment's own genomic
sequence rather than inheriting it. Its output is `shard_<i>_VERIFY.json`. Its controls:

- every row's crRNA must equal its guide plus `NNN`, and its mismatch count must be within budget
- the released locus set must be a subset of the corrected one
- forgery: attributing by the **first** deletion word instead of the terminal one must explain
  none of the corrected-only loci, and it explains zero
- corruption: a single-base corruption of the guide must attribute nothing, and it attributes
  nothing while the real terminal word attributes thousands

One detail is worth recording because the first run of the verifier flagged it. Eleven of the
8.08 million corrected-only loci sit at an IUPAC ambiguity code in the reference. A strict
base-by-base comparison counts the ambiguity as a mismatch and pushes those loci one over
budget; the search itself treats a code as matching any base it denotes. The verifier now scores
ambiguity the way the search does and reports how many rows matched only through a code. The
paper's own chr22 note records the same phenomenon at four substitutions.

## Replay

The two alignment tables are 77 GB and are not redistributed. From the deposited records:

    python3 -I verify_budget_extend.py

re-derives every total above from the per-guide records, asserts that the producer and the
independent verifier agree guide by guide at every budget, asserts that each shard's controls
behaved as stated, and includes a control that a single altered count breaks the agreement test.
60 claims, all reproduced.

To regenerate the tables, stage the reference and both binaries, set three variables and run
on a Slurm cluster:

    export PANEL_ROOT=/path/to/this/directory
    export PANEL_FASTA=/path/to/staged/GRCh38          # 25 records, 3,132,404,952 bytes
    export PANEL_BIN=/path/to/both/pinned/binaries
    sbatch panel_extend.sbatch        # 8 tasks, 16 cores, 120 GB, about 20 minutes each
    sbatch verify_panel.sbatch        # 8 tasks, 2 cores, 120 GB, about 5 minutes each
    python3 aggregate_panel.py

Set `--account` in both sbatch files to your own allocation. `panel_extend.py` refuses to run
unless the two binaries' SHA-256 digests and the reference's record count and byte total match
the pinned values, so a wrong input fails rather than producing a number.

## Provenance of these scripts

The three cluster scripts were run with their roots hard coded to the directory they ran in, and
are deposited with those paths replaced by `PANEL_ROOT`, `PANEL_FASTA` and `PANEL_BIN` so they
are portable. Nothing else changed. Both digests are recorded so the substitution is checkable:

| script | as executed | as deposited |
|---|---|---|
| `panel_extend.py` | `c1bade84095d0d496f5bdc0c9c3ebffa280cac0d90162ca57682668f3490b398` | `7f4f0c4eba9cf3a35138392c5a099ae6c92273d4e561c5b2d0920ed5825e1370` |
| `verify_panel.py` | `415ef3befddff67a6585bb1627407645f12e1e9d1d9cdfefcb1c940f442627a6` | `c6c54d46cdc702bcd89f968415fe7bf5962ed4e128708265e8bda443609ad389` |
| `aggregate_panel.py` | `6cfcf4810f520ad3220e0e6e2f2bc79eb7d88a2589560299694a38d7e27431f2` | `e90055159af3ef65846872a4041d0b479530abd6edbea44be066e367c4ed4d50` |

`verify_budget_extend.py`, the checker a reader runs, was written against the deposited records
and needs no cluster: `353806c0c9b22c71ade7770a246cf92a04982a83537780b7674ca8ef7c9ce34b`.

`panel_guides.txt` is the panel, one guide per line with its identifier.
