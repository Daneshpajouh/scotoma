# What the omitted alignment classes cost a published dataset

    python3 -I verify_published_dataset_hole.py                 # unbound, exits 2
    python3 -I verify_published_dataset_hole.py <dir with the two public CSVs>

The second form is the real one. The first cannot bind the absence statements to the files
they are about, says so, and exits 2 rather than 0.

## The claim, stated precisely

CRISPR-Net predicts off-target activity with indels, and Cas-OFFinder enumerated its hg38
candidates. In its two deposited indel datasets, all 238,705 enumerated RNA-bulge rows on
sixteen guides place the deletion at spacer positions 2 to 19, never at 1 or 20, and the DNA
slots are likewise 1 to 19. Those sixteen guides need 241 distinct deletion words, 21 of which
sit only at position 1 or 20. A scan of the whole GRCh38 primary assembly for all 241, at the
settings the datasets' own rows exhibit, finds that those 21 words reach 35,272 loci and
10,888 off-target sequences no interior word reaches within the budget.

**What this costs the files is an alignment class, not the sites.** Of the 10,888 sequences,
10,854 are in the published files anyway: 6,267 as ungapped twenty-mers at three or four
mismatches, and 4,587 also inside a DNA-bulge twenty-one-mer. Only **34** are absent outright.
A model fitted to these files has therefore seen the sites and has never seen the endpoint
deletion that reaches them.

### A correction, recorded

An earlier version of this directory stated that 10,887 of the 10,888 were absent from the
published data, and the paper said so. That was wrong. It came from testing absence only
against the 19-mers that the files' own *gapped* rows spell out, which calls a site absent
when the files plainly contain it under a different alignment. Under that narrow test all
10,888 come out absent; under the conservative test used here, 34 do. The claim in the paper
now follows the conservative test, and both are computed and printed so the difference cannot
hide again. The reason the error survived is also recorded: the old checker read the absence
flags off a deposited column instead of recomputing them.

## The two definitions

**Unreachable.** A (guide, off-target 19-mer) pair the scan reached, no word of which sits at a
spacer position the generator's loop visits: its position set misses 2..19 entirely, so only a
deletion at position 1 or 20 produces it. Each one is checked here to be within budget at an
endpoint and beyond budget at every interior position, so re-anchoring an endpoint gap inward
cannot report it either; the closest interior word to any of them is at distance four.

**Held by the published files.** A 19-mer occurs in the files for a guide when it is a
contiguous 19-column window of any protospacer those files report for that guide, gap and pad
characters removed. This is deliberately the generous direction: a sequence counts as held even
when the row holding it is ungapped or carries a DNA bulge. The narrow alternative is computed
too and reported as `absent_gapped_only`.

## Replay from the raw output

    tar xzf scan_raw_outputs.tgz -C <scan>          # per-record scan outputs, 25 records
    python3 -I make_reached.py <scan> <scan>        # -> reached_positions.tsv, reached_loci.tsv
    python3 -I derive_hole.py <scan> <csvdir> <out> # -> the two tables below

`make_reached.py` reduces 850,463 alignment rows to the two statistics everything else uses.
`derive_hole.py` builds both tables from those plus the public CSVs. The verifier runs
`derive_hole.py` itself and compares its output to the deposited tables byte for byte, so a
table that was edited rather than derived fails.

## Files

    scan_raw_outputs.tgz        the scan's own per-record output, 25 records
    queries are in guides_meta.json: spacer, deletion-word library, positions per word
    reached_positions.tsv       122,170 (guide, 19-mer) pairs with the positions reaching them
    reached_loci.tsv            769,330 loci, 35,272 of them only at position 1 or 20
    per_guide_hole.tsv          per guide, every count the paper quotes
    unreachable_sequences.tsv   the 10,888, each with whether the files hold it and how
    gap_position_census.json    the deposited census of the two CSVs, recomputed when bound

## What this cannot see

The scan is trusted. Nothing here re-reads GRCh38, so a fabricated alignment row in
`scan_raw_outputs.tgz` would pass: that is the same blind spot the genome checker documents,
and closing it needs the reference. Everything downstream of the scan is recomputed.

## Public sources

Repository `github.com/OrensteinLab/CRISPR-Bulge` at commit
`0843c485b1d476f899bedb93a7dd2df42547c1e1`. The two files, with the digests the verifier pins:

    CIRCLE_seq_10gRNA_wholeDataset.csv  86d88947ef9c0a1b545814c4e4241e4f6103d15a4d056d697a4ea5b6017648a8
    elevation_6gRNA_wholeDataset.csv    56f1c95a0a3d61b873af984f7a2110e27538540b49904b913982d63764dcbb55
