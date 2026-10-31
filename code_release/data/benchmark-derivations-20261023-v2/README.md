# Benchmark derivations, 2026-10-23 (version 2)

Scientific Director seat. A data source and a script for benchmark figures that until now
appeared only in notes, so a whole-token figure gate finds them in a deposited data file.

## Files

- `derived_figures.tsv`: one row per figure, with the value bare and with thousands commas.
- `derive_figures.py`: writes that table. Hashes every input first and stops on any difference.
- `paper_fwd_ATTGAGATAGTGTGGGGAAGNGG.csv`, `paper_keys.tsv`: the paper-guide SWOffinder outputs
  (milestone 3), listed in `paper_inputs.sha256`.
- `expected_crispritz.sha256`: copied unchanged from rerun kit v5. The two CRISPRitz passes
  (about 208 MB) are not here; rerun kit v5 regenerates them and this list pins them.
- `figure_lookup.py`: looks up bare figures as whole tokens in data files and notes.
- `test_boundaries.py`: the token cases. Version 1 of figure_lookup.py failed 12 of 24: it
  counted a figure inside a hex hash as a data hit, and a bare input missed copies printed with
  commas or space groups. Version 2 fails none.

## Rerun

    python3 derive_figures.py PAPER_DIR CRZ_DIR > derived_figures.tsv

PAPER_DIR is this folder. CRZ_DIR is the OUTDIR of rerun_benchmark.sh (kit v5), which holds
crz_dna/ and crz_rna/.

## Figures it supports

- 4,569: distinct PAM loci on both strands for the paper guide (2,357 plus, 2,212 minus).
  It is not the number of distinct PAMs in the 5,110 forward rows; that is 4,386.
- 4,569 + 724 extra minus rows - 183 lost minus PAMs = 5,110 rows (identity row = 1).
- 34,652: rows with bulge type X (no bulge) in each CRISPRitz pass. A row count, not a site count.
- 2,592,182: all CRISPRitz rows, both passes.

## Checks run on 2026-10-23

- Exit 0 on the pinned inputs; every value matches milestone 3 and timing-parity.
- A changed paper_keys.tsv gives DIFFERENT and exit 1. A missing CRISPRitz folder gives MISSING
  and exit 1.

Occurrence is not provenance. A figure found in `derived_figures.tsv` is computed by
`derive_figures.py` from pinned inputs; a figure found only in a note is not.
