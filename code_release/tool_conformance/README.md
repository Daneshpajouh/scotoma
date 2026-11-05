# Conformance of bulge-aware tools against the exact minimum

    python3 -I verify_tool_conformance.py

Thirty-two checks, recomputed from the raw tool outputs here and from the guide spacers. The
exact minimum is never read from a file: it is derived per guide by set construction. No genome,
no executable, no network.

## The criterion

It has two sides, and each is a theorem rather than a convention:

* **complete**, emit no fewer than the r distinct single-deletion words, because r is the exact
  minimum number of equal-radius balls covering the deletion union, with replacement centres
  anywhere in the four letters. The covering theorem is over concrete centres; it does not admit
  wildcards, and an earlier version of this file said it did;
* **minimal**, emit no more, because those r words already cover it. Wildcards are excluded by
  a separate result rather than by assumption: among patterns that accept nothing outside the
  union, no pattern can carry a wildcard at four letters, so the only sound centres are the r
  words themselves. Soundness is doing the work there, and without it one all-wildcard pattern
  would accept everything and the minimum would collapse to one.

Below the minimum a tool loses sites. Above it a tool repeats work and output. The minimum is
the unique point that is neither, which is what makes the criterion decidable.

## The verdict on the audited panel, 110 guides, minimum 1,634 words

| route | emits | lost | output rows | repeats |
|---|---:|---:|---:|---:|
| Cas-OFFinder-bulge v1.2 beta, 2016 wrapper | 1,491 words | 143 | n/a | n/a |
| Cas-OFFinder 3.0.0b3 and develop | 1,561 words | 73 | n/a | n/a |
| CRISPRitz 2.7.1 | 1,634 words | 0 | 2,200 | 566 |
| CRISPRitz 2.9.0 | 1,634 words | 0 | 2,200 | 566 |

**No audited tool emits exactly the minimum.** Two are below it and lose sites. Two reach it and
then report it in 2,200 alignment rows, of which 566, or 25.7%, repeat a locus already reported,
one locus seven times. Every repeat group was checked to be one underlying target sequence whose
members differ only in where the gap sits inside a run; the verifier asserts both.

A third party reports the same redundancy independently and its cost: the CHOPOFF preprint calls
CRISPRitz's output "very redundant", and restricted its published comparison to Levenshtein
distance three because the distance-four output reached 160 GB and was "too large to parse".

## Inputs

    per_guide_counts.tsv              550 rows: five audited routes, 110 guides each, with the
                                      words each emits and its deficit
    crispritz_per_guide.tsv           the same panel for the two CRISPRitz engines
    crispritz_2.9.0_rna1.targets.tsv  raw native output, 2,200 alignment rows
    crispritz_2.7.1_rna1.targets.tsv  raw native output from the older engine

The Cas-OFFinder rows come from executing the unmodified parser and the official wrapper
dry-run; the CRISPRitz rows are that engine's own emitted alignments on the planted fixture.
Counts on a fixture measure what a route emits per guide, not genome-wide site totals.

## A control, because a criterion that nothing can satisfy is not a criterion

The final check builds an emitter that returns exactly the distinct words and confirms it scores
zero lost and zero redundant. The test can be passed; no audited tool passes it.

## What a second hostile pass changed

The completeness axis used to read `distinct_deletion_words_recovered` from the summary table.
A reviewer rotated every raw protospacer A to C to G to T, preserving the 2,200 rows, the 1,634
locus keys, the 566 repeats and the sevenfold maximum, and the verifier still reported that all
1,634 distinct words were recovered although **no row contained a valid deletion word at all**.

The word sets are now derived from the raw alignments: every row's target must be a deletion
word of its own guide's spacer, carry an NGG PAM, and be an exact fixture match, and the derived
per-guide word count must reach r. Re-run, the rotation attack fails six assertions and a
route-specific spacer relabelling fails fourteen, while the unmodified deposit passes all 42.
