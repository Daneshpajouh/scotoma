# The genome-scale end-to-end result, and how to re-check it

This directory holds the complete evidence for the paper's genome-scale comparison, together with a
single command that recomputes every number the paper states from these tables and prints a
claim-by-claim verdict:

~~~sh
python3 -I verify_genome_endtoend.py
~~~

It needs no genome, no network and no cluster. Exit status is 0 only when every claim is reproduced.

## What was run

Two builds of the same released program, identical flags, differing by one line of source: the
released CPU binary and a corrected binary that admits the terminal RNA-bulge placement. Scope is the
entire GRCh38 primary assembly, 25 records, chr1 to chr22, chrX, chrY and chrM, 3,088,286,401
assembled bases. Operating settings throughout: 20-base spacer plus NGG, four mismatches, zero DNA
bulges, one RNA bulge, both strands, CPU mode, the 110-guide panel.

## The result

| quantity | released | corrected |
|---|---:|---:|
| loci, one RNA bulge | 937,934 | 972,875 |
| loci, zero RNA bulges | 34,573 | 34,573 |
| deposited measured-active rows returned | 108 | 237 |

- Corrected only: **34,941** loci. Released only: **zero**, and the zero holds on each of the 25
  records separately, not merely in total, and on each of the 110 guides separately.
- Every one of the 34,941 is the same kind of alignment: class RNA, deletion at one-based spacer
  position 20, the terminal base, at four mismatches. All three fields take a single value across the
  whole export. There is no second class, no other position, no other mismatch count.
- Candidate overhead 34,941/937,934, which is 3.7253 percent.
- The bulge-free control is identical for the two builds, 34,573 loci, in total and on every guide.

## The independent agreement

The paper's headline empirical number, 129, is derived by a pattern-level census over a deposited
third-party assay table. This run never used it: it asked which of the 260 deposited measured-active
rows each binary returns. The released build returns 108 and the corrected build 237, a difference of
exactly 129. The remaining 23 rows are reference-matched but outside the paper's stated NGG scope, and
237 plus 23 is 260, so the accounting closes.

## The one flagged record, and why it does not reach the result

Twenty-four records pass every control. **chr22 is flagged `REVIEW_REQUIRED` for exactly one reason**:
it carries five alignments, two in the released run and three in the corrected run, whose matched PAM
is not literally `GG`. The program accepted them through an IUPAC-ambiguous base in the reference; a
strict NGG reading would exclude them. The flag is retained here rather than suppressed.

It does not touch the difference. **None of the 34,941 corrected-only alignments has a PAM that is not
literally `GG`**, the verifier checks this explicitly, and the count is zero. The five ambiguous
alignments lie in the set both builds find, so they cancel. Every other record reports zero such
alignments.

## What these tables do not establish

- **The ungapped baseline comparison is a magnitude check and must not be read as agreement.** The
  collector reports it as `MAGNITUDE_CHECK_ONLY`: observed bulge-free total 34,573 against a deposited
  34,652, a difference of 79 loci or 0.23 percent, with per-guide differences up to about ten. The
  deposited enumeration counts a different quantity than this tool's locus key, its enumerator is
  unavailable, and it carries no reference hashes or locus keys. Nothing here depends on it: the result
  is a difference between two builds on byte-identical inputs.
- This is sequence-similarity nomination, not cleavage and not clinical consequence. The 237 against
  108 concerns which deposited measured-active rows each build can report, not what happens in a cell.
- Timing and memory measurements are recorded per chromosome in the run directory, not in these
  tables, so this verifier makes no performance claim.

## Files

| file | contents |
|---|---|
| `per_chromosome.tsv` | per-record locus counts, control state, host and job id |
| `per_guide.tsv` | per-guide counts for all 110 guides, both builds, gapped and bulge-free |
| `corrected_only_alignments.part000.tsv` | all 34,941 corrected-only alignments with crRNA, DNA and class |
| `released_only_alignments.part000.tsv` | the released-only alignments: header only, no rows |
| `difference_classes.tsv` | the single difference class over the whole assembly |
| `active_sites.part000.tsv` | all 260 deposited measured-active rows with what each build returned |
| `ungapped_baseline_total.tsv` | the magnitude-only baseline comparison, per guide |
| `verify_genome_endtoend.py` | the claim-by-claim replay described above |

## The census route, run on its own

`census_predicts_genome.py` derives every return bit in `active_sites.part000.tsv` from the
assay sequences alone. It reads exactly two files, `../data/sites129.json` and that table, and
touches no genome, no executable and no network. The rule it applies is the paper's
pattern-generation rule and nothing else: a build reports a locus through a single RNA deletion
exactly when some deletion pattern in its own pattern set lies within the mismatch budget of the
19 nt protospacer anchored at that locus' PAM, and the PAM is NGG; the released pattern set omits
the terminal deletion and the corrected set contains it.

    python3 -I census_predicts_genome.py

Eighteen checks. The load-bearing one is that all 520 return bits agree, row by row, so the
agreement is not an agreement of two totals: a permutation that preserves both totals fails it.
Six falsification controls follow. Five vary the rule, each withholding a different spacer
position or changing the budget or the PAM policy, and each must fail to reproduce the bits; the
sixth exchanges two return bits while preserving the totals and must be caught. A prediction that
survived a wrong rule would establish nothing, so the controls are part of the claim.

`verify_genome_endtoend.py` now carries the same row-by-row prediction, and recomputes the
PAM-side boundary, the deletion placement and the mismatch count of all 34,941 exported
alignments from their own strings and coordinates rather than reading the annotation beside
them. It rejects each of seven deliberately corrupted copies of its own inputs: a renamed
record, a broken per-guide bulge-free control, an off-by-one coordinate, a relocated deletion
that keeps its labels, a permutation of two return bits, a rewritten flag reason, and a changed
total.

## What a second hostile pass changed

A reviewer with execution access built sixteen further corruptions and the 45-check version
accepted all of them: one locus repeated 34,941 times, coordinates translated out of the
assembly by 2^40, invented chromosomes and strands, counts moved between subgroups until one
went negative, destroyed assay identities, and a corrected-only target the released patterns
can in fact reach at three mismatches. The checker is now 66 claims and derives each of those
quantities instead of reading a label. Re-run against twenty-one corruptions, the seven earlier
ones and fourteen rebuilt from the reviewer's list, it rejects all twenty-one and still passes
the unmodified export.

One attack it cannot see: changing a mismatching DNA base to a different mismatching base while
retaining `reference_verified=1`. Detecting that needs the reference sequence, which is not
deposited. The limitation is stated rather than papered over, and `reference_verified` should be
read as the producing run's assertion, not as evidence this checker reproduces.

A cross-check worth naming: the 25 record lengths the checker pins sum to 3,088,286,401, which
is the assembled base count the run itself recorded.

The same reviewer also defeated an inference, not a number. A deliberately truncated rule, one
admitting only four of the nineteen non-terminal deletion positions, reproduces all 520 return
bits and passed every check. The reason is structural and is now asserted in
`census_predicts_genome.py`: these rows are terminal-selected, so all 260 corrected bits equal
the single predicate "the PAM is NGG", and the decisive non-terminal placement is position 19 on
248 of the 260. The 520-bit agreement is therefore site-by-site evidence for the terminal
mechanism, not an identification of the pattern rule. Identification is done separately, at
budget zero, where each of 704 distinct words is matched only by itself: the complete rule
matches all 704 and the truncated rule only 225.
