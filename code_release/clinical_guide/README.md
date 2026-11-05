# The omission on a guide that is in the clinic

This directory holds the raw outputs and a one-command replay for the paper's comparison on
**sg1617**, the BCL11A $+58$ erythroid enhancer guide used in the clinical trials for sickle cell
disease and beta-thalassemia (Cancellieri et al., *Nature Genetics* 55(1):34-43, 2023,
doi:10.1038/s41588-022-01257-y). Spacer: `CTAACAGTTGCTTTTATCAC`.

~~~sh
python3 -I verify_clinical_guide.py
~~~

Recomputes every number the paper states for this guide and prints a claim-by-claim verdict. No
genome, no network, no cluster. Exit status 0 only when all 34 claims are reproduced.

## What was run

Two builds of the same released program, differing by one line of source, on the whole GRCh38
primary assembly: 25 records, 3,088,286,401 assembled bases. The reference hashed to the same
value the build gate recorded (`9cce8b9264...`), the scaffolds were excluded so the universe
matches the panel run exactly, and both binaries hashed to the build gate's values
(`b8991a89...` released, `cee92311...` corrected). Settings: 20-base spacer plus NGG, both
strands, zero DNA bulges, CPU mode; one RNA bulge against a zero-bulge control; and the mismatch
budget swept from 0 to 6.

## The result

Loci under the paper's convention, which is (guide, chromosome, strand, PAM-side forward
boundary), zero-based, with no 3-bp merging.

| mismatch budget | released | corrected | corrected only | released only |
|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 0 | 0 |
| 1 | 5 | 5 | 0 | 0 |
| 2 | 34 | 35 | **1** | 0 |
| 3 | 415 | 464 | **49** | 0 |
| 4 | 4,414 | 4,778 | **364** | 0 |

- **Released only is zero at every budget.** The repair never removes a locus.
- The bulge-free control is identical as a multiset at every budget.
- At the deployed budget of four, the 364 loci are 8.2465 percent of the released build's own
  output, more than twice the 3.7253 percent measured across the 110-guide panel.
- Every one of the 364 is the same class: an RNA deletion at one-based spacer position 20, the
  terminal base, at four mismatches, and none has a PAM that is not literally `GG`.

## The ranking error, which is the part that matters

A missed locus does not vanish at a larger budget; it reappears with the wrong mismatch count.
For **every** locus the released build cannot report, its best available alignment costs
**exactly one extra mismatch**: all 49 of the budget-three misses are reported at budget four
with an inflated count, and all 364 of the budget-four misses are reported at budget five with an
inflated count, the inflation being 1 in every single case.

The clearest instance is `chr3:191,469,717` on the plus strand. Its true single-bulge alignment
is a terminal deletion at **two** mismatches, `CTAATAGTTTCTTTTATCA-AGG`. The released build cannot
place that deletion, so its best alignment for the same locus is a deletion one position earlier
at **three** mismatches. At a two-mismatch budget the site is absent from its output altogether;
at a larger budget it is present but looks one mismatch worse than it is.

Any screen that ranks or thresholds candidates by mismatch count therefore under-ranks precisely
the closest matches.

## Scope, stated plainly

- This is **sequence-similarity nomination, not cleavage** and not a clinical consequence. It
  concerns which candidate loci a tool reports and with what mismatch count, not what happens in
  a cell. Cancellieri et al. make the same distinction: the possibility of a somatic alteration
  does not imply a functional one.
- The guide identity is as published by Cancellieri et al., who give `sgRNA-1617:
  CTAACAGTTGCTTTTATCAC` in their Methods and describe it as the guide in the sickle cell and
  beta-thalassemia trials. No equivalence to any particular product's sgRNA designation is claimed
  here.
- Budgets 5 and 6 were run and their alignment counts and output hashes are recorded in
  `sweep_record.json`, but those raw outputs are large and are not deposited, so the verifier
  makes no claim that depends on them. Every claim it checks is replayable from the files here.
- Timing is recorded in `run_record.json` for provenance only; no performance claim is made from
  a single-guide run.

## Files

| file | contents |
|---|---|
| `released_rna1.tsv`, `corrected_rna1.tsv` | raw tool output at the deployed budget, one RNA bulge |
| `released_rna0.tsv`, `corrected_rna0.tsv` | the bulge-free control at the same budget |
| `sweep_mm{0,1,2,3}_{released,corrected}_rna{0,1}.tsv` | the budget sweep outputs that the verifier uses |
| `run_record.json` | gates, hashes, staged record count and bases, per-run timing |
| `sweep_record.json` | alignment counts and output hashes for every sweep run, budgets 0 to 6 |
| `verify_clinical_guide.py` | the 34-claim replay described above |
