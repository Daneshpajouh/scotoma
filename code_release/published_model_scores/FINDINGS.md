**Encoder convention, pinned before any scoring.** The published encoder uses `_` for a biological gap and `-` for zero-vector left padding. Each deposited 23-character pair was converted by replacing only alignment `-` with `_`; every non-gap character was fixed. The encoder then prepended one padding `-` to each sequence, producing a 24-column, seven-channel tensor. It copies each guide `N` from the DNA symbol in that column; the deposited guide PAM is `NNN`. Real DNA PAMs were retained, but that guide mask makes PAM columns matches, so these scores do not estimate PAM specificity.

The original vendored encoder and the reused forward encoder agreed exactly on all 34,941 treatment and 34,941 control tensors, and all 260 measured-active tensors. Treatment and control each converted 34,941 gap characters and preserved 1,572,345 non-gap characters. Every alignment was representable. Passing a biological `-` through unchanged fails the gap-versus-padding negative control. See [the frozen policy](published_model_aligned_scores/ENCODER_POLICY.md) and [the exhaustive encoder audit](published_model_aligned_scores/results/encoder_convention.json).

**1. Whole-assembly nomination scores, with the matched control beside them.** The unchanged deposited NumPy forward pass and checkpoint were reused, without fitting, training or tuning. The threshold was strictly greater than exact `3/1000` throughout.

| Arm | Scored loci | Above cutoff | Percent above |
|---|---:|---:|---:|
| Corrected only | 34,941 | 1,029 | 2.944964% |
| Released returned, matched | 34,941 | 187 | 0.535188% |

The arms match exactly within each of 73 guides, at four protospacer mismatches, one RNA bulge and literal NGG. The control contains distinct released-returned loci selected by the lowest fixed SHA-256 locus ranks per guide, without inspecting scores. Multiple eligible alignments of one locus were represented by the lexicographically smallest DNA alignment. Treatment uses the deposited terminal alignment; control uses a released interior alignment. Chromosome, sequence context and gap position were not matched, so this is a descriptive comparison rather than a causal effect of the repair.

Raw released outputs were recovered by read-only transfer from the original run and reconciled to the deposited per-guide and per-chromosome counts: 1,831,905 raw alignment rows, 937,934 distinct loci. There were 577,485 eligible control loci for the treatment guides; 542,544 were not selected and were not scored. 2 otherwise eligible raw alignment rows with ambiguous DNA/PAM symbols were excluded and not scored; the published encoder cannot represent those symbols. No neighbouring score was substituted.

Released-only export: 0 data rows, confirmed by parsing it. The enumeration totals remain released 937,934 versus corrected 972,875; the ungapped build control is 34,573 versus 34,573, including agreement per guide. These baseline tables contain counts, not scoreable sequence pairs. Their separate deposited magnitude baseline is 34,652; it is not an alignment control or a claim of exact reference agreement.

**2. Full score distributions and paired guide breakdown.** The absolute difference in above-cutoff rates is 2.409776 percentage points. Most rows in both arms fall below the cutoff. All row scores are deposited; the following distribution summaries retain the control.

| Statistic | Corrected only | Matched released returned |
|---|---:|---:|
| Mean | 0.00108808782 | 9.14623725e-05 |
| Population SD | 0.0165596397 | 0.00213051356 |
| Minimum | 2.22745637e-15 | 7.94545193e-18 |
| 1% | 1.25407031e-12 | 7.2523816e-14 |
| 5% | 3.94486797e-11 | 2.97888037e-12 |
| 25% | 8.04570441e-09 | 4.4124659e-10 |
| Median | 2.95714996e-07 | 1.23697224e-08 |
| 75% | 1.12092756e-05 | 2.90997244e-07 |
| 95% | 0.00112311009 | 2.51551483e-05 |
| 99% | 0.0145008842 | 0.000911719121 |
| Maximum | 0.94206242 | 0.225948925 |

| Score interval | Corrected only | Matched released returned |
|---|---:|---:|
| [0, 0.000001) | 20,542 | 29,051 |
| [0.000001, 0.0001) | 9,554 | 4,909 |
| [0.0001, 0.001) | 2,980 | 649 |
| [0.001, 0.003) | 836 | 145 |
| [0.003, 0.01) | 567 | 126 |
| [0.01, 0.1) | 392 | 56 |
| [0.1, 0.5) | 59 | 5 |
| [0.5, 1] | 11 | 0 |

Quantiles use linear interpolation. Histogram bins include their lower endpoint; only the last includes its upper endpoint. The cutoff comparison remains strict regardless of histogram bin boundaries.

Within the 73 scored guides, treatment has a higher above-cutoff rate on 58, an equal rate on 13, and a lower rate on 2. The complete paired breakdown, including guides without a corrected-only row, is [per_guide_scores.tsv](published_model_aligned_scores/results/per_guide_scores.tsv); zero-row groups are explicitly not scored.

**3. Measured-active rows, split from the two return columns.** The assay sequences were joined by `row_index` from `sites129.json`, with guide, chromosome, strand, coordinate, real PAM, positive read count and all return bits rechecked from the sequences. Scores use the deposited terminal alignment and the same `NNN` guide mask in every arm. The corrected-only count was derived from the bits, not imported from the paper.

| Released / corrected return bits | Meaning | Scored rows | Above cutoff | Percent above |
|---|---|---:|---:|---:|
| 0 / 1 | Corrected only; omitted in-scope active sites | 129 | 19 | 14.728682% |
| 1 / 1 | Both returned; measured-active control | 108 | 67 | 62.037037% |
| 0 / 0 | Neither returned; non-NGG | 23 | 13 | 56.521739% |
| 1 / 0 | Released only | 0 | 0 | Not scored |

Thus 19/129 measured-active sites omitted by the released enumerator exceed the scorer's cutoff, beside 67/108 returned active sites. Those denominators have different mismatch compositions and are not a matched enrichment comparison. All 129 omitted in-scope active sites occur in the scored treatment export, with identical aligned sequence/PAM and replay-consistent scores. Released returns 108 active rows; corrected returns 237. The other 23 are outside literal NGG scope. Their source mismatch annotation includes PAM substitutions; matching here uses protospacer mismatches recomputed from the alignment.

If 'released misses' includes out-of-scope PAMs, the total is 32/152 (21.053%), beside the returned active control 67/108 (62.037%). The non-NGG rows are scored but excluded from the terminal-omission claim.

**4. Measured-active control at the same mismatch count, and its limit.**

| Comparison | Corrected-only active sites | Released-returned active control |
|---|---:|---:|
| All available rows at four mismatches | 19/129 (14.729%) | 3/21 (14.286%) |
| Equal counts within shared guides at four mismatches | 5/17 (29.412%) | 2/17 (11.765%) |

The first row matches mismatch count but not guide composition. The second uses the 7 guides present in both active arms; within each guide it takes the smaller arm size from both arms, ordered by original assay row index before looking at scores. The complete selected identities and paired per-guide summaries are in [summary.json](published_model_aligned_scores/results/summary.json). This small selected subset does not establish general enrichment, cleavage prediction or ranking impact. In particular, the four-mismatch active rates in the first row are very similar; their score distributions are not evidence that omitted active sites are more important than retained ones.

These are in-silico nomination scores on retrospective, selected inputs. The combined-data checkpoint is unchanged; the published cutoff came from a CIRCLE-seq leave-one-guide-out sensitivity analysis, and calibration of that cutoff for this checkpoint or these inputs is not established. Prior deposited validation is reused, including its recorded benchmark discrepancies. No independence from training data, prospective validation, cleavage outcome, clinical risk, top-ranked loss or altered ranking is claimed.

**Proposed body sentence.**

> At the published CRISPR-Net cutoff, $19/129$ missed measured-active sites qualify versus $3/21$ returned sites at four mismatches; genome-wide, guide- and mismatch-matched nomination counts are $1{,}029/34{,}941$ versus $187/34{,}941$.

The sentence occupies 3 printed lines in the current body font, width and margins (limit: four). To pay for it, remove exactly the following sentence from `paper_v7/sections/03_method.tex`, in the paragraph beginning “Twenty-four records pass every control”; the removed sentence occupies 5 lines in the same layout test. Its audit detail can remain in the release documentation. No manuscript file was edited.

```tex
The tables and one command
are at \path{code_release/genome_endtoend/}, where the checker recomputes every number above
from the sequences and coordinates and a deposited script runs it on twenty-three corrupted
copies: it rejects twenty-two, among them one locus repeated $34{,}941$ times and counts moved
between subgroups so that every total still agrees, and accepts one invented reference base,
which needs the reference to catch.
```

The isolated typography check establishes the sentence length and a positive line budget; it does not establish the final manuscript pagination after insertion. See [sentence_layout.pdf](published_model_aligned_scores/layout/sentence_layout.pdf) and [layout_receipt.json](published_model_aligned_scores/layout/layout_receipt.json).

**Replay and integrity.** [The producer](published_model_aligned_scores/produce.py) reuses the byte-identical validated forward pass. [The separately coded verifier](published_model_aligned_scores/verify.py) reconstructs the matched control with SQL rather than the producer’s heap, replays the original encoder and independent recurrent/convolution implementation on every scored row, and rebuilds every summary and the rendered numerical report. It checks the [pinned input SHA-256 table](published_model_aligned_scores/INPUT_SHA256.tsv). The checkpoint SHA-256 is `a7f6aa381520f5c68fa1f099a6ef3ebc3b8ce846709b97dfde2053f26ca62f80`.

Run `published_model_aligned_scores/reproduce.sh` from this directory. All Python execution uses `nice -n 10 /usr/local/bin/python3 -I`, with at most four computational threads and observed resident memory below the required limit. Logs and exact resident-memory receipts are deposited. The original reference assembly was not re-read here: reference verification remains the producing run’s assertion, while this replay verifies the pinned exported strings and their run provenance.

**Verifier exit code: `0`.** The corruption control accepts the unmodified input with exit `0` and rejects a one-byte DNA-input change with exit `1`, before scoring. The receipt is [corruption_control/receipt.json](published_model_aligned_scores/corruption_control/receipt.json). Every treatment and control row score is available in [treatment_scores.tsv](published_model_aligned_scores/results/treatment_scores.tsv) and [control_scores.tsv](published_model_aligned_scores/results/control_scores.tsv); measured rows are in [active_scores.tsv](published_model_aligned_scores/results/active_scores.tsv).
