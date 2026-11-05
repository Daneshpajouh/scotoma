# Published checkpoint audit

## Validation result

**The NumPy forward pass passes the deposited-data behaviour check and the independent arithmetic replay.** Every data row in all four files below was scored, before scoring the target TSV. Row numbers are one-based data rows excluding the header; CSV physical line numbers are one greater. There was no row sampling, refitting, or label-based architecture selection. The CIRCLE-seq check includes training data, so its performance is a reconstruction check, not a new generalization estimate.

| Deposited file | Data rows scored, inclusive | Label 1 | Label 0 | AUROC | Trapezoidal AUPRC |
|---|---:|---:|---:|---:|---:|
| CIRCLE-seq | 1–584,949 | 7,371 | 577,578 | 0.987971998455 | 0.750611518433 |
| Six-guide indel file | 1–213,943 | 60 | 213,883 | 0.982071233650 | 0.355008180746 |
| Kleinstiver, five guides | 1–95,829 | 54 | 95,775 | 0.994556686679 | 0.319802355499 |
| Listgarten, 22 guides | 1–383,463 | 56 | 383,407 | 0.995451355497 | 0.317184732118 |

| Own labelled file | Label 0 mean | Label 1 mean | Label 0 median | Label 1 median |
|---|---:|---:|---:|---:|
| CIRCLE-seq | 0.00333737445645 | 0.366075651865 | 1.07956986797e-07 | 0.332410428612 |
| Six-guide indel file | 0.00181124531804 | 0.475291590314 | 2.6910779946e-08 | 0.572837543061 |

The [published Table 2](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.201903562) reports AUROC/AUPRC of 0.995/0.329 on the five-guide test and 0.995/0.317 on the 22-guide test for the combined CIRCLE/Elevation/SITE model. Both reproduced AUROCs and the 22-guide AUPRC agree at the paper's printed precision. **The five-guide AUPRC does not: 0.319802355499 versus 0.329.** The six README example outputs are also not all reproduced; the largest absolute difference is 0.257397676413, and every example and difference is deposited in `results/readme_examples.json`. These discrepancies remain unexplained; exact reproduction of all published output values is not claimed.

The architecture was recovered from the actual scoring JSON at [upstream commit 0843c485b1d476f899bedb93a7dd2df42547c1e1](https://github.com/JasonLinjc/CRISPR-Net/tree/0843c485b1d476f899bedb93a7dd2df42547c1e1), rather than selected by optimizing validation scores: SAME convolutions of widths 1, 2, 3, and 5 with ten filters each and ReLU, concatenation **including the original seven channels**, a 24-by-47 sequence, bidirectional LSTM with 15 units per direction and **hard-sigmoid recurrent gates**, aligned forward/backward outputs flattened to 720, ReLU dense layers of 80 and 20, and a sigmoid output. Dropout is inactive at inference. Tensor names and shapes are read directly with h5py; the forward pass uses NumPy. No historical inference framework is imported or installed.

The separately coded verifier uses a vectorized encoder, window-matrix convolutions, separate per-gate recurrent products, and a different sigmoid expression. It compares sampled validation encodings and **every target encoding** against the deposited original encoder. It replayed all 1,278,184 validation and 10,888 target scores; the largest absolute difference from the producer was 5.77315972805e-15. Validation AUROC is also recomputed with an independent rank-sum calculation. All exact AUROC fractions, labels, source-row identities and per-row scores are deposited.

## Score result and its scope

**For the literal strings in the supplied TSV, 9,046/10,888 (83.082292%) score strictly above 3/1000 (0.003).** This is the paper's stated cutoff in the sensitivity-analysis section preceding Figure 14; it was selected using CIRCLE-seq leave-one-guide-out fits, not calibrated for this combined-data checkpoint. There are 0 scores equal to that exact rational cutoff. Decisions use integer arithmetic on the unrounded binary scores, not rounded displays.

**The supplied TSV is not a set of 10,888 absent inputs.** Recomputing its flags from both complete source CSVs confirms that 10,854 words already occur as contiguous 19-base windows of an ungapped deposited protospacer for the same guide; only 34 are absent under that conservative test. Of the already represented words, 9,019 exceed the cutoff; among the 34 absent words, 27 exceed it (79.411765%). Presence is tested after removing gaps/padding and the final three PAM bases. The original table's `held_by` categories are also rebuilt, not trusted.

**These are literal-input nomination scores, not cleavage measurements or fully specified endpoint-alignment scores.** Each 20-base `spacer` and 19-base `off_target_19mer` is passed unchanged to the published encoder, which separately left-pads both to 24 columns. No PAM, flanking sequence, explicit `_` bulge, or gap placement is invented. Consequently the last three supplied spacer bases occupy columns used for PAM in the deposited labelled data; separate left padding does not reconstruct the omitted endpoint alignments. The TSV alone does not determine the full PAM-bearing, aligned inputs needed for the usual biological interpretation of the scorer. **A count of additional complete aligned nominations is therefore not established by this run.** The all-row count above must not be presented as that count.

The primary distribution uses float64 arithmetic with the deposited float32 parameters promoted unchanged. A complete float32 rerun gives 9,046 above the cutoff and changes 0 target decisions; its largest absolute score difference is 6.12971682057e-06. The independent verifier also replays that precision check. Float32 accumulation order is checked with a rounding tolerance; the float64 replay uses a much tighter tolerance.

| All literal target scores | Value |
|---|---:|
| Minimum | 4.24285812726e-10 |
| First quartile | 0.00914469474008 |
| Median | 0.260426167391 |
| Third quartile | 0.99999999999999889 |
| Maximum | 1 |
| Mean | 0.497450332814 |
| Population standard deviation | 0.477339815388 |

Quartiles use linear interpolation. The full distribution, additional quantiles, per-guide distributions and presence-group distributions are in `results/summary.json`.

| Score interval | Rows |
|---|---:|
| [0, 0.000001) | 68 |
| [0.000001, 0.0001) | 406 |
| [0.0001, 0.001) | 748 |
| [0.001, 0.003) | 620 |
| [0.003, 0.01) | 967 |
| [0.01, 0.1) | 2,059 |
| [0.1, 0.5) | 831 |
| [0.5, 1] | 5,189 |

## Deposit, replay and corruption control

The deposit is [published_model_scores](published_model_scores/README.md). The producer is `produce.py`; the separate forward implementation and complete replay are in `verify.py`. The producer records completion of all validation files and the validation gate before any target score is produced. All inputs, the architecture, original encoder, evaluation code and source receipts are included; reruns need no network. The two Table 2 test CSVs and the architecture were retrieved at the pinned commit and checked against the Git blob identities already recorded by the earlier audit.

Weight SHA-256:

```text
a7f6aa381520f5c68fa1f099a6ef3ebc3b8ce846709b97dfde2053f26ca62f80
```

The weight's Git blob identity is `ab49d51f3ea99e17478768859ac0fe9f05a657fb`. The fixed input/source SHA-256 table is `pins.json`; its own SHA-256 is embedded in both producer and verifier. Each output is bound by `results/outputs.sha256.json` and independently replayed. All report measurements are rendered from those checked results by `render_findings.py`.

Run from this directory:

```sh
nice -n 10 /usr/local/bin/python3 -I -B published_model_scores/produce.py
nice -n 10 /usr/local/bin/python3 -I -B published_model_scores/verify.py
nice -n 10 /usr/local/bin/python3 -I -B published_model_scores/corruption_control.py
```

Threads are limited to four. Exact package versions are recorded in `environment.json`. The recorded producer/verifier/control logs are in `logs/`.

The corruption control passed for 2 inputs: the target TSV and the weight HDF5. In each case an equal-length copy with **exactly one changed byte** was supplied to the same verifier, which returned exit status 1 and a named SHA-256 rejection. The corrupted copies and full receipts are deposited in `corruption_control/`; `corruption_control.py` recreates and checks them.

## Proposed body sentence and space trade

This sentence reports only the literal-input measurement, with no claim of reconstructed aligned nominations. It is six short lines:

```text
After validation on the two deposited labelled files
(AUROC 0.987972 and 0.982071), literal PAM-free guide/19-mer
inputs to CRISPR-Net exceeded the paper's 0.003 cutoff
for 9,046/10,888 pairs, including 27/34 words absent under
the deposited protospacer substring test; these are in-silico
nomination scores for the supplied inputs.
```

Pay for it by replacing the read-depth comparison span in the existing **Genomic annotations** paragraph of `paper_v7/sections/03_method.tex`, beginning “A single-implementation read-depth comparison” and ending at the four-mismatch-control comparison (`p=0.334`). The manuscript was read only; no manuscript or other directory was modified.
