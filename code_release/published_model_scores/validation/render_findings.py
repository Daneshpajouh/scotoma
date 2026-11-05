#!/usr/bin/env python3
"""Render the findings from the numbers that the verifier replays."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def render():
    v=json.loads((HERE/'results/validation.json').read_text())
    s=json.loads((HERE/'results/summary.json').read_text())
    presence=json.loads((HERE/'results/presence_audit.json').read_text())
    precision=json.loads((HERE/'results/precision_check.json').read_text())
    examples=json.loads((HERE/'results/readme_examples.json').read_text())
    control=json.loads((HERE/'corruption_control/receipt.json').read_text())
    log=(HERE/'logs/verifier.log').read_text().splitlines()
    receipt=next(json.loads(line) for line in reversed(log) if line.startswith('{') and '"status": "PASS"' in line)
    names=[('CIRCLE-seq','CIRCLE_seq_10gRNA_wholeDataset.csv'),
           ('Six-guide indel file','elevation_6gRNA_wholeDataset.csv'),
           ('Kleinstiver, five guides','Kleinstiver_5gRNA_wholeDataset.csv'),
           ('Listgarten, 22 guides','Listgarten_22gRNA_wholeDataset.csv')]
    rows=[]
    for label,name in names:
        d=v[name]
        rows.append(f"| {label} | 1–{d['rows']:,} | {d['positive']:,} | {d['negative']:,} | {d['auroc']:.12f} | {d['auprc_trapezoidal']:.12f} |")
    means=[]
    for label,name in names[:2]:
        d=v[name]
        means.append(f"| {label} | {d['label_0_mean']:.12g} | {d['label_1_mean']:.12g} | {d['label_0_median']:.12g} | {d['label_1_median']:.12g} |")
    a=s['all'];absent=s['by_in_published_files']['no'];held=s['by_in_published_files']['yes']
    q=a['quantiles_linear']
    histogram=[]
    for b in a['histogram']:
        interval=f"[{b['lower']}, {b['upper']}{']' if b['upper_inclusive'] else ')'}"
        histogram.append(f"| {interval} | {b['rows']:,} |")
    cutoff=s['cutoff']['exact']
    proposal=(f"After validation on the two deposited labelled files\n"
              f"(AUROC {v[names[0][1]]['auroc']:.6f} and {v[names[1][1]]['auroc']:.6f}), literal PAM-free guide/19-mer\n"
              f"inputs to CRISPR-Net exceeded the paper's 0.003 cutoff\n"
              f"for {a['above_3_over_1000']:,}/{a['rows']:,} pairs, including {absent['above_3_over_1000']}/{absent['rows']} words absent under\n"
              f"the deposited protospacer substring test; these are in-silico\n"
              f"nomination scores for the supplied inputs.")
    assert len(proposal.splitlines())==6
    largest_example=max(e['absolute_difference'] for e in examples)
    return f"""# Published checkpoint audit

## Validation result

**The NumPy forward pass passes the deposited-data behaviour check and the independent arithmetic replay.** Every data row in all four files below was scored, before scoring the target TSV. Row numbers are one-based data rows excluding the header; CSV physical line numbers are one greater. There was no row sampling, refitting, or label-based architecture selection. The CIRCLE-seq check includes training data, so its performance is a reconstruction check, not a new generalization estimate.

| Deposited file | Data rows scored, inclusive | Label 1 | Label 0 | AUROC | Trapezoidal AUPRC |
|---|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

| Own labelled file | Label 0 mean | Label 1 mean | Label 0 median | Label 1 median |
|---|---:|---:|---:|---:|
{chr(10).join(means)}

The [published Table 2](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.201903562) reports AUROC/AUPRC of 0.995/0.329 on the five-guide test and 0.995/0.317 on the 22-guide test for the combined CIRCLE/Elevation/SITE model. Both reproduced AUROCs and the 22-guide AUPRC agree at the paper's printed precision. **The five-guide AUPRC does not: {v[names[2][1]]['auprc_trapezoidal']:.12f} versus 0.329.** The six README example outputs are also not all reproduced; the largest absolute difference is {largest_example:.12g}, and every example and difference is deposited in `results/readme_examples.json`. These discrepancies remain unexplained; exact reproduction of all published output values is not claimed.

The architecture was recovered from the actual scoring JSON at [upstream commit 0843c485b1d476f899bedb93a7dd2df42547c1e1](https://github.com/JasonLinjc/CRISPR-Net/tree/0843c485b1d476f899bedb93a7dd2df42547c1e1), rather than selected by optimizing validation scores: SAME convolutions of widths 1, 2, 3, and 5 with ten filters each and ReLU, concatenation **including the original seven channels**, a 24-by-47 sequence, bidirectional LSTM with 15 units per direction and **hard-sigmoid recurrent gates**, aligned forward/backward outputs flattened to 720, ReLU dense layers of 80 and 20, and a sigmoid output. Dropout is inactive at inference. Tensor names and shapes are read directly with h5py; the forward pass uses NumPy. No historical inference framework is imported or installed.

The separately coded verifier uses a vectorized encoder, window-matrix convolutions, separate per-gate recurrent products, and a different sigmoid expression. It compares sampled validation encodings and **every target encoding** against the deposited original encoder. It replayed all {receipt['validation_rows_replayed']:,} validation and {receipt['target_rows_replayed']:,} target scores; the largest absolute difference from the producer was {receipt['largest_absolute_independent_score_difference']:.12g}. Validation AUROC is also recomputed with an independent rank-sum calculation. All exact AUROC fractions, labels, source-row identities and per-row scores are deposited.

## Score result and its scope

**For the literal strings in the supplied TSV, {a['above_3_over_1000']:,}/{a['rows']:,} ({a['percent_above']:.6f}%) score strictly above {cutoff} (0.003).** This is the paper's stated cutoff in the sensitivity-analysis section preceding Figure 14; it was selected using CIRCLE-seq leave-one-guide-out fits, not calibrated for this combined-data checkpoint. There are {a['equal_3_over_1000']} scores equal to that exact rational cutoff. Decisions use integer arithmetic on the unrounded binary scores, not rounded displays.

**The supplied TSV is not a set of {a['rows']:,} absent inputs.** Recomputing its flags from both complete source CSVs confirms that {held['rows']:,} words already occur as contiguous 19-base windows of an ungapped deposited protospacer for the same guide; only {absent['rows']} are absent under that conservative test. Of the already represented words, {held['above_3_over_1000']:,} exceed the cutoff; among the {absent['rows']} absent words, {absent['above_3_over_1000']} exceed it ({absent['percent_above']:.6f}%). Presence is tested after removing gaps/padding and the final three PAM bases. The original table's `held_by` categories are also rebuilt, not trusted.

**These are literal-input nomination scores, not cleavage measurements or fully specified endpoint-alignment scores.** Each 20-base `spacer` and 19-base `off_target_19mer` is passed unchanged to the published encoder, which separately left-pads both to 24 columns. No PAM, flanking sequence, explicit `_` bulge, or gap placement is invented. Consequently the last three supplied spacer bases occupy columns used for PAM in the deposited labelled data; separate left padding does not reconstruct the omitted endpoint alignments. The TSV alone does not determine the full PAM-bearing, aligned inputs needed for the usual biological interpretation of the scorer. **A count of additional complete aligned nominations is therefore not established by this run.** The all-row count above must not be presented as that count.

The primary distribution uses float64 arithmetic with the deposited float32 parameters promoted unchanged. A complete float32 rerun gives {precision['float32_above']:,} above the cutoff and changes {precision['changed_cutoff_decisions']} target decisions; its largest absolute score difference is {precision['maximum_absolute_score_difference']:.12g}. The independent verifier also replays that precision check. Float32 accumulation order is checked with a rounding tolerance; the float64 replay uses a much tighter tolerance.

| All literal target scores | Value |
|---|---:|
| Minimum | {q['0']:.12g} |
| First quartile | {q['0.25']:.12g} |
| Median | {q['0.5']:.12g} |
| Third quartile | {q['0.75']:.17g} |
| Maximum | {q['1']:.12g} |
| Mean | {a['mean']:.12g} |
| Population standard deviation | {a['population_sd']:.12g} |

Quartiles use linear interpolation. The full distribution, additional quantiles, per-guide distributions and presence-group distributions are in `results/summary.json`.

| Score interval | Rows |
|---|---:|
{chr(10).join(histogram)}

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

The corruption control passed for {len(control['tests'])} inputs: the target TSV and the weight HDF5. In each case an equal-length copy with **exactly one changed byte** was supplied to the same verifier, which returned exit status {control['tests'][0]['verifier_exit_status']} and a named SHA-256 rejection. The corrupted copies and full receipts are deposited in `corruption_control/`; `corruption_control.py` recreates and checks them.

## Proposed body sentence and space trade

This sentence reports only the literal-input measurement, with no claim of reconstructed aligned nominations. It is six short lines:

```text
{proposal}
```

Pay for it by replacing the read-depth comparison span in the existing **Genomic annotations** paragraph of `paper_v7/sections/03_method.tex`, beginning “A single-implementation read-depth comparison” and ending at the four-mismatch-control comparison (`p=0.334`). The manuscript was read only; no manuscript or other directory was modified.
"""


if __name__=='__main__':
    (HERE.parent/'FINDINGS.md').write_text(render())
