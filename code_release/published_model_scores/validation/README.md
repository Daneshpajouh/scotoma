# Deposited checkpoint scores

This package reconstructs the deposited CRISPR-Net scoring checkpoint using NumPy and h5py. It validates on every row of both supplied labelled CSVs and both published Table 2 test CSVs before scoring the target table. See [FINDINGS.md](../FINDINGS.md) for the measured results and the unresolved published-output discrepancies.

The target result concerns the literal guide/19-mer strings passed to the published encoder. The supplied table contains neither PAMs nor explicit gap alignments. No missing context is filled in. These scores cannot establish a count of additional complete aligned nominations. Its conservative presence flags are independently rebuilt from the complete labelled CSVs.

Run from the parent task directory:

```sh
nice -n 10 /usr/local/bin/python3 -I -B published_model_scores/produce.py
nice -n 10 /usr/local/bin/python3 -I -B published_model_scores/verify.py
nice -n 10 /usr/local/bin/python3 -I -B published_model_scores/corruption_control.py
```

For another producer replay without overwriting the recorded outputs, use `--out published_model_scores/replay_results`; pass that path as `--results` to the verifier. All producer output paths must remain within the parent task directory. The scripts set a four-thread maximum and suppress bytecode writes. There is no network access in the scripts.

Files:

- `model_numpy.py`: primary forward pass, stored parameters promoted to float64.
- `produce.py`: validates first, records all row scores, rebuilds presence flags, then scores targets and runs a float32 check.
- `verify.py`: independently coded forward pass and summaries; replays every validation and target score, all reported derived quantities, the original encoder, source hashes, and the precision check.
- `corruption_control.py`: creates and tests a one-byte-corrupted target file and weight file. `corruption_control/receipt.json` records both rejections.
- `pins.json`: exact SHA-256 pins for the source material and all inputs. Its digest is embedded in producer and verifier.
- `inputs/`: complete immutable copies of the supplied data/weights/target table and the two published benchmark CSVs.
- `sources/`: pinned upstream architecture, encoder, data encoder, inference/evaluation code and source receipts. These upstream files are preserved as evidence; the historical inference code is never imported.
- `results/`: all scores, exact row indices, metrics, distributions, presence audit, precision check and execution gate trace.
- `environment.json`: measured runtime/dependency versions.
- `logs/`: executed producer, full verifier and corruption-control logs.
- `render_findings.py`: deterministic rendering of the parent findings from the verified results.

The producer needs NumPy and h5py. The verifier additionally uses SciPy for its independent AUROC rank calculation. The measured versions are NumPy 2.5.1, h5py 3.15.1 and SciPy 1.17.1. No fitting or historical inference-framework installation is performed.

The decision threshold is the exact rational 3/1000 with a strict comparison. It comes from the paper's CIRCLE-seq cross-validation sensitivity analysis and is not separately calibrated for this combined-data checkpoint. Source-row numbers exclude CSV headers. Every validation file is scored in full; none is sampled.

Weight SHA-256: `a7f6aa381520f5c68fa1f099a6ef3ebc3b8ce846709b97dfde2053f26ca62f80`.

All historical sources are pinned to repository commit `0843c485b1d476f899bedb93a7dd2df42547c1e1`. Published comparison: [Lin et al., Advanced Science 2020, Table 2 and sensitivity-analysis section](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.201903562).
