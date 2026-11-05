# The population law, executed

For single RNA-bulge expansion the 2016 `cas-offinder-bulge` route emits a complete word set
exactly when `g[0] == g[1]` and `g[-1] == g[-2]`. Both comparisons are on the guide alone: no
genome, no scan, no search. On a uniform four-letter guide population they hold together on
`1/16`, so the route is **incomplete on 15/16 of all guides, at every length `n >= 4`**.

The point of this directory is the evidence class. The fraction is a one-line consequence of the
deficit formula in the paper, but a consequence of a formula is a claim about our reading of the
generator. Here the generator itself is run.

## What is executed

`execute_population_law.py <n>` runs the vendored 2016 source, unmodified, in its own dry-run
mode, once per guide, over all `4^n` guides, and reads back the pattern file it writes. For each
guide it compares the emitted word set against the required set computed from the definition
(delete one base, keep distinct results), and against the two-character predicate.

The vendored source is `../conformance/upstream/cas-offinder-bulge`, v1.2 beta (2016-07-06). It is
the generator CRISPRseek reaches through backend 2.4.1.

Deposited records, each written by that script:

| n | guides | short | fraction | executed disagreements with the predicate |
|---|---|---|---|---|
| 4 | 256 | 240 | 15/16 | 0 |
| 5 | 1,024 | 960 | 15/16 | 0 |
| 6 | 4,096 | 3,840 | 15/16 | 0 |

5,376 runs of the shipped generator. Every emitted word lay inside the required set, and no guide
was found where the executed completeness differed from the two-character predicate.

## What is not executed here

The `3.0.0b3`/develop native parser's `3/4` is read from its source loop. It is not run in this
directory and is not claimed as an execution. The paper labels it the same way.

Lengths above 6 are stated in the paper from the formula, not from a run. `n = 7` would be 16,384
invocations of the generator; it was not performed, and the paper does not say it was.

## Replay

    python3 -I execute_population_law.py 4      # 256 runs, seconds
    python3 -I execute_population_law.py 6      # 4,096 runs, minutes
    python3 -I verify_population_law.py         # recompute the predicate count; must exit 0
    python3 -I verify_population_law.py --control

The checker never reads a stated fraction. For each record it re-enumerates the whole guide space,
counts the guides the predicate calls incomplete, and requires the executed count to equal it
exactly, the parts to sum to `4^n`, and the record to report zero disagreements.

`--control` removes a single short guide from one record and requires the check to fail. A checker
that cannot fail is not a check, so the control is part of the replay and not an extra.

    shasum -a 256 -c MANIFEST.sha256         # this directory is self-contained

The manifest covers every file here. Re-running the producer overwrites its own record, so verify
the manifest before a replay if you want to compare against what was deposited.
