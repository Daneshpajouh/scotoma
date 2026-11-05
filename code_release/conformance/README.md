# Single-deletion pattern conformance

For an uppercase ACGT guide `g` of length `n`, deleting each base produces exactly `r = 1 + sum(g[i] != g[i-1] for i in range(1,n))` distinct words, where `r` is the run count and is computable in `O(n)` without building words. At a nonnegative integer Hamming budget `m <= n-3`, the unique minimum sound family is this deletion set, even against wildcard competitors, so a sound emitter with fewer than `r` distinct words cannot be complete. The count is a rejection screen until membership in the sound family is established; this suite checks exact word-set equality and reports repeated rows separately.

## Run the check

From this directory, this one-line check on the patched panel must pass:

```sh
nice -n 10 /usr/local/bin/python3 -I suite.py --guides inputs/guides.tsv --emitted results/patch_trailing_after.emitted.json
```

For your tool, replace `--emitted` with your emitted words. Use a JSON object mapping guide IDs to lists of strings, or a TSV with exactly the header `guide_id<TAB>pattern`; the guide TSV header is `guide_id<TAB>guide`. Missing guide IDs mean empty emission; unknown IDs, duplicate guide IDs and malformed input are rejected. For a single guide, use `--guide SPACER --emitted words.json`, with a JSON list of words. The guide and words contain no PAM, gap symbols or external padding; keep repeated rows in the list if you want redundancy measured.

For an alignment engine, use the bundled `inputs/fixture_rna.fa` and `inputs/guides.tsv` at zero mismatches, zero DNA bulges and one RNA bulge, then strip the fixed PAM and gap symbols from targets while retaining guide ownership. Positive-mismatch genomic targets are not emitted pattern centres and are unsuitable for this projection. The deposited alignment adapters validate gap placement, guide identity, mismatch budget and PAM before comparing the resulting words.

`exact` means every required word appears once; `incomplete` lists every missing word; `redundant` reports repeated rows of a complete family; `invalid` lists every word outside the required family. Invalidity takes priority over incompleteness, and incompleteness over redundancy, while all axes remain in the JSON report. Exit codes are `0` for exact, `1` for a nonexact family, and `2` for malformed input. A count at least `r`, including a raw row count, cannot repair a missing word.

`reference.py` is the standalone producer: `run_count(g)` constructs no deletion words, `required_words(g)` deletes once per run, and `required(g)` returns `(r, word_set)`. Its CLI prints JSON; `--count-only` uses only the linear count. Materializing all words costs `O(n*r)` output space and time; the generate-and-diff suite is not claimed to run in linear time. Both programs use only the standard library and work under isolated Python invocation.

## Scoreboard

The panel has 110 guides requiring 1634 guide-specific words. These are rescored deposited outputs; executed/static labels retain the evidence column of Paper A's Table 1, copied to [the authority snapshot](inputs/authority/03_method.tex). “Executed” does not mean a fresh native search was run here.

| Route | Evidence from Table 1 | Distinct / required | Missing (guides) | Repeated word rows | Verdict |
|---|---|---:|---:|---:|---|
| Cas-OFFinder wrapper 2016 (wrapper) | executed | 1491 / 1634 | 143 (98) | 0 | incomplete |
| Cas-OFFinder 3.0.0b3 (native_b3) | executed | 1561 / 1634 | 73 (73) | 0 | incomplete |
| Cas-OFFinder develop (native_develop) | executed | 1561 / 1634 | 73 (73) | 0 | incomplete |
| Cas-OFFinder 3.0.0b / 3.0.0b2 (native_early_static) | static | 1491 / 1634 | 143 (98) | 0 | incomplete |
| CRISPRitz 2.7.1 (crispritz_2.7.1) | executed | 1634 / 1634 | 0 (0) | 566 | redundant; exact word set |
| CRISPRitz 2.9.0 (crispritz_2.9.0) | executed | 1634 / 1634 | 0 (0) | 566 | redundant; exact word set |
| Cas-OFFinder 2 Rust 1.0.0 (rust_cli) | executed | 1561 / 1634 | 73 (73) | 529 | incomplete |
| Cas-OFFinder 2 Rust 1.0.0 (rust_forward) | executed | 1561 / 1634 | 73 (73) | 529 | incomplete |
| Cas-OFFinder 2 Rust 1.0.0 (rust_reverse) | executed | 1561 / 1634 | 73 (73) | 529 | incomplete |
| CRISPRseek 1.52.0 default bulge route (crisprseek) | executed | 1491 / 1634 | 143 (98) | 0 | incomplete |
| CRISPRme 2.1.14 / CRISPRme+ 2.8.0 (source) | static | N/A | N/A | N/A | not certifiable |
| CHOPOFF 1.0.3 (source) | static | N/A | N/A | N/A | not certifiable |
| SWOffinder 1.0.1 (source) | static | N/A | N/A | N/A | not certifiable |
| CALITAS 7a3456c (source) | static | N/A | N/A | N/A | not certifiable |

All numerical cells above: [`scoreboard-r3`](evidence/scoreboard-r3.stdout), exit `0`. Full per-guide missing words, invalid words, multiplicities and scope are in [the results](results/scoreboard.json). Native parser map keys are deduplicated independently of their placement metadata; repeated alignment rows remain visible for the fixture searches. CHOPOFF's source-component count is checked, but a count summary cannot certify an unavailable native emitted set. CRISPRme pipelines and alignment-only architectures receive no invented word-set verdict.

[The output inventory](inputs/output_inventory.json) accounts for every deposited text table/output in the requested evidence directories. Genomic and clinical outputs are not pattern libraries or exhaustive planted universes, so they are marked inapplicable to this test rather than classified from mismatching genomic targets.

## Apply the patch

[wrapper.patch](wrapper.patch) is a unified diff against [the untouched vendored source](upstream/cas-offinder-bulge). Apply it from an upstream checkout containing that file with `patch -p1 < /absolute/path/to/wrapper.patch`, or run the local validation below. It replaces the RNA deletion bound in each PAM-orientation branch:

```diff
-for i in range(1, len(target)-bulge_size):
+for i in range(len(target)-bulge_size+1):
```

The patch modifies 2 source lines (2 removed, 2 added; 20 physical unified-diff lines with headers and context). Each independent branch needs its endpoint-excluding bound repaired; the existing dictionary already deduplicates words. This is the smallest local loop-bound repair justified here, not a claim of globally minimal diff encoding.

| Tested emission, separately for each PAM orientation | Before | After |
|---|---:|---:|
| Guide-specific words | 1491 | 1634 |
| Missing guide-specific words | 143 | 0 |
| Guides incomplete | 98 | 0 |
| Pooled fixture words | 1490 | 1629 |
| Missing pooled fixture words | 139 | 0 |

Patch measurements: [`validate-patch-r2`](evidence/validate-patch-r2.stdout), exit `0`; detailed [patch validation](results/patch_validation.json). Both orientations are tested at a fixed substitution budget, using the real wrapper's dry-run mode, per guide and pooled. The fixture has 1629 pooled required words; pooled counts must not be substituted for per-guide counts. The patch concerns RNA deletion generation; full wrapper search/postprocessing and DNA bulges are outside this validation.

## Reproduce and test a regression

```sh
nice -n 10 /usr/local/bin/python3 -I reproduce.py
nice -n 10 /usr/local/bin/python3 -I selftest.py
```

The replay needs Python and the system `patch` command, uses only bundled inputs, and writes only here. It rebuilds the scoreboard and diff, applies the diff to a disposable copy, executes the original and patched dry-runs, and runs the suite and tests. Add `--manuscript` to also compile and render local paper copies; that optional check needs TeX and Poppler at the paths in `check_sentence.py`. Original inputs can be recollected from this repository with `collect_inputs.py`; recollection is not needed for replay.

The self-test mutates the actual producer to remove `g[:-1]`, rather than supplying a hand-written bad list. The producer still reports `r=20`, but emits 19 words; the unmodified CLI suite returns `1` and identifies `GAGTAGCGCGAGCACAGCT`. The positive control returns `0`. The self-test exits successfully only when those outcomes occur: [`selftest-r2`](evidence/selftest-r2.stdout), exit `0`. The test suite passes 13 tests, including 5460 independently enumerated guides through length 6, terminal-run collapse, invalid replacements preserving counts, duplicate padding, per-guide ownership, malformed input and corrupted raw alignments: [`unit-tests-r3`](evidence/unit-tests-r3.stdout), exit `0`.

## Claim boundary

This is conformance of a single-deletion word family at a stated Hamming budget. The source theorem's sound-family result applies through `m <= n-3`; Paper A's unrestricted concrete-centre floor is separately scoped to `n >= 20` and `m <= 5`, with applicability flags in the producer. Beyond those hypotheses the program still checks deletion-set equality and makes no minimum-cover claim. An invalid word means outside this specified family; it is not an assertion of unsoundness outside the theorem's range.

Passing does not establish full search-engine completeness at other budgets, PAM policies, strands, references or devices. This deposit makes no claim about cleavage, biological risk, ranking impact, tool safety, upstream acceptance or adoption. Fresh work here is local pattern generation and conformance checking; no genome-scale search is performed.

Exact commands, working directories, script hashes, stdout/stderr paths and exit codes are in [the run ledger](evidence/runs.jsonl). Source/snapshot hashes are in [the manifest](inputs/manifest.json). Earlier diagnostic runs remain in the ledger; use the final successful stage named above for the published measurements.

## A note on the evidence files

`evidence/*.stdout` are the recorded child runs, and the runner writes each command's absolute
argv into them. In this copy the producing machine's prefix is rewritten to `<repo>/` and its home
directory to `<home>`; nothing else in them changed. Re-running any check here will rewrite those
files with **your** absolute paths. That is the log, not the evidence: the evidence is
`results/*.json` and the exit codes, which are path-independent.
