> Paths in this file are relative to the directory it sits in. It was written when this
> package was staged one level deeper and its links carried a `conformance/` prefix; the
> prefix is removed here so every link resolves from the deposit.

# Findings

The completed deposit is [conformance/](README.md). The standalone producer, suite, minimal wrapper diff, planted-regression self-test, evidence snapshots and run logs are all within this directory. No paper or upstream source was edited outside this directory.

The following commands were executed in `conformance/`, under `nice` with thread caps. Every numerical result below comes from these runs or their recorded child commands; [the ledger](evidence/runs.jsonl) contains the exact absolute argv, working directory, stdout/stderr and exit code. Names and version identifiers follow the Table 1 source snapshot.

| Run | Command in deposit | Exit | Output |
|---|---|---:|---|
| `scoreboard-r3` | `/usr/local/bin/python3 -I scoreboard.py` | `0` | [evidence/scoreboard-r3.stdout](evidence/scoreboard-r3.stdout) |
| `validate-patch-r2` | `/usr/local/bin/python3 -I validate_patch.py` | `0` | [evidence/validate-patch-r2.stdout](evidence/validate-patch-r2.stdout) |
| `selftest-r2` | `/usr/local/bin/python3 -I selftest.py` | `0` | [evidence/selftest-r2.stdout](evidence/selftest-r2.stdout) |
| `unit-tests-r3` | `/usr/local/bin/python3 -I test_conformance.py` | `0` | [evidence/unit-tests-r3.stdout](evidence/unit-tests-r3.stdout) |
| `check-sentence-r3` | `/usr/local/bin/python3 -I check_sentence.py` | `0` | [evidence/check-sentence-r3.stdout](evidence/check-sentence-r3.stdout) |

## Scoreboard

The panel requires 1634 words over 110 guides. The suite rescored actual wrapper files, parser map keys and raw fixture alignments, and cross-checked 1160 deposited summary rows against the derived values. Its [inventory](inputs/output_inventory.json) records why each genomic/clinical output or count-only artifact is inapplicable as an emitted pattern set.

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

Executed/static is Table 1's evidence state, not a promotion made by this run. Fresh native engines, R calls and untested pipelines were not run here. In particular, the early native betas remain static loop emulations; CHOPOFF has source-component count evidence but no deposited native emitted word set; CRISPRme delegates to separately tested engines; SWOffinder and CALITAS have no deletion-library count verdict. CRISPRseek's pooled written file is checked for equality with the wrapper union, and per-guide ownership comes from the executed wrapper snapshots, not from crediting a shared word to every possible owner.

CRISPRitz has exact projected word sets and redundant alignment rows. Native parser dumps carry placement metadata, which is collapsed to their unique map keys; it is not counted as repeated emitted patterns. These distinctions are retained in the [machine-readable scoreboard](results/scoreboard.json), with every missing word in per-guide reports.

## Patch and regression test

The [unified diff](wrapper.patch) modifies 2 source lines, one RNA loop bound per PAM orientation: 2 removed and 2 added, or 20 physical diff lines including context and headers. Both independently guarded branches omit endpoints, so each requires a bound repair; existing dictionary deduplication makes further changes unnecessary for word-set exactness. This establishes the smallest local bound repair justified here, not global byte-minimality.

For each PAM orientation, fresh official dry-runs change guide-specific emission from 1491 to 1634, missing words from 143 to 0, and incomplete guides from 98 to 0. On the deposited planted fixture, pooled emission changes from 1490 to 1629 of 1629, with pooled missing words changing from 139 to 0. The patched sets have 0 invalid words and 0 redundant rows. These are pattern-generation measurements, not a new fixture search by the native engine; [patch validation](results/patch_validation.json) records both orientations and every child exit.

The [self-test](results/selftest.json) actually mutates `required_words` to discard the terminal word. Its mutant producer exits `0` with 19 words while retaining count 20; the suite exits `1` and names `GAGTAGCGCGAGCACAGCT` as missing. The positive suite exits `0`, and the self-test harness exits `0` only after requiring both outcomes. The separate test run exits `0` with 13 tests passing, including an independent all-position oracle on 5460 guides through length 6.

## Body sentence and its space payment

Proposed sentence, naming the actual deposit:

```tex
The emitter, suite and patch are in
\path{HOSTILE_REVIEW_2026-11-01/R19_conformance_standard/conformance/}.
```

Replace this existing sentence in `paper_v7/sections/03_method.tex`:

> No audited route emits the minimum family and nothing else: two fall below its word set and two attain it but serialise it in more rows.

The new deposit supplies the scoreboard and exact-emitter control that make the old audit-only summary redundant. In the paper's actual typography, the new sentence prints in 2 lines and the removed sentence in 2; local full-paper copies compile with 10 body pages before and 10 after. Commands and exits are above; [the layout report](results/sentence_check.json), context PDF (a local layout-check build, not deposited) and rendered context are deposited. No manuscript edit has been applied.

Earlier layout attempts are retained: the initial inline driver exited `1` on timeout ([record](evidence/layout_initial_timeout.json)); a subsequent longer sentence was rejected by a layout check exiting `1` ([rejected proposal](results/rejected_sentence.json)). A context-render lookup also exited `1` while still searching for the longer draft's wording; the final lookup uses the final sentence. The final sentence and full-paper build pass the checks above. The initial scorer also exposed repeated parser metadata as if it were emission; final reports correctly use pattern-map keys, with earlier stdout retained for provenance.

## Scope

The linear-time run count is a necessary screen. Exact sound-family conformance additionally requires word membership/equality; generating and comparing all materialized words is not claimed to be linear. The Hamming-budget and soundness hypotheses are explicit in the producer and README. This deposit concerns single-deletion pattern completeness only: no cleavage, risk, ranking-impact, tool-safety, whole-search certification or adoption claim is made.
