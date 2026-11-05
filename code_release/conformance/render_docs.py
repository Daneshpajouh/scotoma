#!/usr/local/bin/python3 -I
"""Render the README and findings directly from executed results and their run records."""
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE / "results" / name).read_text())


def main():
    summary = load("panel_summary.json")
    patch = load("patch_validation.json")
    mutation = load("selftest.json")
    sentence = load("sentence_check.json")
    runs = [json.loads(line) for line in (HERE / "evidence/runs.jsonl").read_text().splitlines()]

    def latest(label):
        return next(row for row in reversed(runs) if row.get("requested_label", row["label"]) == label)

    score_run, patch_run, self_run, tests_run, sentence_run = [latest(label) for label in
        ("scoreboard", "validate-patch", "selftest", "unit-tests", "check-sentence")]
    for row in (score_run, patch_run, self_run, tests_run, sentence_run):
        if row["exit_code"] != 0:
            raise ValueError("cannot document a failed acceptance stage as passing")
    test_stderr = (HERE / tests_run["stderr"]).read_text()
    tests = int(re.search(r"Ran (\d+) tests", test_stderr).group(1))
    stdout = (HERE / tests_run["stdout"]).read_text()
    exhaustive = json.loads(next(line for line in stdout.splitlines() if line.startswith('{"exhaustive_guides"')))
    scoreboard = (HERE / "results/SCOREBOARD.md").read_text().split("\n\nScopes", 1)[0]
    before, after = patch["runs"]["trailing_before"], patch["runs"]["trailing_after"]
    mutant = mutation["cases"]["terminal-dropped"]
    missing = mutant["suite_report"]["per_guide"]["guide"]["missing_words"][0]

    def provenance(row):
        return f"[`{row['label']}`]({row['stdout']}), exit `{row['exit_code']}`"

    readme = f"""# Single-deletion pattern conformance

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

The panel has {summary['guides']} guides requiring {summary['required_guide_specific_words']} guide-specific words. These are rescored deposited outputs; executed/static labels retain the evidence column of Paper A's Table 1, copied to [the authority snapshot](inputs/authority/03_method.tex). “Executed” does not mean a fresh native search was run here.

{scoreboard}

All numerical cells above: {provenance(score_run)}. Full per-guide missing words, invalid words, multiplicities and scope are in [the results](results/scoreboard.json). Native parser map keys are deduplicated independently of their placement metadata; repeated alignment rows remain visible for the fixture searches. CHOPOFF's source-component count is checked, but a count summary cannot certify an unavailable native emitted set. CRISPRme pipelines and alignment-only architectures receive no invented word-set verdict.

[The output inventory](inputs/output_inventory.json) accounts for every deposited text table/output in the requested evidence directories. Genomic and clinical outputs are not pattern libraries or exhaustive planted universes, so they are marked inapplicable to this test rather than classified from mismatching genomic targets.

## Apply the patch

[wrapper.patch](wrapper.patch) is a unified diff against [the untouched vendored source](upstream/cas-offinder-bulge). Apply it from an upstream checkout containing that file with `patch -p1 < /absolute/path/to/wrapper.patch`, or run the local validation below. It replaces the RNA deletion bound in each PAM-orientation branch:

```diff
-for i in range(1, len(target)-bulge_size):
+for i in range(len(target)-bulge_size+1):
```

The patch modifies {patch['modified_source_lines']} source lines ({patch['removed_lines']} removed, {patch['added_lines']} added; {patch['unified_diff_physical_lines']} physical unified-diff lines with headers and context). Each independent branch needs its endpoint-excluding bound repaired; the existing dictionary already deduplicates words. This is the smallest local loop-bound repair justified here, not a claim of globally minimal diff encoding.

| Tested emission, separately for each PAM orientation | Before | After |
|---|---:|---:|
| Guide-specific words | {before['valid_distinct_words']} | {after['valid_distinct_words']} |
| Missing guide-specific words | {before['missing_words']} | {after['missing_words']} |
| Guides incomplete | {before['guides_incomplete']} | {after['guides_incomplete']} |
| Pooled fixture words | {before['fixture_pooled_emitted']} | {after['fixture_pooled_emitted']} |
| Missing pooled fixture words | {before['fixture_pooled_missing']} | {after['fixture_pooled_missing']} |

Patch measurements: {provenance(patch_run)}; detailed [patch validation](results/patch_validation.json). Both orientations are tested at a fixed substitution budget, using the real wrapper's dry-run mode, per guide and pooled. The fixture has {summary['pooled_required_words']} pooled required words; pooled counts must not be substituted for per-guide counts. The patch concerns RNA deletion generation; full wrapper search/postprocessing and DNA bulges are outside this validation.

## Reproduce and test a regression

```sh
nice -n 10 /usr/local/bin/python3 -I reproduce.py
nice -n 10 /usr/local/bin/python3 -I selftest.py
```

The replay needs Python and the system `patch` command, uses only bundled inputs, and writes only here. It rebuilds the scoreboard and diff, applies the diff to a disposable copy, executes the original and patched dry-runs, and runs the suite and tests. Add `--manuscript` to also compile and render local paper copies; that optional check needs TeX and Poppler at the paths in `check_sentence.py`. Original inputs can be recollected from this repository with `collect_inputs.py`; recollection is not needed for replay.

The self-test mutates the actual producer to remove `g[:-1]`, rather than supplying a hand-written bad list. The producer still reports `r={mutant['r']}`, but emits {mutant['emitted_words']} words; the unmodified CLI suite returns `{mutant['suite_exit']}` and identifies `{missing}`. The positive control returns `{mutation['cases']['positive']['suite_exit']}`. The self-test exits successfully only when those outcomes occur: {provenance(self_run)}. The test suite passes {tests} tests, including {exhaustive['exhaustive_guides']} independently enumerated guides through length {exhaustive['maximum_exhaustive_length']}, terminal-run collapse, invalid replacements preserving counts, duplicate padding, per-guide ownership, malformed input and corrupted raw alignments: {provenance(tests_run)}.

## Claim boundary

This is conformance of a single-deletion word family at a stated Hamming budget. The source theorem's sound-family result applies through `m <= n-3`; Paper A's unrestricted concrete-centre floor is separately scoped to `n >= 20` and `m <= 5`, with applicability flags in the producer. Beyond those hypotheses the program still checks deletion-set equality and makes no minimum-cover claim. An invalid word means outside this specified family; it is not an assertion of unsoundness outside the theorem's range.

Passing does not establish full search-engine completeness at other budgets, PAM policies, strands, references or devices. This deposit makes no claim about cleavage, biological risk, ranking impact, tool safety, upstream acceptance or adoption. Fresh work here is local pattern generation and conformance checking; no genome-scale search is performed.

Exact commands, working directories, script hashes, stdout/stderr paths and exit codes are in [the run ledger](evidence/runs.jsonl). Source/snapshot hashes are in [the manifest](inputs/manifest.json). Earlier diagnostic runs remain in the ledger; use the final successful stage named above for the published measurements.
"""
    (HERE / "README.md").write_text(readme)

    command_rows = []
    for row in (score_run, patch_run, self_run, tests_run, sentence_run):
        script = Path(row["argv"][2]).name
        command_rows.append(f"| `{row['label']}` | `/usr/local/bin/python3 -I {script}` | `{row['exit_code']}` | [{row['stdout']}]({row['stdout']}) |")
    findings = f"""# Findings

The completed deposit is [conformance/](conformance/README.md). The standalone producer, suite, minimal wrapper diff, planted-regression self-test, evidence snapshots and run logs are all within this directory. No paper or upstream source was edited outside this directory.

The following commands were executed in `conformance/`, under `nice` with thread caps. Every numerical result below comes from these runs or their recorded child commands; [the ledger](conformance/evidence/runs.jsonl) contains the exact absolute argv, working directory, stdout/stderr and exit code. Names and version identifiers follow the Table 1 source snapshot.

| Run | Command in deposit | Exit | Output |
|---|---|---:|---|
""" + "\n".join(row.replace("](evidence/", "](conformance/evidence/") for row in command_rows) + f"""

## Scoreboard

The panel requires {summary['required_guide_specific_words']} words over {summary['guides']} guides. The suite rescored actual wrapper files, parser map keys and raw fixture alignments, and cross-checked {summary['summary_rows_crosschecked']} deposited summary rows against the derived values. Its [inventory](conformance/inputs/output_inventory.json) records why each genomic/clinical output or count-only artifact is inapplicable as an emitted pattern set.

{scoreboard}

Executed/static is Table 1's evidence state, not a promotion made by this run. Fresh native engines, R calls and untested pipelines were not run here. In particular, the early native betas remain static loop emulations; CHOPOFF has source-component count evidence but no deposited native emitted word set; CRISPRme delegates to separately tested engines; SWOffinder and CALITAS have no deletion-library count verdict. CRISPRseek's pooled written file is checked for equality with the wrapper union, and per-guide ownership comes from the executed wrapper snapshots, not from crediting a shared word to every possible owner.

CRISPRitz has exact projected word sets and redundant alignment rows. Native parser dumps carry placement metadata, which is collapsed to their unique map keys; it is not counted as repeated emitted patterns. These distinctions are retained in the [machine-readable scoreboard](conformance/results/scoreboard.json), with every missing word in per-guide reports.

## Patch and regression test

The [unified diff](conformance/wrapper.patch) modifies {patch['modified_source_lines']} source lines, one RNA loop bound per PAM orientation: {patch['removed_lines']} removed and {patch['added_lines']} added, or {patch['unified_diff_physical_lines']} physical diff lines including context and headers. Both independently guarded branches omit endpoints, so each requires a bound repair; existing dictionary deduplication makes further changes unnecessary for word-set exactness. This establishes the smallest local bound repair justified here, not global byte-minimality.

For each PAM orientation, fresh official dry-runs change guide-specific emission from {before['valid_distinct_words']} to {after['valid_distinct_words']}, missing words from {before['missing_words']} to {after['missing_words']}, and incomplete guides from {before['guides_incomplete']} to {after['guides_incomplete']}. On the deposited planted fixture, pooled emission changes from {before['fixture_pooled_emitted']} to {after['fixture_pooled_emitted']} of {after['fixture_pooled_required']}, with pooled missing words changing from {before['fixture_pooled_missing']} to {after['fixture_pooled_missing']}. The patched sets have {after['invalid_words']} invalid words and {after['redundant_rows']} redundant rows. These are pattern-generation measurements, not a new fixture search by the native engine; [patch validation](conformance/results/patch_validation.json) records both orientations and every child exit.

The [self-test](conformance/results/selftest.json) actually mutates `required_words` to discard the terminal word. Its mutant producer exits `{mutant['producer_exit']}` with {mutant['emitted_words']} words while retaining count {mutant['r']}; the suite exits `{mutant['suite_exit']}` and names `{missing}` as missing. The positive suite exits `{mutation['cases']['positive']['suite_exit']}`, and the self-test harness exits `{self_run['exit_code']}` only after requiring both outcomes. The separate test run exits `{tests_run['exit_code']}` with {tests} tests passing, including an independent all-position oracle on {exhaustive['exhaustive_guides']} guides through length {exhaustive['maximum_exhaustive_length']}.

## Body sentence and its space payment

Proposed sentence, naming the actual deposit:

```tex
{sentence['new_sentence_tex']}
```

Replace this existing sentence in `paper_v7/sections/03_method.tex`:

> {sentence['old_sentence_tex'].replace(chr(10), ' ')}

The new deposit supplies the scoreboard and exact-emitter control that make the old audit-only summary redundant. In the paper's actual typography, the new sentence prints in {sentence['new_printed_lines']} lines and the removed sentence in {sentence['old_printed_lines']}; local full-paper copies compile with {sentence['paper_builds']['before']['body_pages']} body pages before and {sentence['paper_builds']['after']['body_pages']} after. Commands and exits are above; [the layout report](conformance/results/sentence_check.json), [context PDF](conformance/build/manuscript/after/main.pdf) and rendered context are deposited. No manuscript edit has been applied.

Earlier layout attempts are retained: the initial inline driver exited `1` on timeout ([record](conformance/evidence/layout_initial_timeout.json)); a subsequent longer sentence was rejected by a layout check exiting `1` ([rejected proposal](conformance/results/rejected_sentence.json)). A context-render lookup also exited `1` while still searching for the longer draft's wording; the final lookup uses the final sentence. The final sentence and full-paper build pass the checks above. The initial scorer also exposed repeated parser metadata as if it were emission; final reports correctly use pattern-map keys, with earlier stdout retained for provenance.

## Scope

The linear-time run count is a necessary screen. Exact sound-family conformance additionally requires word membership/equality; generating and comparing all materialized words is not claimed to be linear. The Hamming-budget and soundness hypotheses are explicit in the producer and README. This deposit concerns single-deletion pattern completeness only: no cleavage, risk, ranking-impact, tool-safety, whole-search certification or adoption claim is made.
"""
    (HERE.parent / "FINDINGS.md").write_text(findings)
    print(json.dumps({"readme": "conformance/README.md", "findings": "FINDINGS.md", "tests_passed": tests,
                      "new_sentence_lines": sentence["new_printed_lines"]}, sort_keys=True))


if __name__ == "__main__":
    main()
