# T1: field-level single-RNA-bulge completeness audit

**VERIFIED-COMPUTATIONALLY  -  principal result:** **CRISPRitz passes the tested deletion-position check.** Unmodified native CRISPRitz v2.9.0 and v2.7.1 each recovered all **1,634 guide-specific distinct deletions on all 110 real panel guides**, with zero deficits. The Cas-OFFinder 2016 wrapper is short on **98/110** guides; the native **3.0.0b3/develop parser** is short on **73/110**. The evidence supports a specific Cas-OFFinder enumeration defect, not a universal failure of current off-target tools.

**VERIFIED-COMPUTATIONALLY  -  audit scope:** sources were retrieved on 2026-10-29. All downloads, builds, fixtures, reports and scripts are under this output directory. No RECOMB_2027 file was modified, no commits or pushes were made, and no cluster jobs were submitted. The panel input is the actual 110-guide `multiguide_20261017/match_changeseq.json`, not the synthetic panel shipped elsewhere in the repository. Its SHA-256 is `7f24817f9e5f3b37fee8fb1b160a9219b22322ec1d6117df5365fb2aebee5693`.

**PROVED  -  release-status boundary:** the audited native Cas-OFFinder README [line 1](https://github.com/snugel/cas-offinder/blob/0a9ac005859641ca215add6c6dbed147ae80f7c0/README.md#L1) marks version 3 as not production ready. The native source result therefore concerns its named beta/develop snapshots; this audit does not measure deployment prevalence or establish that every current production tool fails.

## Verdicts and deciding source lines

**PROVED  -  interpretation:** the established criterion concerns distinct complete deletion words, after removing fixed PAM notation and external length-padding. Raw deletion positions, repeated alignments, short prefix probes and selected genomic loci are different objects. A demonstrated count below r rejects a deletion library; no count at or above r by itself certifies complete search. The following verdicts keep that distinction explicit.

| Tool and source version | Actual architecture | Distinct single-deletion family and duplicate handling | Verdict and execution scope |
|---|---|---|---|
| CRISPRitz v2.9.0; additionally v2.7.1 | TST traversal over indexed targets, branching on gaps; no standalone deletion-pattern file | n possible RNA-gap placements imply r distinct deletion words. Identical full output lines are deduplicated; equivalent gap placements within a run remain separate alignments. | **VERIFIED-COMPUTATIONALLY: PASS positional fixture**, both native versions, all 110 guides. No general genome-search completeness certificate. |
| CRISPRme 2.1.14; additionally active CRISPRme+ 2.8.0 | Workflow delegating intact guides and bulge budgets to CRISPRitz | No independent deletion-pattern family or deduplication at the wrapper level. Docker engine pins are 2.7.1 and 2.9.0 respectively. | **PROVED:** delegated architecture. **VERIFIED-COMPUTATIONALLY:** forwarding checks passed and both pinned engines tested separately. Whole pipelines **UNTESTED**. |
| SWOffinder v1.0.1 | Smith to Waterman-style dynamic programming and constrained traceback | No deletion-pattern library. Alignment/locus selection is not pattern deduplication. | **PROVED: count test N/A**. Native runtime **UNTESTED**. |
| CHOPOFF tag v1.0.3 (Project.toml says 1.0.2) | Hybrid symbolic alignment probes plus final Levenshtein alignment | RNA-gap graph includes every position; its unpadded single-deletion projection contains all r words. Concrete prefix probes are deduplicated and can merge more than full deletion words. | **PROVED:** all-position source component. **VERIFIED-COMPUTATIONALLY:** source-derived transcription has zero deficits. Native application **UNTESTED**; strict count test **N/A to ordinary truncated prefix probes**. |
| Cas-OFFinder native banner 3.0.0; release tag 3.0.0b3 and develop | Explicit deletion-pattern enumeration into a map | n-1 raw deletion positions; r-e distinct words, where e=1[g[n-1]!=g[n-2]], for n>=2. Keys collapse duplicate words; position metadata is retained. | **PROVED / VERIFIED-COMPUTATIONALLY: FAIL** on 73/110 guides. Complete unmodified parser bodies executed separately. Full search runtime unavailable because OpenCL device discovery failed. |
| Cas-OFFinder-bulge v1.2 beta, 2016 wrapper | Explicit deletion-pattern enumeration | n-2 raw deletion positions; r-s-e distinct words, where s=1[g[0]!=g[1]], for n>=3. Dictionary keys collapse duplicate patterns. | **VERIFIED-COMPUTATIONALLY: FAIL** on 98/110 guides using the official unmodified dry-run. |

**PROVED  -  CRISPRitz deciding code.** At [searchOnTST.cpp:946 to 951](https://github.com/pinellolab/CRISPRitz/blob/0715d27bf6a3c2375ce2c4a236a5100c1cacfd88/sourceCode/CRISPR-Cas-Tree/searchOnTST.cpp#L946-L951), `if (bR > 0)` branches to the next guide character at the same tree node. No interior-only position bound excludes an endpoint. The guide reversal for trailing PAMs is at [1076 to 1079](https://github.com/pinellolab/CRISPRitz/blob/0715d27bf6a3c2375ce2c4a236a5100c1cacfd88/sourceCode/CRISPR-Cas-Tree/searchOnTST.cpp#L1076-L1079). Whole-line deduplication is defined at [122 to 147](https://github.com/pinellolab/CRISPRitz/blob/0715d27bf6a3c2375ce2c4a236a5100c1cacfd88/sourceCode/CRISPR-Cas-Tree/searchOnTST.cpp#L122-L147) and invoked at line 1432. v2.7.1 has the corresponding RNA branch at [939 to 944](https://github.com/pinellolab/CRISPRitz/blob/4e70530031626b7a8377d8e0e32179f390d809db/sourceCode/CRISPR-Cas-Tree/searchOnTST.cpp#L939-L944). Detailed native execution: [CRISPRitz report](crispritz_audit/REPORT.md).

**PROVED  -  CRISPRme deciding code.** The reference route forwards the unchanged guide file and `-bRNA $bRNA` at [submit_job_automated_new_multiple_vcfs.sh:398 to 399](https://github.com/pinellolab/CRISPRme/blob/170aeaf1e2e0e8c2379aa14e49d73c4246150ed1/PostProcess/submit_job_automated_new_multiple_vcfs.sh#L398-L399); [Dockerfile:8](https://github.com/pinellolab/CRISPRme/blob/170aeaf1e2e0e8c2379aa14e49d73c4246150ed1/Dockerfile#L8) pins CRISPRitz 2.7.1. Its active successor forwards the same inputs at [shell:635 to 637](https://github.com/pinellolab/crisprme-plus/blob/6de886182cd68f02c0a528200a1f5f9e4546d1c2/PostProcess/submit_job_automated_new_multiple_vcfs.sh#L635-L637), with CRISPRitz v2.9.0 pinned at [Dockerfile:29](https://github.com/pinellolab/crisprme-plus/blob/6de886182cd68f02c0a528200a1f5f9e4546d1c2/Dockerfile#L29). Executing source function bodies with command capture verified forwarding on all six guides for both wrappers. Those checks are not full pipeline executions. [CRISPRme report](crisprme_audit/REPORT.md).

**PROVED  -  SWOffinder deciding code.** [SmithWatermanOffTargetSearchAlign.java:269 to 283](https://github.com/OrensteinLab/SWOffinder/blob/c20e581a2bcfad841afd62fe667285aadf4394f2/SmithWatermanOffTarget/SmithWatermanOffTargetSearchAlign.java#L269-L283) implements the DP recurrence; [221 to 230](https://github.com/OrensteinLab/SWOffinder/blob/c20e581a2bcfad841afd62fe667285aadf4394f2/SmithWatermanOffTarget/SmithWatermanOffTargetSearchAlign.java#L221-L230) traces a guide base against a target gap. It outputs candidate reference alignments at [357 to 376](https://github.com/OrensteinLab/SWOffinder/blob/c20e581a2bcfad841afd62fe667285aadf4394f2/SmithWatermanOffTarget/SmithWatermanOffTargetSearchAlign.java#L357-L376). It does not emit the deletion-pattern set requested by the count test.

**PROVED  -  CHOPOFF deciding code.** [motif_path_templates.jl:124 to 137](https://github.com/JokingHero/CHOPOFF.jl/blob/67d70c343c58628888fcd0e6ba6446483c1bc934/src/motif_path_templates.jl#L124-L137) iterates `for i in 1:len` and inserts RNA-gap edges, including both endpoints. Gap-path conversion is at [334 to 345](https://github.com/JokingHero/CHOPOFF.jl/blob/67d70c343c58628888fcd0e6ba6446483c1bc934/src/motif_path_templates.jl#L334-L345); symbolic duplicates are removed at [364 to 368](https://github.com/JokingHero/CHOPOFF.jl/blob/67d70c343c58628888fcd0e6ba6446483c1bc934/src/motif_path_templates.jl#L364-L368). Concrete mapped prefix probes are collapsed with `sa = unique(sa)` at [db_prefix_hash.jl:392](https://github.com/JokingHero/CHOPOFF.jl/blob/67d70c343c58628888fcd0e6ba6446483c1bc934/src/db_prefix_hash.jl#L392), and final alignment occurs at [419 to 423](https://github.com/JokingHero/CHOPOFF.jl/blob/67d70c343c58628888fcd0e6ba6446483c1bc934/src/db_prefix_hash.jl#L419-L423). This architecture makes a raw prefix-probe count an invalid substitute for a count of full deletion balls.

**PROVED  -  Cas-OFFinder deciding code.** The native develop loop at [cas-offinder.cpp:673](https://github.com/snugel/cas-offinder/blob/0a9ac005859641ca215add6c6dbed147ae80f7c0/cas-offinder.cpp#L673) uses the strict bound `j < sline[0].find_last_not_of('N') + 1 - i`. With i=1 and an ACGT guide followed by PAM placeholders, this allows zero-based j=0,...,n-2 and omits j=n-1. The release has the same bound at [3.0.0b3:669](https://github.com/snugel/cas-offinder/blob/1b7d803112cbb8ab5228d2b8ccc0589c28764b33/cas-offinder.cpp#L669). The map-key deduplication is explicit in [cas-offinder.cpp:14 to 18](https://github.com/snugel/cas-offinder/blob/0a9ac005859641ca215add6c6dbed147ae80f7c0/cas-offinder.cpp#L14-L18). The 2016 wrapper uses `range(1, len(target)-bulge_size)` at [cas-offinder-bulge:72](https://github.com/hyugel/cas-offinder-bulge/blob/b5d6ab4087814e3ccb7e1b2cc7a280286119ceb1/cas-offinder-bulge#L71-L75), excluding both endpoints. Its dictionary and unique-key output are at [48](https://github.com/hyugel/cas-offinder-bulge/blob/b5d6ab4087814e3ccb7e1b2cc7a280286119ceb1/cas-offinder-bulge#L48) and [87 to 88](https://github.com/hyugel/cas-offinder-bulge/blob/b5d6ab4087814e3ccb7e1b2cc7a280286119ceb1/cas-offinder-bulge#L87-L88). [Baseline report](cas_baselines/REPORT.md).

**PROVED  -  CALITAS deciding code.** [SequentialGuideAligner.scala:210](https://github.com/editasmedicine/calitas/blob/7a3456c47c4f920d7ef9eaa49ed8e61ee6c46bef/calitas/src/main/scala/com/editasmedicine/aligner/SequentialGuideAligner.scala#L210) selects `Mode.Glocal`; [295 to 300](https://github.com/editasmedicine/calitas/blob/7a3456c47c4f920d7ef9eaa49ed8e61ee6c46bef/calitas/src/main/scala/com/editasmedicine/aligner/SequentialGuideAligner.scala#L295-L300) passes complete guide and target sequences to the aligner. [315 to 320](https://github.com/editasmedicine/calitas/blob/7a3456c47c4f920d7ef9eaa49ed8e61ee6c46bef/calitas/src/main/scala/com/editasmedicine/aligner/SequentialGuideAligner.scala#L315-L320) selects noncompeting alignments. No deletion-pattern count is defined. [Alignment-tool report](alignment_audit/REPORT.md).

## Shared six-guide count table

**VERIFIED-COMPUTATIONALLY:** all guides below have n=20. Cells are **distinct count / deficit**. Wrapper cells come from official unmodified dry-run output; native Cas cells come from executing the exact complete parser source; CRISPRitz cells are distinct deletion words in emitted native alignments; CHOPOFF cells are the source-derived unpadded RNA-gap projection and are not native emitted probes. Both Cas snapshots agree, as do both CRISPRitz versions.

| Guide | Sequence | n | r | 2016 wrapper | Native Cas 3.0.0b3/develop | CRISPRitz 2.9.0 and 2.7.1 | CHOPOFF source component |
|---|---|---:|---:|---:|---:|---:|---:|
| g82 | GGAGAAGGTGGGGGGGTTCC | 20 | 9 | 9 / 0 | 9 / 0 | 9 / 0 | 9 / 0 |
| g97 | GAGAAGGTGGGGGGGTTCCA | 20 | 10 | 8 / 2 | 9 / 1 | 10 / 0 | 10 / 0 |
| g101 | GGGGCCACTAGGGACAGGAT | 20 | 13 | 12 / 1 | 12 / 1 | 13 / 0 | 13 / 0 |
| g0 | GATTTCTATGACCTGTATGG | 20 | 16 | 15 / 1 | 16 / 0 | 16 / 0 | 16 / 0 |
| g2 | GTACAGGCTGCACCTGTCAG | 20 | 18 | 16 / 2 | 17 / 1 | 18 / 0 | 18 / 0 |
| g46 | GAGTAGCGCGAGCACAGCTA | 20 | 20 | 18 / 2 | 19 / 1 | 20 / 0 | 20 / 0 |

**PROVED  -  N/A rows:** CRISPRme, CRISPRme+, SWOffinder and CALITAS have no independent deletion-library count for any of these guides; their emitted-pattern count and deficit are recorded as N/A in the normalized [60-row selected table](audit_selected_counts.tsv). Its separate columns distinguish actual pattern enumeration, observed native alignment words, and source-derived words. No untested program is labeled PASS.

## All-guide totals, original producers and bounded builds

**VERIFIED-COMPUTATIONALLY:** the complete panel requires 1,634 per-guide distinct words, or 1,629 pooled words after cross-guide collisions. The wrapper emitted 1,491 per-guide words: 12 guides had no deficit, 53 missed one and 45 missed two. Each native Cas parser emitted 1,561: 37 guides had no deficit and 73 missed one. Each tested CRISPRitz engine recovered all 1,634; its fixture contained all 1,629 pooled words with AGG PAMs. The [550-row native/source-parser table](audit_all110_native_counts.tsv) preserves every numerical row.

**VERIFIED-COMPUTATIONALLY:** byte-identical copies of the repository's `code_release/make_scheme_lemmas.py` and `make_pattern_sets.py --control` ran successfully in `cas_baselines/`, reproducing the published 110-guide run identity, endpoint counts and omission law. The official wrapper dry-run independently reproduced 98/110. No dedicated existing native-v3 producer was found in the searched repository paths; its omission was reconfirmed using newly executed verbatim complete `parseInput()` source bodies from the pinned release and develop sources. Those bodies were byte-identical across the two snapshots. The new native 73/110 panel figure is explicitly a parser-source result.

**VERIFIED-COMPUTATIONALLY  -  execution limits:** native CRISPRitz builds completed cleanly with the already installed GNU compiler and OpenMP. Native Cas-OFFinder also compiled, but OpenCL device discovery failed with `clGetDeviceIDs Failed: -30`; the enumeration parser was then executed independently of OpenCL. Native CHOPOFF was not built: Julia was absent and no release binary asset was available. Full CRISPRme workflows and native SWOffinder/CALITAS searches were not run. Their source mechanisms answer whether this count test applies; they remain UNTESTED for full search behavior.

## Every executed control and its outcome

**VERIFIED-COMPUTATIONALLY:** **54 negative-control instances were rejected and both positive controls were accepted.** Full per-instance records are in [audit_results.json](audit_results.json) and the original sub-audit logs. The control groups are:

**PROVED  -  control interpretation:** the CHOPOFF raw-count control rejects an incorrect measurement of distinctness, not a complete deletion family merely because its raw count exceeds r. All completeness controls use distinct sets; no converse of the one-sided count criterion is assumed.

| Component | Control | Actual outcome |
|---|---|---|
| CRISPRitz, 6 guides  *  2 engines | Remove one distinct observed deletion | 12/12 REJECTED |
| CRISPRitz, 6 guides  *  2 engines | Pad an incomplete word set with duplicate rows until raw count exceeds r | 12/12 REJECTED |
| CRISPRitz, 6 guides  *  2 engines | Execute native search with RNA budget zero | 12/12 REJECTED as a one-RNA search; both engines produced zero RNA rows on all 110 guides |
| Cas baseline checker | Omit final deletion on r=20 guide | REJECTED |
| Cas baseline checker | Substitute a duplicate while preserving raw count 20 | REJECTED |
| Cas baseline checker | Actual native-parser output on r=20 guide | REJECTED |
| Cas baseline checker | Actual wrapper output on r=20 guide | REJECTED |
| Cas baseline checker | Invalid nucleotide | REJECTED |
| Native parser | Malformed guide length | REJECTED, exit 1 |
| Copied repository producer | Invalid guide length | REJECTED, exit 1 |
| Copied repository producer | 109-guide input against 110-guide gate | REJECTED, exit 1 |
| Cas baseline checker | Complete reference deletion set | ACCEPTED, as required |
| CRISPRme source command checks, both versions | Force forwarded RNA budget to zero | 2/2 REJECTED |
| CRISPRme source command checks, both versions | Replace guide-file path | 2/2 REJECTED |
| CRISPRme source command checks, both versions | Request zero-RNA search while checking for one-RNA forwarding | 2/2 REJECTED |
| CHOPOFF source-derived graph | Remove first RNA-gap edge on r=20 guide | REJECTED: 19<20 |
| CHOPOFF source-derived graph | Remove last RNA-gap edge on r=20 guide | REJECTED: 19<20 |
| CHOPOFF source-derived graph | Report raw gap-template count as distinct count on r=9 guide | REJECTED: reported 20, actual 9 |
| CHOPOFF source-derived graph | Suppress all RNA-gap paths on r=20 guide | REJECTED: 0<20 |
| CHOPOFF source-derived graph | Remove final raw position inside terminal CC run on r=9 guide | ACCEPTED: all 9 distinct words remain |

**VERIFIED-COMPUTATIONALLY:** CHOPOFF's source transcription was additionally checked against the direct deletion-set oracle for every panel guide and all 21,844 ACGT guides of lengths 1 to 7. This validates the transcription's deletion component, not execution of native CHOPOFF.

## Abstract-ready wording and claim boundary

**PROVED / VERIFIED-COMPUTATIONALLY  -  three sentences, also in [ABSTRACT_READY.txt](ABSTRACT_READY.txt):**

> A source-level deletion-library audit identified incomplete single-RNA-bulge enumeration in the 2016 Cas-OFFinder wrapper on 98 of 110 CHANGE-seq guides and in the Cas-OFFinder 3.0.0b3/develop parser on 73 of 110 guides. In contrast, unmodified CRISPRitz 2.9.0 recovered every distinct single deletion for all 110 guides in a synthetic-reference positional test, as did the 2.7.1 engine pinned by CRISPRme. Source analysis showed that CRISPRme delegates to CRISPRitz, CHOPOFF includes every deletion position in its symbolic graph, and SWOffinder and CALITAS use alignment methods outside the deletion-pattern count test.

**REFUTED  -  for the audited positional criterion:** "All current bulge-search tools share this deletion-position omission." CRISPRitz's tested output is a direct counterexample. **CONJECTURE / UNTESTED:** a broader claim that every admissible target is found at every mismatch budget by any entire audited search pipeline is not established here.

**VERIFIED-COMPUTATIONALLY  -  artifact map:** `audit_selected_counts.tsv` is the normalized six-guide table; `audit_all110_native_counts.tsv` contains all numerical native/parser rows; `audit_results.json` combines results and every control; `artifact_sha256.json` pins the final report, tables, scripts and representative raw outputs. Each subdirectory retains its official source checkout, immutable commit evidence, executed script and raw logs. `assemble_audit.py` regenerates the combined numerical artifacts. Hub: no (direct source retrieval and isolated source/native execution supplied the required evidence; unpublished repository context was not sent to an external review service).
