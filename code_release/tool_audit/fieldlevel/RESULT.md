# A5: Cas-OFFinder 2 positional audit and wrapper distribution

Executed 2026-10-31, using the supplied offline source snapshots. **Cas-OFFinder 2 v1.0.0 retains the terminal RNA-bulge omission in its tested CPU search.** On each strand it recovered **1,561/1,634 guide-specific distinct deletion words**, missing one unique terminal word for **73/110 guides**. All first and interior deletion positions were recovered. Static inspection confirms three bundled wrapper copies with unchanged position loops, an embedded derivative, and additional call/install routes; the evidence does not establish affected user counts.

## 1. What was built and executed

The supplied `pnucolab/cas-offinder-rust` checkout is tag `v1.0.0`, commit **d7df0f45ec72a70f36510bf4df99d2c00c472105**, dated **2026-08-29T15:09:28+09:00**. The CLI identifies itself as “Cas-OFFinder 2 - v1.0.0”; its library manifest is version 1.0.1. Tracked upstream source files were unchanged (`git diff --exit-code` passed); the vendor directory and Cargo configuration are local additions. The design note was read in full, but the results below come from execution, rather than from its stated goals.

Rust **1.98.1** and Cargo **1.98.1** were used. The first strict offline release CLI build failed on upstream ambiguous `Result` glob imports in `cli_utils.rs` (preserved in [the failure log](logs/build_cli_strict_failed.log)). It then built successfully with **`RUSTFLAGS='-C codegen-units=1 -A ambiguous_glob_imports'`**, offline vendored sources, and one Cargo job. This changes a compiler lint setting, not the source or search logic. See [build_offline.sh](build_offline.sh), [build log](logs/build_cli.log), and [CLI help output](logs/cli_help.log).

The CPU scheduler in upstream `search.rs:1099–1109` uses `thread::available_parallelism()` and spawns that many workers; it does not consult the supplied native-library thread caps. To respect the compute rules, **the CLI was built and invoked for help, but its multithreaded search was not run**. The permitted small Rust harness instead includes the original source files without editing them and directly calls:

1. `precompute_pam_precheck`, with the same forward/reverse-complement pattern preparation as the native CPU path;
2. **`search_chunk_myers`**, including its PAM precheck, Myers candidate sweep, `classify_myers_candidate` and `traceback_all`;
3. **`convert_matches_myers`**, including native coordinate and reverse-complement conversion.

The adapter is [serial_bridge.rs](cpu_harness/src/serial_bridge.rs); its driver is [main.rs](cpu_harness/src/main.rs). It uses one thread, one entire small chromosome buffer at a time, and the original bit4 encoding. It adds guide labels and two explicit per-type bulge-count columns to its TSV output. It does not replace the native candidate generator or traceback with an oracle. **GPU kernels, the public threaded scheduler, streamed chromosome/chunk assembly, and the CLI result writer were not tested.** Zero-bulge controls also call this same Myers function directly, rather than the CLI's separate popcount branch.

### Fixture and budgets

The original [CRISPRitz audit report](../T1_cross_tool_audit/crispritz_audit/REPORT.md), `run_audit.py`, and `run/fixture_manifest.json` were read. The brief's top-level `crispritz_audit/guides.txt` is absent; `run_audit.py` names **`run/guides.txt`**, which is present and was used. The original `run/chrAudit.fa` is absent, so its exact sequence was reconstructed using the supplied method, without a filesystem search. The 110 guides were cross-checked against the named panel and the guide file; the panel SHA-256 is **7f24817f9e5f3b37fee8fb1b160a9219b22322ec1d6117df5365fb2aebee5693**.

The RNA fixture contains all **1,629 pooled distinct** 19-base deletion words, each followed by AGG, separated and bounded by 64 A bases: **140,158 bases**, matching the original manifest. Summing distinct deletion words separately for each guide gives **1,634**, since some words are shared between guides. For each 20-base guide, the independent expected set is `{guide[:i] + guide[i+1:] : i=0,…,19}`, with size `r = 1 + number of adjacent-base changes`. The search pattern is `guide + NNN`, with filter `N×20 + NGG`.

The additional DNA fixture plants every distinct word obtained by inserting each of A/C/G/T at each of the 21 guide boundaries: **64 distinct words per guide**, **7,040 guide-specific / 7,035 pooled words**, and **619,144 bases**. Reverse-complement copies of both complete references were searched separately. The primary searches used:

| Fixture | Mismatches | DNA bulges | RNA bulges |
|---|---:|---:|---:|
| RNA, forward and reverse reference | 0 | 0 | 1 |
| DNA, forward and reverse reference | 0 | 1 | 0 |
| RNA-disabled control | 0 | 0 | 0 |
| Exact/nonmatch control reference | 0 | 1 | 1 |
| Exact/nonmatch/mismatch control reference | 0, then 1 | 0 | 0 |

These budgets cover the planted edit types. Every reported word was checked at its **own planted coordinate and expected strand**, not merely somewhere in the output. Projection removes alignment gaps and the three-base PAM. [fixtures/manifest.json](fixtures/manifest.json) records the guide IDs, counts and input hashes; [audit.py](audit.py) prepares, executes and evaluates the fixtures.

### RNA result

| Measure, independently on each strand | Result |
|---|---:|
| Expected guide-specific distinct words | 1,634 |
| Recovered at planted loci | 1,561 |
| Missing distinct words | 73 |
| Guides with a distinct-word deficit | 73 |
| Guides without a distinct-word deficit | 37 |
| Missing first-position words / interior-position words | 0 / 0 |
| Missing unique terminal-position words | 73 |
| Unexpected RNA words | 0 |
| Raw RNA alignment rows | 2,090 |

**All 110 guides emit exactly the 19 gap placements at guide positions 1–19, and no position-20 gap placement**, in guide 5′→3′ coordinates, on both strands. When the final base belongs to a run of at least two identical bases, the terminal deletion word is also obtainable by an allowed interior deletion. Therefore 37 guides recover all `r` words despite the missing terminal alignment; the other 73 recover `r−1`. Raw alignment counts exceed `r` and would conceal the defect if used alone.

| Guide | r | Recovered, forward | Recovered, reverse | Missing word positions |
|---|---:|---:|---:|---|
| g82 | 9 | 9 | 9 | none |
| g97 | 10 | 9 | 9 | 20 |
| g101 | 13 | 12 | 12 | 20 |
| g0 | 16 | 16 | 16 | none |
| g2 | 18 | 17 | 17 | 20 |
| g46 | 20 | 19 | 19 | 20 |

All 110 per-guide results appear below and in [rna_all110.tsv](rna_all110.tsv) / [rna_reverse_all110.tsv](rna_reverse_all110.tsv). The 73 missing words and their planted coordinates are in [rna_missing.tsv](rna_missing.tsv) and [rna_reverse_missing.tsv](rna_reverse_missing.tsv). Raw outputs are [run/rna.tsv](run/rna.tsv) and [run/rna_reverse.tsv](run/rna_reverse.tsv).

The deciding upstream code is explicit. In `cas-offinder-lib/src/traceback.rs:238–240`:

```rust
let first_non_n: usize = pam_is_n.iter().position(|&x| !x).unwrap_or(0);
let last_non_n: usize = pam_is_n.iter().rposition(|&x| !x).unwrap_or(0);
let forbidden_skip: usize = if is_rc { first_non_n } else { last_non_n };
```

And at `traceback.rs:349–352`:

```rust
if rb < max_rna && max_rna > 0 {
    let pattern_pos_is_n = pam_is_n[i - 1];
    let is_forbidden = (i - 1) == forbidden_skip;
    if !pattern_pos_is_n && !is_forbidden {
```

For `guide + NNN`, the non-N endpoints are the guide endpoints. The RC remapping forbids the same terminal guide base in the opposite orientation. Thus changing to Myers does not remove the omission: traceback deliberately retains it. This is source evidence at [the pinned release](https://github.com/pnucolab/cas-offinder-rust/blob/d7df0f45ec72a70f36510bf4df99d2c00c472105/cas-offinder-lib/src/traceback.rs#L228-L240), corroborated by the executed CPU test.

### DNA result and interpretation

On **both strands**, every guide recovered **61/64 distinct one-DNA-bulge words** at their planted spans: **6,710/7,040** guide-specific words, with **330 absent DNA-bulge words**. All interior insertion boundaries and the guide/PAM insertion boundary were recovered. The three absent words per guide are precisely `b + guide` for the three bases `b != guide[0]`: insertions before the first guide base. The fourth leading-insertion word, with the first base repeated, can be aligned using an allowed internal insertion.

This is a **DNA-bulge alignment/span deficit**, not 330 completely unreported physical loci: each of these words also contains the exact guide as a suffix, and the engine reports that suffix as `X`, shifted one base on the forward reference (with the corresponding reverse coordinate). The independent validator checked all 440 leading-insertion sites, including the repeated-first-base case, for that shifted exact hit on each strand. Native traceback emits/returns immediately at `i == 0` and prohibits a DNA-gap transition past the opposite pattern end (`traceback.rs:279–319, 365–379`), explaining the boundary behavior.

Full records: [dna_all110.tsv](dna_all110.tsv), [dna_reverse_all110.tsv](dna_reverse_all110.tsv), [dna_missing.tsv](dna_missing.tsv), and [dna_reverse_missing.tsv](dna_reverse_missing.tsv). There are 80 typed DNA alignment rows per guide, collapsing to 61 words; each complete DNA run has 9,240 rows including exact suffix hits. No unexpected typed DNA words were observed.

### Executed controls and independent verification

The toy smoke test recovered 19/20 distinct deletion words, including the first, and missed exactly the unique terminal word. The full controls passed:

- **Own exact sites:** all 110 guides' own exact sites were reported on each strand under each of the three control configurations: 220/220 per configuration.
- **Deliberate nonmatch:** the planted `ATATATATATATATATATAT + AGG` site, in its own chromosome, produced **zero rows for all 110 guides**, including at joint DNA/RNA budgets of 1/1. Its minimum guide-to-planted-spacer edit distance was independently calculated as 9.
- **Mismatch budgets:** one substitution was planted for every guide. All 110 were rejected at zero mismatch/zero bulges and accepted at one mismatch/zero bulges.
- **RNA budget zero:** zero RNA rows on the complete deletion fixture for all 110 guides.
- **Evaluation controls:** deleting one distinct word was rejected as incomplete for all 110 guides; duplicating the remaining words beyond `r` did not repair the deficit, also for all 110.

[verify_results.py](verify_results.py), separate from the search and summary code, independently validated **all 35,300 raw full-run rows** against the FASTA and guide: exact physical sequence, coordinates, strand, gap counts, mismatch count, PAM and per-type budgets. Forward/reverse alignment multisets were identical after coordinate/strand normalization. It also checked the per-guide RNA formula `r − I(last two bases differ)` and the three leading DNA insertion deficits. See [verification.json](verification.json), [logs/verification.log](logs/verification.log), [controls.json](controls.json), and [summary.json](summary.json).

**Claim boundary:** this establishes an executed positional counterexample for the pinned CPU algorithm at the listed budgets. It does not certify genome-scale completeness, other mismatch/bulge budgets, another release, GPU parity, or biological activity. The prior 129 measured sites/40 guides were not rerun or reattributed to this Rust release in this audit.

## 2. Who ships or calls the 2016 wrapper

This part is **static inspection only**. No downstream package, Docker image, server, or clinical pipeline was executed. `usage_evidence/SOURCES.txt` gives **file-specific commit/date metadata**, not a verified common checkout commit for every project. Retrieval was **2026-10-31T21:37:51Z**. Immutable source links and local deciding lines are supplied below; the complete line-numbered excerpts are in [usage_deciding_lines.md](usage_deciding_lines.md), and AST/byte comparisons and file hashes are in [usage_comparison.json](usage_comparison.json).

The reference is `hyugel/cas-offinder-bulge`, file commit **1762e98aa9a9a41727bd1e702e1e6e8f4b50bb6b**, **2019-06-08T00:58:50Z**. Its banner is `Cas-OFFinder-bulge v1.2 beta (2016-07-06)` (line 9). Its DNA loops use `for i in range(1, len(target))` (lines 55, 68), and its RNA loops use `for i in range(1, len(target)-bulge_size)` (lines 59, 72). For one RNA bulge, these emit zero-based guide indices 1–18, excluding both endpoints; for one DNA bulge they emit internal boundaries 1–19. The comparison used parsed Python ASTs, not only matching banner strings.

| Project | What the fetched source establishes | Position loops compared with upstream | Role/default evidence |
|---|---|---|---|
| CRISPRseek 1.53.1 | Bundles an adapted callable wrapper in `inst/extdata/Cas-OFFinder/cas-offinder-bulge`; helper loads it with `source_python()` for the 2.4.1 backend. | Both orientations' DNA/RNA inner loops, including bodies, are AST-identical. Whole file differs. | Version-specific bulge backend; 3.0.0b3 helper branch invokes the native binary directly. Exported-caller/default trace has the source gap explained below. |
| OffRisk | Local runner invokes `python` plus `const.CAS_OFFINDER_BULGE_PATH`; API runner calls `/v1/cas-offinder-bulge/`. | Configured wrapper body was not supplied, so identity cannot be checked. Its separate CRISPRitz-emulation code also has the same excluded-endpoint intervals, with renamed variables. | Optional search tool; guide documents FlashFry as the API and UI default, with default RNA/DNA budgets zero. |
| IGU_CHANGEseq | `changeseq/dev/run_casoffinder.py` embeds an attributed derivative of the wrapper, not an external `cas-offinder-bulge` invocation. | Forward DNA/RNA **position-range ASTs** are identical; full loop bodies differ because guide-ID handling was added. No reverse-PAM branch is implemented there. | Bulge mode selected when either budget is nonzero. This development script's shown top-level call sets RNA=1, DNA=1; default use by the whole CHANGE-seq pipeline is not established. |
| aryeelab/cas-designer | Docker recipe downloads `cas-offinder-bulge` from the unversioned rgenome URL and moves it into `/bin`. | No downloaded wrapper body or built image supplied; cannot compare. | Docker installation route. Recipe does not establish runtime/default invocation or a built image's bytes. |
| pinellolab/CRISPR_WGS_Detection | Bundles `supervised/wgs_tools/cas-offinder-bulge/cas-offinder-bulge`; README documents its invocation. | Both orientations' full inner loops are AST-identical; whole-file changes are the two `cas-offinder` → `./cas-offinder` command paths. | Bundled supervised-tool utility; no top-level pipeline default established by these two files. |
| BE-Designer | Bundles the wrapper and calls `./cas-offinder-bulge` from `run_cas_offinder()` in its normal `main()` flow. | **Whole file byte-identical** to upstream; thus loops are identical. | Default shown call uses GPU, but **both bulge budgets are fixed to zero**, so the supplied default flow does not exercise the RNA positional omission. |
| dhslab docker image recipe | Recipe for `ghcr.io/dhslab/docker-casoffinder-bulge` clones upstream wrapper and installs it in `/usr/local/bin`, alongside native 3.0.0b3. | Unpinned clone; no built image or installed wrapper supplied, so installed bytes cannot be compared. | Dedicated image recipe; comments identify Nextflow use, but those workflow files and a running/published image were not inspected. Default CMD is `/bin/bash`. |

### CRISPRseek: exported API, wrapper branch, and the exact missing edge

The fetched `NAMESPACE` **line 16** is `export(getOfftargetWithBulge)`; line 20 exports `offTargetAnalysis`. The utilities file **line 10** says `# Helper functions for getOfftargetWithBulge()` (with backticks in the original), and lines 13–17 define a helper that takes results from `offTargetAnalysis`. The executable wrapper-selection evidence is:

```text
getOfftargetWithBulgeUtils.R:152  prepCasOFFinderBulge <- function(version = c("2.4.1", "3.0.0b3")) {
getOfftargetWithBulgeUtils.R:153    if (version == "2.4.1") {
getOfftargetWithBulgeUtils.R:154      wrapper_path <- system.file("extdata", "Cas-OFFinder/cas-offinder-bulge", package = "CRISPRseek")
getOfftargetWithBulgeUtils.R:172      source_python(wrapper_path)
getOfftargetWithBulgeUtils.R:173      return(cas_offinder_bulge)
```

The bundled Python file defines `cas_offinder_bulge(argv1, argv2, argv3, argv4)` at **line 13**, has the unchanged RNA loops at **64 and 77**, and invokes `Popen([argv4, ...])` at **line 110**. The alternate utilities branch at **174–183** constructs and returns a direct native-binary runner for version 3.0.0b3. `getCasOFFinder` at **87–88** defaults its own `match.arg` selection to 2.4.1.

**The definition of exported `getOfftargetWithBulge()` and its actual call to `prepCasOFFinderBulge()` are not present in the supplied files.** Neither is the `offTargetAnalysis()` definition. Thus the provable trace is exported API name → documented helper association → version-specific Python loader → bundled callable wrapper → native-binary call; the **public-function-to-helper call edge and public default version remain unverified**. In particular, the utilities' vector default alone is insufficient to establish the exported function's default behavior. Completing that specific trace requires the pinned `R/getOfftargetWithBulge.R` caller source; it was not searched for or fetched, following the no-search/no-network instructions. This gap does not affect the verified bundled-copy identity.

### Deciding call/install lines for the other projects

- **OffRisk**, `app/off_target.py:330–334`: `args = ["python", const.CAS_OFFINDER_BULGE_PATH, path_in, c_g_option, ...]`, followed by `run_external_proc(args)`. Its API route is at **290**. The user guide's **109** is `search_tools: List[str] = ["flashfry"]`, and **233** says `Select search tool. The default is FlashFry.` Its separate CRISPRitz emulation at **588, 594** uses `range(1, len(sequence))` and `range(1, len(sequence) - bulge_size)`; this is a property of that adapter, not of native CRISPRitz.
- **IGU_CHANGEseq**, `run_casoffinder.py:42` attributes `https://github.com/hyugel/cas-offinder-bulge`; **62, 68** contain the original forward position bounds. **253–256** select bulge mode unless both budgets are zero; **266** calls the embedded function, **281** sets `mm, RNAbb, DNAbb, PU = 4, 1, 1, "C"`, and **287–295** execute that configured example. **102** invokes the native `cas-offinder` binary.
- **aryeelab/cas-designer**, `Dockerfile:27`: `wget http://www.rgenome.net/static/files/cas-offinder-bulge`; **30** makes it executable, and **32** moves `cas-*` into `/bin`.
- **CRISPR_WGS_Detection**, wrapper **59, 72** retain the RNA loops; **96** calls `Popen(['./cas-offinder', ...])`. Its README **7** documents `python3 cas-offinder-bulge input.txt G output.txt`.
- **BE-Designer**, `be_designer.py:31–32`: `self.bulge_dna = 0`, `self.bulge_rna = 0`; **188–190** put these budgets into the input; **208** calls `Popen(('./cas-offinder-bulge', ... 'G', ...))`; **260–261** call target finding and off-target search in `main()`.
- **dhslab**, `Dockerfile:28–29`: `git clone --depth 1 https://github.com/hyugel/cas-offinder-bulge.git /opt/cas-offinder-bulge` followed by `install -m755 ... /usr/local/bin/cas-offinder-bulge`. Native version 3.0.0b3 is selected at **20**; **36** sets `CMD ["/bin/bash"]`.

### File commits and dates

| Fetched project/file | Commit | Date (UTC) |
|---|---|---|
| [CRISPRseek NAMESPACE and DESCRIPTION](https://github.com/bioc/CRISPRseek/blob/604cea08aede432b36cfc09b141c7b602a88dd24/NAMESPACE) | `604cea08aede432b36cfc09b141c7b602a88dd24` | 2026-10-13 15:39:19 |
| [CRISPRseek R utilities](https://github.com/bioc/CRISPRseek/blob/ef901b7bb8623be2663e02883a02f597056d62f0/R/getOfftargetWithBulgeUtils.R) | `ef901b7bb8623be2663e02883a02f597056d62f0` | 2025-02-25 18:55:48 |
| [CRISPRseek wrapper](https://github.com/bioc/CRISPRseek/blob/b71c045e4c08295077eb832bcedbe3f483f0de61/inst/extdata/Cas-OFFinder/cas-offinder-bulge) | `b71c045e4c08295077eb832bcedbe3f483f0de61` | 2024-11-21 00:39:07 |
| [OffRisk Python caller](https://github.com/IsanaVekslerLublinsky/OffRisk/blob/223f91e3ee6554475eba1c1c2c477a2ebda8ef3f/app/off_target.py) | `223f91e3ee6554475eba1c1c2c477a2ebda8ef3f` | 2023-09-13 17:53:05 |
| [OffRisk guide](https://github.com/IsanaVekslerLublinsky/OffRisk/blob/b9744361b0aa1a6fc01c042d874f4eb5a7acdf0e/OffRisk_user_guide.md) | `b9744361b0aa1a6fc01c042d874f4eb5a7acdf0e` | 2023-08-23 19:33:10 |
| [IGU_CHANGEseq development script](https://github.com/Interventional-Genomics-Unit/IGU_CHANGEseq/blob/3c9d06b813e365ef1184ff4e4559fdb3b353bca8/changeseq/dev/run_casoffinder.py) | `3c9d06b813e365ef1184ff4e4559fdb3b353bca8` | 2024-09-04 19:58:13 |
| [aryeelab/cas-designer Dockerfile](https://github.com/aryeelab/cas-designer/blob/b9820f3b7ba8e38abae30a0893f279398c1f8c2d/Dockerfile) | `b9820f3b7ba8e38abae30a0893f279398c1f8c2d` | 2021-06-26 19:40:40 |
| [CRISPR_WGS_Detection wrapper and README](https://github.com/pinellolab/CRISPR_WGS_Detection/blob/bbf6b6454531ebf398315551d4c2fe820caf5797/supervised/wgs_tools/cas-offinder-bulge/cas-offinder-bulge) | `bbf6b6454531ebf398315551d4c2fe820caf5797` | 2021-05-27 09:09:49 |
| [BE-Designer caller](https://github.com/hyugel/BE-Designer/blob/9669b7ba0d580e19c6d957827858f3e9945610f6/be_designer.py) | `9669b7ba0d580e19c6d957827858f3e9945610f6` | 2020-07-19 23:54:41 |
| [BE-Designer wrapper](https://github.com/hyugel/BE-Designer/blob/13e920db175db8addee1713279ed81f5a331b5e9/cas-offinder-bulge) | `13e920db175db8addee1713279ed81f5a331b5e9` | 2019-11-11 07:55:34 |
| [dhslab Dockerfile](https://github.com/dhslab/dhslab-docker-images/blob/56b06f34107046c1b651306328d7b6e993abfda9/docker-casoffinder-bulge/Dockerfile) | `56b06f34107046c1b651306328d7b6e993abfda9` | 2026-07-27 15:08:43 |
| [Upstream reference wrapper](https://github.com/hyugel/cas-offinder-bulge/blob/1762e98aa9a9a41727bd1e702e1e6e8f4b50bb6b/cas-offinder-bulge) | `1762e98aa9a9a41727bd1e702e1e6e8f4b50bb6b` | 2019-06-08 00:58:50 |

The 2026 CRISPRseek package metadata and dhslab recipe show recent distribution evidence. They do not establish downloads, deployed service versions, default wrapper use across the packages, or how many analyses encountered a missing active site.

## 3. Exact commands, compute limits and artifacts

All computation was serial, with no network/cluster access. The environment used `OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4`, `RAYON_NUM_THREADS=1`, `CARGO_BUILD_JOBS=1`, and single-unit Rust code generation. Each build had an external 1,200-second `timeout`; each search had a 240-second `timeout` inside the bounded Python driver. There were no process pools. The harness calls no worker-spawning scheduler. One-minute loads recorded before builds/runs were below 40. Search processes ran in approximately 0.01–1.3 seconds; progress was emitted every ten guides. The recorded maximum child RSS was approximately 5.3 MB (the child metric is cumulative over a driver's children, not total build/orchestrator memory).

**Priority/platform detail:** every computational command was invoked with `nice -n 15`. The sandbox denied its `setpriority` syscall, leaving the **already inherited niceness of 20**, as recorded in build and run logs. That is lower priority than 15; no attempt was made to raise it or escape the sandbox. An initial attempt to add a Darwin `ulimit -d` cap was rejected as unsupported and aborted before compilation; the final scripts do not claim an OS-enforced data-segment cap. Fixtures and outputs are small, and no file over 1 GB was copied. Native threaded CLI search was avoided specifically to prevent uncontrolled local parallelism.

From this directory, the final complete reproduction was executed as:

```sh
sh reproduce.sh > logs/reproduce.log 2>&1
```

[reproduce.sh](reproduce.sh) contains the exact ordered commands and environments. Its substantive commands are:

```sh
sh build_offline.sh > logs/build_cli.log 2>&1
sh build_harness.sh > logs/build_harness.log 2>&1
timeout 60s nice -n 15 python3 audit.py prepare > logs/prepare.log 2>&1
timeout 120s nice -n 15 python3 audit.py smoke > logs/smoke_driver.log 2>&1
timeout 1200s nice -n 15 python3 audit.py full > logs/full_driver.log 2>&1
timeout 60s nice -n 15 python3 verify_results.py > logs/verification.log 2>&1
timeout 60s nice -n 15 python3 usage_audit.py > logs/usage_audit.log 2>&1
timeout 15s nice -n 15 build/target/release/cas-offinder > logs/cli_help.log 2>&1
```

Both build scripts ultimately execute `timeout 1200s nice -n 15 cargo build --release --offline --jobs 1` in the relevant crate directory, using generated `.cargo/config.toml` entries pointing to the supplied vendor directory. Each exact harness invocation, budget, exit code, elapsed time and child RSS is in [logs/commands.jsonl](logs/commands.jsonl); first-pass and final-reproduction entries are retained. One example is:

```sh
timeout 240s nice -n 15 build/target/release/positional-audit \
  fixtures/rna.fa fixtures/guides.tsv 0 0 1 run/rna.tsv 10
```

The run artifacts are in [run/](run/), logs in [logs/](logs/), and input/output/source checksums in [SHA256SUMS](SHA256SUMS). [cluster_bundle.tar.gz](cluster_bundle.tar.gz) contains a self-contained source/vendor/input bundle and separate [smoke](cluster/smoke.sbatch) and [full](cluster/full.sbatch) Slurm scripts. **No cluster submission is needed to obtain the reported results: the full audit already executed locally.** The Slurm scripts are prepared and syntax-checked, not cluster-tested. If the lead wants a cluster reproduction, [cluster/README.md](cluster/README.md) gives the exact smoke-then-full submission commands and expected outputs; no remote execution is claimed.

The archive was also extracted into a fresh temporary directory, with no preexisting build or run outputs, and **`sh reproduce.sh --staged` completed a fresh offline build and the full audit**. All 22 result files compared were byte-identical to the primary run ([bundle_validation.json](bundle_validation.json)); the fresh CLI/harness builds took 14.38/8.90 seconds. Their logs are [bundle_build_cli.log](logs/bundle_build_cli.log) and [bundle_build_harness.log](logs/bundle_build_harness.log). This is a local bundle validation, not a Slurm execution.

Durable work products were written only in this output directory; the isolated bundle check also used the permitted temporary directory. Neither paper nor `code_release/` was edited.

## 4. Five-sentence paragraph for Paper A

Using a serial harness calling the unchanged CPU search and traceback functions of Cas-OFFinder 2 v1.0.0 (`d7df0f45`), we recovered 1,561 of 1,634 guide-specific distinct single-deletion words across 110 guides at zero mismatches and one RNA bulge, with identical results on the reverse-complement fixture ([execution](summary.json), [validation](verification.json)). All 73 missing words were unique terminal deletions, consistent with the release's explicit `forbidden_skip` rule; first and interior deletion positions were recovered, while the terminal gap placement was absent for every guide ([source](https://github.com/pnucolab/cas-offinder-rust/blob/d7df0f45ec72a70f36510bf4df99d2c00c472105/cas-offinder-lib/src/traceback.rs#L228-L240), [per-guide results](rna_all110.tsv)). Exact-site, deliberate nonmatch, and edit-budget controls passed ([controls](controls.json)). Static source snapshots show copies of the 2016 wrapper with unchanged position loops in CRISPRseek, BE-Designer and CRISPR_WGS_Detection, an embedded derivative in IGU_CHANGEseq, and wrapper call or installation routes in OffRisk and two Docker recipes ([source excerpts](usage_deciding_lines.md), [comparisons](usage_comparison.json)). These observations establish distribution and available execution routes, while OffRisk's FlashFry default, BE-Designer's zero-bulge configuration and the absent CRISPRseek exported-caller definition prevent inferring universal default wrapper use or affected user counts.

## Appendix: all 110 guides

Positions are guide 5′→3′, 1-based for RNA deletions. Forward and reverse recovery counts agree for every row. Each guide recovered RNA alignment placements 1–19; placement 20 was absent even when its word was represented by an interior placement. DNA recovery is **61/64** for every guide on each strand, with the leading-boundary qualification described above.

| Guide | Sequence | r | Recovered, each strand | Missing RNA word position | DNA words, each strand |
|---|---|---:|---:|---|---:|
| g0 | GATTTCTATGACCTGTATGG | 16 | 16 | none | 61/64 |
| g1 | GGGAACCCAGCGAGTGAAGA | 14 | 13 | 20 | 61/64 |
| g2 | GTACAGGCTGCACCTGTCAG | 18 | 17 | 20 | 61/64 |
| g3 | GTCAGGGTTCTGGATATCTG | 16 | 15 | 20 | 61/64 |
| g4 | GCGTGACTTCCACATGAGCG | 18 | 17 | 20 | 61/64 |
| g5 | GCTTGTCCGTCTGGTTGCTG | 16 | 15 | 20 | 61/64 |
| g6 | GCTGCAGAAACAGCAAGCCC | 15 | 15 | none | 61/64 |
| g7 | GTCCCTAGTGGCCCCACTGT | 14 | 13 | 20 | 61/64 |
| g8 | GAAGGTGGCGTTGTCCCCTT | 12 | 12 | none | 61/64 |
| g9 | GGGGGGTTCCAGGGCCTGTC | 10 | 9 | 20 | 61/64 |
| g10 | GGCATACTCATCTTTTTCAG | 15 | 14 | 20 | 61/64 |
| g11 | GGACAGTAAGAAGGAAAAAC | 12 | 11 | 20 | 61/64 |
| g12 | GGTGGTGTGGGCCCAGGAGG | 12 | 12 | none | 61/64 |
| g13 | GAAGTTGACTTACTGAAGAA | 15 | 15 | none | 61/64 |
| g14 | GTCGAGAAAAGCTTTGAAAC | 13 | 12 | 20 | 61/64 |
| g15 | GGGGGTTCCAGGGCCTGTCT | 11 | 10 | 20 | 61/64 |
| g16 | GGCCGAGATGTCTCGCTCCG | 17 | 16 | 20 | 61/64 |
| g17 | GGCAGAAACCCTGGTGGTCG | 13 | 12 | 20 | 61/64 |
| g18 | GAGAATCAAAATCGGTGAAT | 14 | 13 | 20 | 61/64 |
| g19 | GCCCTGGCCAGTCGTCTGGG | 14 | 14 | none | 61/64 |
| g20 | GAAGATGATGGAGTAGATGG | 17 | 17 | none | 61/64 |
| g21 | GGTGGATGATGGTGCCGTCG | 16 | 15 | 20 | 61/64 |
| g22 | GCAGATGGAATCATCTAGGA | 17 | 16 | 20 | 61/64 |
| g23 | GGCCAGTACCACAGCAGGCT | 16 | 15 | 20 | 61/64 |
| g24 | GAAGAAACTGAGAAGCATGA | 16 | 15 | 20 | 61/64 |
| g25 | GAGCAGGGCTGGGGAGAAGG | 13 | 13 | none | 61/64 |
| g26 | GAAGCATGACGGACAAGTAC | 17 | 16 | 20 | 61/64 |
| g27 | GGTTTCACCGAGACCTCAGT | 15 | 14 | 20 | 61/64 |
| g28 | GCTGTCCTGAAGTGGACATA | 17 | 16 | 20 | 61/64 |
| g29 | GGGCTGCAGGGGAGCTGGGC | 13 | 12 | 20 | 61/64 |
| g30 | GCTGGTACACGGCAGGGTCA | 16 | 15 | 20 | 61/64 |
| g31 | GCTGACCCCGCTGGGCAGGC | 14 | 13 | 20 | 61/64 |
| g32 | GTTTGCGACTCTGACAGAGC | 18 | 17 | 20 | 61/64 |
| g33 | GGGACTCTACATCTGCAAGG | 16 | 16 | none | 61/64 |
| g34 | GGGGATTTCTATGACCTGTA | 14 | 13 | 20 | 61/64 |
| g35 | GCACGTGGCCCAGCCTGCTG | 16 | 15 | 20 | 61/64 |
| g36 | GACACCTTCTTCCCCAGCCC | 12 | 12 | none | 61/64 |
| g37 | GGGGGATTGTGGGGCTGCAG | 12 | 11 | 20 | 61/64 |
| g38 | GCATACTCATCTTTTTCAGT | 16 | 15 | 20 | 61/64 |
| g39 | GCTGCCGCCCAGTGGGACTT | 14 | 14 | none | 61/64 |
| g40 | GGTCACCACGAGCAGGGCTG | 16 | 15 | 20 | 61/64 |
| g41 | GTGGAGCGCAGTGGTCTCCG | 17 | 16 | 20 | 61/64 |
| g42 | GGAAACTTGGCCACTCTATG | 14 | 13 | 20 | 61/64 |
| g43 | GGGCGGTGCTACAACTGGGC | 14 | 13 | 20 | 61/64 |
| g44 | GTGGTACTGGCCAGCAGCCG | 16 | 15 | 20 | 61/64 |
| g45 | GCTGTTTCTGCAGCCGCTTT | 15 | 15 | none | 61/64 |
| g46 | GAGTAGCGCGAGCACAGCTA | 20 | 19 | 20 | 61/64 |
| g47 | GCCGTGGCAAACTGGTACTT | 14 | 14 | none | 61/64 |
| g48 | GTATGGAAAATGAGAGCTGC | 16 | 15 | 20 | 61/64 |
| g49 | GGCACCAACTGGATGGATCA | 15 | 14 | 20 | 61/64 |
| g50 | GCTGGCGATGCCTCGGCTGC | 17 | 16 | 20 | 61/64 |
| g51 | GGGTATTATTGATGCTATTC | 15 | 14 | 20 | 61/64 |
| g52 | GGTGAGGGAGGAGAGATGCC | 15 | 15 | none | 61/64 |
| g53 | GTCCCCTGAGCCCATTTCCT | 12 | 11 | 20 | 61/64 |
| g54 | GCTACTCTCTCTTTCTGGCC | 16 | 16 | none | 61/64 |
| g55 | GCACAAGGCTCAGCTGAACC | 16 | 16 | none | 61/64 |
| g56 | GGCCCAGCCTGCTGTGGTAC | 15 | 14 | 20 | 61/64 |
| g57 | GAGGTTCACTTGATTTCCAC | 14 | 13 | 20 | 61/64 |
| g58 | GTGCGGCAACCTACATGATG | 17 | 16 | 20 | 61/64 |
| g59 | GTAGAGCGGAGGCAGGAGGC | 16 | 15 | 20 | 61/64 |
| g60 | GATGCTATTCAGGATGCAGT | 18 | 17 | 20 | 61/64 |
| g61 | GGTACCTATCGATTGTCAGG | 16 | 16 | none | 61/64 |
| g62 | GTCACCAATCCTGTCCCTAG | 15 | 14 | 20 | 61/64 |
| g63 | GAAGGCTGAGATCCTGGAGG | 15 | 15 | none | 61/64 |
| g64 | GAACAAGGTGTTCCCACCCG | 12 | 11 | 20 | 61/64 |
| g65 | GTTGGAGAAGCTGCAGGTGA | 16 | 15 | 20 | 61/64 |
| g66 | GTCTCCCTGATCCATCCAGT | 16 | 15 | 20 | 61/64 |
| g67 | GCTGCCCCGCCTGCCCAGCG | 14 | 13 | 20 | 61/64 |
| g68 | GGTGACAAGTGTGATCACTT | 17 | 17 | none | 61/64 |
| g69 | GGCAGCATAGTGAGCCCAGA | 17 | 16 | 20 | 61/64 |
| g70 | GGGATCAGGTGACCCATATT | 14 | 14 | none | 61/64 |
| g71 | GTCTGGGCGGTGCTACAACT | 16 | 15 | 20 | 61/64 |
| g72 | GGTGCACAGTGGGGTCAGCA | 16 | 15 | 20 | 61/64 |
| g73 | GTGAGTAGAGCGGAGGCAGG | 17 | 17 | none | 61/64 |
| g74 | GGATTTCCTCCTCGACCACC | 13 | 13 | none | 61/64 |
| g75 | GCCTCTCCAGCCAGGGGCTG | 14 | 13 | 20 | 61/64 |
| g76 | GAGGGCTCACCAGAGGTAGG | 15 | 15 | none | 61/64 |
| g77 | GATTTCCTCCTCGACCACCA | 14 | 13 | 20 | 61/64 |
| g78 | GGCCCCACTGTGGGGTGGAG | 12 | 11 | 20 | 61/64 |
| g79 | GATAACTACACCGAGGAAAT | 15 | 14 | 20 | 61/64 |
| g80 | GCGCCGAGAAGGAAGTGCTC | 16 | 15 | 20 | 61/64 |
| g81 | GCCCTGCTCGTGGTGACCGA | 16 | 15 | 20 | 61/64 |
| g82 | GGAGAAGGTGGGGGGGTTCC | 9 | 9 | none | 61/64 |
| g83 | GGAATCTGCCTAACAGGAGG | 14 | 14 | none | 61/64 |
| g84 | GGACCGCAGCCAGCCCGGCC | 13 | 13 | none | 61/64 |
| g85 | GCTGGGCAGGAGCCCCCTCC | 12 | 12 | none | 61/64 |
| g86 | GCAAACATGCTGTCCTGAAG | 16 | 15 | 20 | 61/64 |
| g87 | GATGTGGGAGGCTCAGTTCC | 15 | 15 | none | 61/64 |
| g88 | GGACTGAGGGCCATGGACAC | 15 | 14 | 20 | 61/64 |
| g89 | GAGACCCTGCTCAAGGGCCG | 14 | 13 | 20 | 61/64 |
| g90 | GAGCCACATTAACCGGCCCT | 13 | 12 | 20 | 61/64 |
| g91 | GTCCCCTCCACCCCACAGTG | 13 | 12 | 20 | 61/64 |
| g92 | GGGGCAGCTCCGGCGCTCCT | 14 | 13 | 20 | 61/64 |
| g93 | GCTGTGTTTGCGTCTCTCCC | 16 | 16 | none | 61/64 |
| g94 | GGGCAATGGATTGGTCATCC | 13 | 13 | none | 61/64 |
| g95 | GCTTCGGCAGGCTGACAGCC | 16 | 16 | none | 61/64 |
| g96 | GATGCAGAGACCCTGCTCAA | 17 | 17 | none | 61/64 |
| g97 | GAGAAGGTGGGGGGGTTCCA | 10 | 9 | 20 | 61/64 |
| g98 | GGCGCCCTGGCCAGTCGTCT | 15 | 14 | 20 | 61/64 |
| g99 | GACATTAAAGATAGTCATCT | 17 | 16 | 20 | 61/64 |
| g100 | GGTCCCGGTGGTGTGGGCCC | 11 | 11 | none | 61/64 |
| g101 | GGGGCCACTAGGGACAGGAT | 13 | 12 | 20 | 61/64 |
| g102 | GCTCGGGGACACAGGATCCC | 14 | 14 | none | 61/64 |
| g103 | GGCCACGGAGCGAGACATCT | 17 | 16 | 20 | 61/64 |
| g104 | GTGGGATCGGAGCAGTTCAG | 16 | 15 | 20 | 61/64 |
| g105 | GAAGCGTGATGACAAAGAGG | 16 | 16 | none | 61/64 |
| g106 | GCATTTTCTTCACGGAAACA | 13 | 12 | 20 | 61/64 |
| g107 | GGAAGGAGGAGGCCTAAGGA | 12 | 11 | 20 | 61/64 |
| g108 | GCAGCATAGTGAGCCCAGAA | 17 | 17 | none | 61/64 |
| g109 | GTGTGTGAGTATGCATCTCC | 19 | 19 | none | 61/64 |
