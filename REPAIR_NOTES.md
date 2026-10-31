# What changed in the local deposit candidate

This task did not edit the source repositories, papers or abstracts. All changes
below exist only under `staging/scotoma/`. The final audit detected concurrent
changes to paper/build files and the source release README and manifest.
`SOURCE_PRESERVATION.json` records that drift without inferring a writer or
reverting any file. Scientific source code and certificate payloads remained
unchanged; the final claim ledger was checked against the revised references
and bound to `PAPER_SOURCE_AT_AUDIT.json`. The original local Scotoma base
matched every blob in the freshly retrieved public Git tree
`0de5b8d3d7a30fcece66c8ebc99899ade73dfa7e`. The public tree had no
`code_release/`. `sources/BASE_DIFFERS_FROM_PUBLIC.json` records zero differences.
No remote write, commit, push or pull request was performed.

## Reader entry point and dependency closure

- The root README now names the two current registered papers and the correct
  `code_release/` entry point. The detailed enumerator manual is retained under
  `docs/SCOTOMA_USAGE.md`.
- `reproduce.py` verifies the complete inventory and works in a disposable copy.
  It clears scientific input overrides, bounds numerical-library threads and
  creates writable temporary matplotlib state. The original package commands
  remain visible in its command map and execution log.
- `requirements.txt` pins the tested Python dependencies. The pinned NumPy
  version needs Python 3.12 or later; the validation interpreter was 3.13.9.
  Individual standard-library checkers may support older Python versions.
- `upper_rows/extend.py` imports python-flint only inside the producer's exact
  reconstruction routine. The verification route never calls that routine or
  its optional highspy optimizer. All seven upper-row certificate replays and
  the transport audit ran successfully without either package installed.
- `data/panel110.txt` and `data/panel_match.json` are lossless views of the
  already released 110 spacers. Two structure scripts now read this included
  panel rather than an absent sibling `multiguide_20261017` directory.
- `reproduce_external.py` checks the exact byte count and SHA-256 of the
  reader's own input before an assay or CFD-table analysis. It does not fetch
  or redistribute those inputs. Missing input and wrong input are tested failures.

## Exact budget-ten native replay

The original native precomputation repeated a window/far-arm convolution for
every shell. The candidate groups shells at an anchor by `(t,u,fl,fr)` and sums
their coefficient-weighted window polynomials before convolving with the shared
far-arm distance distribution. Distributivity gives exactly the same anchor
table. The positive certificate sums and the existing signed-128-bit final
accumulator remain in use; the far convolution also uses signed-128-bit integers.
The dynamic program, support definitions, normalizers and certificate vectors
are unchanged.

The original precomputation is retained in the same source, selected with
`DP_REFERENCE_PRECOMPUTE=1`. The candidate's startup calibration compares both
paths against literal four-letter enumeration on nine instances, 265728 centres
in total, including free guards, binary flips, stars and far arms. It also
checks the maximizing witnesses. All nine complete deposited certificates
were replayed by both implementations, with identical exact maxima, witnesses,
masses and class counts. Both rejected the overloaded and outside-support
controls. The original full replay took about 492 seconds and the grouped
candidate about 98 seconds, including calibration, in the observed workstation
runs. These are single-run observations under the recorded load.

The complete original transcript is in the preparation audit at
`verification/reference_m10/`; the candidate transcript is
`verification/m10_closure.log`. This is a computation optimization, not a new
mathematical result.

## Failed checks and their disposition

- macOS `/usr/bin/time -l` failed at the sandbox-blocked `sysctl kern.clockrate`
  after successful verifier commands. The wrapper failure is retained under
  `verification/timer_failures/`. Command exits and RSS were recorded inside
  the reader runner using `RUSAGE_CHILDREN`. Later runs use `wait4` directly.
  A timing-wrapper failure is never presented as a failed scientific result.
- The legacy index gate still expected budget-nine reach 15 and 61 panel
  guides, although its included individual-cell index reaches 16 and 89.
  Its expectations and labels now describe that historical index subset.
  This subset does not supply the current universal budget-nine theorem;
  the full-weight base and the five closure certificates do.
- The initial dials command list incorrectly invoked the full-source RNA
  positional producer without its stdin CSV. It failed. That command is now
  in the explicitly external CHANGE-seq route, while `dials` contains only
  the three self-contained derived-summary producers. This does not credit
  the positional measurement as reproduced from the tree.
- The initial path-structure check exposed an absent sibling panel path.
  It now reads the included lossless panel view and was rerun successfully.
- The direct grid wrapper now rejects an empty or incomplete 432-row inventory,
  requires 412 direct pairs in a full replay, and rejects an empty selection.
  The gap wrapper requires exactly the thirteen declared certificate pairs
  plus its separate calibration fixture. Empty inventories cannot report success.
- A new standard-library checker independently exhausts the ACGTA soundness
  example. It distinguishes the unrestricted four-ball cover from a minimum
  five-ball sound cover over all 256 ambient words.

The initial invocation failures remain under `verification/initial_failures/`.
The negative-input test deletes a scientific input from each of the eight
mandatory packages, empties the direct grid and gap inventories, and tests
missing and wrong external files. Each case must fail as required.

## Packaging and prose

Git metadata, backups, caches, machine executables and an old historical claims
ledger are excluded. The five large threshold polynomial outputs are regenerated
from the included family definitions and are not deposited as oversized files.
All referenced finite boundary certificates remain included.

Eleven existing proof/report Markdown files were transcribed to ASCII in the
candidate. Greek names use their spelled-out forms; `Q` denotes the rationals;
floor, ceiling, relations and subscripts use ASCII notation. No certificate
number, sequence or numerical data file was changed by this transcription.
`ASCII_TRANSCRIPTIONS.json` records the files and original character classes.
The tool-audit report's manifest was updated for this typography change; that
manifest check remains a check of deposited evidence, not a native rerun.

The root SHA-256 inventory covers every proposed deposited byte. Nested
inventories were regenerated to match this candidate. `.gitattributes` preserves
all deposited bytes across line-ending configurations. Fresh display logs, when
included, normalize only local workspace paths and prose typography; the
preparation audit retains the original logs. No full third-party assay table,
raw reference assembly, ancestry input or external CFD pickle is in staging.
