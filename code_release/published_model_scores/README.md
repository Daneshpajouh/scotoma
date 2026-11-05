# What a published scorer does with the omitted alignments

This directory backs one paragraph of the companion CRISPR paper: the published CRISPR-Net
checkpoint, run as published, scored on the alignments a released bulge-aware enumerator cannot
return, with a matched control from the alignments it does return.

    validation/   the checkpoint reimplemented from its deposited weights, and the check that it
                  reproduces published behaviour on the authors' own labelled data BEFORE use
    aligned/      the scoring run on the deposited genomic alignments, with its matched control,
                  its encoder-convention audit, and its corruption control

## To check it, run one thing

    ./check.sh

That runs the seal check, then the corruption control, then the full replay, in the order they
depend on. A reviewer does not need to know the order; the reason it matters is at the end of this
file.

## Read this order

1. `validation/FINDINGS.md`. The forward pass is NumPy, read from the deposited weights with
   h5py; no historical inference framework is installed. It is not used until it reproduces the
   checkpoint's behaviour on the checkpoint's own data, and the one place it does not reproduce a
   published figure is stated there rather than omitted.
2. `aligned/ENCODER_POLICY.md`. The published encoder distinguishes a biological gap from
   zero-vector padding. That convention is pinned before anything is scored, audited against the
   original encoder on every tensor, and a negative control confirms that passing the wrong gap
   character through fails.
3. `aligned/FINDINGS.md`. The scores, every one of them beside its control.

## The public inputs are named and hashed, not re-hosted

Four of the authors' own labelled files, 86 MB together, are not copied here. They are published
by the model's authors and are pinned by SHA-256 in `validation/pins.json`, whose own digest is
embedded in both the producer and the verifier. Download them from the model's own repository at
the commit named in `validation/sources/` and the verifier checks each digest before reading it.

| file | SHA-256 |
|---|---|
| `CIRCLE_seq_10gRNA_wholeDataset.csv` | `86d88947ef9c0a1b545814c4e4241e4f6103d15a4d056d697a4ea5b6017648a8` |
| `Listgarten_22gRNA_wholeDataset.csv` | `74fc9151fc175cf48aaf9712eeb1216ac42306a72453f041b17e8d6b9389097b` |
| `elevation_6gRNA_wholeDataset.csv` | `56f1c95a0a3d61b873af984f7a2110e27538540b49904b913982d63764dcbb55` |
| `Kleinstiver_5gRNA_wholeDataset.csv` | `2bd3860b7a78af42b1c1a29c13c16263d519e54b5b7829921d9e4fc8bba81049` |

The checkpoint itself is pinned the same way, at
`a7f6aa381520f5c68fa1f099a6ef3ebc3b8ce846709b97dfde2053f26ca62f80`. The per-row score dumps for
those four validation files are also not copied; their digests are pinned and the verifier
regenerates them from the public inputs. Everything the paper's paragraph quotes is in
`aligned/results/`, which is complete.

## What is claimed, and what is not

These are in-silico nomination scores on retrospective inputs. The deposited guide column carries
`NNN` where a PAM would be, and the published encoder copies the guide's PAM letters from the DNA
side, so PAM columns are matches by construction and these scores do not estimate PAM
specificity. Treatment and control are matched by guide and by mismatch count and by nothing
else, so the genome-wide comparison is descriptive. On the measured-active rows at four
mismatches the two arms are close, and that is reported as the null it is. No cleavage outcome,
clinical risk, ranking impact or prospective validation is claimed.

`FINDINGS.md` at this level is the aligned scoring report, where `aligned/verify.py` expects to
find it; `aligned/FINDINGS.md` is the same file, kept so that directory reads on its own.

## One change made to this copy before it was deposited

Five log files recorded absolute paths from the machine that produced them, and a 124 MB SQLite
scratch file from a verifier run had been left in the tree. The scratch was moved out, and in
those five logs the producing machine's directory prefix was replaced by `<repo>/` and its home
directory by `<home>`; nothing else in them was altered. `aligned/SHA256SUMS` was therefore
resealed with `aligned/seal.py --write` after those edits, so the manifest describes this copy and
`aligned/seal.py --verify` passes on it. The score files, inputs, producer and verifier were not
touched.

The repository's top-level `SHA256SUMS` does not list this subtree. It has no writer script, and
hashes are not entered by hand here, so this directory seals itself instead: `aligned/SHA256SUMS`
is written and checked by `aligned/seal.py`, and `validation/SHA256SUMS` covers the
reimplementation. Run `aligned/seal.py --verify` and `aligned/verify.py` to check both.

## Changes made to this copy before it was deposited, with the reason

The producing run recorded the machine's own absolute paths in five places: `local_sources.json`,
which is the provenance map of where each input was copied from, and four logs. In this copy the
project prefix is rewritten to `<repo>/` and the home directory to `<home>`; nothing else in those
files changed, and no scored input, score, score table or weight was touched. Rewriting
`local_sources.json` changes its digest, so it was re-pinned by computation in `aligned/pins.json`
and `aligned/INPUT_SHA256.tsv`, and `pins.json`'s own digest was updated where `aligned/integrity.py`
embeds it. `local_sources.json` went from `a75d777664b5...` to `e332d217d8af...` and `pins.json`
from `2db76bd46dec...` to `58fd64f18484...`; every other pinned digest is unchanged, including the
checkpoint's.

`aligned/corruption_control/changed_input/` is generated by `aligned/corruption_control.py`, which
writes it with absolute symlinks into the producing directory. Those would dangle for anyone who
cloned this repository, so each was replaced by a real copy of the deposit's own corresponding
file. After that, `corruption_control.py` reports `status PASS` with the corrupted copy rejected at
exit 1 before scoring and the unmodified copy accepted at exit 0, `seal.py --verify` passes, and
`verify.py` passes at exit 0 having replayed 70,142 scores to 1.6e-15. Running the control again
locally will rewrite its own log with your absolute paths; that is the log, not the evidence.

## Replay order, and why the corrupted copy is not shipped

`aligned/corruption_control.py` builds `aligned/corruption_control/changed_input/` itself, as
symlinks into the directory it runs in. Symlinks in a published tree would dangle for anyone who
cloned, and copying them in instead makes the control mutate a clean file and stop rejecting
anything, so the generated tree is not deposited. The generator and its receipt are, and
`aligned/reproduce.sh` already runs them in the right order. Run either

    aligned/reproduce.sh

or, for the checks alone,

    python3 -I aligned/corruption_control.py && python3 -I aligned/verify.py

`verify.py` on its own will stop with a missing `changed_input` until the control has built it;
that is the order, not a fault. Measured here before deposit: the control reports `status PASS`
with the corrupted copy rejected at exit 1 before scoring and the unmodified copy accepted at exit
0, and `verify.py` then exits 0 having replayed 70,142 scores to 1.6e-15.
