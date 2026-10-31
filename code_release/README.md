# Exact pattern-cover certificates and derived tables

This is code_release/ inside the Scotoma reproducibility snapshot.
Start one directory above, with README.md and reproduce.py. The central checks are:

```sh
python3 -B reproduce.py m9_full_weight
python3 -B reproduce.py m8_closure
python3 -B reproduce.py m9_closure
python3 -B reproduce.py m10_closure
python3 -B reproduce.py m10_r11
python3 -B reproduce.py infinite_families
python3 -B reproduce.py thresholds
python3 -B reproduce.py private_sphere
```

Run these from the snapshot root, not from this directory. Individual package
READMEs also give their direct verifier commands relative to code_release/.
The wrapper checks inventory hashes and runs a disposable copy so generated
reports do not change this deposit. It requests the complete measured threshold
grid. Core checkers do not run an optimizer.

CLAIM_COVERAGE.tsv at the root states which result groups can be checked from
this tree and which need external inputs. VERIFICATION.md distinguishes fresh
replay from stored evidence. The source-table, genome and cohort limitations
are in EXTERNAL_INPUTS.md. Producer and legacy gate files are retained for
inspection, but the supported fresh replay targets are exactly those named
in reproduce.py. Stored transcripts are not new runs.

data/panel110.txt and data/panel_match.json are lossless guide-list views of
data/producers_2026-10-20/PANELREF_GRCh38_NGG.json. No assay rows were added to
make them. data/LICENSE-DATA.md governs the derived data.
