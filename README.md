# Scotoma and the exact pattern-cover certificates

This snapshot accompanies:

- "Exact minimum pattern sets for bulge-aware CRISPR off-target search".
- "The exact reach of fractional certificates for single-deletion pattern covers".

Authors: Amirhossein Daneshpajouh, Megan Fowler and Kay C. Wiese.

Start in the directory containing this README. `scotoma/` contains the
enumerator; `code_release/` contains the mathematical certificates and derived
tables. The complete scientific certificate inputs are included. Large
third-party assay tables, reference assemblies and ancestry inputs are not.
See `EXTERNAL_INPUTS.md` and `CLAIM_COVERAGE.tsv` for the precise limits.

## Verify the snapshot and replay certificates

Use Python 3.12 or later for the pinned full suite; validation used Python
3.13.9. NumPy and a C++17 compiler cover the central certificates.
The broader suite also uses SciPy, matplotlib, pytest, Node 20 or later and
Boost headers. `requirements.txt` pins the Python versions used for this
validation. No optimization solver is needed for the exact checkers. The
upper-row producer's optional python-flint/highspy dependencies are not needed
by its verification route.

Native sources use `__int128` and POSIX threads in addition to C++17. Use a
compiler and environment that provide those extensions; the validated compiler
was Apple clang 21.

```sh
python3 -m venv ../scotoma-env
../scotoma-env/bin/python -m pip install -r requirements.txt
../scotoma-env/bin/python -B reproduce.py integrity
../scotoma-env/bin/python -B reproduce.py m9_full_weight
../scotoma-env/bin/python -B reproduce.py m8_closure
../scotoma-env/bin/python -B reproduce.py m9_closure
../scotoma-env/bin/python -B reproduce.py m10_closure
../scotoma-env/bin/python -B reproduce.py m10_r11
../scotoma-env/bin/python -B reproduce.py m10_r11_kernel
../scotoma-env/bin/python -B reproduce.py infinite_families
../scotoma-env/bin/python -B reproduce.py thresholds
../scotoma-env/bin/python -B reproduce.py private_sphere
```

The runner checks the complete SHA-256 inventory first and executes in a
disposable copy. Missing or changed inputs fail before a verifier starts.
Reports, native binaries and regenerated polynomial records go into that copy,
which is removed on exit. The deposit is left unchanged. Set `TMPDIR` to a
writable location with sufficient free space. The full suite is:

```sh
../scotoma-env/bin/python -B reproduce.py all
```

`all` runs every self-contained target named in `reproduce.py`. The threshold
target requests all 412 direct boundary pairs, as well as the infinite-family
and general-structure algebra and controls. The 20 formula-derived table rows
are distinct from those 412 direct measurements. The five large generated
`thresholds/infinite_r*.certificate.json` records are deliberately omitted:
`prove_family.py` regenerates them from the supplied small family inputs in
temporary space. Every required finite boundary input remains included.

The mass-twenty radius-nine replay was measured over all 1162261467 classes;
it needs about 4.7 GB of peak process RSS and about 4.9 GB of sampled temporary
disk including the runner's copy. Allow at least 6 GB for each. The upper-row
suite was observed near 8.1 GB peak process RSS; allow at least 10 GB for the
complete serial suite. These are workstation observations, not cluster sizing
measurements or bounds on arbitrary future runs. `VERIFICATION.md` records
every measured command, verdict, wall time and memory scope.

The grouped shell precomputation in the budget-ten native checker retains the
original path under `DP_REFERENCE_PRECOMPUTE=1`. Its startup calibration compares
both paths with literal four-letter enumeration. The complete certificate
replay must still recover the deposited maximum in all nine cells and reject
the planted invalid inputs. `REPAIR_NOTES.md` describes the local changes.

For native upper-row compilation, `CXX` may name a C++17 compiler and
`BOOST_INCLUDE` the directory containing `boost/`. The runner recognizes the
usual Homebrew include directory when present. On this validation workstation,
SciPy work used `/usr/local/bin/python3`.

## Enumerate a reader-supplied assembly

```sh
../scotoma-env/bin/python -m scotoma.cli --help
../scotoma-env/bin/python -B reproduce.py enumerator
../scotoma-env/bin/python -m scotoma.cli universe \
  --fasta GRCh38.fa.gz --label GRCh38 --out-npz GS_GRCh38.npz
../scotoma-env/bin/python -m scotoma.cli scan \
  --assembly haplotype.fa.gz --universe GS_GRCh38.npz \
  --panel code_release/data/panel110.txt --out-json scan.json
```

The synthetic suite has 142 tests. Passing it checks the supplied software;
it does not reproduce the paper's large-genome measurements without their
external inputs. `docs/SCOTOMA_USAGE.md` retains the detailed enumerator manual.

## Reproduce external-data analyses

After acquiring and checking the exact inputs in `EXTERNAL_INPUTS.md`:

```sh
../scotoma-env/bin/python -B reproduce_external.py changeseq /path/to/CHANGEseq_CR_Lazzarotto_2020_dataset.csv
../scotoma-env/bin/python -B reproduce_external.py guideseq /path/to/GUIDEseq_CR_Lazzarotto_2020_dataset.csv
../scotoma-env/bin/python -B reproduce_external.py newguideseq /path/to/NewGUIDEseq_CR_Lazzarotto_2020_dataset.csv
../scotoma-env/bin/python -B reproduce_external.py cfd-table /path/to/mismatch_score.pkl
```

These commands reject missing inputs and wrong sizes or hashes. The three
large assay CSVs were not available during this snapshot's local validation;
their end-to-end analyses are not represented as fresh passes. The read-depth
comparison and full source-table selection require CHANGE-seq. The ancestry
analyses additionally lack complete historical cohort specifications. No
command here claims to regenerate those missing specifications.

## Licence

Code: MIT, under the supplied `LICENSE` files. Derived tables and annotations:
the terms in `code_release/data/LICENSE-DATA.md`. Third-party underlying data
retain their source terms; this snapshot grants no new rights over them.
