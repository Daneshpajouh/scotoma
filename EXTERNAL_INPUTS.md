# External inputs and the limits of this deposit

Checked on 2026-10-30 in Vancouver (2026-10-31 UTC). The mathematical
certificate packages have no external scientific-data inputs. Python packages,
a C++17 compiler, Boost headers and Node are software prerequisites, as listed
in the staging README. They are not hidden scientific inputs.

## Three assay tables

Fetch the archive from the authors' [CRISPR-Bulge repository](https://github.com/OrensteinLab/CRISPR-Bulge).
The repository README explicitly directs readers there for the LFS files.
The [datasets.zip pointer](https://github.com/OrensteinLab/CRISPR-Bulge/blob/5e21042abac4f4e3c0ff125d9ee7c8078d98b6fd/files/datasets.zip)
was checked through the read-only GitHub connector. It identifies the object
below. A normal source ZIP can contain only this pointer. Fetching raw reads
from BioProject PRJNA625995 does not supply these reprocessed aligned CSVs.

```sh
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/OrensteinLab/CRISPR-Bulge.git
git -C CRISPR-Bulge checkout 5e21042abac4f4e3c0ff125d9ee7c8078d98b6fd
git -C CRISPR-Bulge lfs pull --include=files/datasets.zip
shasum -a 256 CRISPR-Bulge/files/datasets.zip
unzip CRISPR-Bulge/files/datasets.zip \
  'datasets/CHANGEseq/include_on_targets/CHANGEseq_CR_Lazzarotto_2020_dataset.csv' \
  'datasets/GUIDEseq/include_on_targets/GUIDEseq_CR_Lazzarotto_2020_dataset.csv' \
  'datasets/NewGUIDEseq/include_on_targets/NewGUIDEseq_CR_Lazzarotto_2020_dataset.csv'
```

These are reader instructions. No clone, checkout or download of the large
archive was performed by this task, and no existing repository was changed.

| Input | Bytes | SHA-256 |
|---|---:|---|
| `files/datasets.zip` | 524400344 | `f892f70ba4ac3b05b03b2171b4ad38746630de08ad630e650f355dd61203eab0` |
| `datasets/CHANGEseq/include_on_targets/CHANGEseq_CR_Lazzarotto_2020_dataset.csv` | 878760836 | `68f3e1d7eb9e9e18e3ddc859752be8d156378a3e681b0fff329fe679366adc11` |
| `datasets/GUIDEseq/include_on_targets/GUIDEseq_CR_Lazzarotto_2020_dataset.csv` | 381626484 | `a3398cba2ef452ec37ebebda38a44aa62615c66165a43699691ca210f467bcd2` |
| `datasets/NewGUIDEseq/include_on_targets/NewGUIDEseq_CR_Lazzarotto_2020_dataset.csv` | 207094458 | `ac9ca55a384a9af5ba37d667766e6723a7c211b5fa1d54875e450db23a13a23b` |

The CSV sizes and digests come from the existing pinned, assay-specific
`verification/bulge_shares/README.md` and checker. The archive object ID and
size were freshly confirmed. The large CSV bytes were not fetched or freshly
hashed in this task. Do not substitute the similarly named tables outside
`include_on_targets/`.

With these files, `verification/bulge_shares/share_strata.py` checks each
assay's shares. The CHANGE-seq file additionally supplies
`code_release/make_sites129.py --check`, `verify_missed_activity.py`,
`validate_terminal_law.py`, `make_rna_position_series.py`,
`read_depth/matched_bin.py` and `verification/cfd_gap_asymmetry/cfd_gap_arms.py`.
`make_sites129.py`, `verify_missed_activity.py`, `validate_terminal_law.py`
and `make_rna_position_series.py` read the CSV on stdin. The read-depth
script reads `CHANGESEQ_CSV`. None can reconstruct the missing source table
from `sites129.json`: that is a selected, annotated subset.

`verify_fpbudget_positional_propensity.py` additionally reads
`datasets/CrisprNet.csv` from this same archive, given as its path argument
or `CFD_DATASETS_ZIP`. It is a retained legacy producer rather than a result
counted in the current paper ledger, and is not one of the three tables above.

The complete archive and CSVs are excluded. The existing data licence covers
our derivation and annotation; it does not grant new rights over the underlying
third-party measurements. This task does not infer a redistribution licence
from public availability.

## CFD lookup tables

`code_release/verify_cfd_gap_table.py` requires `mismatch_score.pkl`, given as
its first argument or `CFD_MM_PKL`. Fetch the tables from the authors' pinned
CRISPRme source:

| Input | Source | Bytes | SHA-256 |
|---|---|---:|---|
| `mismatch_score.pkl` | [PostProcess/mismatch_score.pkl](https://github.com/pinellolab/CRISPRme/blob/170aeaf1e2e0e8c2379aa14e49d73c4246150ed1/PostProcess/mismatch_score.pkl) | 11199 | `eba85c300e484207e5cfe67b8717146cb49e2ee24a4b1b2a4af4960baaa891a5` |
| `PAM_scores.pkl` | [PostProcess/PAM_scores.pkl](https://github.com/pinellolab/CRISPRme/blob/170aeaf1e2e0e8c2379aa14e49d73c4246150ed1/PostProcess/PAM_scores.pkl) | 345 | `48bd31f0159382836d0648987ac4f2444d7509f2d04513b9c6f040b3e56049de` |

Both were freshly retrieved and hashed into this task's private audit area,
`sources/`. Neither is in staging. The structural checker uses its restricted
unpickler. The released four-arm CFD script embeds its scoring constants and
instead needs the external CHANGE-seq CSV. Verifying the table's structure
does not reproduce an AUROC comparison.

## References and assemblies

Fresh source checks used the official
[HPRC Release 2 assembly documentation](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/README.md)
and [assembly index](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/2a7e3a5b1ca10454762699176b1a0a1027ccbf21/data_tables/assemblies_release2_v1.0.index.csv).
The index contains 466 rows: 464 donor haplotypes from 232 samples, plus
GRCh38 and CHM13. Fetch the S3 objects named in the `assembly` column and
check their `assembly_md5` sidecars. An example is:

```sh
aws s3 --no-sign-request cp \
  s3://human-pangenomics/working/HPRC/HG00235/assemblies/release2/HG00235_hap1_hprc_r2_v1.0.1.fa.gz .
```

The source also lists the original reference files:

- GRCh38: `s3://human-pangenomics/working/HPRC_PLUS/GRCh38/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`.
- CHM13: `s3://human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz`.

The recorded GRCh38 census uses 25 primary chromosomes, including mitochondria,
with 3088286401 bases and 226176633 distinct canonical NGG windows. Filter by
the actual chromosome records, not by a filename substring. A download of a
different primary/no-alt/decoy scope is not an identical input. The CHM13
producer metadata is incomplete (`meta_unreadable`); 227157495 universe codes
are recorded, but the exact historical FASTA digest is not deposited. The
official reference URLs are acquisition routes, not a claim that a fresh
download was byte-matched to that historical CHM13 run.

`python3 -m scotoma.cli universe` and `scan` operate on reader-supplied FASTAs.
The primary references and assemblies are omitted because of size and because
the task calls for external acquisition, not because every HPRC FASTA is
necessarily access-controlled or forbidden to redistribute. The deposited
per-guide summaries permit reanalysis of their counts; they do not rerun a
whole-genome enumeration.

## The ancestry cohorts are not specified well enough for an identical replay

Paper A distinguishes 252 haplotypes (181 donors) from a separate 232-donor,
single-copy-autosomal, diploid analysis. They must not be interchanged.
The public acquisition route for the assemblies is the pinned HPRC index above;
sample labels are in the official
[hprc_release2_sample_metadata.csv](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/sample/hprc_release2_sample_metadata.csv).

For an identical 252-haplotype replay the reader additionally needs the original
252-row inclusion manifest, its donor/haplotype/file/label mapping, the assembly
and reference checksums, and the exact analysis input set. Only 97 local
`ABS_*.json` records were found in `results_2026-10-20/hapabs/`; the original
`hap_manifest.tsv` is absent here. Neither the missing 155 records nor the
full inclusion manifest is invented or deposited. `analyse_hapabs.py` is local
but not part of the existing public package; its historical decision-rule
wording is corrected in the project notes, so it is not silently promoted
into a new public analysis entry point.

For the 232-donor result the reader needs all original per-donor projection
inputs, the exact 206133695-site single-copy-autosomal universe, its construction
and checksum, and the original analysis and label joins. Those inputs are not
in this tree. Downloading every Release 2 assembly is insufficient to reproduce
the same predicates and analysis without those missing specifications.
No exact fetch list for the missing cohort subsets can honestly be supplied.
Both ancestry results remain unreproducible from this deposit, even after
following its currently available acquisition links.

## Other external inputs and retained evidence

- Re-running native comparator pipelines requires the pinned upstream tools,
  reference inputs and command harnesses named in `code_release/tool_audit/REPORT.md`.
  Five deposited audit files can be checked by hash; that is not a fresh native
  execution. Several linked audit subdirectories and the combined internal
  audit JSON are absent from the release.
- The benchmark derivation README names `CRZ_DIR` from `rerun_benchmark.sh`
  kit v5. Its two complete CRISPRitz passes (about 208 MB) and the kit are not
  included. `expected_crispritz.sha256` pins their expected files but cannot
  regenerate them. The shipped SWOffinder fixtures and derived tables support
  only the stated partial replay.
- Maize B73 AGPv4 is the external reference behind the maize scan. The emitted
  and exclusive locus lists are shipped and their set arithmetic can be rerun;
  rerunning the original genome scan also needs the reference and original
  scan inputs. A list-level recount is not a new alignment or genome scan.
- The five Tracking-seq2 sequence windows are embedded in the cellular package.
  Their sequence-distance check is self-contained. Rechecking biological
  validation requires the authors' supplementary source, especially
  `41467_2026_76778_MOESM4_ESM.xlsx` for
  [Tracking-seq2](https://doi.org/10.1038/s41467-026-76778-9), plus the
  independent targeted-amplicon evidence. Those experiments are not rerun here.
- `code_release/verify_tolerance_bound.py` also names a JSON of 459
  per-haplotype shares. It is not shipped and no exact acquisition manifest
  is recorded. Its method-only check cannot reproduce the five empirical rows.

No raw genome, ancestry dataset or full third-party assay table is in staging.
