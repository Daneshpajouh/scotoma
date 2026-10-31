# Data manifest and provenance

Written 2026-10-22 so that every stage of the empirical chain is checkable **without** the 878 MB
source table. A referee reaching `sites129.json` can otherwise verify 260 to 129 to 40 and nothing
upstream of the 260.

## Files, with SHA-256

| file | bytes | sha256 |
|---|---|---|
| `bulge_census_full_idx_dna.json` | 139810 | `34b448e31f17b23dbc176bbfbed2c552567610c4d1d6aec03687754aa57c6eb8` |
| `exacel_enumeration_summary.json` | 4561 | `ce4a22d8d5da39312e01a8ef36581368f111440354aa992ae40afcb6757203bd` |
| `exacel_windows.tsv.gz` | 344909 | `a0a251cc338963af8d24481f8b7bac91039d3aad8330e0b990c4a1c479c28072` |
| `minimality_cert.json` | 81818 | `c6fddb15bc52798113b013e3763c4947aac76fb7553c279ad50bd27895a997cd` |
| `sites129.json` | 89301 | `12db5cb98cddc3ef2f5125cbeb121457869398909234bc07f6d77ddc6a64bc9a` |

## `sites129.json`: what it is and where each stage comes from

| stage | value | how it is obtained |
|---|---|---|
| source deposit | CHANGE-seq tables distributed with CRISPR-Bulge | third party, cited in the paper, **not redistributed here** |
| source table | `CHANGEseq_CR_Lazzarotto_2020_dataset.csv`, 878,760,836 bytes | cluster path in the capture script |
| rows in the source | 4,936,278 | one pass, counted |
| measured-active rows (`reads > 0`) | 202,040 | filter `reads > 0` |
| active **off-target** rows | 201,931 | the 202,040 less 109 on-target rows |
| bulge-bearing measured actives | 23,692 | `Align.#Bulges >= 1` with the gap count agreeing |
| **records in this file** | **260** | the single-bulge measured actives at the four-mismatch budget that the complete pattern set nominates |
| exclusively missed by the 2016 wrapper's position set | **129** | of those 260 |
| distinct guides carrying them | **40** | of the 110-guide panel |

**Record-level invariants, all holding on all 260 records**, and any of them failing means this file is
not the file the paper describes:

- `deleted_position` is **19 for every record**, so every site is at the RNA terminal and none at the
  DNA slot.
- `aligned_target` with its gap removed equals `target`.
- The gap index in `aligned_target` equals `deleted_position`.
- The target occurs inside its own 73-base window with the orientation the strand field gives.

**THE CHAIN IS RE-EXECUTABLE, AND THE SELECTION STEP IS NAMED.** `../make_sites129.py` writes this
file from the deposit. Run end to end against the copy held here it reproduces it **set-identically**,
by a SHA-256 over the sorted record keys (chrom, aligned sgRNA, aligned target, mismatch count,
reads) rather than by agreement of counts:
`dcc09aa424db55c216e356619f89a29120c737fc6f9fcd9b2bf1fb3aa9b0e479`, 260 records on both sides.

The selection is **two** requirements on a single-RNA-bulge measured-active row within the
four-mismatch budget, and both are needed:

1. the deleted guide base sits at the **terminal spacer position** -- this alone selects
   **583 rows over 76 guides**; and
2. the guide is **at risk under the terminal law**, its last two spacer bases differing -- this is
   what takes those 583 to the **260 over 48 guides** in this file.

The 323 rows that requirement 2 removes are the guides the law rescues, and $323 + 260 = 583$
exactly. Requirement 2 was already shipped, inside the one-line conjunction at
`../verify_missed_activity.py` line 26; what was missing was a producer and a manifest line pointing
at it, both of which now exist. Applying requirement 1 without requirement 2 yields 583 and looks
like a reproduction failure, which is why they are separated here instead of left in a conjunction a
reader has to split. History in
`RESOLVED_2026-10-22_the_260_reproduces_and_the_filter_was_shipped_all_along.md`.

So a referee can now verify the whole chain from this folder: deposit to 260 by re-execution, and
260 to 129 to 40 by the three scripts that read this file. The stage table above and the per-file
hashes make it auditable; `make_sites129.py` makes it re-executable.

## Licence, by file

Every file listed in the table above, and the generated tables and figures in this repository, are
licensed **CC BY 4.0**, https://creativecommons.org/licenses/by/4.0/, with the exception noted next.

Copyright (c) 2026 Amir Daneshpajouh, Megan Fowler and Kay C. Wiese. Please cite the paper.

**Derived-data exception, stated explicitly.** `sites129.json` is **derived from a third-party
deposit** (the CHANGE-seq tables distributed with CRISPR-Bulge). We license our derivation and
annotation of it under CC BY 4.0; **the underlying measurements remain under their own deposit's
terms**, the deposit is cited in the paper, and none of it is redistributed here. `exacel_windows.tsv.gz`
and the enumeration summary are derived from the public GRCh38 assembly, which is likewise not
redistributed.

## Producer inputs shipped with the release, added 2026-10-22

Until 2026-10-22 five scripts in this folder could not run from it: their data paths were bare
relative directories resolved against the working directory, and the inputs lived in the paper tree
rather than here, so a reader who received only this folder could not run them from anywhere. The
inputs below now ship with the release, 2.9 MB in total, and every path resolves against the script.
They are **derived tables produced by our own enumeration**, carry no third-party records, and fall
under the same derived-data terms as the rest of this folder: our derivation is CC BY 4.0.

| file | bytes | SHA-256 |
| --- | --- | --- |
| `producers_2026-10-20/BULGE_CENSUS_size1_FULL_mm4.json` | 25,952 | `558c73addedcdc7648ab9660c2fe5bb02e47315f9bec5f44f3c1340091fb38d6` |
| `producers_2026-10-20/BULGE_CENSUS_size1_mm4.json` | 25,953 | `f82b8b81703fea05b23e3209d3be8b30cc3ae1d5ea73173e24bc723374d36911` |
| `producers_2026-10-20/BULGE_CENSUS_size2_mm4.json` | 26,423 | `75ffda6c49aee7376474fefe2f0c7472b2e50e0ff8969c6498d13b96b3febcdd` |
| `producers_2026-10-20/PAM_PANEL_FACTOR.json` | 731 | `74ce148eb80eb6c51501bc43f4cfc256fb367db2136c95dcf0e0379ecd3284d4` |
| `producers_2026-10-20/PAM_PANEL_FACTOR_PERGUIDE.json` | 53,946 | `12000f2d9aeb68dce8eabba363db3e6fb3fa07207bb50d3b6f2fe4b5a66e7273` |
| `producers_2026-10-20/PANELREF_CHM13_NGG.json` | 41,477 | `cf2490226dc9fadaca2659e2e41706aff8f1982437827f55b48d8968c336e266` |
| `producers_2026-10-20/PANELREF_GRCh38_NGG.json` | 41,819 | `6f0d396789d575959b1907bf710ce2fa80e3352c5593636d69620773ff82074b` |
| `producers_2026-10-20/TAXONOMY_UNION_mm4.json` | 2,401 | `37ffc524b3f5ce3234d9474070407c5eba87deae6c72afcb33ba88de2a022216` |
| `producers_2026-10-20/bulge2_census_idx_dna2.json` | 142,010 | `b346945342d1b4a40abb61b675fa42ccfaab15b5f998d9c3d0efbde3ff43d8a8` |
| `producers_2026-10-20/bulge2_census_idx_rna2.json` | 94,997 | `13b1ab2c9fa7e5e46ba92b96aa0e9d535a4b9aac6cf5c62728b4d88ae0cf9b5e` |
| `producers_2026-10-20/bulge_census_full_idx_dna.json` | 139,810 | `34b448e31f17b23dbc176bbfbed2c552567610c4d1d6aec03687754aa57c6eb8` |
| `producers_2026-10-20/bulge_census_full_idx_rna.json` | 100,361 | `8f8c03c4971656cba9dbf34ab8f37ec38955872a52d777d36f48ae6bc2264d23` |
| `producers_2026-10-20/bulge_census_idx_rna.json` | 91,673 | `43b9d296f96a03e64fcdf9b69557b2f15d8f111ef8d3e11162944bc582c1d9d0` |
| `producers_2026-10-20/bulgemix_census_idx_mixed.json` | 1,996,822 | `97a9e72d155676ad8a805f15c2d2742f21ceaf7496337781895f6a5b5a190070` |
| `producers_2026-10-20/idx_dna.json` | 132,770 | `9169477ab549201b29bddb35e845b81ac9ef6808969d71855e05f966c7862ef1` |
| `mm6_pending/BULGE_CENSUS_mm3.json` | 24,591 | `39c6272b18e2bcc44d1cb68c8b1a0ef442c94bd9bb9b824c337bd37281fd270d` |
| `mm6_pending/BULGE_CENSUS_mm4.json` | 25,952 | `558c73addedcdc7648ab9660c2fe5bb02e47315f9bec5f44f3c1340091fb38d6` |
| `mm6_pending/BULGE_CENSUS_mm5.json` | 26,634 | `5976d9e902f5b1e82fde3dff2ebbe61b7cea979e38f9e2b26b92dd7549199d68` |
| `mm6_pending/BULGE_CENSUS_mm6.json` | 27,162 | `81a2946f5423a1413c233d0c796113305140b7318180a815ec0f114bfaa09761` |

**Two files appear twice under different names, stated rather than left for a reader to notice.**

- `bulge_census_full_idx_dna.json` and `producers_2026-10-20/bulge_census_full_idx_dna.json` are byte-identical (`34b448e31f17b23dbc176bbfbed2c552567610c4d1d6aec03687754aa57c6eb8`).
- `mm6_pending/BULGE_CENSUS_mm4.json` and `producers_2026-10-20/BULGE_CENSUS_size1_FULL_mm4.json` are byte-identical (`558c73addedcdc7648ab9660c2fe5bb02e47315f9bec5f44f3c1340091fb38d6`).

Both duplications are the paper tree's own history, not a packaging error: the four-point dial reads
its four budgets from a directory the cluster jobs delivered into, and one of those budgets is the
same table the single-bulge census produced under its own name. We ship both paths rather than
rewrite the scripts to share one, because changing which file a gate reads is a larger change than
shipping 26 KB twice, and a referee comparing hashes should find the coincidence documented.

## Benchmark deposit, added 2026-10-23

Two packages the appendix's benchmark figures are read from. The loci package carries the
per-guide comparison the scope A and scope B counts come out of; the derivations package
recomputes the figures that had no deposited source from pinned inputs, and carries the
script that does it so a reader can rerun rather than trust the table.

| file | bytes | SHA-256 |
| --- | --- | --- |
| `benchmark-loci-deposit-20261023/README.md` | 2,354 | `1bf573fe5e3260e2eb89093eb2d9e54ae786f10c799e798b16d90606730da426` |
| `benchmark-loci-deposit-20261023/SHA256SUMS` | 820 | `e704dcd40a6841b8bdfeb6e887e39849c7e2750e150f3ce115aa2eab00e75f57` |
| `benchmark-loci-deposit-20261023/compare_panel_idx1.py` | 3,966 | `a26dfaf77d097e24e104379f9450ae45f70fc5f49990b0add4ece537b59c069f` |
| `benchmark-loci-deposit-20261023/idx1_panel_compare_out.txt` | 758 | `bad0aa86cee3330b202d5c2f6703613e02b939e30f0b09fa6c8cca5c8d85758f` |
| `benchmark-loci-deposit-20261023/idx1_panel_either_crz_only.tsv` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `benchmark-loci-deposit-20261023/idx1_panel_either_swo_only.tsv` | 56 | `8608661b2946f34f7f00243e16ef64095b4361b448ab22e7bebedff8139ea5b5` |
| `benchmark-loci-deposit-20261023/idx1_panel_per_guide.tsv` | 3,839 | `3fa5f578555015295469b827e01e3c74f53c6f7c959beaab545160c01f1cad37` |
| `benchmark-loci-deposit-20261023/idx1_panel_rna_only_crz_only.tsv` | 264,694 | `4a8ec7025f29dd93f0a2a937bb560a4f4d18bffcb5f0f6992ed78557335d50d3` |
| `benchmark-loci-deposit-20261023/idx1_panel_rna_only_swo_only.tsv` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `benchmark-loci-deposit-20261023/inputs.sha256` | 5,412 | `8f680e4712b358f2dd2db37d26fdd1a6c8d47bdf9e65812dd6c6eebee649fdf0` |
| `benchmark-derivations-20261023-v2/README.md` | 2,083 | `d260cf085158fd163e31bdc17a3d9e19949512215e14625a82b3c1a7a9f274b7` |
| `benchmark-derivations-20261023-v2/SHA256SUMS` | 777 | `773766efc73790562cf906d1b44e1e8b2cfd0836755b9dae24c9611bdd2cf9bc` |
| `benchmark-derivations-20261023-v2/derive_figures.py` | 5,208 | `f4743b87a01345066c7347b88c673d782d421dea33d55fbd0fe13c161a7921a2` |
| `benchmark-derivations-20261023-v2/derived_figures.tsv` | 1,585 | `86250737f3b3280cdc31ffb2218f049560661d48f299447d99b20c5a26a09d27` |
| `benchmark-derivations-20261023-v2/expected_crispritz.sha256` | 4,576 | `19f629548b2829e6830fe1f45abeb079918fcb3f2501744aebedd2b679a6750a` |
| `benchmark-derivations-20261023-v2/figure_lookup.py` | 4,620 | `b43695629e504743ef5d2375ba6c02260cc549fedde996248c2264d1b57d6ed5` |
| `benchmark-derivations-20261023-v2/paper_fwd_ATTGAGATAGTGTGGGGAAGNGG.csv` | 507,984 | `193c47d11f2c53037777c6084de3d3b0bfe0ada0de3410d2771b765da2dcfc35` |
| `benchmark-derivations-20261023-v2/paper_inputs.sha256` | 185 | `3a0b3bdc2f9fbb82728acda46781730008858c3dbf43c5c220525e5e70baeafd` |
| `benchmark-derivations-20261023-v2/paper_keys.tsv` | 317,034 | `a2a17eb62fb937e46d33830fd7625ae5ea4b811662b924d891916d12895d5180` |
| `benchmark-derivations-20261023-v2/test_boundaries.py` | 1,914 | `7ddbbcb3c1e2519a5eba7b83513080a5d67c2f1485c3c5b445d3c6adf1f054e8` |
