# Fresh verification and measured resources

All 30 self-contained replay targets passed. The separate external CFD-table target also passed using a freshly fetched, hash-checked input excluded from the deposit.
The three large assay tables and the ancestry analyses were not rerun. They are explicitly excluded from the fresh-pass count.

Validation used Python 3.13.9, NumPy 2.5.1, SciPy 1.17.1, matplotlib 3.10.7, pytest 9.0.2, Node 22.15.1 and Apple clang 21 with C++17. See ENVIRONMENT.json and requirements.txt.

| Replay target | Verifier result | Wall seconds | Peak process RSS bytes | Sampled scratch allocation bytes |
|---|---|---:|---:|---:|
| ceiling | PASS | 2.068 | 154533888 | 229486592 |
| cellular | PASS | 1.041 | 19513344 | 143478784 |
| cover_families | PASS | 1.576 | 28409856 | 229330944 |
| dials | PASS | 21.051 | 156336128 | 229613568 |
| empirical_gates | PASS | 14.895 | 39796736 | 229330944 |
| enumerator | PASS | 4.732 | 155795456 | 233443328 |
| external_cfd | PASS | 4.116 | 158924800 | 221372416 |
| five_run | PASS | 4.191 | 35635200 | 229330944 |
| floor_rows | PASS | 19.793 | 155828224 | 229330944 |
| gap_census | PASS | 35.936 | 86638592 | 229330944 |
| infinite_families | PASS | 338.199 | 151535616 | 321626112 |
| insertion_nofloor | PASS | 168.143 | 1486143488 | 229330944 |
| m10_closure | PASS | 97.930 | 2180022272 | 229801984 |
| m10_r11 | PASS | 11.777 | 803848192 | 229486592 |
| m10_r11_kernel | PASS | 6.116 | 72269824 | 229486592 |
| m8_closure | PASS | 91.177 | 153452544 | 229486592 |
| m9_closure | PASS | 12.301 | 895385600 | 229650432 |
| m9_full_weight | PASS | 70.854 | 4684251136 | 4896038912 |
| maize_loci | PASS | 2.098 | 218218496 | 229486592 |
| monotonicity | PASS | 5.278 | 72679424 | 229330944 |
| multibulge | PASS | 1.024 | 18104320 | 139276288 |
| pooled | PASS | 8.942 | 26771456 | 229330944 |
| private_sphere | PASS | 18.549 | 308527104 | 229330944 |
| r20_m8 | PASS | 29.456 | 4873994240 | 229486592 |
| scope | PASS | 9.974 | 74727424 | 229330944 |
| selector | PASS | 3.155 | 23969792 | 229330944 |
| structure | PASS | 15.061 | 155222016 | 229527552 |
| thresholds | PASS | 1425.968 | 845217792 | 521986048 |
| tool_evidence | PASS | 4.156 | 157630464 | 227827712 |
| upper_rows | PASS | 113.710 | 8029552640 | 1028313088 |
| validation_budget | PASS | 42.754 | 154992640 | 229527552 |

RSS is the maximum reported by wait4 or RUSAGE_CHILDREN for waited command subtrees. It is not the sum of simultaneously resident workers. Disk was sampled every 0.5 seconds and includes the disposable deposit copy. Wall time includes integrity checks, copying, compilation where applicable, replay and cleanup. Runs were observed on this workstation, with some independent work concurrent; this is not an isolated performance benchmark.

The original macOS time wrapper failed after successful commands because its clock-rate system query was blocked. The raw timing failures and successful command exit records are retained in the preparation audit. The table uses actual verifier exits and the separately recorded rusage, not the failing timing-wrapper exit.

For m9_full_weight the measured RSS is 4684251136 bytes and sampled allocated scratch is 4896038912 bytes. The native residue array alone has 4 * 3^19 = 4649045868 bytes. The README estimate of roughly 5 GB for memory and temporary disk is supported by this run; reserve at least 6 GB for each. The upper-row route peaked at 8029552640 bytes of process RSS; allow at least 10 GB for the complete serial suite. No local measurement substitutes for the missing N=22 cluster sizing probe.

Every row is run through python3 -B reproduce.py TARGET, except external_cfd, which is python3 -B reproduce_external.py cfd-table /path/to/mismatch_score.pkl. Validation used /usr/local/bin/python3. reproduce.py contains the exact subordinate commands. The threshold target uses verify_all.py --full-grid and checks all 412 direct boundary pairs.

Display transcripts are in replay_logs/. They preserve numerical results while normalizing local task paths and prose typography. Their original unsanitized transcripts and command/resource JSONs are retained outside the proposed deposit in verification/.

The eight mandatory certificate packages all passed. The broader result ledger has 73 explicit claim groups: 50 are reproducible from the staged tree with the stated software prerequisites, and 23 require external or missing inputs. Counting is one result group per TSV row, not one numeral or sentence. Infinite conclusions use their analytic proofs in addition to finite exact checks. Conditional reanalysis of deposited records is distinguished from replaying the complete source experiment. The final ledger was reconciled with the current paper references after concurrent workspace changes were detected.
