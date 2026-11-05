# Independent endpoint audit

Execution: COMPLETED.

PAM policy: `paper_ngg_core_only`.

| Quantity | Paper | Recomputed or [lower, upper] | Verdict |
|---|---:|---:|---|
| terminal_rna_rows | 583 | 583 | SURVIVED |
| exposed_terminal_rna_rows | 260 | 260 | SURVIVED |
| exclusive_terminal_losses | 129 | 129 | SURVIVED |
| exclusive_terminal_guides | 40 | 40 | SURVIVED |
| first_rna_losses | 0 | 0 | SURVIVED |
| dna_losses | 0 | 0 | SURVIVED |

Data SHA-256: `68f3e1d7eb9e9e18e3ddc859752be8d156378a3e681b0fff329fe679366adc11`.
Bytes: 878760836.
Runtime seconds: 31.262851.
Peak RSS bytes: 28905472.

Active means reads > 0, parsed with Decimal; labels and training thresholds are unused.
Rows retain their multiplicity and own-guide attribution. Exact zero-substitution,
zero-bulge rows are excluded from off-target counts. Admission uses the deposited
substitution count and checks it against the alignment, including the two fixed
G bases of NGG. A gap does not consume a substitution. This annotation census
retains all deposited PAMs, as stated in the paper appendix. Strict-NGG annotation
subcounts and out-of-scope counts are also reported in summary.json.

Search policy: {"annotation_admission": "deposited_mismatches_at_most_budget_including_fixed_PAM_Gs", "loss_scope": "exact_NGG_only_non_NGG_explicitly_out_of_scope", "same_register_distance": "spacer_substitutions_only", "window_registers": "exact_NGG_both_strands"}.

ALTERNATE_POLICY denotes a deliberately different scope, not a paper refutation.

Exclusivity requires failure of every emitted same-PAM pattern and every bulge-free
exact-NGG register on both strands of the 73-base window. The all-window bulged
sensitivity result is separate in summary.json. See README.md for definitions.

Controls: {"budget_changes_loss_states": true, "budget_changes_terminal_admission": true, "complete_restores_all": true, "first_only_restores_terminal": true, "interior_annotations_reached": true, "last_only_cannot_increase_loss": true}.

Anchor policy: schema3-all-occurrences-no-window-genomic-origin.
All-row anchor census: {"invalid_window": 81, "multiple": 10353, "unique": 4925844}.
Admitted active bulge anchor census: {"multiple": 15, "unique": 23677}.
Multiple occurrences remain placement-unresolved: no window genomic origin is
supplied. Each occurrence is evaluated; incomplete context permits either outcome.
Agreement across occurrences identifies a classification, not a genomic placement.
Unresolved-anchor contributions and count bounds are in summary.json;
unresolved_anchors.jsonl retains every affected admitted row and all placements.

Endpoint candidates, losses and changed-control rows are in remote JSONL files.
Exclusive identity hashes and multiplicity are in exclusive_identities.json.
The parent dual-policy result compares these with the data-only reference manifest.


## Both native policies

| PAM policy | Final RNA rows | Exposed rows | Exclusive losses | Guides |
|---|---:|---:|---:|---:|
| pam_inclusive_deposit | 583 | 260 | 148 | 43 |
| paper_ngg_core_only | 583 | 260 | 129 | 40 |

A null value requires consulting the lower/upper bounds in summary.json.
Reference data contain unlabelled candidates, not a supplied list of exclusive sites.
identity_comparison.json lists every set/multiplicity difference both against the
complete candidate list and against its independently classified loss subset.
Raw counts were computed before this comparison; reference identities did not filter them.

Combined runtime seconds: 62.059258.
