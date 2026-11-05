# Lead verification of A5's Cas-OFFinder 2 result (31 Oct 2026, ~16:35 PDT)

Program: the real Cas-OFFinder 2 command-line program, "Cas-OFFinder 2 - v1.0.0", built by A5 from the
unmodified tag v1.0.0 (commit d7df0f45, git diff clean; only RUSTFLAGS -A ambiguous_glob_imports and
vendored crates), binary sha256 prefix 985c9f4d874c64ee. Invoked in CPU mode (C) by the lead, not via
A5's library harness. Input: A5's planted RNA fixture (fixtures/rna.fa: all 1,629 pooled distinct 19-base
deletion words of the 110 CHANGE-seq panel guides, each followed by AGG), pattern NNNNNNNNNNNNNNNNNNNNNGG,
0 mismatches, 0 DNA bulges, 1 RNA bulge; one crRNA line per guide (input_rna_local.txt).
Run: exit 0, 0.22 s, 2,090 output rows, all of bulge type RNA.
Scoring: the lead's own scorer (score2.py), which never imports A5 code: for each guide, the distinct
single-deletion words of its spacer versus the gap-removed DNA field (PAM stripped) of that guide's rows.
Result: 1,561 of 1,634 guide-specific distinct deletion words recovered; 73 missing; every one of the 73
guides with a miss is missing exactly its own terminal-deletion word and nothing else; 73 is exactly the
number of panel guides whose last two bases differ (the guides whose terminal-deletion word is unique).
This matches A5's harness result (1,561/1,634; 73/110) and the source reading (traceback.rs L228-240).
Scope: planted fixture at these budgets; not a genome-scale run; not GPU; says nothing about measured
activity of any site.
