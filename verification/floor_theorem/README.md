# Minimum covering number for a single-deletion pattern library

Scripts that check, by exhaustive computation, every quantitative claim in the
minimality result and the counterexamples that bound it. Each script needs only
the standard library and prints what it checked. Three of them read data, as
"Reading the output" below says; one of those needs a third-party table that is
not released here.

Run them with any Python 3.9 or later:

    python3 verification/floor_theorem/verify_construction.py

## The statement being checked

For an alphabet of size at least four, a string of length at least twenty, and a
mismatch budget of at most five, the distinct single-base deletion patterns of
that string cannot be covered at that budget by fewer centres than there are
patterns, where the competing centres may be any words at all. The minimum is
attained by the patterns themselves.

## What each script establishes

| Script | What it checks |
|---|---|
| `canonical_relabelling_lemma.py` | Every coordinate of a deletion library carries at most two letters and switches between them at most once, so the canonical relabelling is a distance-preserving bijection. Exhaustive for lengths 2 to 8, plus 20,000 random 20-mers and adversarial shapes. |
| `verify_construction.py` | The witness family at budget five: distance to its own centre, distance to every other centre, and all pairwise distances, against the eighteen predicted values. |
| `verify_m4.py` | The same three conditions at every budget from one to five, obtained by restoring coordinates from the budget-five family. |
| `scope.py` | How far the budget-four family reaches in length and run count, and where one extension of it stops. The m sweep extends the support with coordinates 17 and 18; that extension holds to m = 5 and fails at m = 6 and 7. That is one support failing, not a statement that no support works at those budgets. The expected YES and NO pattern is written into the script, which exits 1 if the run departs from it. |
| `pooled_counterexample.py` | The panel behind the body's 1,628, checked end to end from its released guide list `data/synthetic_panel_110.txt` (the original list, over the relabelled alphabet). The 1,628 belongs to this SYNTHETIC panel, not to the real CHANGE-seq panel of the next row, which has the same pattern count and needs all 1,629 balls by the pooled theorem. Checked here: 110 guides of length 20, pairwise disjoint libraries, 1,629 patterns, run counts 9 to 20, the redundant pattern and its 16 witnesses present in the library, and its radius-4 ball inside their union over all 341,716 words. The control, dropping one witness, runs every time. |
| `pooled_volume_bound.py` | The volume lower bound on covering the REAL pooled panel, computed from its 1,629 patterns in `data/pooled_patterns_1629.txt`, not from the synthetic panel of the row above: exact union 381,985,082 words, ball 341,716, so at least 1,118 centres, 68.6%. The body prints 839 (51.5%), which rests on a weaker union bound of 286,613,932; both are valid lower bounds and 1,118 is the stronger. Takes about 5 minutes and 11 GB, stated in the script. The 1,629 patterns in `data/pooled_patterns_1629.txt` are the distinct single-base deletions of the 110 guides of the CHANGE-seq measurement table of Lazzarotto et al. (2020): 1,634 counted guide by guide, 5 of them shared by guides offset by one base. The synthetic panel of the row above is written over 0123, matches the real one in guide count and pattern count, and has pairwise disjoint libraries, which the real panel does not. `pooled_patterns_producer.py` rebuilds the file from the CHANGE-seq table. |
| `pooled_patterns_producer.py` | Where the 1,629 patterns come from. Reads the CHANGE-seq table (not redistributed here; its source and SHA-256 are in the script), rebuilds `data/pooled_patterns_1629.txt` from its 110 guides and checks the result byte for byte: 110 guides over ACGT, 1,634 patterns guide by guide, 1,629 distinct, exactly 5 shared, each by two guides offset by one base and listed, SHA-256 b749f0869f5d4a97d1e6c6c48f750694dbec54460846ad4730b5cba76136fa02. The control, dropping one guide, runs every time and must change the digest. Prints the input table's SHA-256 so a run log shows which copy was read. |
| `pooled_redundancy.py` | That a pattern of a pooled library can be redundant: B_4(p) lies inside the union of sixteen other patterns' radius-4 balls, so such a library needs one fewer ball than it has patterns. Checks the MECHANISM behind the body's 1,628 figure, not the 110-guide construction, which `pooled_counterexample.py` checks on the released list `data/synthetic_panel_110.txt`. Its control, dropping one witness, is part of the run and must leave the ball uncovered. |
| `sharp_threshold.py` | The shortest length at which a cyclic witness template exists at budget four, over every four-coordinate support and letter assignment of the template. It searches cyclic templates only, as its own docstring says, so it does not rule out a non-cyclic family at a shorter length. |
| `floor_counterexample.py` | That the statement is false outside its range: an exhaustive search showing four centres covering the union of five deletion balls. |
| `ball_overlaps.py` | Exact intersection sizes of two Hamming balls at every distance, and the volume lower bound they give, for comparison with the minimality result. |
| `cylinder_decomposition.py` | That a Hamming ball is exactly a union of cylinders, how many are needed, and how close that is to the information-theoretic minimum. |

## Reading the output

Every script prints the quantity it computed beside the value it was checking
against, and a control where one exists. `floor_counterexample.py` prints the
covering set it found, so the refutation can be checked by hand.

Three scripts read data. `pooled_counterexample.py` reads the synthetic guide
list `data/synthetic_panel_110.txt`, and `pooled_volume_bound.py` reads the
pattern list `data/pooled_patterns_1629.txt`; both lists are released here.
`pooled_patterns_producer.py` reads the CHANGE-seq table, a third-party dataset
that is not redistributed, and rebuilds the pattern list from it. Every other
script operates on the alphabet directly and reads no file.
