# Minimum covering number for a single-deletion pattern library

Scripts that check, by exhaustive computation, every quantitative claim in the
minimality result and the counterexamples that bound it. Each script is
self-contained, needs only the standard library, and prints what it checked.

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
| `pooled_counterexample.py` | The panel behind the body's 1,628, rebuilt and checked end to end from its released guide list: 110 guides of length 20, pairwise disjoint libraries, 1,629 patterns, run counts 9 to 20, the redundant pattern and its 16 witnesses present in the library, and its radius-4 ball inside their union over all 341,716 words. The control, dropping one witness, runs every time. |
| `pooled_volume_bound.py` | The volume lower bound on covering the pooled panel, from the released 1,629 patterns: exact union 381,985,082 words, ball 341,716, so at least 1,118 centres, 68.6%. The body prints 839 (51.5%), which rests on a weaker union bound of 286,613,932; both are valid lower bounds and 1,118 is the stronger. Takes about 5 minutes and 11 GB, stated in the script. |
| `pooled_redundancy.py` | That a pattern of a pooled library can be redundant: B_4(p) lies inside the union of sixteen other patterns' radius-4 balls, so such a library needs one fewer ball than it has patterns. Checks the MECHANISM behind the body's 1,628 figure, not the 110-guide construction, whose synthetic guide list is not in the release. Its control, dropping one witness, is part of the run and must leave the ball uncovered. |
| `sharp_threshold.py` | The shortest length at which a cyclic witness template exists at budget four, over every four-coordinate support and letter assignment of the template. It searches cyclic templates only, as its own docstring says, so it does not rule out a non-cyclic family at a shorter length. |
| `floor_counterexample.py` | That the statement is false outside its range: an exhaustive search showing four centres covering the union of five deletion balls. |
| `ball_overlaps.py` | Exact intersection sizes of two Hamming balls at every distance, and the volume lower bound they give, for comparison with the minimality result. |
| `cylinder_decomposition.py` | That a Hamming ball is exactly a union of cylinders, how many are needed, and how close that is to the information-theoretic minimum. |

## Reading the output

Every script prints the quantity it computed beside the value it was checking
against, and a control where one exists. `floor_counterexample.py` prints the
covering set it found, so the refutation can be checked by hand. No script
depends on any dataset; all of them operate on the alphabet directly.
