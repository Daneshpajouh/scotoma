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
| `scope.py` | The parameter range: pattern lengths to 30, run counts to 31, and the first budget at which the construction fails. |
| `sharp_threshold.py` | The shortest length admitting a separated witness family at budget four, searched over every four-coordinate support and letter assignment. |
| `floor_counterexample.py` | That the statement is false outside its range: an exhaustive search showing four centres covering the union of five deletion balls. |
| `ball_overlaps.py` | Exact intersection sizes of two Hamming balls at every distance, and the volume lower bound they give, for comparison with the minimality result. |
| `cylinder_decomposition.py` | That a Hamming ball is exactly a union of cylinders, how many are needed, and how close that is to the information-theoretic minimum. |

## Reading the output

Every script prints the quantity it computed beside the value it was checking
against, and a control where one exists. `floor_counterexample.py` prints the
covering set it found, so the refutation can be checked by hand. No script
depends on any dataset; all of them operate on the alphabet directly.
