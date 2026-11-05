# Scope fixtures for the fractional-certificate paper

From the release root, with Python 3.10 or later:

```
python3 -I certificates/paper_B_scope/verify_scope.py
python3 -I certificates/paper_B_scope/pooled_counterexample.py
```

The first command checks the two separated target words at (N,r,m)=(19,9,9),
all 192 three-letter column/assigned-binary-letter cases used in the counting
obstruction, and all 2,310 insertion coordinates of the released 110-guide
panel. Its input is the panel reference under data/producers_2026-10-20/.

The second command is the existing exact synthetic-panel verifier, copied
with its guide list. It checks all 110 guide lengths, run counts nine through
twenty, disjoint deletion libraries and all 1,629 distinct patterns. It then
checks containment of the removed pattern's entire radius-four ball in the
union of sixteen other library balls, all 341,716 words. Dropping one neighbour
is a mandatory failing control. Thus the remaining 1,628 balls cover the
same union. This synthetic panel is distinct from the measured panel whose
fractional and integral optima are both 1,629.

All mathematical decisions use integer arithmetic and the standard library.
These commands print their outcomes and write no files. SHA256SUMS covers
the local package; the measured panel remains at its existing release path.
