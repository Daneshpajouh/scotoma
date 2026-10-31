# Full fractional weight at guide length twenty

PROVED, with exact finite verification: the canonical library at word length
19, run count 20 and radius 9 has fractional optimum 20. The JSON measure has
one unit in each of its twenty private regions. The library is a matching cover.

Run from the release root:

```
python3 -I certificates/m9_full_weight/verify.py
```

Python 3.10 or later and a C++17 compiler are required. The Python code uses
only the standard library. The replay compiles the supplied integer transform
and checks all 1,162,261,467 ternary centre classes. Independent swaps of unused
letters account for all four-letter centres. Library centres have load one;
the maximum among other centres is 39995753666499/40000000000000.
CRT uniqueness and intermediate integer bounds are checked before use.
The two transform passes require approximately 5 GB of memory and a similar
temporary disk allocation. Four worker threads are used by default.

The copied producer and transform are retained without mathematical changes.
The wrapper verifies the SHA-256 manifest, uses a temporary working directory,
recomputes privacy and every regional mass, and compares a fresh exact transform
with the deposited record. It writes no files into the release. The optional
`--metadata-only` checks support, mass and the deposited record but does not
replay the transform; its output makes that distinction explicit.
