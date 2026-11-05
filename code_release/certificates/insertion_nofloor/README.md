# DNA insertion counts and matching-floor refutations

**PROVED.** For a nonempty length-`n` guide over `{A,C,G,T}`, the complete
one-letter insertion family has `3n+4` distinct concrete words. The `n+1`
one-wildcard insertion patterns have the same radius-`m` target union, with
the wildcard ignored in the mismatch count. Word balls and wildcard patterns
have different costs. Common PAM wildcards do not enter this spacer model.

**REFUTED.** Neither family size is a universal covering floor. For guide `A`
at budget one the seven concrete words need exactly four balls, with centres
`AA, AC, AG, AT`. The single pattern `A*` suffices instead of both original
patterns. The target is all 16 two-letter words. At each coordinate a letter
absent from three proposed centres produces a target missed by all three.
The checker also exhausts all 560 distinct ambient triples literally.

Run from this directory with Python 3.10 or later:

```sh
nice -n 10 /usr/local/bin/python3 -I verify.py verify
```

The verifier uses only the standard library, uses integer arithmetic for all
mathematical decisions, runs in one process, reads paths relative to its own
file, writes no files, and exits zero only when every check passes. Omitting
`verify` runs the same checks. A different Python 3.10+ executable may be
substituted. No solver, network, cluster, or files outside this folder are
needed. Keep the scheduling request on a shared machine. If the operating
system rejects the priority change, inspect the inherited nice value before
running on a shared host.

**VERIFIED-COMPUTATIONALLY.** All 14 exact cases have a written upper cover
and a finite branching-DAG proof that one fewer centre cannot cover. See
`table.json`, `exact_table.csv`, and `certificates/*.json.gz`. In all cases
competitors include every four-letter centre, including those outside the
insertion family. The exact radius-one values for the representatives are:

| Guide | Family size | Target count | Exact word-ball minimum |
|:---|---:|---:|---:|
| A | 7 | 16 | 4 |
| AA | 10 | 37 | 6 |
| AC | 10 | 46 | 7 |
| AAA | 13 | 67 | 9 |
| AAC | 13 | 85 | 10 |
| ACA | 13 | 94 | 12 |
| AAAA | 16 | 106 | 12 |
| AAAC | 16 | 133 | 13 |
| AACA | 16 | 151 | 15 |
| AACC | 16 | 142 | 13 |
| ACAC | 16 | 160 | 16 |
| ACCA | 16 | 151 | 15 |

**VERIFIED-COMPUTATIONALLY.** At radius two, both `AA` and `AC` have target
size 64 and exact minimum four, attained by `AAA, CCC, GGG, TTT`. Run-length
isometries are checked on all 340 guides of lengths one through four. Counts,
positive-slot counts, and full-coordinate saturation are checked on every
guide of lengths zero through seven, 21,845 guides in total. The empty guide
occurs only in these count checks; the smallest-witness statement assumes a
nonempty guide.

**PROVED.** For `A^n`, the family is `B_1(A^(n+1))` and its radius-`m` union
is `B_(m+1)(A^(n+1))` when `1 <= m <= n`. Select `n+1-m` coordinates. The
all-`A` word and the three non-`A` single-coordinate changes at those positions
give a sound cover of size `1+3(n+1-m)`. One wildcard at each selected position
gives a sound cover of size `n+1-m`. Every boundary target has `m+1` non-`A`
positions, so it meets the selected set whose complement has size `m`.

**VERIFIED-COMPUTATIONALLY.** `long_witness.json` selects one-based positions
2 through 18 at `n=20,m=4`. Its 52 word centres and 17 wildcard patterns are
written in `long_concrete_cover.txt` and `long_wildcard_cover.txt` and checked
against the certificate. Their union is exactly `B_5(A^21)`, of size 5,467,456.
The complete family has 64 words and 21 wildcard patterns. Under the
positive-slot convention it has 61 words and 20 patterns, with the same target
union and the same smaller covers. The 52 and 17 are upper bounds, with no
claim of optimality.

**PROVED.** The long verification quotients by the set of non-`A` positions:
independent permutations of `C,G,T` at each coordinate preserve both the
target and these particular written covers. A support of size `w` therefore
represents exactly `3^w` words. The checker exhausts all `2^21` supports,
representing all `4^21` ambient words, and checks literal word membership on
5,008 small homopolymer word/budget instances to calibrate the formulas.

**PROVED / VERIFIED-COMPUTATIONALLY.** An additional witness contains all four
letters: for `g=A^17CGT`, remove `A^18CGT` from the 64-word family. Its 63
remaining words are in `run_extension_cover.txt` and `run_extension_witness.json`.
The checker exhausts the removed radius-four ball's 522,649 members and finds
a retained neighbour covering each; every other original ball stays present.
This verifies the union containment without enumerating the full ambient space.

**PROVED (certificate semantics).** For each finite case the target oracle
deletes each possible target coordinate and compares the result with the
guide. This independently reconstructs insertion-family distances. Every
ambient centre is retained, empty on the target, or supplied with a retained
centre whose target intersection contains it. Replacing a discarded centre
by its checked superset cannot increase a cover's cardinality.

At a DAG branch, select an uncovered target. Every possible covering centre
available for it must appear as a child. In child order, choose that centre
and exclude it and all earlier choices from the remainder. Every possible
cover belongs to a branch according to its first choice for that target.
Leaves check an exhausted centre budget, an uncovered target with no choice,
an integer volume bound, or too many uncovered targets with disjoint sets of
candidate centres. The checker requires exact state equality on DAG reuse,
rejects cycles and unreachable nodes, and checks every leaf inequality.
The upper cover and the infeasibility of `K-1` prove the exact value.

**VERIFIED-COMPUTATIONALLY (rejection controls).** The default replay rejects
all seven primary planted-invalid claims: deleting `AA` from the four-centre
cover; using that cover at radius zero; omitting an ambient centre; omitting
a DAG branch; asserting a false volume bound; changing the family count to
six; and claiming that unit weights on the seven insertion words give a
valid radius-one packing. It also rejects a tampered 63-word construction.
The two saved cover controls can be checked individually; both must exit 1:

```sh
nice -n 10 /usr/local/bin/python3 -I verify.py cover controls/missing_AA.json
nice -n 10 /usr/local/bin/python3 -I verify.py cover controls/wrong_radius.json
```

Each prints `REFUTED: checker rejected claim: uncovered target CA` on stderr.
A single finite certificate can be replayed with
`verify.py certificate certificates/A_m1.json.gz` using the same executable.

For file integrity, run from this directory:

```sh
shasum -a 256 -c SHA256SUMS
```

**CONJECTURE (open questions).** General exact covering numbers and a
classification for arbitrary insertion guides remain open beyond the stated
cases. The wildcard matching floor itself is refuted, including for sound
one-wildcard competitors and the positive-slot convention at length 20 and
budget four. Neither 52 nor 17 is certified optimal. These witnesses do not
determine covering minima for the measured panel, or for mixed or multi-base
edit classes, and do not settle the pooled deletion-pattern problem.
