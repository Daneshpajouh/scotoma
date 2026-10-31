# Exact m=8 closure at guide lengths 21 and 22

**PROVED (computer-assisted).** The three certificates close the short-length
cells needed for `K = r` at every guide length `n >= 20`, every run count and
every mismatch budget `m <= 8`. The appendix theorem `thm:eightall` gives the
case split with padding, step-down, the length selector and budget monotonicity.
The JSON filenames use library word length `N = n - 1`.

| Guide length n | Runs r | Exact mass | Exact maximum radius-8 load |
|---:|---:|---|---|
| 21 | 21 | 1099560000000000/53458540628729 > 20 | 1 |
| 22 | 21 | 1596291840000000/76178068400927 > 20 | 9519495622400899/9522258550115875 < 1 |
| 22 | 22 | 11706140160000000/533246478806489 > 21 | 1 |

Run from the release root with Python 3.9 or later:

```
python3 -I certificates/m8_closure/verify.py verify --workers 4
```

The checker uses only the standard library. It reads certificates relative to
its own file, prints its results, writes no files, and exits zero only when
all checks pass. `--workers` accepts 1 to 8. It directly checks every support
skeleton, computes the exact total mass, and maximizes the load over every
four-letter centre using ternary symmetry classes and an exact dynamic program.
The maximum represents all `3^20` or `3^21` centre classes, including nonbinary
centres. Literal expansion checks the maximizing centre. A small full-space
exhaustion and deliberately overloaded/outside-support controls are included.
This proves the integral covering number, without asserting a fractional optimum.

For integrity, run from this directory:

```
shasum -a 256 -c SHA256SUMS
```

A row `(g,j,t,f_left,f_right,v)` denotes a uniform guarded shell at
`P_j = 1^j 0^(N-j)`. Its size is
`S = binom(L,t) binom(t,f_left) binom(R,8-t) binom(8-t,f_right) 2^(8-f_left-f_right)`,
where `L=max(0,j-g)` and `R=max(0,N-j-g)`. Its per-word weight is
`v H / (S A)`, with `A` the JSON normalizer and `H` the least common multiple
of every admissible shell size over all listed guards and all anchors `0..N`.
The lcm is over all admissible shells, not only shells listed with positive
weight. Overlapping shells add their weights. The checker recomputes this scale.

The certificates were built and checked by one exact implementation and
re-checked by `verify.py`, which was written independently from the
definitions; the two agree on every fraction above. The `verify` subcommand
runs the self-test, mass, support, dynamic-program and negative-control
checks, and needs no files outside this folder.
