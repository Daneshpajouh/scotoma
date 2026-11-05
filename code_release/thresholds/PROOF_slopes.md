# PROVED: the asymptotic slope at every measured run count r=6,...,10

Date: 2026-10-29. The supplied five-run theorem has asymptotic slope 3/5.

**PROVED.** For each fixed r in {6,7,8,9,10},

\[
\lim_{n\to\infty}\frac{M(n-1,r)}{n}=\frac35.
\]

**PROVED: the upper bound.** Restriction E4 and E21 give

\[
M(n-1,r)\le M(n-1,5)=\left\lfloor\frac{3n-11}{5}\right\rfloor.
\]

**PROVED: matching lower slopes.** The exact constructions in this package give

| r | Proved lower bound | Length range | Certificate input |
|---:|---|---|---|
| 6 | floor((3n-12)/5) | n>=13 | family_r6.json |
| 7 | floor((3n-14)/5) | n>=15 | family_r7_lower14.json |
| 8 | floor((3n-15)/5) | n>=15 | family_r8.json |
| 9 | floor(3(n-9)/5) | n>=27 | family_limit_r9.json |
| 10 | floor(3(n-10)/5) | n>=30 | family_limit_r10.json |

The first three use the same exact polynomial and large-tail-distance proof
as `PROOF_r6.md`. The r=8 construction has 36 private orbits, with masses
having denominator 379. Its certificate has 257 maximal prefix profiles,
84,810 polynomials, and 2,841,135 nonnegative Newton coefficients. The
analytic tail bound starts at p=66 and is

\[
\frac{4116374}{909979}\left(\frac{171}{175}\right)^{66}<\frac{99}{100}.
\]

**PROVED: simple r=9 and r=10 lower-bound constructions.** Give total mass
one to each listed prefix orbit, with a uniform tail of weight b-d, where
b=floor(3C/5) and C=n-r. The displayed d values are verified unique-nearest
prefix distances.

| r | d | Prefixes, each of mass one |
|---:|---:|---|
| 9 | 1 | 0000011X, X0011111, 000X0001, 0111X111 |
| 9 | 2 | 001111XX, XX000011, 0X0011X1, 1111101X, X0100000 |
| 10 | 1 | 000000110, 100111111, 000011101, 010001111 |
| 10 | 2 | 00111101X, X01000011, 0111101X1, 0X0100001, 1000X0000, 1111X1110 |

**VERIFIED-COMPUTATIONALLY.** The nine-run construction has 73 maximal
profiles, 15,695 polynomials, and 345,290 nonnegative Newton coefficients.
Its tail bound, from p=43 onward, is

\[
\frac{1031}{243}\left(\frac{29}{30}\right)^{43}<\frac{99}{100}.
\]

The ten-run construction has 69 maximal profiles, 15,525 polynomials,
and 357,075 nonnegative coefficients. Its tail bound from p=45 is

\[
\frac{1079}{243}\left(\frac{29}{30}\right)^{45}<\frac{99}{100}.
\]

**PROVED.** All these constructions have a constant distance and mass
description, exact mass one per private region, and exact load at most one
at every centre. Their finite polynomial certificates cover arbitrarily
long tails, rather than sampling a maximum length. Dividing the upper and
lower bounds by n and applying the squeeze theorem proves the stated limits.

**CONJECTURE / unresolved.** These slope results do not determine the exact
finite thresholds for r=9,10 or for arbitrary r>=11. They do not prove a
two-variable closed form. The r=9 and r=10 lower bounds are deliberately
conservative and are not asserted to be exact thresholds.

```sh
/usr/local/bin/python3 -B prove_family.py family_r8.json
/usr/local/bin/python3 -B prove_family.py family_limit_r9.json
/usr/local/bin/python3 -B prove_family.py family_limit_r10.json
```
