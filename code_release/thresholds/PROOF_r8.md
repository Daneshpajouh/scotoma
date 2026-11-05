# PROVED: the eight-run threshold, with its longer exceptional window

Date: 2026-10-29. Alphabet size four; N=n-1.

**PROVED, by exact computer-assisted certificates and the seven-run theorem.**

\[
\boxed{M(n-1,8)=\left\lfloor\frac{3n-15}{5}\right\rfloor+
\mathbf1_{\{n\in E_8\}},\qquad n\ge10,}
\]

\[
E_8=\{11,16,21,26\}\ \cup\ \{13+5j:0\le j\le13\}.
\]

Thus the last exceptional length is 78, and the first failure of the
previously plausible floor((3n-14)/5) boundary in the remaining residue is n=83.

## 1. PROVED: a 36-orbit lower construction

For n>=15 put C=n-8 and b=floor((3C+9)/5)=floor((3n-15)/5). Each listed prefix has
a uniform tail of weight b-d. The displayed numerator is divided by 379
and is the total mass **per prefix**, with X expanded uniformly to {2,3}.

| Prefix pair | d | Mass numerator |
|---|---:|---:|
| 0000000, 1111111 | 0 | 68 |
| 000011X, X001111 | 1 | 22 |
| 0001110, 1000111 | 1 | 128 |
| 000X000, 111X111 | 1 | 311 |
| 0011110, 1000011 | 1 | 136 |
| 001111X, X000011 | 1 | 14 |
| 00111X1, 0X00011 | 1 | 34 |
| 00111XX, XX00011 | 2 | 26 |
| 00X0001, 0111X11 | 1 | 85 |
| 0111X0X, X1X0001 | 3 | 150 |
| 0111X10, 10X0001 | 2 | 18 |
| 0111X1X, X0X0001 | 2 | 55 |
| 0111XX1, 0XX0001 | 2 | 71 |
| 0X0001X, X0111X1 | 2 | 112 |
| 0X0011X, X0011X1 | 2 | 142 |
| 1X0001X, X0111X0 | 3 | 20 |
| X0011XX, XX0011X | 3 | 87 |
| X0111XX, XX0001X | 3 | 37 |

**PROVED by exact prefix distances.** These are private sphere orbits,
with mass one in each of the eight regions. Their tail weights are
admissible for C>=7. The gain formula, polynomial degree bound, cumulative
profile domination, Newton expansion, and moment argument are those
proved in `PROOF_r6.md`, with prefix length seven and beta=9.

**VERIFIED-COMPUTATIONALLY.** All 2,187 centre prefixes give 1,107 distinct
profiles and 257 maximal profiles. At p=0,...,65 and all five residues,
84,810 polynomials have 2,841,135 nonnegative Newton coefficients. No
shifted starting points or skipped finite cases are needed. For p>=66,

\[
\ell(c,p)\le\frac{4116374}{909979}
 \left(\frac{171}{175}\right)^{66}<\frac{99}{100}<1.
\]

**PROVED.** This gives M(n-1,8)>=b for all n>=15. The source is
`family_r8.json`; `prove_family.py` regenerates the infinite certificate.

## 2. PROVED: the final missing upper residue

Restriction gives M(n-1,8)<=M(n-1,7). Beyond the finite exceptions in the
seven-run formula, this upper bound is floor((3n-14)/5). It equals b except
when n congruent_to 3 mod 5.

**PROVED by a parametric private-dual certificate.** At every integer z>=0,

\[
n=83+5z,\qquad r=8,\qquad m=47+3z,
\]

full weight is impossible. The certificate is `upper_r8_n83.json`.
It uses 31 groups of centres of tail weight at most two, polynomial
numerators, and a degree-twelve common denominator. For z>=51, exact
polynomial inequalities establish nonnegative multipliers, coverage of
every private orbit, and total multiplier strictly below eight. The 51
initial cases z=0,...,50 have separate exact private-dual certificates; the
last is n=333,m=197.

The private-dual implication is proved in `PROOF_r7.md`: every packing of
mass eight must be entirely private. This is not a claimed cover of all U,
and the multiplier sum is not claimed to equal tau*.

**VERIFIED-COMPUTATIONALLY.** The candidate basis was reconstructed by exact
integer evaluations and interpolation with a degree bound of fifteen.
The full polynomial residual was then checked identically, so acceptance
does not rely on interpolation alone. The independent standard-library
checker `verify_additional.py` ignores those derivation expressions: it
regenerates the degree-at-most-two incidence numerators from actual tail
counts, checks every polynomial inequality after z=u+51, and checks all
51 finite bridges.

## 3. PROVED: finish the finite exceptions

For n>=46 outside n congruent_to 3 mod 5, the seven-run upper bound equals b. The new
upper family handles every n congruent_to 3 mod 5 with n>=83. Exact boundary pairs in
the measured table cover n=10,...,50 and the remaining lengths
53,58,63,68,73,78 in that residue. They give an excess of one precisely
on E_8. The few seven-run exceptions below 46 which do not persist at
eight runs are also covered by those finite pairs.

Together with the uniform lower construction, these arguments establish
the boxed all-length eight-run formula.

## 4. VERIFIED-COMPUTATIONALLY: execution and control

```sh
/usr/local/bin/python3 -B prove_family.py family_r8.json
/usr/local/bin/python3 -B verify_additional.py
/usr/local/bin/python3 -B extra_controls.py
```

**REFUTED / executed control.** Using this construction at n=83,m=47
instead of b=46 preserves privacy and regional masses, but has maximum
load 9513874/9465525>1. The exact verifier rejects it.
