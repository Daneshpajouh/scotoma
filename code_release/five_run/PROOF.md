# PROVED: the five-pattern full-weight conjecture, for every n >= 10

Date: 2026-10-29.

**PROVED (exact computer-assisted proof).** For every integer n >= 10, with N=n-1,

\[
\boxed{\tau^*(N,5,m)=5\quad\text{for every integer }
0\le m\le\left\lfloor\frac{3n-11}{5}\right\rfloor.}
\]

**PROVED, combining this result with the supplied E13 and E16.**

\[
\boxed{M(n-1,5)=\left\lfloor\frac{3n-11}{5}\right\rfloor
\qquad(n\ge10).}
\]

The proof below constructs a packing of mass five at the boundary radius. Its
load proof consists of a finite, exact polynomial certificate for small tail
distances and an analytic bound for every larger distance. It does not infer a
packing from the cover values in the question. The established upper direction
and budget monotonicity are used as premises and are not recomputed here.

## 1. PROVED: an explicit packing family

Put

\[
C=n-5\ge5,\qquad b=\left\lfloor\frac{3C+4}{5}\right\rfloor.
\]

Write X for an independently chosen letter in {2,3}. An orbit (t,w) consists of
all words with prefix t, with X expanded uniformly, and a length-C tail having
exactly w nonzero entries, also uniformly. Its size is

\[
2^{\#X(t)}\binom Cw3^w.
\]

Assign the following **total mass to each listed orbit**, not to each word:

| Prefixes (each prefix is a separate orbit) | Tail weight | Mass of each orbit |
|---|---:|---:|
| 1111, 0111, 0011, 0001, 0000 | b | 1/5 |
| 001X, X011 | b-1 | 2/5 |
| 011X, X001, 0X00, 11X1 | b-1 | 4/5 |

Thus there are eleven orbits. Their weights do not need to be fitted as
functions of n. For clarity, if A=binom(C,b)3^b and
B=binom(C,b-1)3^(b-1), the individual-word weights are respectively
1/(5A), 1/(5B), and 2/(5B) in the three rows.

Every library prefix is uniquely nearest to itself. The six other prefix
orbits have the following unique nearest library prefixes, each at distance one:

| Prefix | Unique nearest library prefix |
|---|---|
| 001X, X011 | 0011 |
| 011X | 0111 |
| X001 | 0001 |
| 0X00 | 0000 |
| 11X1 | 1111 |

Consequently every supported word has distance exactly b to its assigned
library word and distance at least b+1 to the other four. All support is
private and lies in U. Each private region receives mass one: the middle region
receives 1/5+2/5+2/5, and each other region receives 1/5+4/5. The total mass is
exactly five. Tail weights are admissible because C>=5 and 3<=b<=C.

It remains to prove that every radius-b ball has load at most one.

## 2. PROVED: exact load formula and four sufficient prefix profiles

The packing is invariant under swapping 2 and 3 independently in each prefix
coordinate, permuting tail coordinates, and permuting the three nonzero letters
independently in each tail coordinate. These transformations fix the library.
Therefore every centre is represented by (c,p), where c is in {0,1,X}^4 and p is
the number of nonzero tail coordinates. There are 81(C+1) centre classes,
including all centres outside U.
For centre prefixes, X denotes the representative letter 2; invariance makes
the representative with letter 3 equivalent.

For a uniform weight-w tail, compare its distance to a fixed weight-p centre
tail with its distance w to the zero tail. Let S be the distance reduction.
At each of the p changed centre positions the reduction is -1 if the supported
tail letter is zero, 0 if it is a different nonzero letter, and +1 if it matches.
Define

\[
T_{C,w,p}(k)=\Pr(S\ge k).
\]

With binomial coefficients outside their integer range interpreted as zero,

\[
T_{C,w,p}(k)=
\frac{R(C,w,p,k)}{\binom Cp3^p},
\tag{1}
\]

\[
R(C,w,p,k)=
\sum_{a=0}^{p}\binom wa\binom{C-w}{p-a}3^{p-a}
\sum_{e=\max(0,p+k-a)}^{a}\binom ae2^{a-e}.
\tag{2}
\]

Indeed, a is the number of occupied tail positions among the p centre positions;
conditional on a, the number e of matching letters is binomial(a,1/3), and
S=a+e-p. This establishes (1) to (2) by counting, with no approximation.

For each fixed prefix c, let A_d/5 be the mass of the five library-prefix
orbits at prefix distance d, and let B_d/5 be the corresponding mass of the
six free-letter orbits. Free-letter assignments are averaged here. The exact
load is

\[
\ell(c,p)=\frac15\sum_{d=0}^4
\left[A_dT_{C,b,p}(d)+B_dT_{C,b-1,p}(d-1)\right].
\tag{3}
\]

There are 44 distinct pairs (A,B). Each is dominated, separately in the
cumulative sums of A and B, by at least one of the following four actual
prefix profiles:

| Row | Representative c | A_0,...,A_4 | B_0,...,B_4 |
|---:|---|---|---|
| 0 | 0011 | 1,2,2,0,0 | 0,4,8,8,0 |
| 1 | 001X | 0,1,3,1,0 | 1,3,4,8,4 |
| 2 | 0101 | 0,2,3,0,0 | 0,0,16,4,0 |
| 3 | 011X | 0,1,2,2,0 | 2,3,1,10,4 |

Here domination means sum_{d=0}^j A_d(upper) >= sum_{d=0}^j A_d(lower),
and the same for B, for every j=0,...,4. All A totals are five and all B totals
are twenty. The two tail probabilities in (3) decrease with d, so summation by
parts proves that this domination also dominates the load. It suffices to bound
the four rows in the table.

**VERIFIED-COMPUTATIONALLY.** The verifier constructs the profiles directly from
the eleven orbits by literal free-letter expansion, and checks the cumulative
inequalities for every one of the 81 centre prefixes. The complete domination
map is retained in `proof_certificate.json`.

## 3. PROVED using a finite exact certificate: every tail distance 0 <= p <= 45

This step is a polynomial proof over all lengths, not a sweep up to a maximum
length. Write

\[
C=5q+s,\qquad s\in\{0,1,2,3,4\},\qquad
b=3q+\left\lfloor\frac{3s+4}{5}\right\rfloor.
\]

For a fixed residue s, a fixed p, and a fixed profile (A,B) in the table, define

\[
F_{s,p,A,B}(q)=
5\binom{5q+s}{p}3^p
-\sum_{d=0}^4\left[
A_dR(5q+s,b,p,d)+B_dR(5q+s,b-1,p,d-1)\right].
\tag{4}
\]

**PROVED.** This is a rational polynomial in q of degree at most p. In each
summand of (2), binom(w,a) has degree a and binom(C-w,p-a) has degree p-a,
and all other factors are constant with respect to q. This remains true when
an integer binomial coefficient vanishes because its lower argument exceeds
its nonnegative upper argument: its falling-factorial polynomial vanishes there
as well.

The domain C>=5 and C>=p is precisely q>=q_0, where

\[
q_0=\max\left(1,\left\lceil\frac{p-s}{5}\right\rceil\right).
\]

For any polynomial F of degree at most p, Newton's identity gives

\[
F(q_0+z)=\sum_{k=0}^{p}\Delta^k F(q_0)\binom zk.
\tag{5}
\]

For integers z>=0 all binom(z,k) are nonnegative. Hence the finite set of exact
conditions

\[
\Delta^kF(q_0)\ge0\qquad(k=0,\ldots,p)
\tag{6}
\]

proves F(q)>=0 for every admissible integer q, with no upper bound on q.

**VERIFIED-COMPUTATIONALLY.** All 920 polynomials in (4), corresponding to
five residues  *  46 distances  *  four profiles, satisfy (6). All 21,620 Newton
coefficients are nonnegative integers. There are no exceptional lengths or
additional shifts of q_0. The verifier regenerates the coefficients from (2)
and (4); it does not trust stored coefficients. It also checks the Newton
reconstruction at an extra evaluation point. Every coefficient is provided in
`proof_certificate.json`.

The complete finite check is the following algorithm, using integer arithmetic:

```text
for s = 0,...,4:
  for p = 0,...,45:
    q0 = max(1, ceil((p-s)/5))
    for each of the four profiles:
      v = [F(q0), F(q0+1), ..., F(q0+p)]
      while v is nonempty:
        assert v[0] >= 0
        v = [v[i+1]-v[i] for i = 0,...,len(v)-2]
```

**PROVED.** Combining the verified finite conditions (6) with the degree proof
and (5) establishes every four-row load inequality for every C>=5 and every
0<=p<=min(C,45).

As a useful exact boundary check, the row-0 load at p=1 simplifies to

\[
\ell_0(1)=\frac{5b-4}{3C}\le1.
\tag{7}
\]

The construction therefore detects the same integer boundary as the established
upper direction, rather than relying on a rounded numerical margin.

## 4. PROVED: every tail distance p >= 46

Take t=7/5. Conditional on an occupied tail position, the expectation of
t raised to its gain is (2+t)/3; for an unoccupied position it is 1/t.
The p positions sample without replacement from C positive values: w copies
of (2+t)/3 and C-w copies of 1/t. The elementary symmetric mean inequality
therefore gives

\[
\mathbb E[t^S]\le
\left[\frac{1-w/C}{t}+\frac wC\frac{2+t}{3}\right]^p.
\tag{8}
\]

For completeness, the elementary symmetric mean inequality follows by
equalizing any two positive population entries while preserving their sum.
The elementary symmetric sum of order p has the form
E_p(others)+(x+y)E_{p-1}(others)+xy E_{p-2}(others), so equalization increases
it. Its maximum at fixed sum is attained when all entries are equal, yielding
the right side of (8) after division by binom(C,p). Thus this step does not
assume independence of occupied positions.

For a supported orbit with prefix distance delta from the centre and assigned
library-prefix distance d in {0,1}, acceptance requires S>=delta-d. By Markov's
inequality and (8), its accepted fraction is at most

\[
t^{d-\delta}
\left[\frac{1-w/C}{t}+\frac wC\frac{2+t}{3}\right]^p.
\tag{9}
\]

Both supported tail weights w satisfy

\[
\frac wC\le\frac bC\le\frac35+\frac4{5C}.
\]

Since C>=p>=46, the bracket in (9) is at most

\[
h=\frac{1-\rho}{t}+\rho\frac{2+t}{3}
=\frac{11749}{12075}<1,
\qquad \rho=\frac35+\frac4{5\cdot46}.
\tag{10}
\]

The prefix multiplier, summed over the eleven weighted orbits and their
free-letter assignments, is

\[
V(c)=\frac15\sum_{\delta=0}^4
\left[A_\delta t^{-\delta}+B_\delta t^{1-\delta}\right]
\le\frac{169}{49}.
\tag{11}
\]

**VERIFIED-COMPUTATIONALLY, with exact fractions.** Inequality (11) holds for
all 81 prefixes; equality occurs for row 0. It can also be checked just on the
four dominating rows. The exact scalar inequality checked by the verifier is

\[
\frac{169}{49}\left(\frac{11749}{12075}\right)^{46}
<\frac{99}{100}<1.
\tag{12}
\]

**PROVED.** Equations (9) to (12) give the uniform bound

\[
\ell(c,p)\le V(c)h^p
\le\frac{169}{49}h^{46}<1
\qquad(p\ge46).
\]

This covers arbitrarily long tails and every possible number of changed centre
positions. Together with Section 3, every centre has load at most one.

## 5. PROVED: conclude the full stated range

At radius b the explicit packing is feasible and has mass five, so
tau*(N,5,b)>=5. The five library balls supply a cover of mass five, proving
tau*(N,5,b)=5. The supplied budget monotonicity E16 then gives tau*(N,5,m)>=5
for every 0<=m<=b, while the library cover still gives tau*<=5.

The conjectured full-weight direction is therefore proved for **every integer
n>=10, r=5, and every integer 0<=m<=floor((3n-11)/5)**. Together with the supplied
strict upper bound beyond this budget, this gives the exact formula for M
stated at the start. This proof makes no claim for r>=6 or n<10.

## 6. VERIFIED-COMPUTATIONALLY: executable checks and failing controls

Run from this directory:

```sh
/usr/local/bin/python3 -B verify_k3.py --write-results
/usr/local/bin/python3 -B independent_literal_check.py
```

`verify_k3.py` is standalone and uses only the Python standard library. All
loads, masses and tight/short decisions use `fractions.Fraction`; polynomial
coefficients and combinatorial counts use exact integers. It regenerates the
infinite-range polynomial certificate and exponential bound on every run.
`--n` can replace the finite length list.

**VERIFIED-COMPUTATIONALLY.** The new eleven-orbit packing was checked at
n=10,...,22, 27, 32, 52, 102, 103, 137, 202. Every one of the 63,423 centre
classes across these twenty cells has load at most one. Each cell has mass
five, mass one in every private region, and maximum load exactly one. For these
finite checks the tail incidence is recomputed by selecting actual tail words
and their matches, with denominator binom(C,w)3^w. This differs from the
sampling-of-centre-positions formula (1), whose denominator is binom(C,p)3^p.

**VERIFIED-COMPUTATIONALLY.** The independent literal fixture enumerates all
4^6=4,096 centres and all supported words at N=6,m=2. It compares actual Hamming
ball sums against the eleven-orbit load calculation, including all free-letter
assignments. This is a check of the incidence calculation; the theorem's stated
range remains n>=10.

**REFUTED / rejected as required  -  control 1.** At n=17, increase the radius
from b=8 to 9 and regenerate the same eleven-orbit construction. Its mass remains
five and its support remains private, but the centre with prefix 0011 and one
nonzero tail coordinate has exact load 41/36>1. The complete checker rejects
the packing; its maximum load is 9002/7425 at prefix 0011 and tail weight 3.

**REFUTED / rejected as required  -  control 2.** At n=17,b=8, give the free-letter
orbits tail weight b instead of b-1. The first such orbit, 001X with tail weight
8, has distance 9 from its closest library word and lies outside U. The checker
rejects it at the privacy/support test. Both controls were executed and neither
was accepted.

The machine-readable outcomes are in `verification_results.json`,
`proof_certificate.json`, and `independent_literal_result.json`. The earlier
scripts and logs in this directory retain the unsuccessful fixed-support
candidate and the route to the constant eleven-orbit family. They are
exploratory artifacts, not additional premises of this proof.
