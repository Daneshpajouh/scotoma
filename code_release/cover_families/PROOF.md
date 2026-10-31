# The proposed general-length equivalence is refuted

Date: 2026-10-29. All words, balls, and libraries below use the four-letter alphabet in the question. The admissible parameter range is $5\le r\le N+1$, with integer radii. Statements marked **PROVED** have mathematical proofs below, sometimes with an exhaustive finite polynomial certificate. **VERIFIED-COMPUTATIONALLY** denotes exact finite certificate verification. No novelty or external acceptance claim is made.

## 1. Decisive result

**REFUTED.** The assertion

\[
\tau^*(N,r,m)=r\quad\Longleftrightarrow\quad
m\le\left\lfloor\frac{3(N+1)-11}{5}\right\rfloor
\qquad\text{for every }r\ge5
\]

fails within the requested length range. At $N=16,n=17,r=6,m=8$, the asserted threshold is exactly 8, but

\[
\boxed{\tau^*(16,6,8)=\frac{199683}{33307}
       =6-\frac{159}{33307}<6.}
\tag{1}
\]

**VERIFIED-COMPUTATIONALLY.** Equation (1) has a feasible fractional cover and a feasible packing with exactly the same mass. The cover uses 37 symmetry orbits and covers every one of the 1,631 point orbits in $U$, with minimum coverage 1. The packing uses 37 orbits and has maximum load 1 over all 2,916 centre orbits, including centres outside $U$. The complete coefficients are in `counterexample_N16_r6_m8_cover.tsv` and `counterexample_N16_r6_m8_packing.tsv`, and are embedded in `verify_j1.py`.

**VERIFIED-COMPUTATIONALLY.** A mass-6 packing at $(16,6,7)$, together with (1) and the supplied monotonicity E2, gives $M(16,6)=7$. The $r=5$ certificates below give $M(16,5)=8$. Thus the threshold already depends on $r$ at this length.

**PROVED.** The cover extension in Section 4 gives

\[
\tau^*(16,r,8)\le r-\frac{159}{33307}<r,
\qquad 6\le r\le17.
\tag{2}
\]

**PROVED, by the family in Section 6.** There are infinitely many counterexamples. For every integer $z\ge0$, put

\[
N=16+5z,\qquad n=17+5z,\qquad m=8+3z.
\]

Then $m=\lfloor(3n-11)/5\rfloor$, while $\tau^*(N,r,m)<r$ for every $6\le r\le N+1$. The analytic formula covers $z\ge3$; three explicit finite certificates cover $z=0,1,2$.

**PROVED.** This is a refutation of the fractional full-weight assertion. At the decisive cell the integral covering number is nevertheless $K=6$: (1) is greater than 5, so $K\ge\lceil\tau^*\rceil=6$, and the library gives $K\le6$. No conclusion here identifies $K$, $\tau^*$, or the sound-pattern covering number with one another.

## 2. Exact recomputation for every $n=10,\ldots,22$, $r=5$

**VERIFIED-COMPUTATIONALLY.** Write $b(n)=\lfloor(3n-11)/5\rfloor$. Every entry in both value columns below has a packing and a cover of matching exact mass, with the library supplying the mass-5 cover in the middle column. All centre orbits and all point orbits are checked. E2 then determines the full-weight threshold at each listed length.

| $n$ | $N$ | $b(n)$ | $\tau^*(N,5,b(n))$ | $\tau^*(N,5,b(n)+1)$ |
|---:|---:|---:|---:|---:|
| 10 | 9 | 3 | 5 | $359/72$ |
| 11 | 10 | 4 | 5 | $68413/14670$ |
| 12 | 11 | 5 | 5 | $10922951/2485520$ |
| 13 | 12 | 5 | 5 | $566/117$ |
| 14 | 13 | 6 | 5 | $2820675/610127$ |
| 15 | 14 | 6 | 5 | $119/24$ |
| 16 | 15 | 7 | 5 | $19889/4158$ |
| 17 | 16 | 8 | 5 | $73973/16152$ |
| 18 | 17 | 8 | 5 | $2905/594$ |
| 19 | 18 | 9 | 5 | $4089074/862755$ |
| 20 | 19 | 9 | 5 | $3727/750$ |
| 21 | 20 | 10 | 5 | $2353607/486024$ |
| 22 | 21 | 11 | 5 | $18889661/4027368$ |

**VERIFIED-COMPUTATIONALLY.** The same two boundary checks were also completed for every $n=23,\ldots,36$, and for $n=40,45,50,60$, with agreement with the proposed $r=5$ threshold in every case. The per-cell JSON files and `recompute_extended.log` retain the exact certificates. The 26 certificates in the table above are embedded in the standalone verifier; the extended sweep is not substituted for a general proof.

**CONJECTURE.** For $r=5$ and $n\ge10$, the matching full-weight direction may hold at every length. This work does not provide an all-length packing family for that assertion. The upper, non-full-weight direction is proved in Section 5.

**REFUTED, at a smaller length as well.** Without any lower-length qualification, even the $r=5$ assertion fails at $N=4,n=5,m=1$: its proposed threshold is 0, but $\tau^*(4,5,1)=5$. A mass-5 packing is uniform mass 1 on each of the five prefix orbits `001X`, `01X1`, `0X00`, `1X11`, `X001`. This small exception is additional to the $n=17,r=6$ counterexample and does not explain it away.

## 3. Orbit notation and exact incidence calculation

**PROVED.** Put $L=r-1,C=N-L$. The following symmetries preserve the library and $U$: interchange letters 2 and 3 independently in any of the first $L$ coordinates; permute tail coordinates; and permute the three nonzero letters independently in each tail coordinate. Consequently an orbit is described by

\[
(t,w),\qquad t\in\{0,1,X\}^L,\quad 0\le w\le C.
\]

Its size is $2^{x(t)}\binom Cw3^w$. Here (X) abbreviates choices 2 and 3 in an orbit of ordinary words. It is not a wildcard in a pattern. A listed orbit mass $v$ means weight $v/|O|$ on every word of the orbit; for a cover, every such word is an ordinary radius-$m$ ball centre.

**PROVED.** Averaging either a packing or a cover over this finite symmetry group preserves its mass and feasibility. Conversely, verifying the load or coverage at one representative of every fixed-word orbit suffices. Thus no binary-prefix restriction is used: every prefix in ${0,1,X}^L$ is present.

**PROVED.** For a fixed word $c$ and a counted orbit $O$, the incidence coefficient is $|O\cap B_m(c)|/|O|$. Prefix distance counts are obtained by multiplying, coordinate by coordinate, the possible distance contributions of the counted letter. For example an (X) coordinate contributes (1+q) against a fixed 2, and (2q) against a fixed 0 or 1. In the independent verifier, the tail is counted by the recurrence

\[
H_{i+1}(w+[a\ne0],d+[a\ne c_i])
\mathrel{+}=H_i(w,d),\qquad a\in\{0,1,2,3\},
\]

starting from $H_0(0,0)=1$. The counts with the required tail weight are combined with the prefix distance counts and summed through distance $m$. All counts are integers and every resulting coefficient is a `Fraction`.

**VERIFIED-COMPUTATIONALLY.** Before checking the certificates, `verify_j1.py` compares this recurrence with direct enumeration of the four-letter cube at $N=4,r=3$, for all orbit pairs and every radius, totaling 3,645 exact incidence comparisons.

**PROVED.** Membership in $U$ is $w+d_0(t)\le m$, where $d_0(t)=\min_j d_H(t,0^j1^{L-j})$. Packing checks use every centre orbit, including orbits outside $U$. Cover checks use every orbit in $U$, including nonbinary prefixes and points in overlapping library balls.

## 4. The extension to larger libraries, and the direction of E4

**PROVED.** Suppose a radius-$m$ fractional cover of the canonical $s$-pattern library at length $N$ has mass $s-\epsilon$, with $\epsilon>0$. For any $s\le r\le N+1$, the first $s$ patterns of the $r$-pattern library form an isometric copy of the $s$-pattern library: interchange 0 and 1 in tail coordinates $s-1,\ldots,r-2$, using zero-based indices. Apply that isometry to every centre of the cover. Add coefficient 1 on each of the remaining $r-s$ library centres. Every point of the new union lies either in the covered sub-library union or in one of the added balls. The resulting cover has mass $r-\epsilon$.

**PROVED.** Taking $s=5$ extends the non-full-weight theorem of Section 5 to every admissible $r\ge5$. Taking $s=6$ extends the refutations of Sections 1 and 6 to every admissible $r\ge6$.

**PROVED from the supplied E4.** Full weight for a larger library implies full weight for a consecutive sub-library. The reverse implication is unavailable. The pair $M(16,5)=8$, $M(16,6)=7$ explicitly disproves that reverse implication. Therefore a proof for $r=5$ alone cannot establish the positive direction for all $r\ge5$.

**PROVED from the supplied E6.** The five-pattern requirement is essential to the upper theorem: a library with $r\le4$ does not contain the five-pattern sub-library. At $N=19,m=10$, the proposed threshold is 9, but E6 gives full weight for $r=3,4$. Thus the non-full-weight conclusion itself cannot extend to all smaller $r$.

## 5. An explicit general cover above the proposed threshold

### 5.1 Family coefficients

**PROVED.** For $r=5$, let $C=N-4$ and assume $C\ge m\ge6$. Define $D_5=6m(3m-5)$. Each orbit in a row receives the stated numerator divided by $D_5$:

| Prefix or prefixes | Tail weight | Numerator, per orbit |
|---|---:|---:|
| `0000`, `1111` | 0 | $12m^2-15m-6$ |
| `0001`, `0111` | 1 | (18C(m-2)) |
| `0011` | 0 | $6m^2-21m-6$ |
| `0011` | 1 | (18C) |
| `001X`, `X011` | 0 | (20) |
| `0100`, `1010`, `1101` | 0 | (3(m-2)) |
| `011X`, `X001` | 0 | (18m-40) |

**PROVED.** These are 13 orbit masses, all nonnegative. Their sum is

\[
T_5(C,m)=\frac{18Cm-27C+15m^2-3m-38}{3m(3m-5)}.
\tag{3}
\]

**VERIFIED-COMPUTATIONALLY.** This support shape was extracted from exact covers at several lengths. Its validity below is established by exact symbolic inequalities rather than interpolation. The family is not asserted to be optimal; the separate matching certificates in Section 2 establish the finite optimum values.

### 5.2 Coverage inequality for all lengths in the domain

**PROVED.** For a fixed prefix $t$, let $A_s(t,h)$ be the fraction of prefix words in the counted centre-prefix orbit $s$ at distance at most (h) from $t$. Define it as 0 for $h<0$ and 1 for $h\ge L$. A point of tail weight $w=m-k$ is covered by the orbit with prefix $s$ and tail weight 0 with probability $A_s(t,k)$. For tail weight 1 the probability is

\[
A_s(t,k-1)+\frac{m-k}{C}
\left[-A_s(t,k-1)+\frac23A_s(t,k)+\frac13A_s(t,k+1)\right].
\tag{4}
\]

Indeed, the changed tail coordinate is a zero of the point with probability $1-w/C$, a different nonzero letter with probability $2w/(3C)$, and the same nonzero letter with probability $w/(3C)$. These cases change tail distance from $w$ to $w+1,w,w-1$, respectively. Prefix and tail choices are independent.

**PROVED.** Since all counted centres have tail weight at most 1, a point with $k\ge L+1$ lies in every counted ball. For the remaining points it suffices to consider, for each $t\in\{0,1,X\}^L$, every integer $d_0(t)\le k\le L$. The assumptions $C\ge m\ge6$ make each of these radial cases admissible for $L=4$. There are 262 such cases, plus the common interior case, totaling 263.

**PROVED by an exhaustive finite algebra certificate.** Insert the seven numerator formulas in (4), subtract 1 from the resulting total coverage, and multiply by $D_5$. Set $a=m-6$ and $b=C-m$. Every resulting polynomial has nonnegative rational coefficients in $1,a,b,a^2,ab,b^2$. The complete 263 rows are listed in `coverage_positivity_r5.tsv`; there are 82 distinct polynomial vectors. `r5_symbolic_proof()` independently constructs all these rows from the displayed formulas and checks every coefficient. For example:

* `0001`, $k=0$: the numerator is 0.
* `0000`, $k=0$: the numerator is $3a+12$.
* `0110`, $k=1$: the numerator is $3a+2$.

Because $a,b\ge0$ and $D_5>0$, this finite coefficient certificate proves coverage at every $C\ge m\ge6$, not merely at sampled lengths. The first example also proves that the minimum coverage is exactly 1.

### 5.3 Strict shortness and the stated range

**PROVED.** Put $\delta=5m-3C-4$. Directly from (3),

\[
5-T_5(C,m)
=\frac{(6m-9)\delta-3m+2}{3m(3m-5)}.
\tag{5}
\]

If $\delta\ge1$ and $m\ge6$, its numerator is at least $3m-7>0$. Thus the cover is strictly short whenever $C\ge m\ge6$ and $5m>3C+4$.

**PROVED.** For $N\ge11$, take

\[
m_0=\left\lfloor\frac{3N-8}{5}\right\rfloor+1.
\]

Then $m_0\ge6$, $1\le5m_0-3(N-4)-4\le5$, and $m_0\le(3N-3)/5\le N-4=C$. Therefore (3) to (5) give $\tau^*(N,5,m_0)<5$. For $N=9,10$, the explicit covers of masses $359/72$ and $68413/14670$ handle $m_0=4,5$, respectively. They are embedded in the verifier.

**PROVED, using E2 and Section 4.** For every $n\ge10$, every admissible $r\ge5$, and every integer $m>\lfloor(3n-11)/5\rfloor$,

\[
\boxed{\tau^*(n-1,r,m)<r.}
\tag{6}
\]

This is the non-full-weight direction of the proposed theorem, with its exact range stated. It remains valid despite the refutation of the converse.

## 6. An explicit infinite counterexample family

### 6.1 Coefficients

**PROVED.** Let $z\ge3$ be an integer, $N=16+5z$, $r=6$, $C=11+5z$, and $m=8+3z$. Put

\[
Q(z)=3024z^3+28935z^2+84681z+77345.
\]

Each of the two prefix orbits in a row receives the displayed numerator divided by $Q(z)$. There are 28 centre orbits in total.

| Prefix pair | Tail weight | Numerator, per orbit |
|---|---:|---|
| `00000`, `11111` | 0 | $2808z^3+26442z^2+76599z+69500$ |
| `00000`, `11111` | 1 | (3(5z+11)(66z+161)) |
| `00001`, `01111` | 0 | $1674z^3+16524z^2+49671z+46411$ |
| `00001`, `01111` | 1 | $3(5z+11)(72z^2+441z+524)$ |
| `00011`, `00111` | 1 | $3(5z+11)(234z^2+1626z+2267)$ |
| `0001X`, `X0111` | 0 | $2(216z^2+1749z+2438)$ |
| `00100`, `11011` | 0 | $36z^2+255z+274$ |
| `00110`, `10011` | 1 | (9(4z+29)(5z+11)) |
| `0011X`, `X0011` | 0 | $12(57z^2+353z+476)$ |
| `01001`, `01101` | 0 | $72z^2+639z+1007$ |
| `0111X`, `X0001` | 0 | $2(36z^2+282z+311)$ |
| `011X1`, `0X001` | 0 | $2(72z^2+639z+1007)$ |
| `10000`, `11110` | 0 | $2(18z^2+159z+286)$ |
| `1001X`, `X0110` | 0 | (2(51z+211)) |

### 6.2 General coverage and gap proofs

**PROVED.** Every weight and $Q(z)$ is positive on $z\ge3$. Use the same coverage expression (4) with $L=5$. Since $C\ge m\ge17$, the cases $d_0(t)\le k\le5$ are all admissible; the case $k\ge6$ again collects the full mass. There are 902 prefix/radial cases plus the interior case, totaling 903.

**PROVED by an exhaustive finite algebra certificate.** For every one of these cases, multiply coverage minus 1 by $Q(z)$, and substitute $z=u+3$. Every coefficient in the resulting polynomial in $1,u,u^2,u^3$ is nonnegative. The complete 903 rows appear in `coverage_positivity_r6.tsv`; the verifier constructs and checks them using exact polynomial arithmetic. There are 300 distinct polynomial vectors. In particular the numerator vanishes at prefix `0XX00`, $k=2$, so the minimum coverage is exactly 1. This proves coverage for every $z\ge3$.

**PROVED.** Summing the 28 orbit masses gives

\[
\boxed{T_6(z)=6-
\frac{2(90z^2+513z+338)}{3024z^3+28935z^2+84681z+77345}<6.}
\tag{7}
\]

The coefficient tables establish coverage for all these lengths, while (7) establishes the strict deficit. Also

\[
\left\lfloor\frac{3(N+1)-11}{5}\right\rfloor
=\frac{3N-8}{5}=8+3z=m.
\]

Thus every member refutes the claimed full-weight direction. Section 4 extends every such member from $r=6$ to every $6\le r\le N+1$.

**VERIFIED-COMPUTATIONALLY.** Three explicit covers bridge the initial parameters $z=0,1,2$:

| $z$ | $N$ | $m$ | Cover mass |
|---:|---:|---:|---:|
| 0 | 16 | 8 | $199683/33307<6$ |
| 1 | 21 | 11 | $5991123/1000000<6$ |
| 2 | 26 | 14 | $297552589/49684499<6$ |

For $z=0,2$, a matching packing proves the mass is optimal. The $z=1$ cover is rounded upward coefficient by coefficient to denominator $10^6$, then checked exactly; optimality is not asserted for that rounded cover. These three finite cases and (7) prove the infinite obstruction for every integer $z\ge0$.

## 7. Verification, controls, and the exact scope of the result

**VERIFIED-COMPUTATIONALLY.** Run the standalone standard-library verifier from any working directory:

```sh
/usr/local/bin/python3 -B /absolute/path/to/j1/verify_j1.py
```

It imports no repository modules and uses no optimization packages. Its embedded certificates include both sides of all 13 $r=5$ boundary pairs, both sides of the decisive $r=6$ optimum, the $m=7$ full packing, the other two initial counterexample covers, and the literal $n=5$ exception. It additionally checks the general formulas numerically at the listed $n,r,z$, after proving their symbolic inequalities. `--n`, `--r`, and `--z` change those lists. Larger-$r$ cases use the proved cover-extension reduction; the program does not claim to enumerate their exponentially larger quotients.

**VERIFIED-COMPUTATIONALLY.** All five negative controls were executed, and every invalid certificate exited with code 1:

| Control | Exact rejected value |
|---|---|
| Multiply the decisive cover by $999/1000$ | Minimum coverage $999/1000<1$ |
| Multiply the decisive packing by $1001/1000$ | Maximum load $1001/1000>1$ |
| Use the analytic $r=5$ formula at $N=9,m=4$, outside $m\ge6$ | Minimum coverage $13/14<1$ |
| Use the analytic $r=6$ formula at $z=2$, outside $z\ge3$ | Minimum coverage $382641/386639<1$ |
| Replace $Q(z)$ by the uncorrected larger denominator $Q(z)+12z+87$, at $z=3$ | Minimum coverage $673451/673574<1$ |

**VERIFIED-COMPUTATIONALLY.** The last control is material: an initially fitted rational continuation was infeasible. The corrected denominator in Section 6 repairs its four deficient prefix/radial rows and still gives the strictly positive deficit (7). The failed raw candidate is retained as a candidate artifact, not used as evidence of coverage.

```sh
/usr/local/bin/python3 -B /absolute/path/to/j1/run_controls.py
```

**PROVED / REFUTED / CONJECTURE  -  final scope.**

* **REFUTED:** the requested equivalence for every $r\ge5$, already at $n=17,r=6$, and at infinitely many larger lengths.
* **PROVED:** its non-full-weight direction for every $n\ge10$ and every admissible $r\ge5$, using an explicit rational cover family and two finite exceptions.
* **VERIFIED-COMPUTATIONALLY:** the matching $r=5$ values at all $n=10,\ldots,22$, with further finite checks through $n=60$ as stated above.
* **CONJECTURE:** the all-length $r=5$ full-weight direction on $n\ge10$. No all-length full-mass packing family is claimed.

**VERIFIED-COMPUTATIONALLY, provenance.** Numerical linear programming was used only to propose supports and bases. Candidate vectors were reconstructed as rationals, sometimes by solving a proposed active basis over the rationals Q. A decision was accepted only after exact load, coverage, and mass checks. The standalone verifier uses a different incidence implementation from the search scripts. All produced files are confined to this `j1` directory; no repository edits, commits, pushes, or cluster jobs were used.
