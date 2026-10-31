# PROVED: every fixed run count has slope 3/5 and eventual period five

Date: 2026-10-29. The alphabet has four letters. Write N=n-1, L=r-1,
C=n-r, and use the canonical staircase from the question.

**PROVED (uniform packing, all admissible lengths and run counts).**

\[
\boxed{\tau^*(n-1,r,m)=r\quad\text{whenever}\quad
0\le m\le\left\lfloor\frac{3(n-r)}5\right\rfloor.}
\tag{1}
\]

**PROVED (all fixed r>=5).**

\[
\boxed{\lim_{n\to\infty}\frac{M(n-1,r)}n=\frac35.}
\tag{2}
\]

**PROVED (eventual exact structural form, all fixed r>=5).** There are integers
kappa[r,0],...,kappa[r,4] and a finite Q[r] such that

\[
\boxed{M(5q+s-1,r)=3q+\kappa[r,s]
\quad(s=0,\ldots,4;\ q\ge Q[r]).}
\tag{3}
\]

The constants satisfy the explicit finite bounds

\[
\left\lfloor\frac{3(s-r)}5\right\rfloor
\le\kappa[r,s]\le
\left\lfloor\frac{3s-11}5\right\rfloor.
\tag{4}
\]

**CONJECTURE / unresolved.** This proof does not compute kappa[r,s] or the last
exceptional length for arbitrary r. Equation (3) is an existence and
eventual-form theorem, not an explicit numerical F(n,r). The complete
all-length two-variable formula remains unresolved.

## 1. PROVED: a uniform r-orbit packing

At every library prefix put total mass one on its uniform weight-m tail
orbit. Thus the r supported prefix orbits are exactly the library prefixes,
and every supported word has distance m from its own pattern and
m+|i-j| from any other pattern. The support is private, and the total mass
is r. Assume m<=3C/5, so these tails are admissible. C=0 means m=0 and is
immediate.

Consider any centre, with prefix c and p nonzero tail coordinates. For a
supported tail let A be the number of its occupied positions among those p
centre positions, and E the number of matching nonzero letters. Relative
to distance m from the zero tail, its distance reduction is

\[
S=A+E-p,\qquad A\sim\operatorname{Hypergeom}(C,m,p),\quad
E\mid A\sim\operatorname{Bin}(A,1/3).
\]

For a fixed integer S>=0, the number of library prefixes within prefix
distance S of c is at most 2S+1. Indeed, their smallest and largest indices
i,j satisfy |i-j|=d(P_i,P_j)<=2S, so there are at most 2S+1 indices.
For S<0 none are accepted. Therefore the exact load is at most

\[
\mathbb E[g(S)],\qquad g(s)=\max(0,2s+1).
\tag{5}
\]

## 2. PROVED: replacing the hypergeometric count by a binomial count

For fixed p define

\[
h(a)=\mathbb E\left[g\!\left(Y_1+\cdots+Y_a-p\right)\right],
\qquad\Pr(Y_i=1)=2/3,\quad\Pr(Y_i=2)=1/3.
\]

This is the conditional expectation in (5) given A=a. The function h is
discretely convex. To see this, g(x+y)-g(x) is nondecreasing in x for each
y>0. Averaging over Y gives a nondecreasing increment function. Adding
another positive Y to the partial sum can only increase the expected
increment, which proves h(a+2)-h(a+1)>=h(a+1)-h(a). Linear interpolation
makes h convex on the real interval [0,p].

**PROVED: the needed sampling convex order.** For any fixed population of
C real numbers and p<=C, a sum sampled without replacement is dominated in
convex order by a sum of p independent samples with replacement. Here is
an explicit coupling, so no outside theorem is needed.

First sample p distinct population indices I_1,...,I_p uniformly. Independently
simulate the multiplicity counts from p draws with replacement from C
labels, retain their positive counts, pad them with zeros to p slots, and
permute those p counts uniformly. Call the resulting counts V_1,...,V_p.
Their sum is p and, by symmetry, each has expectation one. Conditional on
the number of positive slots, the selected population labels and their
multiplicities have exactly the with-replacement sampling distribution.
Consequently X=sumV_ia[I_i] has that distribution, while

\[
\mathbb E[X\mid I_1,\ldots,I_p]=\sum_i a[I_i]=Y.
\]

Jensen's inequality gives Ephi(Y)<=Ephi(X) for every convex phi. Apply this to
the population with m ones and C-m zeros, and to phi=h. This proves

\[
\mathbb E[h(A)]\le\mathbb E[h(B)],\qquad
B\sim\operatorname{Bin}(p,m/C).
\tag{6}
\]

The right side is the expectation of g of a sum of p independent gains
with probabilities (1-rho,2rho/3,rho/3) on (-1,0,+1), where rho=m/C.
Increasing rho increases each gain in a coupling. Since g is increasing,
rho<=3/5 lets us replace it by 3/5.

## 3. PROVED: the exact biased-walk inequality

Let Z_p be the sum of p independent increments with probabilities
(2/5,2/5,1/5) on (-1,0,+1). Swapping positive and negative steps gives

\[
\Pr(Z_p=-k)=2^k\Pr(Z_p=k),\qquad k\ge1.
\tag{7}
\]

For g(s)=max(0,2s+1), the one-step expected change is -2/5 at s>=1,
zero at s=0, 1/5 at s=-1, and zero at s<=-2. Hence (7) gives

\[
\begin{aligned}
\mathbb E[g(Z_{p+1})]-\mathbb E[g(Z_p)]
&=-\tfrac25\Pr(Z_p\ge1)+\tfrac15\Pr(Z_p=-1)\\
&=-\tfrac25\Pr(Z_p\ge2)\le0.
\end{aligned}
\tag{8}
\]

Initially g(Z_0)=1. Equations (5) to (8) bound every packing load by one.
The packing mass is r and the library cover mass is r, proving (1).

**PROVED.** Combining (1) with supplied E4 and E21, for every r>=5,

\[
\left\lfloor\frac{3(n-r)}5\right\rfloor
\le M(n-1,r)\le
\left\lfloor\frac{3n-11}5\right\rfloor.
\tag{9}
\]

For each fixed r, division by n proves (2). This is a fixed-r limit; no
claim takes r proportional to n.

## 4. PROVED: a finite centre-tail cutoff depending only on r

Any mass-r packing must be private: summing its r library-ball loads counts
every point at least once, and any positively weighted overlapping point
would make that sum exceed r. The supplied sphere fact E3 then puts its
support in orbits (a,m-d(a)), where a ranges over the finitely many uniquely
nearest prefixes in {0,1,X}^L and 0<=d(a)<=L. Averaging over the library's
symmetries loses no feasibility or mass. Each private region has mass one.

For n>=13r, C=n-r>=12r. Every radius in the interval (9) has

\[
\frac mC\le\frac{3(C+r)}{5C}\le\frac{13}{20}.
\]

Take t=6/5. The sampling-without-replacement elementary symmetric mean
bound, proved by pairwise equalization as in the other proof files, gives
for any supported tail of weight w<=m

\[
\mathbb E[t^S]\le
\left[\frac{1-w/C}{t}+\frac wC\frac{2+t}{3}\right]^p
\le\left(\frac{197}{200}\right)^p.
\]

If its prefix distance from the centre is delta, acceptance requires S>=delta-d.
Markov's inequality bounds its share by t^(d-delta)(197/200)^p. Summing a
total mass r and using d<=L,delta>=0 gives

\[
\ell(c,p)\le r(6/5)^{r-1}(197/200)^p.
\tag{10}
\]

For integer r>=1, r<=2^(r-1). The exact scalar inequality

\[
\frac{12}{5}\left(\frac{197}{200}\right)^{60}<1
\]

therefore implies that (10) is strictly below one whenever p>=60r.
All such centre constraints are redundant for every mass-r candidate.
Thus for n>=61r the full-weight question is an LP on a fixed set of private
prefix variables and centre types (c,p), 0<=p<60r. Their number depends on
r, but not on n.

## 5. PROVED: eventual feasibility is constant in each residue and offset

Fix r>=5, s in {0,...,4}, and write n=5q+s. By (9) only the finitely many
integer offsets k in (4) can occur in m=3q+k. For all sufficiently large q,
all private-prefix tail weights m-d(a) lie in [0,C].

Each incidence numerator in the finite LP of Section 4 is a polynomial in
q, of degree at most p. Its denominator is a positive multiple of
binom(C,p)3^p. This follows directly from the gain count R in the other
proof files: each term has factors binom(m-d,a)binom(C-m+d,p-a).
Clearing positive denominators gives a finite LP whose coefficients are
polynomials in q. Private-region masses equal one, so this polytope is
bounded, lying in [0,1]^D for a fixed finite D.

If feasible, a bounded nonempty polytope has a vertex. Its coordinates
solve one of finitely many active square bases. Every non-identically-zero
basis determinant is a polynomial in q, and Cramer's rule makes all vertex
coordinates rational functions of q. Feasibility of that basis is determined
by signs of finitely many polynomials. Each nonzero polynomial has only
finitely many roots and has a constant sign beyond its largest real root.
Identically zero determinants cannot define a basis; exceptional roots
contribute only finitely many q values.

Consequently each possible basis is either feasible at every sufficiently
large q or infeasible at every sufficiently large q. The finite union of
bases has the same property. For each fixed (r,s,k), full-weight feasibility
therefore stabilizes. Maximizing over the finite nonempty offset interval
gives an integer kappa[r,s] which stabilizes as well. Taking the largest of the
five stabilization points proves (3), with bounds (4).

**PROVED: scope.** The eventual period-five shape is universal for every
fixed r>=5. It does not assert a short expression for its five constants,
an explicit general exception window, or a uniform threshold when r grows
with n. The already proved r=6,7 formulas compute those missing quantities
in those cases; the general explicit F(n,r) remains open.

## 6. VERIFIED-COMPUTATIONALLY: exact supporting checks and control

`verify_general.py` checks every p<=C at C=0,...,100, with m=floor(3C/5),
including the exact hypergeometric-to-binomial convex-order comparison.
It checks the biased-walk reflection and drift identities through p=100,
the chain count for all prefix classes through r=10, and the exact scalar
inequality used in the 60r cutoff. These support the analytic proof;
the all-parameter conclusions above do not extrapolate from that finite grid.

**REFUTED / executed control.** At C=5,r=3,m=4, raising the uniform-orbit
construction above its allowed tail density preserves privacy and all
region masses but gives load 4/3 at the middle library prefix with one
nonzero centre-tail position. The verifier must reject it.

```sh
/usr/local/bin/python3 -B verify_general.py
```
