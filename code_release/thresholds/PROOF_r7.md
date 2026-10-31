# PROVED: the seven-run threshold, including its finite exceptions

Date: 2026-10-29. Alphabet size four; N=n-1.

**PROVED, by exact computer-assisted certificates and the six-run theorem.**

\[
\boxed{M(n-1,7)=\left\lfloor\frac{3n-14}{5}\right\rfloor+
 \mathbf 1_{\{n\in E\}},\qquad n\ge10,}
\]

\[
E=\{11,14,16,21,26,31,36,41\}.
\]

**REFUTED.** The initially plausible continuation floor((3n-13)/5) (even after
allowing the n=14 exception) is wrong. At n=46 its value is 25, while the exact
threshold is 24. Section 3 proves infinitely many such failures at n=46+5z.

## 1. PROVED: a uniform full-weight lower bound

For n>=15 put C=n-7 and b=floor((3C+7)/5)=floor((3n-14)/5). Each prefix below has a
uniform weight-(b-d) tail of length C; X expands independently to {2,3}.

| Prefix pair | d | Total mass per prefix |
|---|---:|---:|
| 000000, 111111 | 0 | 2/13 |
| 000011, 001111 | 0 | 3/13 |
| 00011X, X00111 | 1 | 7/26 |
| 0011X1, 0X0011 | 1 | 3/13 |
| 00X000, 111X11 | 1 | 11/13 |
| 01110X, X10001 | 2 | 3/13 |
| 011110, 100001 | 1 | 5/13 |
| 0111X0, 1X0001 | 2 | 5/13 |
| 0X001X, X011X1 | 2 | 7/13 |
| 10011X, X00110 | 2 | 3/13 |

**PROVED by prefix distances.** All twenty orbits are private sphere orbits;
each of the seven regions receives mass one. Tail weights are admissible
for C>=8. The load formula and degree argument are exactly those proved in
`PROOF_r6.md`, with prefix length six, scale 52, and beta=7 in b=floor((3C+beta)/5).

**VERIFIED-COMPUTATIONALLY.** All 729 centre prefixes are represented by 378
profiles and 54 maximal profiles. For all five residues and p=0,...,58,
15,930 polynomials have 477,900 nonnegative Newton coefficients. No shifts
or finite exceptions are needed. `prove_family.py` regenerates them from
`family_r7_lower14.json`.

**PROVED.** For p>=59 the same sampling-without-replacement argument gives

\[
\ell(c,p)\le\frac{1433}{343}
 \left(\frac{30221}{30975}\right)^{59}<\frac{99}{100}<1.
\]

The tail-weight proportion used here is at most
3/5+7/(5 * 59)=184/295. Thus the packing is feasible at b at every n>=15,
and E16 establishes every smaller radius. In particular M(n-1,7)>=b.

## 2. PROVED: the first upper-bound residue

**PROVED by a parametric private-dual certificate.** For every z>=0, at

\[
n=19+5z,\qquad r=7,\qquad m=9+3z,
\]

full weight is impossible. `upper_r7.json` contains the complete certificate.
For z>=27 its nonnegative dual multipliers are polynomials divided by a
positive polynomial Q(z). For each of the 339 private prefix orbits, their
incidence-weighted sum is at least one, while their sum is strictly below
seven. All relevant polynomial inequalities have nonnegative coefficients
after z=u+27. Exact finite certificates cover z=0,...,26.

**PROVED: why a private dual suffices.** A packing of mass seven must saturate
the seven library-ball constraints. Since the sum of their indicators is at
least one everywhere on U, any packing mass on an overlapping point would
make that summed load exceed seven. Consequently every mass-seven packing
is supported entirely on private points. A feasible dual on all private
points, of mass strictly below seven, rules out such a packing. Finiteness of
the LP and the mass-seven library cover then imply tau*<7. The private-dual
mass is not asserted to be the covering optimum or a cover of all of U.

**PROVED using restriction and the six-run formula.** Outside n congruent_to 4 mod 5,
floor((3n-12)/5)=floor((3n-13)/5), so M(n-1,7)<=floor((3n-13)/5) follows from r=6.
The displayed parametric family handles every n congruent_to 4 mod 5 with n>=19.
Together they give this upper bound for all n>=15.

## 3. PROVED: the later transition at n=46

**PROVED by a second parametric private-dual certificate.** Full weight is
impossible for every integer z>=0 at

\[
n=46+5z,\qquad r=7,\qquad m=25+3z.
\]

This is the remaining residue n congruent_to 1 mod 5, at the radius one above b. The
certificate `upper_r7_n46.json` has a degree-eight common denominator and
21 centre groups. Positivity is established after z=u+59, with 59 exact
finite certificates for z=0,...,58. The last bridge is n=336,m=199.

**VERIFIED-COMPUTATIONALLY.** `verify_additional.py` independently regenerates
the incidence polynomials from coordinate-wise tail counts at z=0,1,2.
Each centre has tail weight at most two, so the incidence numerator has
degree at most two. Multiplying by a common positive tail denominator
9C(C-1) produces polynomial inequalities. It checks nonnegative coefficients
after the declared shift, strict denominator positivity, strict shortness,
and all finite bridges, without SymPy, SciPy, or stored slack expressions.

**PROVED.** The bounds floor((3n-13)/5) and b differ only when n congruent_to 1 mod 5.
The second family removes that remaining excess at every such n>=46.
Thus M=b for all n>=15 except possibly n in {16,21,26,31,36,41}.

## 4. VERIFIED-COMPUTATIONALLY: the finite exceptions

Both boundary certificates were verified at every n=10,...,50. They give

| n | M(n-1,7) | b |
|---:|---:|---:|
| 11 | 4 | 3 |
| 14 | 6 | 5 |
| 16 | 7 | 6 |
| 21 | 10 | 9 |
| 26 | 13 | 12 |
| 31 | 16 | 15 |
| 36 | 19 | 18 |
| 41 | 22 | 21 |

At n=10,12,13,15 the value equals b. The lower bound in Section 1,
the two upper families, the six-run theorem, and these finite certificates
establish the boxed all-length formula.

## 5. PROVED: an obstruction to the naive construction

**PROVED.** There is no single constant vector of private-prefix orbit
masses that realizes radius seven at n=16 and also stays feasible in the
limiting tail model with nonzero proportion 3/5. All 339 private prefixes
are allowed in this statement, not just a selected support.

The exact Farkas certificate in `constant_mass_obstruction_r7.json` combines
finite-centre constraints and limiting-centre constraints. Its minimum
coefficient on every private-prefix mass is one, but its total multiplier is

\[
\frac{3922162576875}{562499999968}
=7-\frac{15337422901}{562499999968}<7.
\]

Any constant mass-seven vector satisfying all those constraints would give
the contradiction 7<=that strict upper bound. A construction valid at every
length at asymptotic density 3/5 must satisfy the limiting constraints by
continuity. This explains why a single K3-style vector cannot realize all
the higher finite exceptions. It does not rule out closed forms: the boxed
formula itself is one, using a separate family for its eventual boundary.

## 6. VERIFIED-COMPUTATIONALLY: execution

```sh
/usr/local/bin/python3 -B prove_family.py family_r7_lower14.json
/usr/local/bin/python3 -B verify_additional.py
/usr/local/bin/python3 -B extra_controls.py
```

**REFUTED / executed controls.** A finite-only constant candidate fails at
n=107 with maximum load 412780571/411334200. The first inverse-C candidate
fails at n=68 with maximum load 1514758/1507005. Scaling down the symbolic
upper multipliers or the constant-family obstruction multipliers is also
rejected. Every control outcome is recorded in the two control-result files.
