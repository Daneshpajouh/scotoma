# Independent exact maximum-load verification

## Independence and coordinate convention

This verifier extends the permitted ZX independent tensor lineage. The old
`ZX_independent_maxload/METHOD.md` and `independent_maxload.cpp` supplied the
three-state tensor transform. The constant-block compression, parser, exact
single-centre kernel, arithmetic selection, tests and handoff here were written
for this brief. No ZR3 implementation, upper_rows implementation, gap_census
implementation, or dp_resolved implementation was read, imported or invoked.
Only the permitted ZR3 README, decision certificates and RESULT data, and the
two named calibration data files, were mathematical/format inputs.

The certificate convention is Q_j = 0^j 1^(L-j) 0^C. Complementing 0 and 1
in the first L positions takes Q_j to the brief's P_j = 1^j 0^(N-j).
It fixes X and the constant block and preserves every Hamming distance.
The parser performs this conversion. Output argmax prefixes use the brief's
P convention; final results also give the certificate convention. A literal
test checks this isometry at every one of the 4^5 words of a fixture.

## Exact measure and support

For a point class (p,w), write x(p) for the number of X symbols. Its size is

    s(p,w) = 2^x(p) binom(C,w) 3^w.

A certificate mass a(p,w) is total class mass, so its per-concrete-word density
is a/s. Duplicate classes are merged as Fractions. Floating masses, negative
masses, invalid symbols and weights are rejected. Zero masses are discarded.
The common denominator D is the least common multiple of the reduced density
denominators, and q(p,w) = D a(p,w)/s(p,w) is an exact positive integer.
The total mass is recomputed with Fraction, independently of the producer.

Every word of a point class has the same distance to each P_j:

    w + #{i: p_i differs from 1 if i<j, and from 0 otherwise}.

An X always differs from either binary letter. The minimum over all j must
be at most m for each positive class. Python and C++ both perform this check.

## Why the scan covers every unrestricted centre

Swapping letters 2 and 3 independently at an informative position preserves
the measure and Hamming distance. In the constant block, independent
permutations of the three nonzero letters and permutations of coordinates
also preserve the measure. The complete centre space therefore has precisely

    (C+1) 3^L

classes. A class is represented by its ternary prefix and constant-block
weight v, using suffix 1^v 0^(C-v). Centres outside U are included.
No support restriction or reflection restriction is applied to centres.
The final scan counts both classes and their concrete multiplicities,
checking that the latter sum to 4^N.

## Tensor convolution with a compressed constant block

The informative-coordinate distance operator on per-word densities is

    A(z) = [[1,z,2z], [z,1,2z], [z,z,1+z]].

The inherited integer transform is

    T(a,b,c) = (a+b+2c, a-b, a+b-2c),

with inverse

    ((s+2d+e)/4, (s-2d+e)/4, (s-e)/4).

Direct multiplication diagonalizes A(z) with diagonal
(1+3z, 1-z, 1-z). A nonzero transformed ternary digit contributes one
nontrivial frequency coordinate.

For the constant block, a radial density f_w transforms by the integer matrix

    K[s,w] = sum_j (-1)^j binom(s,j) binom(C-s,w-j) 3^(w-j).

One derivation is to transform the density z^weight: at frequency weight s its
generating polynomial is (1+3z)^(C-s)(1-z)^s. Taking coefficients gives K.
The inverse is K/4^C. The code verifies every entry of K K = 4^C I on every run.
This compression stores C+1 values instead of 3^C ternary suffix values.

For total frequency weight k, the radius-m ball eigenvalue is

    lambda_k = sum_{d=0}^m sum_{j=0}^d
                   (-1)^j binom(k,j) binom(N-k,d-j) 3^(d-j).

The engine populates q, applies T on every prefix axis and K on the suffix
axis, multiplies by lambda_(prefix_frequency_weight + suffix_frequency_weight),
then inverts K and every T. It scans every resulting exact load numerator.
There is no numerical optimization, screening, floating error band, CRT,
sample-based feasibility decision, or producer executable in this path.

The principal array has (C+1)3^L entries. Its layout is v*3^L + prefix_id,
where prefix_id = sum_i digit(p_i)3^i and digit(X)=2. Ties select the smallest
array index. Each thread owns disjoint triples/rows or scan intervals.

## Integer safety and checks

Let M = sum q(p,w)s(p,w), and V = sum_{d=0}^m binom(N,d)3^d. Partial forward
coefficients have absolute value at most M, since absolute transform
coefficients are bounded by concrete class multiplicities. Each |lambda|
is at most V. Every inverse-prefix row has absolute row sum at most one;
each inverse-radial row also does, since sum_w |K[s,w]| <= 4^C.
Thus stored inverse values are bounded by MV; undivided inverse numerators
are bounded by max(4,4^C)MV. The weighted positive mass/load reductions are
bounded by M and MV respectively.

C++ first computes the bound with cpp_int. Signed 128-bit arithmetic is
selected only when it is strictly below 2^120, leaving spare bits for signed
expression temporaries. Otherwise every array entry and arithmetic operation
uses cpp_int. The large-r target cells fit the signed-128 bound; long rational
denominators select cpp_int for the other cells as needed. No fraction is
rounded, and every inverse division checks its exact remainder.

Additional checks are:

1. Recount M from the populated dense array using concrete multiplicities.
2. Require every final load numerator to be nonnegative.
3. Scan exactly (C+1)3^L entries and count exactly 4^N concrete centres.
4. Sum loads with centre multiplicities and require exact equality to M V.
5. Recompute the returned argmax load using the separate direct kernel below.
6. Require the maximum numerator to equal D exactly for a passing certificate.

## Separate single-centre kernel

The Python counting kernel never uses the tensor transform. For a point
prefix at a fixed centre prefix, a binary mismatch contributes z, an X at a
binary centre contributes 2z, and an X at an X centre contributes 1+z.
For the suffix it iterates concrete alphabet choices, retaining nonzero
point weight and Hamming distance. At a zero centre position the choices
are (zero, distance 0, multiplicity 1) and (nonzero, distance 1, multiplicity 3).
At a nonzero centre position they are (zero, distance 1, multiplicity 1),
(equal nonzero, distance 0, multiplicity 1), and (other nonzero, distance 1,
multiplicity 2). The exact cumulative counts give the ball intersection.
Dotting these counts with q and dividing by D gives an exact Fraction load.

## Calibration and adversarial controls

Literal four-letter enumeration tests every concrete point and centre in
15 small cases, including m=0, m=N, C=0 and L=0. Every tensor load and every
direct-kernel load agrees with enumeration. Both arithmetic backends are
tested on the same objects. All suffix histograms through C=5 are also
compared to literal enumeration. A two-thread/single-thread comparison is
additionally required by the smoke test on Linux Slurm; it is not run on
the Mac because of the local process limit.

The N=19, r=3, m=11 gap certificate has an independently scanned exact
maximum 1 and mass 32/11. All 162 tensor entries agree with the direct kernel.
The separate N=19, r=12, m=10 calibration has mass
10.6014609356056... and was also completely scanned locally, after its
estimated allocation was shown to fit the binding resource limits. All
1,594,323 centre classes were checked; the exact maximum was 1 at canonical
prefix 00000000000 and suffix weight 0. It is repeated in the cluster
smoke before a target verdict is accepted.

For each target, two actual wrong certificate JSON files are saved and
read back. One class mass is multiplied by 1000000001/1000000000. Its load at
an independently counted unit-load centre exceeds one by an exact positive
Fraction, proving rejection without another global scan. Another class is
moved to (X^L,C), whose distance from every P_j is N; support validation
rejects it. The controls are repeated at the global argmax after each full
scan. A unit-load seed alone establishes a lower bound, never feasibility.

## Execution states and handoff

Local audits distinguish mass/support/control checks from global maximum
verification. The complete smoke includes literal tests, both N=19
calibrations and the smallest target. Only a successful smoke writes its
hash-bound gate. Full array tasks check that gate against the exact current
source/input identity. Each target writes H3 only after its complete exact
scan, direct argmax check and both controls pass. Missing output, a timeout,
an allocation refusal, an arithmetic discrepancy or a failed control
cannot produce H3. No dense array is written to disk.
