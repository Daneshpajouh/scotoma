# J2 result: (S) is refuted; universal (W) and the threshold remain open

**REFUTED.** (S) is false. At \((m,N)=(19,41)\),
\(c=1^{19}2^40^{18}\) has

\[
\Lambda(c)=\frac{2571333812653}{1837953717600}
>\frac{380}{273}=\max(1,L),\qquad
\Lambda(c)-L=\frac{13009956653}{1837953717600}>0.
\]

**PROVED.** There is also an infinite family of counterexamples in the fully
interior regime: for every integer \(m\ge35\), set
\(N=2m+6\) and \(c=1^{m+1}2^40^{m+1}\). Then \(\Lambda(c)>L>1\).
Thus the correction note's universal reading of a "one length wide" failure
band is false, and failure does not require a contributing end pattern.

**CONJECTURE  -  UNSETTLED.** Universal (W), quantifying over every centre with
\(u\ge2\) and every \(N\ge N^*(m)\), has not been proved or refuted here.
Consequently this deliverable does **not** make the threshold unconditional.
The negative result is accompanied by exhaustive exact search, a proof for the
interior four-free-letter family, and failing controls.

## Exact object and scope

**PROVED  -  definitions.** This report concerns the weighting in E11: the full
chain \(P_j=1^j0^{N-j}\), \(j=0,\ldots,N\), on four letters. Complementing the
two binary letters gives the user's canonical full staircase. Set
\(A_j=\{1,\ldots,N\}\setminus\{j,j+1\}\), ignoring excluded indices outside
the coordinate range, and let \(S_j\) consist of words with exactly \(m\)
coordinates of \(A_j\) replaced by letters 2 or 3. Put unit mass uniformly on
each \(S_j\). Write

\[
\Lambda(c)=\sum_{j=0}^N\frac{|B_m(c)\cap S_j|}{\binom{|A_j|}{m}2^m},
\qquad V(m,N)=\max_c\Lambda(c).
\]

**PROVED  -  distinction.** \(V\) is a maximum load of this specified weighting.
It is neither \(\tau^*\), nor \(K\), nor \(K_{\rm sound}\), nor the user's
budget threshold \(M(N,r)\). Failure of this weighting alone implies no strict
upper bound on any of those covering quantities.

## Exhaustive verification

**VERIFIED-COMPUTATIONALLY.** The standalone `verify_j2.py` implements the exact
count from its elementary selection condition, without importing repository
code. For a pattern, let \(u\) be the global number of free coordinates,
\(e\) the global number of binary mismatches, \(a\) the number of free
coordinates in \(A_j\), \(b\) the number of binary mismatches in \(A_j\), and
\(M_j=|A_j|\). Selecting \(A\) free coordinates and \(B\) binary mismatches,
with \(s\) wrong free letters among the selected coordinates, gives

\[
d(x,c)=m+u+e-2A-B+s.
\]

**PROVED.** Therefore the exact count is

\[
\sum_{A,B}\binom aA\binom bB\binom{M_j-a-b}{m-A-B}2^{m-A}
\sum_{s=0}^{\min(A,\,2A+B-u-e)}\binom As,
\]

with an empty inner sum interpreted as zero. The excluded free coordinates
are already included in global \(u\); counting them a second time is wrong.

**PROVED.** To exhaust centres, fix their global totals \((Z,O,u)\). After a
prefix of length \(p\), retain its zero and one counts \((z,o)\) and last
letter. The next letter completes the contribution for pattern \(j=p\), since
its binary mismatch count is \(e=z+O-o\). Every word is a path; states with
equal counts and last letter have identical possible future scores. Max-plus
state merging thus preserves the exact maximum. Endpoint contributions use
\(N-1\) allowed coordinates, interior contributions \(N-2\). Common integer
scaling avoids rounding. All final comparisons use `fractions.Fraction`.

**PROVED.** Independently swapping letters 2 and 3 at any coordinate preserves
the spheres and loads. Thus searching all \(3^N\) words on \(\{0,1,X\}\)
covers all \(4^N\) actual centres. A class with \(u\) free coordinates has
\(2^u\) representatives. If \(u>2m\), every sphere word has distance at least
\(u-m>m\), so its load is zero; these classes are discharged by proof.

**VERIFIED-COMPUTATIONALLY.** Exact exhaustive outputs are:

| \(m\) | \(N\) | searched domain | maximum | maximum for \(u\ge2\) |
|---|---|---|---|---|
| 19 | 41 | all \(3^{41}\) classes | \(2571333812653/1837953717600\) | same, at \(u=4\) |
| 19 | 51 | all \(3^{51}\) classes | \(1520/1519=L(19,51)\) | \(1273/1519\), at \(u=2\) |
| 19 | 52 | all \(3^{52}\) classes | \(1\) | \(323/400\), at \(u=2\) |
| 35 | 76 | every \(u=4\) class | \(349741/259296\) | only \(u=4\) searched |

**VERIFIED-COMPUTATIONALLY.** The respective evaluated DP-state counts were
3,664,284; 10,245,729; 11,206,104; and 918,340. These are compressed exhaustive
searches, not literal visits to every word. All per-\(u\) maxima and attaining
words are retained in `exhaustive_N*_m*.jsonl`. `summary.json` checks class
coverage, including \(\sum_u2^u\#\text{classes}_u=4^N\) for each full search,
and independently checks 770 local contributions with the released counter.

## A proved formula for the interior four-free-letter family

**PROVED  -  theorem.** Let \(m\ge4\), \(a,b\ge m+1\),
\(N=a+b+4\), \(c=1^a2^40^b\), and \(M=N-2\). Then

\[
F_4(m,N):=\Lambda(c)
=\frac{9\binom{M-3}{m-2}+30\binom{M-3}{m-3}
             +16\binom{M-3}{m-4}}{4\binom Mm}
\]

\[
=\frac{m(m-1)(9N^2+12Nm-87N-5m^2-23m+174)}
 {4(N-2)(N-3)(N-4)(N-m-1)}.
\]

**PROVED  -  proof.** The five block-adjacent patterns are \(j=a+t\),
\(t=0,\ldots,4\). There are no binary mismatches there. At \(t=0,4\), three
free coordinates are allowed; acceptance selects either two, both matching,
or three, with at most two wrong letters. Each endpoint contribution is

\[
\frac{\frac34\binom{M-3}{m-2}+\frac78\binom{M-3}{m-3}}{\binom Mm}.
\]

**PROVED.** At \(t=1,2,3\), exactly two free coordinates are allowed, both must
be selected and matched, giving \(\binom{M-2}{m-2}/(4\binom Mm)\) each.

**PROVED.** Off the block, put \(e=E+1\ge1\), the number of binary mismatches.
Exactly one is excluded, and all four free coordinates are allowed. The
contribution on either side, multiplied by \(\binom Mm\), is

\[
\begin{split}
T_E={}&2\binom{M-4-E}{m-3-E}
+\frac E2\binom{M-4-E}{m-2-E}
+\frac{15}{16}\binom{M-4-E}{m-4-E}\\
&+\frac{11E}{16}\binom{M-4-E}{m-3-E}
+\frac5{16}\binom E2\binom{M-4-E}{m-2-E}
+\frac1{16}\binom E3\binom{M-4-E}{m-1-E}.
\end{split}
\]

**PROVED.** This follows by selecting three free coordinates (zero or one
unselected binary mismatch), or all four (zero through three unselected binary
mismatches), in the preceding exact count. The letter probabilities for four
selected free coordinates are \(15,11,5,1\), divided by 16. All other choices
are impossible. Since \(a,b\ge m+1\), each side includes every potentially
contributing level, and the actual end patterns have zero contribution.

**PROVED.** Generalized Vandermonde gives

\[
\sum_{E\ge0}\binom Eh\binom{M-4-E}{m-k-E}
=\binom{M-3}{m-k-h}.
\]

**PROVED.** To see it combinatorially, count subsets of size
\((M-4)-(m-k)+h+1\) in an ordered set of size \(M-3\) by the position \(E+1\)
of their \((h+1)\)-st element. Consequently each side's tail sums to

\[
\sum_{E\ge0}T_E
=\frac52\binom{M-3}{m-3}+2\binom{M-3}{m-4}.
\]

**PROVED.** Add both tails, both block endpoints and the three interior block
terms, and use Pascal's identity. The resulting coefficient vector is
\((9,30,16)/4\), proving the first formula. Cancelling factorials proves the
rational form. This is a derivation of the previously reconstructed formula;
it does not rely on a fit. \(\square\)

## Exact failure region for this family

**PROVED.** Define

\[
H_m(N)=-8N^3+(13m+59)N^2+(12m^2-127m-93)N
       -5m^3-18m^2+245m-30.
\]

Then

\[
F_4-L=\frac{mH_m(N)}{4(N-2)(N-3)(N-4)(N-m-1)}.
\]

**PROVED.** For \(t\ge0\),

\[
H_m(2m+6+t)=(7m+3)(m^2-33m-64)
 -(32m^2+311m+249)t-(35m+85)t^2-8t^3.
\]

Thus \(H_m\) strictly decreases in \(N\ge2m+6\). For \(4\le m\le34\),
\(m^2-33m-64=(m-34)(m+1)-30<0\). For \(m=35+r\), it equals
\(r^2+37r+6>0\).

**PROVED.** Moreover \(H_m(N)<0\) for \(m\ge2\) and \(N\ge5m/2\). Substitute
\(N=5m/2+t\):

\[
H_m=-8t^3+(59-47m)t^2+(-73m^2+168m-93)t
 -\frac{75m^3-133m^2-50m+120}{4}.
\]

All coefficients are negative for \(m\ge2\): set \(m=2+r\); the last
numerator becomes \(75r^3+317r^2+318r+88\), and the linear coefficient becomes
\(-73r^2-124r-49\).

**PROVED  -  exact classification, family scope.** For \(4\le m\le34\), (S)
holds for every member of this interior family at every \(N\ge2m+6\). For
\(m\ge35\), let \(R_m\) be the unique real root of \(H_m\) in
\((2m+6,5m/2)\). The family violates strict (S) precisely at integer lengths

\[
2m+6\le N\le R_m,
\]

equivalently precisely when \(H_m(N)\ge0\). Equality also violates (S), whose
comparison is strict. Throughout this failure interval \(L>1\): the larger
root of its defining quadratic lies above \(5m/2+1\), as shown below. Outside
the interval, \(F_4<L\) if \(L\ge1\), and \(F_4<L<1\) if \(L<1\).
This is an exact characterization **of this family**, not a classification of
all possible centres.

**PROVED.** This interval grows linearly: \(R_m/m\) tends to the unique root
in \((2,5/2)\) of \(8x^3-13x^2-12x+5=0\). Divide \(H_m(mx)\) by \(m^3\);
its coefficients converge, its leading polynomial is strictly decreasing on
that interval before changing sign, and the bounded roots have that unique
limit. No fitted intercept is needed.

**VERIFIED-COMPUTATIONALLY.** Examples of integer failure intervals for this
interior family are \(m=35:[76,76]\), \(m=50:[106,108]\),
\(m=100:[206,217]\), and \(m=200:[406,435]\). At the first interior onset,

\[
F_4(35,76)=\frac{349741}{259296},\quad L(35,76)=\frac{399}{296},\quad
F_4-L=\frac{217}{259296}>0.
\]

## What survives for the threshold

**PROVED  -  family separation.** Let
\(Q(m,N)=N^2-3(m+1)N+m^2+5m+2\). Since
\(Q(m,5m/2+1)=-m(m+2)/4<0\), its larger root \(A(m)\) exceeds
\(5m/2+1\). Therefore

\[
N^*(m)-1=\lceil A(m)\rceil-1>5m/2.
\]

Every interior-family failure just proved lies strictly below the cell read
by the selector. For this family, (W) holds whenever the family exists and
\(N\ge N^*\), because \(F_4<L\le1\).

**PROVED  -  stronger family bound.** In fact \(F_4<1\) for \(m\ge4\) and
\(N\ge\max(2m+6,5m/2)\). To check it without approximation, subtract the
numerator of \(F_4\) from its positive denominator. At \(m=4+r\),
\(N=5m/2+t\), this polynomial is

\[
\begin{split}
4t^4&+(36r+104)t^3+(111r^2+633r+896)t^2\\
&+(118r^3+990r^2+2765r+2572)t\\
&+\frac{25}{2}r^4+\frac{487}{4}r^3+\frac{923}{2}r^2+818r+576>0.
\end{split}
\]

For \(m\ge8\), \(Q(m,2m+6)=-m^2+5m+20<0\) also guarantees that the family
exists at \(N=N^*-1\). Thus its load is strictly below 1 even at the
selector's premise cell. This closes the relevant bound on this family only.

**PROVED  -  dependency repair.** The correction note's proposed single-cell
identity \(P:V(m,k-1)=L(m,k-1)\), \(k=N^*\), does not alone prove a threshold
that passes at every subsequent length. A sufficient additional hypothesis is

\[
T:\quad V(m,N)<1+1/N\quad\hbox{for every }N\ge k.
\]

Together with the explicit lower exclusion of earlier lengths in
`threshold_audit/audit.md`, these imply the selector
\(n_t=k+[k(L(m,k-1)-1)\ge L(m,k-1)]\). Universal (W) supplies only the
\(u\ge2\) part of a sufficient tail bound. A universal upper bound for
\(u=1\) is still needed; the five-cell enumeration is not a proof. If only
the first passing cell is meant, rather than a persistent threshold, an upper
bound at \(N=k\) suffices in place of the whole tail.

**PROVED.** P is stronger than necessary and is not equivalent to end-pair
maximality among \(u\ge2\): those centres can exceed the end pair while
remaining below \(L\). A failure of P also need not change the selector if
the true maximum still passes its strict rescaling criterion. These are
logical dependency statements, not alleged additional Hamming witnesses.

**VERIFIED-COMPUTATIONALLY.** At \(m=19\), the exhaustive tests give P at
\(N=51\) and (W) at \(N=52\). The selector test is
\(52(L-1)=52/1519<1520/1519=L\), so that cell supports \(n_t=52\).
The (S) refutation at \(N=41\) does not spoil these tests. No uncomputed
infinite tail is inferred from the two finite cells.

## Controls, files and reruns

**VERIFIED-COMPUTATIONALLY.** All accepted claims have the following executed
checks; expected failures were actually observed:

| control | outcome |
|---|---|
| Count formula versus actual sphere words, all \(3^7\) centres, every pattern, \(m=2\) | 17,496 exact comparisons passed |
| All per-\(u\) DP maxima versus that brute-force experiment | all eight strata passed |
| 2/3 symmetry, all \(4^5\) centres, every pattern | 6,144 checks passed |
| Double-count excluded free letters in the formula | **failed as required**, 1,342 disagreements; process exit 1 |
| Interior closed form versus direct summary counts | 588 checks passed |
| Independent coordinate coefficient DP versus released counter on six stated cells | all exact, including \(m=34,35\) onset pair |
| Assert (S) on the three audited failing centres | all three assertions **failed as required** |
| Assert the \(m=35\) interior onset already holds at \(m=34\) | **failed as required** |
| Apply end-pair plateau identity below its scope, \(m=3,4\) | both equalities **failed as required** |
| Extend the clipped one-star lower-exclusion result to \(m=2\) | **failed as required** |
| Treat equality at the rescaling boundary as passing | **failed as required** |
| Three false logical implications concerning P | all three **failed as required**, on explicitly abstract logical assignments |
| Change coefficient 30 to 29 in the binomial proof | polynomial identity **failed as required** |
| Fraction polynomial identities and sign expansions | all passed |
| Remove the length hypothesis from the auxiliary adjacent-swap lemma | **failed as required**, load decreases from \(4/3\) to \(1\) |
| Extend that adjacent swap across a free letter, at a (W) cell with \(u=2\) | **failed as required**, \(71/770\) decreases to \(83/924\) |
| Apply the middle-four formula with a short binary run | **failed as required**, actual \(2549/10296\) differs from \(61/198\) |

**VERIFIED-COMPUTATIONALLY.** The independent audit additionally checked 75
binary lower-exclusion cells, 11 clipped single-star formula cells, 27,239
earlier-length exclusions, and all eight supplied selector arithmetic controls.
Its initial run exposed an incorrectly positioned one-free-letter positive
control; the corrected constructor passes. The failed run is retained and no
mathematical claim relies on it. Its eleven negative controls were rerun from
the main context, all failing as required.

**PROVED  -  auxiliary result.** `w_proof/PARTIAL_PROOF.md` proves that an
adjacent interior swap \(01\to10\) cannot decrease \(\Lambda\) for
\(N\ge2m+3\), for any \(u\). A generating-function proof uses a binomial
midpoint inequality. The result does not permit sorting across a free letter,
as the stated failing control demonstrates, so it does not prove (W).

**VERIFIED-COMPUTATIONALLY.** Its verifier checked 14,094 such swaps, 972
independent polynomial delta identities and 48 direct sphere comparisons.
It also checked the middle-four formula independently at 15 cells and 75
block-local contributions. All three auxiliary negative controls above failed
as required. See `w_proof/verification.json` and `w_proof/root_rerun.log`.

**VERIFIED-COMPUTATIONALLY.** Reproduce from this directory:

```sh
PYTHONDONTWRITEBYTECODE=1 /usr/local/bin/python3 verify_j2.py verify
PYTHONDONTWRITEBYTECODE=1 /usr/local/bin/python3 verify_j2.py mutant-control
PYTHONDONTWRITEBYTECODE=1 /usr/local/bin/python3 proof_checks.py
PYTHONDONTWRITEBYTECODE=1 /usr/local/bin/python3 verify_j2.py cell 41 19
PYTHONDONTWRITEBYTECODE=1 /usr/local/bin/python3 verify_j2.py cell 51 19
PYTHONDONTWRITEBYTECODE=1 /usr/local/bin/python3 verify_j2.py cell 52 19
PYTHONDONTWRITEBYTECODE=1 /usr/local/bin/python3 verify_j2.py cell 76 35 --u 4
PYTHONDONTWRITEBYTECODE=1 /usr/local/bin/python3 summarize.py
PYTHONDONTWRITEBYTECODE=1 /usr/local/bin/python3 threshold_audit/verify_dependencies.py
```

**VERIFIED-COMPUTATIONALLY.** `mutant-control` deliberately exits 1. The other
commands exit 0. `summarize.py` reads the saved exhaustive JSONL logs; regenerate
them with the corresponding `cell` command and a redirect to its existing log
path if needed. The main verifier uses only the standard library. The
independent cross-checks read `code_release/privsphere.py`, SHA-256
`b637fed0c627dea8258670059a1d40465845510d5a802669c28abf4f46281eb3`, with bytecode
disabled. Everything produced is under the requested `j2` directory. No
repository file was edited, committed or pushed; no cluster job was submitted.

**REFUTED:** universal (S), and a universal one-length failure band.
**PROVED:** the interior four-free-letter formula, its exact failure criterion,
its separation from the selector, and (W) on that family.
**CONJECTURE  -  UNSETTLED:** universal (W), universal single-free-letter
optimality, and the unconditional all-budget threshold.
