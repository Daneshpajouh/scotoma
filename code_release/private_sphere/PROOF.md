# (W) is proved by an exact finite certificate

**PROVED, with the finite rational inequalities verified computationally.** For the private-sphere weighting defined in `LEMMA_2026-10-23_the_straddling_load_is_two_terms.md`,

\[
\boxed{\quad u(c)\ge2,\quad N\ge N^*(m)\quad\Longrightarrow\quad \Lambda(c)<1.\quad}
\]

A stronger, shape-independent statement is established:

\[
 m\ge2,\quad \frac m{N-2}\le\frac25,\quad u(c)\ge2
 \quad\Longrightarrow\quad \Lambda(c)<\frac{99}{100}.
\tag{1}
\]

No assertion identifying the maximizing centre is used. The proof combines elementary sampling identities, a finite directed grid, and a geometric tail. The independent verifier uses only the Python standard library and `fractions.Fraction`; no optimizer, repository import, or floating-point decision enters verification.

**VERIFIED-COMPUTATIONALLY.** The largest finite certificate bound is
`98389493/100000000`, at `u=4`. This is an upper bound from the certificate, not the exact maximum load. All 34 strata `u=2,...,35` clear `99/100`. The bound for all `u>=36` is proved analytically below.

## 1. The local event and its positive polynomial bounds  -  PROVED

Use the lemma's orientation \(P_j=1^j0^{N-j}\), \(j=0,\ldots,N\). Put
\(A_j=[N]\setminus\{j,j+1\}\), ignoring nonexistent boundary coordinates, and
\(n_A=|A_j|\). A random point of \(S_j\) chooses a uniform \(m\)-subset \(T\subseteq A_j\) and independently overwrites its coordinates by 2 or 3, each with probability \(1/2\).

Fix a centre. Let \(u\) be its total number of free letters. At this boundary let

* \(s\) be the number of excluded free coordinates; \(a=u-s\) are allowed;
* \(e\) be the number of allowed binary disagreements;
* \(f\) be the number of excluded binary disagreements.

Necessarily \(0\le s\le2\) and \(0\le f\le2-s\). If \(t\) allowed free coordinates and \(b\) allowed binary disagreements belong to \(T\), then the conditional acceptance probability is

\[
 g_{u,e,f}(t,b)=2^{-t}\sum_{h=0}^{\min(t,\,2t-u-e+b-f)}\binom th,
\tag{2}
\]

with an empty sum interpreted as zero. Indeed, if \(h\) overwritten free coordinates disagree with the centre, subtracting \(m\) from the Hamming distance gives
\(u+e+f-2t-b+h\). Equation (2) is exactly its nonpositive event.

For nonnegative rational coefficients \(q_{kl}\), define

\[
 Q(t,b)=\sum_{k=0}^{a}\sum_{l=0}^{e}
 q_{kl}\frac{\binom tk\binom bl}{\binom ak\binom el}.
\tag{3}
\]

If \(Q(t,b)\ge g_{u,e,f}(t,b)\) for every integer \(0\le t\le a,0\le b\le e\), then

\[
 w_j(c)\le \sum_{k,l}q_{kl}\frac{(m)_{k+l}}{(n_A)_{k+l}}
 \le \sum_{k,l}q_{kl}p^{k+l},\qquad p=2/5,
\tag{4}
\]

provided \(m/n_A\le p\). The first equality for the expectation of each normalized monomial follows by counting specified selected coordinates. The second inequality follows from
\((m-i)/(n_A-i)\le m/n_A\le p\); moments of order exceeding \(m\) are zero. Thus (4) is valid for sampling without replacement, without asserting a Bernoulli tail domination.

The certificate supplies either a polynomial (3) or the following exponential bound for each local configuration. Every coefficient, every inequality in (3), and the resulting bound (4) is checked exactly.

## 2. Exponential bounds  -  PROVED

Let \(t_0\ge1\) be rational and set

\[
 A(t_0)=1-p+\frac p2(t_0+t_0^2),\qquad B(t_0)=1-p+pt_0.
\]

A selected free coordinate contributes a gain of 1 or 2, equiprobably; a selected binary disagreement contributes gain 1. The event in (2) requires total gain at least \(u+e+f\). Therefore

\[
 w_j(c)\le t_0^{-u-e-f}A(t_0)^{u-s}B(t_0)^e.
\tag{5}
\]

To justify the product bound under sampling without replacement, expand its generating product in the selection indicators. All coefficients are nonnegative because \(t_0\ge1\), so each degree-\(d\) moment is at most \(p^d\), exactly as in (4). Markov's inequality then gives (5).

Write \(E_j=e+f\) for the total binary disagreement count. Since \(A,B\ge1\), (5) implies

\[
 w_j(c)\le \left(\frac{A(t_0)}{t_0}\right)^u
             \left(\frac{B(t_0)}{t_0}\right)^{E_j}.
\tag{6}
\]

## 3. Walk geometry and the infinite tail  -  PROVED

At cut \(j\), let \(z_j\) count zeros in the prefix of the centre and \(o_j\) count ones in its suffix. Then
\(E_j=z_j+o_j\). Passing a 0 increases \(z\); passing a 1 decreases \(o\); passing a free letter leaves both unchanged.

For cuts \(i<j\),

\[
 E_i+E_j\ge (j-i)-\#\{\text{free coordinates between the cuts}\}.
\]

Every binary coordinate between the cuts disagrees with at least one of the two staircase patterns. Consequently the number of cuts with \(E_j\le h\) is at most

\[
 C_h\le 2h+u+1.
\tag{7}
\]

For any \(0<q<1\), Abel summation gives

\[
 \sum_jq^{E_j}=(1-q)\sum_{h\ge0}q^h C_h
 \le u+1+\frac{2q}{1-q}.
\tag{8}
\]

Choose \(t_0=2\). Then \(A/t_0=9/10\), \(B/t_0=7/10\), so

\[
 \Lambda(c)\le H(u):=(u+17/3)(9/10)^u.
\tag{9}
\]

Exactly,

\[
 H(36)=\frac{7509466514979724803946715958257547}
 {8000000000000000000000000000000000}<99/100.
\]

Moreover \(H(u+1)/H(u)=9(3u+20)/[10(3u+17)]<1\) for every \(u\ge36\), because \(3u>10\). This proves (1) for every \(u\ge36\).

For the finite strata take \(K=32\). Using (6) with \(t_0=4\), and the weaker level-count consequence \(\#\{j:E_j=h\}\le2h+u+1\) of (7), the entire omitted tail is bounded by

\[
 T_u=\left(\frac{23}{20}\right)^u\left(\frac{11}{20}\right)^K
 \left[\frac{2K+u+1}{1-11/20}+
       \frac{2(11/20)}{(1-11/20)^2}\right].
\tag{10}
\]

This bound holds at every length, including arbitrarily long centres.

## 4. A finite grid covers every centre  -  PROVED

For each fixed \(u=2,\ldots,35\), the certificate provides upper bounds
\(F_u(s,e,f)\) for \(0\le e\le K\), obtained by (3) to (5). At a cut with previous letter \(a_0\) and next letter \(b_0\), use sentinel 1 on the left and sentinel 0 on the right, and put

\[
 s=[a_0=X]+[b_0=X],\qquad f=[a_0=0]+[b_0=1],\qquad e=z+o-f.
\tag{11}
\]

A state is \((z,o,x,a_0)\): zero prefix count, one suffix count, remaining free letters, previous letter. Start at \((0,K,u,1)\). The allowed transitions append

* 0: \((z,o,x)\mapsto(z+1,o,x)\);
* 1: \((z,o,x)\mapsto(z,o-1,x)\);
* \(X\): \((z,o,x)\mapsto(z,o,x-1)\).

Finish at \((K,0,0,\cdot)\). A cut earns reward \(F_u(s,z+o-f,f)\) when \(z+o<K\), and zero otherwise. The graph is acyclic; every path has exactly \(2K+u\) real letters. Its maximum reward \(B_u\) is computed by exact finite dynamic programming.

Here is why this finite graph bounds all lengths. If a centre has active cuts \(E_j<K\), take its first and last active cuts. Throughout the segment between them, \(0\le z,o<K\), by monotonicity of \(z,o\). Preserve that segment in the grid. Extend its prefix from \((0,K)\) by zeros and then ones, so the last added prefix letter is 1; extend its suffix by zeros and then ones to \((K,0)\), so the first added suffix letter is 0. At a first active cut, the original previous letter is necessarily 1, or the left sentinel; at a last active cut the next letter is necessarily 0, or the right sentinel. Thus both boundary rewards are preserved. Put all free letters outside the preserved segment at the initial level \(E=K\), where their rewards are zero. The result is a grid path with all the original active rewards and possibly additional nonnegative rewards.

If there are no active cuts the active sum is zero. Thus, in every case,

\[
 \Lambda(c)\le B_u+T_u.
\tag{12}
\]

**VERIFIED-COMPUTATIONALLY.** `verify_w.py` checks all local inequalities, recomputes \(B_u\) by a backward recursion independent of the producer's forward recursion, checks (10), and verifies \(B_u+T_u<99/100\) in every stratum. This is a finite verification of a universal bound, not enumeration of finitely many original \((m,N)\) cells.

## 5. The parameter range and small budgets  -  PROVED

For \(m\ge6\), at \(n_A=N-2=5m/2\),

\[
 L=2/5+2m/(3m+2)\ge1.
\]

Because \(L\) strictly decreases with \(N\), \(L\le1\) implies \(m/(N-2)\le2/5\). For \(m=2,3,4,5\), the exact values \(N^*=7,10,12,15\) give the same ratio bound. Endpoint spheres have \(n_A=N-1\), so they also satisfy it. Equations (1) to (12) therefore prove (W) for every \(m\ge2\).

For \(m=1\), \(N^*=4\). If \(u>2m\), no sphere word can be within distance \(m\) of the centre: from the local event, a contributing selection requires \(u\le2t\le2m\). For \(u=2\), its sole selected coordinate must be free and match the centre; all binary coordinates must agree with the pattern. The contributing cuts form the zero-disagreement plateau in the binary subsequence. It contains zero, one, or two free coordinates; the sums of the allowed free counts at those cuts are respectively 2, \(1+1\), or \(1+0+1\), or smaller at an endpoint. Therefore

\[
 \Lambda(c)\le 1/(N-2)\le1/2<1.
\]

At \(m=0\), every centre with a free letter has zero load. This closes all budgets.

## 6. Companion lemma for one free letter  -  PROVED

Finite enumeration at five cells alone would not establish the all-length floor consequence. The following argument supplies the missing general upper bound.

Let \(N\ge2m+2\), \(m\ge2\), and put

\[
 \rho=\frac m{N-2},\quad n=N-3,\quad k=m-1,\quad q=k/n,
 \quad h_d=\frac{(k)_d}{(n)_d},\quad h_0=1,\quad h_d=0\ (d>k).
\]

Assume

\[
 F(N,m):=N^2-(3m+2)N+m^2+6m-1\ge0.
\tag{13}
\]

A centre with one free letter contributes only when that coordinate is allowed and selected. Delete it. The two cuts adjacent to it contribute zero; the other cuts inject into the cuts of the resulting binary word with the same boundary types and total binary disagreement. Use a common sphere size parameter \(N-2\). This is exact at interior cuts and dominates endpoint contributions: an endpoint has one extra allowed matching binary coordinate, and deleting that coordinate increases the probability of any increasing event in the selected disagreements. If a row is nonzero, its binary disagreements are at most \(m\), so such a matching coordinate exists.

Conditional on selection of the free coordinate, the load divided by \(\rho\) at a descent of height \(d\) is

\[
 A_d=h_d+\frac d2(h_{d-1}-h_d),\qquad A_0=1;
\]

at a monotone cut of height \(d\ge1\) it is \(h_{d-1}/2\); at an ascent it is zero. These follow directly from (2): a matching free overwrite permits one uncovered binary disagreement; a nonmatching one permits none.

Write \(D\) for the binary word's number of descents, including sentinel boundaries. Every descent height is at least \(D-1\), as in the established binary proof. Assign each monotone stretch to its neighbouring descent; the two arms incident to a descent of height \(d\) contribute at most \(\sum_{i\ge d}h_i\). Consequently the divided load is at most

\[
 \sum_{i=1}^{D}V(d_i),\qquad
 V(d)=A_d+\sum_{i\ge d}h_i.
\tag{14}
\]

The hockey stick gives
\(T:=\sum_{i\ge0}h_i=(n+1)/(n-k+1)\), and hence

\[
 V(0)=1+T,\quad V(1)=T-1/2+q/2,\quad V(2)=T-1.
\]

For \(d\ge2\), both \(h_d\) and its tails decrease by a factor at most \(q\); so does \(h_{d-1}-h_d\), since its next-term ratio is
\((k-d+1)/(n-d)\le k/n\) when \(n\ge2k\). It follows that
\(V(d+1)\le q(d+1)V(d)/d\). The sequence \(V\) is nonincreasing. Furthermore \(3V(2)\le2V(1)\) because \(T<2\), and \(D V(D-1)\) decreases for \(D\ge3\), its next-term ratio being at most \(q(D+1)/(D-1)\le2q<1\). Finally,

\[
 V(0)-2V(1)=2-T-q
 =\frac{F(N,m)}{(N-m-1)(N-3)}\ge0.
\]

Thus (14) is at most \(V(0)\) for every \(D\), and

\[
 \boxed{\ \Lambda(c)\le\rho V(0)
 =\frac m{N-2}+\frac m{N-m-1}=L(m,N).\ }
\tag{15}
\]

For \(N\ge N^*\), condition (13) follows from
\(F=Q+N+m-3\), where
\(Q=N^2-3(m+1)N+m^2+5m+2\ge0\). For \(m=1,N\ge4\), the same calculation has \(k=0\), \(V(0)=2,V(1)=1/2,V(d)=0\) for \(d\ge2\); hence \(\Lambda\le2/(N-2)=L\le1\).

## 7. Consequences and boundaries  -  PROVED

Combining (W), (15), and the established binary theorem shows that **the full-chain private-sphere weighting is feasible without rescaling for every \(N\ge N^*(m)\)**. Each sphere carries mass one, so its mass is \(N+1\). Restriction to consecutive private regions, using E4, gives

\[
 \boxed{\quad N\ge N^*(m)\quad\Longrightarrow\quad
 \tau^*(N,r,m)=K(N,r,m)=r\quad (1\le r\le N+1).\quad}
\]

For \(m=1,N=4\), the binary bound is still \(D/2^{D-1}\le1\), sufficient for feasibility, although its strict uniqueness clause requires the stronger range in the established theorem.

This is a sufficient-length theorem. It does not identify the true \(M(N,r)\), does not restore the refuted all-length law, and does not assert the end pair maximizes \(\Lambda\). The (S) counterexamples remain valid. Sound-pattern minimality and arbitrary-ball covering remain distinct objects.

The distinction between upper and lower bounds is retained: feasibility supplies a **lower** bound \(\tau^*\ge r\), while the library balls give the matching upper bound. In the rescaled certificate below, failure only concerns this particular certificate; it is not an upper bound on \(K\).

## 8. Reproduction and controls  -  VERIFIED-COMPUTATIONALLY

From this directory:

```
/usr/local/bin/python3 -B verify_w.py --controls
/usr/local/bin/python3 -B verify_w.py --control zero-local
```

The first command must exit 0. The second must exit 1. The other planted-invalid controls are `empty-polynomial`, `missing-local`, `lower-path`, `zero-tail`, and `wrong-p`; every one is required to be rejected. `verification_result.json` and `verification.log` record the actual outcomes. The certificate search used a numerical optimizer only to propose nonnegative polynomial coefficients; every retained coefficient, local bound, path decision, and acceptance test is rechecked exactly by the standalone verifier.

**VERIFIED-COMPUTATIONALLY, initial stress tests.** At \((m,N)=(35,76)\), \(1^{36}2^4 0^{36}\) has load \(349741/259296>1\), but \(N^*(35)=93\), outside (W). At \((35,93)\), \(1^{36}2^4 0^{53}\) has load \(289153/395694<1\). The corresponding \(N^*\) tests at \(m=19,100,500\) also lie below one; see `initial_tests.txt`. These are stress tests, not premises of the universal proof.

**PROVED/REFUTED scope note.** One-letter deletion and sorting binary letters across a free block were considered as candidate proof moves; neither is used. The latter fails, for example, at \(N=7,m=2\): `0021000` has load `7/20`, whereas swapping its binary inversion across the free letter gives `0120000`, with load `1/4`. See `explore_sort.txt`.

Hub: no (purely local exact mathematics and local verification; no external transmission).

## 9. The threshold selector survives, with its premise now proved  -  PROVED

Here \(V(m,N)=\max_c\Lambda(c)\), and \(n_t(m)\) denotes the least guide length in \(N\ge2m+3\) at which **this uniform private-sphere weighting, with its single global rescaling**, certifies \(K=r\) simultaneously for every \(r\le N+1\). This is a threshold of a particular certificate, not an assertion about the true covering number or \(M(N,r)\) below that threshold.

For every integer \(m\ge3\), the original selector is now unconditional for this certificate:

\[
 \boxed{\ n_t(m)=N^*(m)+
 \left[\,N^*(m)\bigl(L(m,N^*(m)-1)-1\bigr)
              \ge L(m,N^*(m)-1)\,\right].\ }
\tag{16}
\]

The route does not restore (S). It proves only the two regions the selector needs.

**PROVED: one free letter at \(N^*-1\).** Let

\[
 R=\frac{3m+3+\sqrt{5m^2-2m+1}}2,\qquad N^*=\lceil R\rceil.
\]

At \(N=R-1\), the polynomial in (13) is
\([5m-5-\sqrt{5m^2-2m+1}]/2>0\) for \(m\ge2\). Squaring is valid, and reduces the strict inequality to \(5m^2-12m+6>0\). Also \(F\) increases for \(N\ge R-1\), since \(F'(R-1)=\sqrt{5m^2-2m+1}-1>0\). Therefore (13) holds at \(N=N^*-1\). The companion lemma proves \(\Lambda\le L\) for \(u=1\). The explicit straddling centre in the two-term lemma attains \(L\) because \(N^*-1\ge2m+3\) for \(m\ge3\).

**PROVED: two or more free letters at \(N^*-1\), except finitely many small budgets.** For even \(m\ge8\),

\[
 Q(m,5m/2+2)=m(6-m)/4<0.
\]

For odd \(m\ge11\),

\[
 Q(m,5m/2+5/2)=(-m^2+10m+3)/4<0.
\]

Since these arguments are above the smaller root, the larger root exceeds them. Taking ceilings gives
\(m/(N^*-3)\le2/5\). Thus the universal bound (1) applies at \(N=N^*-1\) for \(m=8\) and every \(m\ge10\).

**VERIFIED-COMPUTATIONALLY: the remaining six budgets.** `finite_cells.py` exhaustively maximizes over all centre classes with \(u\ge2\), by a finite exact grid with the actual endpoint sphere denominators. All \(u>2m\) have zero load by the local event. For each remaining \(u\), it checks every binary composition and every ordering. Its two literal-enumeration controls agree with its DP; deleting the free-letter transition fails both controls. Maximizing words are also evaluated by the direct local sphere-count formula.

| \(m\) | \(N=N^*-1\) | exact maximum over \(u\ge2\) | \(L(m,N)\) |
|---:|---:|---:|---:|
| 3 | 9 | 429/560 | 36/35 |
| 4 | 11 | 1117/1260 | 10/9 |
| 5 | 14 | 5/6 | 25/24 |
| 6 | 16 | 19/21 | 23/21 |
| 7 | 19 | 161/187 | 196/187 |
| 9 | 24 | 135/154 | 81/77 |

Every entry in the third column is less than one. Together with the binary theorem and (15), this proves

\[
 \boxed{V(m,N^*(m)-1)=L(m,N^*(m)-1)\quad\text{for every }m\ge3.}
\tag{17}
\]

**PROVED: no earlier length in the stated range certifies by this rescaling.** The explicit straddling centre has load \(L(m,N)\) whenever \(N\ge2m+3\). Define \(G(N)=N(L(m,N)-1)\). Then

\[
 G'(N)=-1-\frac{2m}{(N-2)^2}
         -\frac{m(m+1)}{(N-m-1)^2}<0.
\]

At \(N=R-1\),

\[
 G(R-1)-1=
 \frac{(m+2)R-m^2-5m-4}{(R-3)(R-m-2)}>0,
\]

because \(R\ge m+3\), making the numerator at least 2. Since \(N^*-2<R-1\), every \(N\le N^*-2\) in \(N\ge2m+3\) has \(L\ge1+1/N\), so this particular rescaled mass cannot exceed \(N\). At \(N\ge N^*\), the unscaled weighting is feasible. At the remaining cell \(N=N^*-1\), (17) supplies the exact maximum.

A rescaled mass \((N+1)/V\) pins \(K=N+1\) precisely when it exceeds \(N\), namely \(V<1+1/N\). Substituting \(N=N^*-1\) and \(V=L\) gives exactly (16), including the strict inequality and its equality case. Restricting the weighting gives all smaller run counts as well.

**PROVED limitation.** Equation (16) is a closed form for the certification threshold just defined. It does not give the true all-length \(M(N,r)\), does not identify the minimum length where some other certificate or argument might prove the covering floor, and does not revive the refuted \(\lfloor(3n-11)/5\rfloor\) law. The arbitrary-ball and sound-pattern floors remain separate.
