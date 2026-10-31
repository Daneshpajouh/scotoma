# Exact integer families for the four-letter staircase

For every integer `h >= 0`:

| N | r | m | K |
|---|---:|---|---:|
| `19+h` | 4 | `13+h` | 4 |
| `19+2h` | 5 | `12+h` | 5 |
| `20+2h` | 6 | `12+h` | 6 |
| `21+2h` | 6 | `13+h` | 6 |
| `19+2h` | 7 | `11+h` | 7 |

The four-run family has `tau* < 3` at every `h`. The other four rows are
integer families; a strict fractional deficit is verified at their seeds.
Such a deficit does not follow from guard transport. In particular,
`tau*(41,5,23)=5` at `h=11` in the five-run row. The verifier explicitly
checks the mass-five packing and maximum load one at that counterexample.

From the release root, with Python 3.10 or later, NumPy and a C++17 compiler:

```sh
python3 -I certificates/infinite_families/verify.py verify
```

Every input is in this folder. The verifier resolves paths relative to itself,
checks the local SHA-256 manifest, prints its results, exits nonzero on failure,
and keeps all generated matrices, eligible lists and binaries in a temporary
directory that is removed on exit. It writes no files in the release. Native
searches run sequentially with at most four worker threads; numerical-library
threads are set to one. No optimizer is needed.

The replay checks every entry of the three finite seed kernels with an
independent support-overlap formula, regenerates all 3,325 singleton-feasible
triples at `(19,4,13)`, and recounts every retained pair witness by independent
coordinate-polynomial multiplication. Literal small-cube tests check the pair
minima. The 426 eligible classes at `(19,5,12)` and the 371 at `(20,6,12)`
are derived again from exact packing loads, then searched exhaustively with
repetition allowed. Both searches must return zero solutions.

The new packings certify `(18,7,10)`, `(19,7,11)` and `(21,6,13)`.
The latter two use constant-coordinate saturation: any hypothetical `r-1`
cover must carry all four letters at every constant coordinate. For `C=N-r+1`,
its total nonzero suffix weight is therefore at least `3C`. The exact caps give
`v+3C*beta > (r-1)*M`, with conservative positive margins `2171/12500` and
`3/5000`. A second polynomial kernel verifies the packing maxima. Reducing
the last adjusted cap by `1/50000` must fail. Overloaded packing and underweight
fractional-cover controls are checked by exact extrema.

The seven seeds for the appendix window are `(19,4,13)`, `(19,5,12)`,
`(20,6,12)`, `(19,11,10)`, `(18,7,10)`, `(19,7,11)` and `(21,6,13)`.
The eleven-run seed has its separate existing replay in
`certificates/m10_r11/`; this package does not repeat that larger certificate.
For any seed `(N0,r0,m0)`, any `2 <= r <= r0`, and `d=ceil((r-1)/3)`,
guard transport proves `K(N,r,m)=r` whenever
`N >= N0` and `0 <= m <= m0+floor((N-N0)/d)`.
The arithmetic guard is `1+3d > r-1`. General transports require a guard set
inside `B_s(0^t)` that no `r-1` radius-`s-1` balls cover; the sufficient volume
test is `V_4(t,s) > (r-1)*V_4(t,s-1)`. These are sufficient conditions.

The overlap inequality is elementary inclusion-exclusion, with a forest
extension proved pointwise in the appendix. Its one-pair case is obtainable
from one Sherali-Adams multiplication round in the concrete binary cover
formulation. No novelty is claimed for that inequality or the hierarchy.
The contribution is their application to symmetry classes, the saturation
condition and the transport step for these staircase libraries.
`K(19,12,10)` remains open in `{11,12}`.

For file integrity, run from this folder:

```sh
shasum -a 256 -c SHA256SUMS
```
