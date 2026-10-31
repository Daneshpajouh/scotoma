# Exact integer closure at N=19, r=11, m=10

**PROVED (computer-assisted).** K(19,11,10)=11 for unrestricted four-letter
centres. The eleven single-deletion pattern balls form a minimum cover of
their radius-ten union U. The appendix theorem `thm:m10r11` states the proof.
The matching fractional optimum is strictly between 9 and 10, so its ceiling
is 10. The integer case split proves the extra unit beyond fractional rounding.
The cells r=12 through 20 at N=19, m=10 remain open; K(19,12,10) is in {11,12}.

Run from the release root with Python 3.11 or later and NumPy:

```sh
python3 -I certificates/m10_r11/verify.py verify
```

Or run from this folder:

```sh
python3 -I verify.py verify
python3 -I verify_kernel_certs.py
shasum -a 256 -c SHA256SUMS
```

Both verifiers resolve inputs relative to their own files and write no files.
They require no inputs outside this folder and do not run an optimization
solver. The full verifier uses one numerical-library thread. The second
verifier independently accumulates every certificate product in Python integers.
A successful full run ends with `M10_R11_VERIFIED K(19,11,10)=11
ceil(tau_star)=10`; the second ends with `KERNEL_CERTS_OK`.

## Exact checks

The informative prefix has length 10, and the constant block has length 9.
Swaps of the two unused informative letters, permutations of constant
coordinates and independent permutations of their three nonzero letters give
3^10 * 10 = 590,490 centre classes. The kernel averages incidence over a
packing orbit. A concrete cover necessarily satisfies every orbit inequality
by the union bound, with integer class multiplicities and no symmetry
assumption on that cover.

`cell_19_11_10.certificate.json` contains the 420-orbit packing and the
412-class fractional cover, both uniform within each orbit. It uses the
staircase orientation 0^j 1^(10-j) followed by nine zeros. `threshold.json`
records the exact packing mass v and the eligibility threshold v-9. Every
centre has packing load at most one. Thus any cover of at most ten balls can
use only the 12,052 classes whose load is at least v-9; 578,438 are excluded.

`verify.py` reconstructs the kernel by support-intersection counts, checks
support and exact packing mass, and scans all classes with rigorous integer
intervals and exact rational fallbacks. A proved accumulation bound guards
every vectorized signed-integer product against overflow. It checks the
entire supplied kernel and reduced matrices after an explicit column
permutation, then verifies all five integer certificates. It separately
checks the matching fractional cover on all 351,740 orbits of U. Small
literal four-letter enumerations calibrate the formula; a wrong-radius
kernel and an inflated dual are rejected as negative controls.

`filter_instance.txt` stores the original interval and kernel table input;
`filter_out_candidates.txt` and `elig_index.npy` fix the eligible-column
order. `pack_matrix.npy` gives all 5,061,840 kernel numerators over
1,693,052,928. `A_reduced.npy` and `rhs_reduced.npy` divide each row and its
right hand side by their gcd. The reconstruction checks every exact entry.

## Five integer certificates

The four named columns all have constant-block weight zero:

| Column | Informative prefix |
|---:|---|
| 0 | 0000000000 |
| 7074 | 1111111111 |
| 1513 | 0000000001 |
| 2065 | 0111111111 |

The files `cert_<name>_z.npy` are indexed by `cert_index.json`. Write A_j
for column j and D=10^12. Every z is a nonnegative integer vector and every
usable column product is at most D. A banned column is omitted only in a
branch that forces its multiplicity to zero.

| Certificate name | Subtracted columns | Banned columns | Signed residual dot | Strict threshold |
|---|---|---|---:|---:|
| omit_both_ends | none | 0, 7074 | 10020000594807 | 10 D |
| end0_omit_1513 | 0 | 1513 | 9020000146095 | 9 D |
| end0_and_1513 | 0, 1513 | none | 8020001389392 | 8 D |
| end7074_omit_2065 | 7074 | 2065 | 9020000146095 | 9 D |
| end7074_and_2065 | 7074, 2065 | none | 8020001389392 | 8 D |

The branches are: neither end occurs; end 0 occurs with next 1513 absent or
present; or end 0 is absent and end 7074 occurs with next 2065 absent or
present. Subtract the indicated one or two unit columns. The remaining
integer sum exceeds 10, 9 or 8 respectively, so the total sum is at least
11 in every case. Subtracted columns remain capped if used again; arbitrary
multiplicities are permitted. The five branches exhaust every nonnegative
integer vector. The eleven library balls give the upper bound by the
definition of U.

## Independent checking and integrity

The original coordinate-polynomial implementation and an independently
written support-intersection checker agree on every exact value used here,
including all 5,061,840 retained matrix entries and all five certificates.
The full standalone verifier packages the latter reconstruction; the
separate Python-integer verifier checks the original matrix and vectors.
All original certificate vectors, the index, threshold, eligible list and
matrices are copied byte for byte. SHA256SUMS binds every input and checker.
