# Exact budget-9 closure certificates

**PROVED (computer-assisted).** For every guide length `n >= 20`, every run count and every integer budget
`0 <= m <= 9`, the unrestricted covering number is `K = r`.

The filenames use word length `N = n - 1`. Each exact measure is supported in
the target union and has mass strictly greater than `r - 1`. The exact maximum
radius-9 load is one in every row.

| File | Guide length n | Runs r | Exact mass | Centre classes |
| --- | ---: | ---: | --- | ---: |
| `packing_N20_r21_m9.json` | 21 | 21 | 25728861133195119900/1285411466839188713 | 3486784401 |
| `packing_N21_r22_m9.json` | 22 | 22 | 1244606387850518652960/59159674329805194787 | 10460353203 |
| `packing_N22_r23_m9.json` | 23 | 23 | 222870000000048450/10022082107013869 | 31381059609 |
| `packing_N23_r24_m9.json` | 24 | 24 | 11174042880000000000000/479696779442938946033 | 94143178827 |
| `packing_N24_r25_m9.json` | 25 | 25 | 18623404799997765191424/771843880811584148951 | 282429536481 |

Run from the release root with Python 3.10 or later and a C++17 compiler:

```sh
python3 -I certificates/m9_closure/verify.py verify
```

The Python code uses only the standard library. The verifier locates every input
relative to its own file and needs no repository files outside this folder. It
compiles the included C++ sources with `c++ -O3 -std=c++17 -pthread`; `CXX` may
name another compatible compiler. All generated binaries and working inputs live
in a temporary subdirectory of this folder, removed when the command finishes.
No precompiled binary is required.

The checker recomputes shell sizes, support and exact total mass, then maximizes
the load over all ternary centre classes using an exact dynamic program. Swapping
the two unused letters at each coordinate proves that this is also the maximum
over all four-letter centres. It compares the computed maximum with the stored
normalizer and rejects planted overload and outside-support forgeries. It prints
all results and exits zero only if every documented certificate and both controls
pass. These are integral covering certificates; no fractional optimum is claimed.

The command also runs the supplied small-space literal enumerations for both
private and tagged shells, including four-letter symmetry and large coefficients,
and checks each maximum witness with a separate Python convolution.
The existing full-mass packing at guide length 20, the five diagonal certificates,
step-down, padding and the sufficient-length bound from guide length 26 give the
universal budget-nine case split. Budget monotonicity supplies smaller budgets.

For integrity, run from this directory:

```sh
shasum -a 256 -c SHA256SUMS
```
