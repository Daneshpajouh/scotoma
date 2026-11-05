# Exact budget-10 closure certificates

**PROVED (computer-assisted).** At `m = 10`, `K = r` for every run count at every `n >= 24`, and for
`r <= 16` at `n = 23`. The open cells are `12 <= r <= n` at `n = 20,21,22`,
and `17 <= r <= 23` at `n = 23`. The previously proved threshold handles `n >= 28`.

The filenames use word length `N = n - 1`. Each exact measure is supported in
the target union and has mass strictly greater than `r - 1`. The exact maximum
radius-10 load is one in every row.

| File | Guide length n | Runs r | Exact mass | Centre classes |
| --- | ---: | ---: | --- | ---: |
| `packing_N22_r12_m10_d2_c2_all.json` | 23 | 12 | 113692351711201623360/10082005511554229881 | 31381059609 |
| `packing_N22_r13_m10_d2_c2_all.json` | 23 | 13 | 122943944145541559040/10129313980325173261 | 31381059609 |
| `packing_N22_r14_m10_d2_c2_all.json` | 23 | 14 | 132160725848045774880/10142931530083131979 | 31381059609 |
| `packing_N22_r15_m10_d2_c2_all.json` | 23 | 15 | 19108261846499760/1358250091297951 | 31381059609 |
| `packing_N22_r16_m10_d2_c2_all.json` | 23 | 16 | 194601605895496800/12803107753430537 | 31381059609 |
| `packing_N23_r24_m10_d2_c1_all.json` | 24 | 24 | 7425758183451928320/322698929582077393 | 94143178827 |
| `packing_N24_r25_m10_d2_c1_all.json` | 25 | 25 | 762537851/31250000 | 282429536481 |
| `packing_N25_r26_m10_d2_c0_all.json` | 26 | 26 | 25831142183/1000000002 | 847288609443 |
| `packing_N26_r27_m10_d2_c0_all.json` | 27 | 27 | 895134238309190880/33969511411405147 | 2541865828329 |

Run from the release root with Python 3.10 or later and a C++17 compiler:

```sh
python3 -I certificates/m10_closure/verify.py verify
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

Step-down from the four diagonal files handles all run counts at lengths 24-27,
and step-down from `(n,r)=(23,16)` handles smaller run counts at that length. The
four smaller length-23 files are direct certificates of already implied cells.
The older integer certificate for `(n,r,m)=(20,11,10)` is a separate package.

For integrity, run from this directory:

```sh
shasum -a 256 -c SHA256SUMS
```
