# The floor at every length: checkers for the weighting argument

These files check the main-text theorem "The floor at every length, by a weighting rather than a packing" and the appendix subsection "The weighting at every length" (its proposition, the tail table, the threshold table and the table of eight (m, n) pairs). The m <= 5 case of that theorem is checked separately, by the scripts in `verification/floor_theorem/`.

Keep these files together in one folder. They use only the Python standard library and were tested with Python 3.9.6. Each file finds its sibling files in its own folder, so you can run them from any working directory.

## What each file checks

| File | What it computes | Paper claim it backs |
|---|---|---|
| `privsphere.py` | The private-sphere weighting: each pattern's sphere, its weight, privacy, and the closed-form share of a ball. Run it alone for its four controls (C1 to C4). | Steps 1 and 2: privacy and the share formula. |
| `step3_dp.py` | Step 3: the exact maximum load M that any radius-m ball, centred anywhere, collects from the full chain at r = n, by a dynamic program over centres in {0,1,X}. Run it alone to compare the DP with brute force at small sizes and print the m = 6, n = 20 maximum. | Every value of M in the finite part, the threshold table and the table of eight pairs. |
| `frac_alln_tail_check.py` | The tail bound B_N and N_0 for m = 6 to 12, tests of the length and neighbour lemmas, a per-centre test of the bound, and M and R at chosen lengths (through `step3_dp.py`). | The tail table, the two lemmas, and every M and R printed in the subsection. |
| `topr.py` | The run-count-aware (top-r) bound on the maximum, using the same dynamic program. It is a library with no command line. `check_dp(N, m)` compares its lambda = 0 value with an independent per-pattern sum. | The per-run-count rows at m = 8 in the appendix (for example the bracket [101/85, 811/680] at n = 20, r = 6, and 349/306 at n = 21, r = 7). The two scripts that print those rows, `topr_run.py` and `close_r6.py`, are in this folder. |
| `topr_run.py` | Prints the per-run-count rows at m = 8: the largest r the bound proves at each length, with the lambda = 0 controls. | The m = 8 rows at n = 20, 21 and 23, and m = 7 at n = 20. |
| `close_r6.py` | Closes r = 6 at n = 20: 811/680 at lambda = 5/136 against 6/5. | The statement that r = 6 is proved at m = 8, n = 20. |
| `verify_ceiling_certificate.py` | A **sampled** consistency check of a solver certificate: arithmetic, invariance, the ceiling on about a hundred centres, and the stored witness. It does **not** establish the maximum load, and says so; `verify_ceiling_exhaustive.py` does that. Takes an unnormalised `cert_19_*.json`. | Consistency of the solver's own output. |
| `cert_19_9_9.json`, `cert_19_10_7.json`, `cert_19_11_6.json` | The three certificates themselves, as data. | The same three rows. |
| `verify_ceiling_exhaustive.py` | The independent sweep over the C_m duals: recomputes each certificate's mass and sweeps EVERY centre class, not a sample. Standard library only. Takes one or more `data/floor_certificates/cert_m*.json`. | The C_9 to C_12 dual certificates, exhaustively. |
| `data/floor_certificates/cert_m9_r9.json` and the m = 10, 11, 12 files | The four duals normalised to maximum ball load at most 1. | The same four rows. |
| `certify_cell.py`, `dump_floor_certificates.py` | The solvers that produced the certificates, with their own sweeps. **These two need SciPy**; everything else here is standard library. You do not need them to check the certificates, only to regenerate them. | Provenance of the certificate files. |
| `emit_ceiling_certificate.py` | Writes a certificate file from a solved cell. It **refuses to write** unless the sampled checker, `verify_F_independent.py` and the exhaustive sweep all pass on the record it is about to write, and unless the total mass is above r-1. Standard library only. | Provenance of the certificate files. |
| `verify_F_independent.py` | Checks the closed form the sweep funnels through, over its whole domain. Called by the emitter; runs standalone too. | The share formula the orbit reduction rests on. |

## What to run, and what a pass looks like

```
python3 frac_alln_tail_check.py --lemmas
```
Prints the tail table, one row per m from 6 to 12: N_0, n_0 = N_0 + 1, B at N_0 (below 1), B at N_0 - 1 (above 1), the worst pair (2, 0), `rises=0` (B_N never rises from N = 2m + 2 to 100) and `True` (B_N < 1 from N_0 to 100). Then it prints the length lemma (`0 increases` at each m) and the neighbour lemma (`0 violations` at N = 2 to 8). A pass means the table equals the appendix's tail table and every count is 0. About one minute.

```
python3 frac_alln_tail_check.py --centres 6 27 3
```
Tests every centre over {0,1,X} within distance 3 of a chain pattern at m = 6, N = 27 against its own bound. A pass prints `552790 non-pattern centres; 0 above B(N,s,E0)`. The same line reports a deliberately false bound (`control (count E+1): 72279 above`), which shows the test can fail. About five minutes. The appendix also quotes `--centres 7 31 3` and `--centres 8 35 3`, which take longer.

```
python3 frac_alln_tail_check.py --finite 12 29 30 31 32
```
Takes N = n - 1. For each length it prints the exact maximum M and R, the largest run count that M certifies (r(M - 1) < M, or every r when M <= 1). `every_r=True` means that length is certified for every run count. This command prints the four m = 12 rows of the table of eight pairs (n = 30 to 33). Any row of the subsection can be reprinted the same way. Larger m and n take longer: this one takes about four minutes.

```
python3 step3_dp.py
python3 privsphere.py
```
Self-tests. `step3_dp.py` must print `ok` at its three brute-force sizes before it reports anything. `privsphere.py` runs its controls C1 to C4, and C4 must show the weighting failing on the paper's own counterexample.

## Reference outputs

`logs/` holds the outputs these files gave in a clean-room run, for comparison line by line, plus the full list of 143 finite maxima (m = 6 to 12, every n from 20 to n_0 - 1). See `logs/README.md` for how each log was made.

## Checking the dual certificates in one command

```
python3 verify_ceiling_exhaustive.py data/floor_certificates/cert_m9_r9.json \
    data/floor_certificates/cert_m10_r7.json data/floor_certificates/cert_m11_r6.json \
    data/floor_certificates/cert_m12_r4.json
```

A pass prints one line per certificate, each reporting the mass against the bound it must clear and the number of centre classes swept against the number that exist, then `ALL PASS` and exit 0. The sweep is exhaustive: 78,732 of 78,732 classes at m = 9, and 10,206, 3,645 and 459 at m = 10, 11 and 12.

It fails loudly. Multiplying one block's mass by ten in `cert_m11_r6.json` gives `FAIL: 416 classes carry load above 1`, `maximum load 1259402/582027, certificate says 1`, and exit 1.

`VERIFY_CEILING.md` is the step-by-step guide to checking the dual certificates.
