# The budget-ten extensions, with everything needed to check them

The companion paper's budget-ten frontier certifies `K = r` through run counts 10, 15, 19, 23 and
24 at guide lengths 20, 21, 22, 23 and 24. Three of those rows rest on ten cells closed after the
first certificate package was deposited, and this directory holds their evidence. Indexed by word
length `N = n - 1`, the ten are `(20,11)`, `(20,12)`, `(21,12)`, `(21,13)`, `(22,20)`,
`(22,21)`, `(20,14)`, `(20,15)`, `(21,19)` and `(22,23)`, all at budget `m = 10`.

Four cells were closed on 3 November, after `independent_scan/H3_summary.json` had been
returned. Each of their independent records is therefore deposited beside that summary as
`independent_scan/cell_N_r_10_result.json` rather than by rewriting a returned artefact, which
is why the summary still reports six and still hashes as it did.

`(20,15)` carries the narrowest margin of the ten: its mass `8193217224/585141539` sits only
`1235678/585141539` above `r - 1 = 14`, about `0.0021`, so the exact arithmetic is the whole
argument there and a float would not settle it. `(22,23)` is the full diagonal `r = n` at guide
length 23, so step-down from it closes every run count at that length and no cell is left open
there. With `(21,19)` each of guide lengths 21, 22 and 23 now rests on its own certificate rather
than on padding from a shorter length.

## Why ten cells settle three table rows

For each cell a fractional packing has total mass strictly above `r - 1`, which forces `K >= r`,
and the library of `r` balls gives `K <= r`, so `K = r`. Step-down in run count and padding then
carry each closure to every smaller run count and every greater length, which is what turns ten
cells into the certified prefixes the table prints.

## What is here

    cells/cell_N_r_10/decision.certificate.json   the packing, one entry per support orbit
    cells/cell_N_r_10/best.certificate.json       the best packing found before the decision
    cells/cell_N_r_10/RESULT.json                 that cell's own record: tau bounds, K, status
    cells/cell_N_r_10/final_verification.json     the producer's closing receipt
    cells/cell_N_r_10/run.log                     the run's own log
    independent_scan/H3_summary.json              the independent exhaustive maximum-load scan
    independent_scan/cell_N_r_10_result.json      the same scan's record for each of the last four
    independent_scan/exact_tensor.cpp             that scan's implementation
    independent_scan/run_verify.py                its driver
    independent_scan/local_checks.py              its local validation
    independent_scan/METHOD.md                    what it computes and how
    budget_ten_frontier.json                      the 120-cell ledger the table and figure derive from
    generate_frontier.py                          the producer that turns the ledger into both
    verify_m10_extensions.py                      the check below

Inside the producer records, one bookkeeping field names the run directory the figure was first
reported from. The account component of those paths reads `YOUR_ACCOUNT` in every deposited copy;
nothing is computed from them and no check reads them.

The two certificate implementations are deliberately separate. The packing was produced by one
program; the maximum load was then re-derived over the complete centre-class space by an
independently written exhaustive scan, which also rejects two forged inputs per cell, a rescaled
class and a moved class. Agreement between them is the point, so both are here.

## Checking it, with no cluster and no sweep

    python3 -I verify_m10_extensions.py

In exact rational arithmetic it re-derives each cell's total mass by **summing the deposited
packing**, rather than reading any field that claims a total, then requires: the mass exceeds
`r - 1`; the cell's own record agrees with that sum to the last digit; the independent scan agrees
with it too; the scan reports an exhaustive exact pass; both of its forgery controls were rejected;
`K_lower = K_upper = r` with the cell marked closed; and the ledger marks every one of the ten
closed, under `n = N + 1`.

    python3 -I verify_m10_extensions.py --control

removes one unit from one cell's packing at a time, ten times over, and requires the check to
fail every time. Each run fails on the cell it perturbed and on no other, which is what makes
the passing run mean something. Running it on every cell, rather than on one, exercises both
routes by which a mass reaches this check: the returned six-cell summary and the four
single-cell records.

## What is not claimed

The exhaustive scan was run on a cluster and its returned summary is deposited, not re-run here;
re-running it needs the resources named in `independent_scan/METHOD.md`. This package lets a
reader confirm that the deposited packings, the cells' own records, the independent scan and the
printed table all agree, and that a perturbed packing is caught. It does not re-execute the scan.
