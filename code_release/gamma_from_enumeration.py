#!/usr/bin/env python3
"""Detectability gap and validation budget for a finite enumerated candidate set.

For a guide with a finite candidate set of N off-target sites, assume (A_k): if the guide
is not clean, at least k enumerated sites are cleaved at a detectable level. Then every
unclean guide puts false-alarm mass at least k/N on the sites, so the detectability gap is
at least k/N, with no assumption on scores, calibration, ROC shape or continuity.

Two sampling models are tabulated.
  * with replacement (i.i.d. draws): the sample size of Theorem 2(i) at gamma = k/N,
    m0 = min{m : (1 - k/N)^m <= delta1}. At small k it exceeds N, which signals that the
    model is wrong for a finite enumerated set.
  * without replacement (the experiment): the hypergeometric zero-defect plan,
    m = min{m : C(N-k, m) / C(N, m) <= delta1}, bounded by N by construction.

The unit is the site (a genomic window), not the distinct 23-mer: cleavage happens at a site,
and the site count is the conservative denominator, since gamma = k/N is smaller with the
larger N and the requirement therefore larger. N excludes the on-target window.

Input: data/exacel_enumeration_summary.json, the output of enumerate_sites.py for the
exa-cel spacer. Run: python3 gamma_from_enumeration.py
"""
import json
import math
import pathlib
from fractions import Fraction as F

from m0_search import m0_hypergeom, m0_iid

HERE = pathlib.Path(__file__).resolve().parent
S = json.loads((HERE / "data" / "exacel_enumeration_summary.json").read_text())
N_SITES = S["off_target_windows"]
N_SEQ = S["distinct_23mers"] - 1          # the on-target 23-mer is one of the distinct sequences
assert N_SITES == S["windows_total"] - S["on_target_windows"]
assert N_SITES > N_SEQ

print(f"spacer {S['spacer']}, {S['reference']}, <= {S['max_mismatches']} mismatches, NGG and NAG PAMs")
print(f"  genomic windows returned by the scan       : {S['windows_total']:,}")
print(f"  on-target windows (zero mismatches)        : {S['on_target_windows']}")
print(f"  candidate off-target sites N (windows)     : {N_SITES:,}   <- denominator used")
print(f"  distinct off-target 23-mers (for contrast) : {N_SEQ:,}")
print("  mismatch histogram of windows: " + ", ".join(
    f"{k}mm:{v:,}" for k, v in S["mismatch_histogram_windows"].items()))
print()

print("WITH REPLACEMENT (i.i.d. draws): m0 = min{m : (1 - k/N)^m <= delta1}")
for delta1 in ("0.05", "0.025", "0.01"):
    print(f"delta1 = {delta1}")
    print(f"  {'population':>10} {'N':>8} {'k':>3} {'gamma = k/N':>13} {'m0':>9} {'m0/N':>7}")
    for tag, N in (("sites", N_SITES), ("sequences", N_SEQ)):
        for k in (1, 2, 5, 10):
            m = m0_iid(F(k, N), delta1)
            print(f"  {tag:>10} {N:>8,} {k:>3} {k / N:>13.3e} {m:>9,} {m / N:>7.3f}")
print()

print("WITHOUT REPLACEMENT (hypergeometric zero-defect plan), N = candidate off-target sites")
print("  'closed form' is ceil(N * (1 - delta1^(1/k))), the leading-order expression; it is")
print("  never smaller than the exact requirement.")
print(f"  {'delta1':>8} {'k':>4} {'exact sites':>12} {'closed form':>12} {'exact / N':>10}")
for delta1 in ("0.05", "0.025", "0.01"):
    for k in (1, 2, 5, 10, 25, 50):
        ex = m0_hypergeom(N_SITES, k, delta1)
        cf = math.ceil((1 - float(delta1) ** (1 / k)) * N_SITES)
        print(f"  {delta1:>8} {k:>4} {ex:>12,} {cf:>12,} {ex / N_SITES:>10.4f}")
print()
print("READING. At k = 1 the certificate costs 95 per cent of the candidate set at delta1 = 0.05:")
print("a census in all but name, which is the price of a zero-false-alarm claim when a single")
print("missed site would break it. The price falls fast in k. The operative question for a trial is")
print("not 'is this guide clean' but 'if it is not, at how many sites at minimum', which is a")
print("mechanism question about the editing enzyme and the assay, not something the validation")
print("data can answer (Proposition 2 of the paper).")
