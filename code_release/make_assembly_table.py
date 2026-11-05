#!/usr/bin/env python3 -I
"""Producer for Table 2, the paired per-guide assembly effect, and section 5's derived figures.

WHY THIS EXISTS. Table 2 was the only table in the paper with no producer. On 2026-10-21 that cost
us two errors at once: Q1 at six mismatches was printed 1.011 when the data give 1.010463, and the
sign column reported only rising/falling, so at three mismatches 58 of 110 unchanged guides were
simply absent and 21 + 31 did not reach 110. Both are the kind of defect a gate catches for free and
a reader catches expensively.

Quantiles use the LINEAR convention (numpy 'linear', type 7): h = (n-1)p, interpolating between
neighbours. That convention is not a free choice here, it is the one that reproduces the published
table, and the nearest-rank alternative disagrees in four cells.

WHAT THIS DOES NOT COVER, stated so the gate is not mistaken for complete. The three-term
decomposition in section 5 (gain 2.371%, GRCh38-only loss 2.128%, copy-number 0.000% / 0.464%
absolute at four mismatches) needs per-site canonical-code membership in the opposite assembly's
universe, which the held per-guide count tables do not carry. Those figures come from cluster jobs
66001436 (forward) and 66012449 (reverse) and are checked there, not here. Note for anyone editing
that passage: 1.867% is the TWO-term loss L = only38 - copynumber, NOT the GRCh38-only term. The
paper reports copy-number separately, so its loss term is 2.128%; using 1.867% beside a separate
copy-number term double-counts copy-number. check_superseded.py guards this.

Run from the RECOMB_2027 directory:  python3 -I code_release/make_assembly_table.py
"""
import json
import sys

# 2026-10-22: this was a bare relative directory, resolved against the WORKING DIRECTORY,
# so the script passed from the paper directory and raised FileNotFoundError when run from
# its own folder, which is where the README runs everything else. Four scripts shared the
# defect and the release did not carry the inputs at all, so a referee could not run them
# from anywhere. The release now ships them under data/ and the path resolves against this
# file, preferring the in-release copy and falling back to the paper tree for development.
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_REL = os.path.join(_HERE, "data", "producers_2026-10-20")
D = _REL if os.path.isdir(_REL) else os.path.join(os.path.dirname(_HERE), "producers_2026-10-20")
TOL_Q = 5e-4            # published to three decimals
TOL_PCT = 0.01          # 1 per cent relative, matching the other producers

# published Table 2: threshold -> (min, Q1, median, Q3, max, rising, unchanged, falling)
EXPECT = {
    "3": (0.667, 0.978, 1.000, 1.000, 1.462,  21, 58, 31),
    "4": (0.913, 0.993, 1.006, 1.020, 1.559,  62, 11, 37),
    "5": (0.911, 1.002, 1.012, 1.030, 1.239,  90,  1, 19),
    "6": (0.995, 1.010, 1.018, 1.051, 2.929, 107,  0,  3),
}
# section 5 prose, recomputable from the same two tables
EXPECT_PROSE = {
    "median_pct_at_4": 0.57,      # "+0.57% at four mismatches"
    "median_pct_at_6": 1.83,      # "+1.83% at six"
    "max_gain_pct_at_6": 192.9,   # "one guide gains 192.9% at six mismatches"
    "moved_at_3": 52,             # "from 52 of 110 at three mismatches"
    "moved_at_6": 110,            # "to all 110 at six"
}


def quantile(sorted_vals, p):
    n = len(sorted_vals)
    h = (n - 1) * p
    lo = int(h)
    if lo + 1 >= n:
        return sorted_vals[lo]
    return sorted_vals[lo] + (h - lo) * (sorted_vals[lo + 1] - sorted_vals[lo])


def load(name):
    with open(f"{D}/{name}") as fh:
        return {r["guide_index"]: r["sites"] for r in json.load(fh)["rows"]}


def main():
    grch38, chm13 = load("PANELREF_GRCh38_NGG.json"), load("PANELREF_CHM13_NGG.json")
    guides = sorted(set(grch38) & set(chm13))
    if len(guides) != 110:
        print(f"GATE FAIL: expected 110 paired guides, got {len(guides)}", file=sys.stderr)
        return 1

    bad, rows = 0, []
    for t in ("3", "4", "5", "6"):
        ratios, up, eq, dn = [], 0, 0, 0
        for i in guides:
            g, c = grch38[i][t], chm13[i][t]
            if not g:                      # no GRCh38 sites: ratio undefined, excluded
                continue
            ratios.append(c / g)
            up, eq, dn = (up + 1, eq, dn) if c > g else \
                         (up, eq, dn + 1) if c < g else (up, eq + 1, dn)
        ratios.sort()
        got = (ratios[0], quantile(ratios, .25), quantile(ratios, .5),
               quantile(ratios, .75), ratios[-1], up, eq, dn)
        exp = EXPECT[t]
        okq = all(abs(a - b) <= TOL_Q for a, b in zip(got[:5], exp[:5]))
        okn = got[5:] == exp[5:]
        total_ok = (up + eq + dn) == len(ratios)
        if not (okq and okn and total_ok):
            bad += 1
        rows.append((t, got, exp, okq, okn, total_ok, len(ratios)))

    # section 5's derived prose figures
    def med_at(t):
        r = sorted(chm13[i][t] / grch38[i][t] for i in guides if grch38[i][t])
        return quantile(r, .5)

    def moved_at(t):
        return sum(1 for i in guides if grch38[i][t] and chm13[i][t] != grch38[i][t])

    prose = {
        "median_pct_at_4": (med_at("4") - 1) * 100,
        "median_pct_at_6": (med_at("6") - 1) * 100,
        "max_gain_pct_at_6": (max(chm13[i]["6"] / grch38[i]["6"]
                                  for i in guides if grch38[i]["6"]) - 1) * 100,
        "moved_at_3": moved_at("3"),
        "moved_at_6": moved_at("6"),
    }
    for k, want in EXPECT_PROSE.items():
        got = prose[k]
        ok = (got == want) if isinstance(want, int) else \
             abs(got - want) <= TOL_PCT * max(abs(want), 1e-9) + 5e-2
        if not ok:
            bad += 1
        prose[k] = (got, want, ok)

    if bad:
        print(f"GATE FAIL: {bad} published figure(s) disagree with the producers. "
              f"No table written.", file=sys.stderr)
        for t, got, exp, okq, okn, tot, n in rows:
            if not (okq and okn and tot):
                print(f"  mm<={t}: got {tuple(round(x, 6) for x in got)} "
                      f"expected {exp}", file=sys.stderr)
        for k, (got, want, ok) in prose.items():
            if not ok:
                print(f"  {k}: got {got} expected {want}", file=sys.stderr)
        return 1

    print(f"GATE PASS: Table 2 and 5 derived figures reproduce from {D} "
          f"({len(guides)} paired guides).")
    print()
    print(r"% generated by code_release/make_assembly_table.py -- do not hand-edit")
    print(r"mm $\le$ & min & Q1 & median & Q3 & max & guides rising / unchanged / falling \\")
    print(r"\midrule")
    for t, got, exp, *_ in rows:
        mn, q1, md, q3, mx, up, eq, dn = got
        bold = r"\textbf{%.3f}" % md if t in ("4", "6") else "%.3f" % md
        print(f"{t} & {mn:.3f} & {q1:.3f} & {bold} & {q3:.3f} & {mx:.3f} "
              f"& {up} / {eq} / {dn} \\\\")
    print()
    for k, (got, want, _) in prose.items():
        print(f"  {k:20s} {got!r:>12} (published {want})")
    print()
    print("NOT GATED HERE: the three-term decomposition (2.371 / 2.128 / 0.000 at mm<=4) needs "
          "per-site\n  universe membership and is gated by cluster jobs 66001436 and 66012449.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
