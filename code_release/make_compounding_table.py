#!/usr/bin/env python3 -I
"""Producer for the dial-compounding result: rank correlations and the compounded per-guide factor.

Every input is a held per-guide producer in producers_2026-10-20/, the same five tables
make_dials_figure.py uses. The script recomputes each published figure and EXITS NON-ZERO without
printing a table if any disagrees, so a silent data swap cannot yield a plausible-looking result.

Only FOUR dials enter the product. The two bulge entries are one dial at two settings, nested, and
measured here at rho = 0.772; multiplying both would double-count the same choice.

Run from the RECOMB_2027 directory:  python3 -I code_release/make_compounding_table.py
"""
import json, sys

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
EXPECT = {                      # published figures, checked to the precision stated in the paper
    ("mismatch", "bulge<=2"): 0.716,
    ("PAM", "bulge<=2"):      0.391,
    ("PAM", "mismatch"):      0.253,
    ("assembly", "bulge<=2"): -0.152,
    ("bulge=1", "bulge<=2"):  0.772,
}
EXPECT_PROD = {"median_of_products": 16166.1, "product_of_medians": 13787.9,
               "p90": 46610.6, "max": 132955.7, "min": 666.0}
TOL = 0.01                      # 1 per cent, matching the dials-figure gate

def med(xs):
    s = sorted(xs); n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2

def rank(xs):
    order = sorted(range(len(xs)), key=lambda k: xs[k]); r = [0.0] * len(xs); i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]: j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1): r[order[k]] = avg
        i = j + 1
    return r

def pearson(a, b):
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = sum((x - ma) ** 2 for x in a) ** .5
    db = sum((y - mb) ** 2 for y in b) ** .5
    return num / (da * db) if da and db else 0.0

def spearman(a, b): return pearson(rank(a), rank(b))

def load():
    G = {r["guide_index"]: r["sites"] for r in json.load(open(f"{D}/PANELREF_GRCh38_NGG.json"))["rows"]}
    C = {r["guide_index"]: r["sites"] for r in json.load(open(f"{D}/PANELREF_CHM13_NGG.json"))["rows"]}
    S1 = {r["guide_index"]: r for r in json.load(open(f"{D}/BULGE_CENSUS_size1_FULL_mm4.json"))["per_guide"]}
    S2 = {r["guide_index"]: r for r in json.load(open(f"{D}/BULGE_CENSUS_size2_mm4.json"))["per_guide"]}
    PM = {r["guide_index"]: r for r in json.load(open(f"{D}/PAM_PANEL_FACTOR_PERGUIDE.json"))["4"]["per_guide"]}
    gs = [i for i in sorted(set(G) & set(C) & set(S1) & set(S2) & set(PM)) if G[i]["4"]]
    return gs, {
        "assembly": [C[i]["4"] / G[i]["4"] for i in gs],
        "PAM":      [PM[i]["factor"] for i in gs],
        "mismatch": [G[i]["6"] / G[i]["4"] for i in gs],
        "bulge<=2": [S2[i]["factor_merged"] for i in gs],
        "bulge=1":  [S1[i]["factor_merged"] for i in gs],
    }

def main():
    gs, d = load()
    if len(gs) != 110:
        print(f"GATE FAILED: expected 110 guides with all five dials, got {len(gs)}"); return 2
    bad = []
    for (a, b), exp in EXPECT.items():
        got = spearman(d[a], d[b])
        if abs(got - exp) > 0.001:
            bad.append(f"rho({a},{b}) = {got:.3f}, published {exp:.3f}")
    keys = ["assembly", "PAM", "mismatch", "bulge<=2"]
    prod = [d["assembly"][k] * d["PAM"][k] * d["mismatch"][k] * d["bulge<=2"][k] for k in range(len(gs))]
    ps = sorted(prod)
    pm = 1.0
    for k in keys: pm *= med(d[k])
    got_prod = {"median_of_products": med(prod), "product_of_medians": pm,
                "p90": ps[int(.9 * len(ps))], "max": ps[-1], "min": ps[0]}
    for k, exp in EXPECT_PROD.items():
        if abs(got_prod[k] - exp) / exp > TOL:
            bad.append(f"{k} = {got_prod[k]:,.1f}, published {exp:,.1f}")
    if bad:
        print("GATE FAILED, published figures not reproduced:")
        for b in bad: print("  " + b)
        return 3
    print(f"guides with all five dials: {len(gs)}")
    print("GATE PASSED: all five published correlations and all five product figures reproduced\n")
    print("           " + "".join(f"{k:>11}" for k in keys))
    for a in keys:
        print(f"{a:>10} " + "".join(f"{spearman(d[a], d[b]):>11.3f}" for b in keys))
    print(f"\nnested check  bulge=1 vs bulge<=2 : {spearman(d['bulge=1'], d['bulge<=2']):.3f}")
    print(f"\nproduct of medians {pm:>12,.1f}x | median of products {med(prod):>12,.1f}x")
    print(f"min {ps[0]:>10,.1f}x | p90 {ps[int(.9*len(ps))]:>12,.1f}x | max {ps[-1]:>12,.1f}x")
    return 0

if __name__ == "__main__":
    sys.exit(main())
