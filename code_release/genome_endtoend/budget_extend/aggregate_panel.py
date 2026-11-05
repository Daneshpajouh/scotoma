#!/usr/bin/env python3
"""Aggregate the PANEL_EXTEND shard records.  Producer-independent of panel_extend.py:
it only reads the per-guide JSON records and re-derives every aggregate from them."""
import os
import json, os, sys

ROOT = os.environ.get("PANEL_ROOT", os.path.dirname(os.path.abspath(__file__)))
shards = sorted(d for d in os.listdir(ROOT) if d.startswith("shard_"))
per_guide = {}
for s in shards:
    p = os.path.join(ROOT, s, "record.json")
    with open(p) as fh:
        rec = json.load(fh)
    for g, v in rec["guides"].items():
        if g in per_guide:
            print("DUPLICATE GUIDE", g, file=sys.stderr); sys.exit(2)
        per_guide[g] = v

print("guides aggregated:", len(per_guide))
budgets = ["mm5", "mm6", "mm7"]
bad_relonly = []
bad_attrib = []
rows = []
for b in budgets:
    tot_c = tot_r = tot_co = tot_ta = 0
    nguides_def = 0
    for g, v in per_guide.items():
        d = v[b]
        tot_c += d["corrected"]; tot_r += d["released"]
        tot_co += d["corrected_only"]; tot_ta += d["terminal_attributable"]
        if d["released_only"] != 0:
            bad_relonly.append((g, b, d["released_only"]))
        if d["terminal_attributable"] != d["corrected_only"]:
            bad_attrib.append((g, b, d["corrected_only"], d["terminal_attributable"]))
        if d["corrected_only"] > 0:
            nguides_def += 1
    rows.append((b, tot_r, tot_c, tot_co, tot_ta, nguides_def,
                 100.0 * tot_co / tot_c if tot_c else 0.0))

print()
print(f"{'budget':>7} {'released':>14} {'corrected':>14} {'corr-only':>12} {'term-attr':>12} {'guides>0':>9} {'share%':>8}")
for b, r, c, co, ta, ng, sh in rows:
    print(f"{b:>7} {r:>14,} {c:>14,} {co:>12,} {ta:>12,} {ng:>9} {sh:>8.3f}")

print()
print("released_only != 0 anywhere :", "NONE" if not bad_relonly else bad_relonly[:10])
print("terminal_attributable != corrected_only :", "NONE" if not bad_attrib else bad_attrib[:10])

# guides with zero deficit at every budget
zero_all = sorted(g for g, v in per_guide.items() if all(v[b]["corrected_only"] == 0 for b in budgets))
print()
print("guides with zero deficit at every budget:", len(zero_all))
print(zero_all)
# consistency: the deficit set should be identical across budgets
sets = {b: set(g for g, v in per_guide.items() if v[b]["corrected_only"] > 0) for b in budgets}
print("deficit guide set identical across budgets:", sets["mm5"] == sets["mm6"] == sets["mm7"])

out = {"n_guides": len(per_guide),
       "rows": [{"budget": b, "released": r, "corrected": c, "corrected_only": co,
                 "terminal_attributable": ta, "guides_with_deficit": ng, "share_pct": sh}
                for b, r, c, co, ta, ng, sh in rows],
       "zero_deficit_guides": zero_all,
       "released_only_violations": bad_relonly,
       "attribution_violations": bad_attrib}
with open(os.path.join(ROOT, "PANEL_AGGREGATE.json"), "w") as fh:
    json.dump(out, fh, indent=2, sort_keys=True)
print()
print("AGGREGATE_OK" if (not bad_relonly and not bad_attrib) else "AGGREGATE_FLAGGED")
