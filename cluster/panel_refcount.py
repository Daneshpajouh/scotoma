#!/usr/bin/env python3
"""Per-guide, per-threshold reference candidate counts for a whole panel.

Generalises the single-guide class-factor measurement: the budget multiplier from widening the
search class is a per-guide quantity, and N at <=6 spans more than twenty-fold across this panel,
so one guide cannot stand for it.
"""
import argparse, json, sys, time

import numpy as np

from scotoma.panel import load_panel
from scotoma.coding import spacer_of
from scotoma.universe import load_universe, panel_reference_counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--changeseq", action="store_true")
    ap.add_argument("--thresholds", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5, 6])
    ap.add_argument("--pam", default="NGG")
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    t0 = time.time()
    codes, mult, meta = load_universe(a.universe, with_mult=True)
    print(f"universe {meta.get('label')} {codes.size:,} codes "
          f"(recorded PAM {meta.get('pam', '<unrecorded>')}), {time.time()-t0:.0f}s", flush=True)

    guides, report = load_panel(a.panel, changeseq=a.changeseq)
    print(f"panel: {report['distinct_guides']} guides, {report['rejected']} rejected", flush=True)
    if guides.size == 0:
        raise SystemExit("GATE FAIL empty panel")

    def progress(i, n):
        if i % 10 == 0 or i == n:
            print(f"  chunk {i}/{n}  {time.time()-t0:.0f}s", flush=True)

    rows = panel_reference_counts(codes, mult, guides,
                                  thresholds=tuple(a.thresholds), progress=progress)
    for r in rows:
        r["spacer"] = spacer_of(int(guides[r["guide_index"]]))

    # Every guide's own on-target must be found; that is the gate.
    found = sum(1 for r in rows if r["sites"][0] >= 1)
    print(f"\nGATE on-targets found: {found} of {len(rows)}", flush=True)

    payload = {"universe": meta, "universe_codes": int(codes.size), "pam": a.pam,
               "panel_report": report, "thresholds": a.thresholds,
               "on_targets_found": found, "guides": len(rows), "rows": rows,
               "wall_seconds": round(time.time() - t0, 1)}
    with open(a.out_json, "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"wrote {a.out_json}  total {time.time()-t0:.0f}s", flush=True)

    sites4 = sorted(r["sites"][4] for r in rows)
    sites6 = sorted(r["sites"][6] for r in rows)
    def med(v): return v[len(v)//2]
    print(f"\nN at <=4: min {sites4[0]:,}  median {med(sites4):,}  max {sites4[-1]:,}")
    print(f"N at <=6: min {sites6[0]:,}  median {med(sites6):,}  max {sites6[-1]:,}")
    factors = sorted(r["sites"][6]/r["sites"][4] for r in rows if r["sites"][4])
    print(f"per-guide <=4 -> <=6 class factor: min {factors[0]:.1f}x  "
          f"median {med(factors):.1f}x  max {factors[-1]:.1f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
