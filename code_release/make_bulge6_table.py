#!/usr/bin/env python3 -I
"""Price the bulge dial at six mismatches, score the registered prediction, emit the appendix rows.

Job 66532070 (`BULGEMM6`, Fir) runs the same launcher that produced every published bulge figure in
this paper, with `MM=6` instead of `MM=4` and a separate work directory. This consumes its output.

WHY THE mm4 ASSERTIONS ARE HERE. The mm4 census from that same launcher reproduces the paper exactly:
base 34,652, added 1,175,323, panel factor 34.9179x, per-guide median 40.135x, RNA over DNA 3.41x.
Those are asserted below, so if the baseline ever drifts this refuses to emit an mm6 comparison
against a baseline that no longer matches the paper.

THE PREDICTION, registered 2026-10-21 in REGISTERED_2026-10-21_bulge_factor_at_six_mismatches.md
BEFORE any output existed: the one-base panel factor at six mismatches is LOWER than 34.9x, in the
range 10x to 30x. KILL CRITERION: >= 34.9x falsifies the diminishing-returns model and is to be
reported as such, not reinterpreted. A value below 34.9x but outside 10x to 30x is a RANGE MISS and
is to be reported as a miss.

Run from RECOMB_2027/. Exit 0 if the mm4 calibration holds (and prints the mm6 comparison when the
file is present), 1 if the calibration fails, 2 if the mm6 file is malformed.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
# 2026-10-22: the release now carries its own copy of the census files, so every gate in this
# folder runs from the folder a referee receives. The parent-tree path is kept as a fallback
# for the development tree, where mm6_pending is where the cluster jobs deliver it.
_IN_RELEASE = os.path.join(HERE, 'data', 'mm6_pending')
PEND = _IN_RELEASE if os.path.isdir(_IN_RELEASE) else os.path.join(os.path.dirname(HERE), 'mm6_pending')
MM4 = os.path.join(PEND, 'BULGE_CENSUS_mm4.json')
# MM6_PATH lets the three prediction branches be exercised on synthetic inputs before the
# real file lands; a scorer that has only ever reported a hit is decoration.
MM6 = os.environ.get('MM6_PATH', os.path.join(PEND, 'BULGE_CENSUS_mm6.json'))
MM5 = os.environ.get('MM5_PATH', os.path.join(PEND, 'BULGE_CENSUS_mm5.json'))

# Round 45: the third point. Registered BEFORE job 66565396 was submitted, in
# REGISTERED_2026-10-21_bulge_factor_at_five_mismatches.md: 29.1x, band 26x to 32x.
# DECISIVE KILL: outside [23.849, 34.918] means the sequence is NOT monotone decreasing
# and the appendix's two-point reading is falsified.
PRED5 = dict(point=29.1, lo=26.0, hi=32.0, mono_lo=23.849, mono_hi=34.918)

PAPER_MM4 = dict(ungapped=34652, add_merged=1175323, panel_factor=34.9, per_guide_median=40.1,
                 rna_over_dna=3.4)
PRED = dict(direction='lower', lo=10.0, hi=30.0, kill_at=34.9)

# The figures the appendix now quotes at six mismatches. Asserted so the paper and the artifact
# cannot drift apart: job 66532070, scored 2026-10-21.
PAPER_MM6 = dict(ungapped=1854173, add_dna=12415044, add_rna=33291743, add_merged=42366572,
                 panel_factor=23.8, per_guide_median=27.5, rna_over_dna=2.68)

def load(p):
    if not os.path.exists(p): return None
    try: return json.load(open(p))
    except Exception as e:
        print(f"MALFORMED {p}: {e}"); sys.exit(2)

def main():
    d4 = load(MM4)
    if d4 is None:
        # same fix as verify_missed_actives.py, 2026-10-21: absence is not a pass.
        print(f"GATE FAILED: mm4 baseline absent at {MM4}; the dial comparison has no anchor.")
        return 2
    t4 = d4['totals']
    bad = []
    if t4['ungapped'] != PAPER_MM4['ungapped']:
        bad.append(f"mm4 base {t4['ungapped']} != the paper's {PAPER_MM4['ungapped']}")
    if t4['add_merged'] != PAPER_MM4['add_merged']:
        bad.append(f"mm4 added {t4['add_merged']} != the paper's {PAPER_MM4['add_merged']}")
    if abs(t4['panel_factor_merged'] - PAPER_MM4['panel_factor']) > 0.05:
        bad.append(f"mm4 panel factor {t4['panel_factor_merged']:.3f} != the paper's "
                   f"{PAPER_MM4['panel_factor']}")
    if abs(t4['per_guide_factor_merged'][0] - PAPER_MM4['per_guide_median']) > 0.05:
        bad.append(f"mm4 per-guide median {t4['per_guide_factor_merged'][0]:.3f} != the paper's "
                   f"{PAPER_MM4['per_guide_median']}")
    rd = t4['add_rna'] / t4['add_dna']
    if abs(rd - PAPER_MM4['rna_over_dna']) > 0.05:
        bad.append(f"mm4 RNA/DNA {rd:.3f} != the paper's {PAPER_MM4['rna_over_dna']}")
    if bad:
        print("CALIBRATION FAILED against the paper; refusing to compare mm6 to a drifted baseline:")
        for b in bad: print("  " + b)
        return 1
    print(f"mm4 baseline calibrates to the paper: base {t4['ungapped']:,}, "
          f"added {t4['add_merged']:,}, panel {t4['panel_factor_merged']:.3f}x, "
          f"per-guide median {t4['per_guide_factor_merged'][0]:.3f}x, RNA/DNA {rd:.2f}x")

    d6 = load(MM6)
    if d6 is None:
        print(f"\nPENDING: job 66532070 has not delivered {os.path.basename(MM6)} yet.")
        print("The prediction stands as registered: LOWER than 34.9x, in 10x to 30x, "
              "killed at >= 34.9x.")
        return 0

    t6 = d6['totals']
    if d6.get('threshold') != 6:
        print(f"MALFORMED: mm6 file reports threshold {d6.get('threshold')}, expected 6"); return 2
    f6 = t6['panel_factor_merged']; f4 = t4['panel_factor_merged']
    print(f"\n==== the bulge dial at six mismatches ====")
    print(f"  ungapped base        {t6['ungapped']:,}   (mm4: {t4['ungapped']:,})")
    print(f"  added by DNA bulge   {t6['add_dna']:,}   (mm4: {t4['add_dna']:,})")
    print(f"  added by RNA bulge   {t6['add_rna']:,}   (mm4: {t4['add_rna']:,})")
    print(f"  added, merged        {t6['add_merged']:,}   (mm4: {t4['add_merged']:,})")
    print(f"  PANEL FACTOR         {f6:.3f}x   (mm4: {f4:.3f}x)")
    pg6 = t6['per_guide_factor_merged']
    print(f"  per-guide median     {pg6[0]:.3f}x (min {pg6[1]:.3f}, max {pg6[2]:.3f})")
    print(f"  RNA over DNA         {t6['add_rna']/t6['add_dna']:.2f}x   (mm4: {rd:.2f}x)")

    bad6 = []
    for k, want in (('ungapped', 'ungapped'), ('add_dna', 'add_dna'),
                    ('add_rna', 'add_rna'), ('add_merged', 'add_merged')):
        if t6[want] != PAPER_MM6[k]:
            bad6.append(f"{k} {t6[want]:,} != the appendix's {PAPER_MM6[k]:,}")
    if abs(f6 - PAPER_MM6['panel_factor']) > 0.05:
        bad6.append(f"panel factor {f6:.3f} != the appendix's {PAPER_MM6['panel_factor']}")
    if abs(pg6[0] - PAPER_MM6['per_guide_median']) > 0.05:
        bad6.append(f"per-guide median {pg6[0]:.3f} != the appendix's {PAPER_MM6['per_guide_median']}")
    if abs(t6['add_rna']/t6['add_dna'] - PAPER_MM6['rna_over_dna']) > 0.01:
        bad6.append(f"RNA/DNA {t6['add_rna']/t6['add_dna']:.3f} != the appendix's {PAPER_MM6['rna_over_dna']}")
    if bad6:
        print("\nGATE FAILED: the mm6 artifact no longer matches the figures the appendix quotes:")
        for b in bad6:
            print("  " + b)
        return 2
    print("\n  all six mm6 figures quoted in the appendix match the artifact")

    print(f"\n==== scoring the registered prediction ====")
    lower = f6 < PRED['kill_at']
    inrange = PRED['lo'] <= f6 <= PRED['hi']
    if not lower:
        print(f"  PREDICTION KILLED. {f6:.3f}x >= {PRED['kill_at']}x. The diminishing-returns model is")
        print(f"  WRONG: bulge policy becomes relatively DEARER as the mismatch budget loosens.")
        print(f"  Report as a falsified prediction, do not reinterpret.")
    elif inrange:
        print(f"  PREDICTION HIT. {f6:.3f}x is below {PRED['kill_at']}x and inside "
              f"{PRED['lo']}x to {PRED['hi']}x.")
    else:
        print(f"  RANGE MISS. {f6:.3f}x is below {PRED['kill_at']}x, so the DIRECTION was right, but")
        print(f"  outside the registered band {PRED['lo']}x to {PRED['hi']}x. Report as a miss.")

    print(f"\n==== appendix rows, at six mismatches ====")
    print(f"one bulge of one base, merged & ${t6['add_merged']:,} & ${f6:.1f}\\times$ \\\\".replace(',', '{,}'))
    d5 = load(MM5)
    if d5 is None:
        print(f"\n  mm5 (job 66565396) pending: {os.path.basename(MM5)} not here yet.")
        print(f"  registered: {PRED5['point']}x, band {PRED5['lo']}x to {PRED5['hi']}x; "
              f"KILLED if outside [{PRED5['mono_lo']}, {PRED5['mono_hi']}] (non-monotone)")
    else:
        if d5.get('threshold') != 5:
            print(f"MALFORMED: mm5 file reports threshold {d5.get('threshold')}"); return 2
        t5 = d5['totals']; f5 = t5['panel_factor_merged']
        print(f"\n==== the bulge dial at FIVE mismatches ====")
        print(f"  ungapped base   {t5['ungapped']:,}")
        print(f"  added, merged   {t5['add_merged']:,}")
        print(f"  PANEL FACTOR    {f5:.3f}x")
        print(f"  per-guide median {t5['per_guide_factor_merged'][0]:.3f}x")
        print(f"\n==== the three-point trend ====")
        print(f"  <=4  {f4:.3f}x      <=5  {f5:.3f}x      <=6  {f6:.3f}x")
        mono = f4 > f5 > f6
        print(f"  monotone decreasing: {'YES' if mono else 'NO'}")
        print(f"\n==== scoring the mm5 prediction ====")
        if not (PRED5['mono_lo'] <= f5 <= PRED5['mono_hi']):
            print(f"  PREDICTION KILLED, and this is the interesting outcome. {f5:.3f}x lies outside")
            print(f"  [{PRED5['mono_lo']}, {PRED5['mono_hi']}], so the bulge dial is NOT monotone in the")
            print(f"  budget. The appendix's two-point decline reading is falsified and must be")
            print(f"  rewritten to say the dependence is not a smooth decline. Do not reinterpret.")
        elif PRED5['lo'] <= f5 <= PRED5['hi']:
            print(f"  PREDICTION HIT. {f5:.3f}x is inside {PRED5['lo']}x to {PRED5['hi']}x, and the")
            print(f"  sequence is monotone decreasing across all three budgets.")
        else:
            print(f"  RANGE MISS. {f5:.3f}x keeps monotonicity but falls outside the registered band")
            print(f"  {PRED5['lo']}x to {PRED5['hi']}x, so the log-base interpolation model is wrong.")
            print(f"  Report as a miss.")

    print(f"\n==== the ordering claim at this budget ====")
    MISMATCH_DIAL = 70.0
    print(f"  bulge, one base   {f6:.1f}x   vs the mismatch dial's {MISMATCH_DIAL}x  -> "
          f"{'bulge larger' if f6 > MISMATCH_DIAL else 'MISMATCH DIAL LARGER'}")
    print(f"  the paper's ordering claim is scoped to four mismatches (round 32), so this is")
    print(f"  additional information and contradicts nothing as written.")
    return 0

# entry point moved to the end of the file so the four-point check is reachable (2026-10-21)

# ---- the four-point dial, added 2026-10-21 when the mm3 census landed ----
# Registered before mm3 was submitted: 43x, band 38x-50x, monotonicity kill at or below 34.918x.
PRED3 = dict(point=43.0, lo=38.0, hi=50.0, mono_floor=34.918)
FOURPOINT = {3: 38.229, 4: 34.918, 5: 28.264, 6: 23.849}
FOURPOINT_PG = {3: 42.656, 4: 40.135, 5: 33.990, 6: 27.495}

def four_point():
    import json as _json, os as _os
    env = {3: "MM3_PATH", 4: "MM4_PATH", 5: "MM5_PATH", 6: "MM6_PATH"}
    got, bad = {}, []
    for b, ev in env.items():
        # Was _os.path.join("mm6_pending", ...), which resolves against the WORKING DIRECTORY
        # while lines 23-29 resolve against __file__. The script therefore passed when run
        # from the parent and failed when run from its own folder, which is where the README
        # runs everything else. One file, two conventions, and only one of them travelled.
        path = _os.environ.get(ev, _os.path.join(PEND, f"BULGE_CENSUS_mm{b}.json"))
        if not _os.path.exists(path):
            print(f"GATE FAILED: census for <={b} absent at {path}; the four-point trend cannot be checked")
            return 2
        t = _json.load(open(path))["totals"]
        pf = t["panel_factor_merged"]
        # the stored factor must equal its own definition
        defn = (t["ungapped"] + t["add_merged"]) / t["ungapped"]
        if abs(pf - defn) > 1e-9:
            bad.append(f"<={b}: stored panel factor {pf} != (ungapped+add_merged)/ungapped {defn}")
        got[b] = (pf, t["per_guide_factor_merged"][0])
    print()
    print("==== the bulge dial at four budgets ====")
    for b in sorted(got):
        pf, pg = got[b]
        exp, expg = FOURPOINT[b], FOURPOINT_PG[b]
        ok = abs(pf - exp) < 5e-4 and abs(pg - expg) < 5e-4
        print(f"  <={b}  panel {pf:7.3f}x  per-guide {pg:7.3f}x  {'ok' if ok else 'DRIFT'}")
        if not ok:
            bad.append(f"<={b}: measured ({pf:.3f}, {pg:.3f}) vs armed ({exp}, {expg})")
    order = sorted(got)
    mono_panel = all(got[order[i]][0] > got[order[i+1]][0] for i in range(len(order)-1))
    mono_pg = all(got[order[i]][1] > got[order[i+1]][1] for i in range(len(order)-1))
    print(f"  monotone decreasing: panel {'YES' if mono_panel else 'NO'}, "
          f"per-guide {'YES' if mono_pg else 'NO'}")
    if not (mono_panel and mono_pg):
        bad.append("the dial is NOT monotone across the four budgets, so the appendix sentence is false")
    pf3 = got[3][0]
    if pf3 <= PRED3["mono_floor"]:
        bad.append(f"<=3 factor {pf3:.3f} is at or below the <=4 factor, monotonicity killed")
    verdict = "HIT" if PRED3["lo"] <= pf3 <= PRED3["hi"] else "RANGE MISS"
    print(f"  mm3 registered 43x in [38x, 50x]; measured {pf3:.3f}x -> PREDICTION {verdict}")
    if bad:
        print("\nGATE FAILED:")
        for b in bad: print("  " + b)
        return 2
    print("  CLEAN: four-point dial matches the armed values and is monotone in both measures")
    return 0

if __name__ == "__main__":
    _rc_main = main()
    _rc_fp = four_point()
    sys.exit(_rc_main or _rc_fp)

