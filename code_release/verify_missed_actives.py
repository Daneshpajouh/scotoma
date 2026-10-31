#!/usr/bin/env python3 -I
"""Gate for the 129 measured active off-target sites the production position set fails to nominate.

Holds the whole chain, not just the headline, because every link was a place the earlier version of
this analysis went wrong: an under-powered denominator (27 sites when 23,692 were held), "sits at an
omitted position" mistaken for "would be missed", and no check that a bulge-free route exists.

Input is the captured site records (spacer, target, 73-base deposit window, reads), produced by
streaming the Lazzarotto 2020 CHANGE-seq table from the cluster. Exits non-zero on any disagreement.
"""
import json, sys, os

EXPECT = dict(at_omitted=260, exclusively_missed=129, guides=40,
              bulgefree_within_4=0, bulgefree_within_6=2,
              reads_min=6, reads_max=724)

# Added 2026-10-21 (round 18). Two facts measured from the CHANGE-seq deposit that the paper now
# states and that a rerun must reproduce. They are recorded here as constants because the capture
# step that produced them streams an 878 MB table from the cluster and is not run on every gate
# pass; the companion script dnarna_split.py regenerates them.
#   1. The DNA omission is provably lossy as a statement about patterns but costs NOTHING measured:
#      of 6,425 DNA-bulge measured actives the wildcard slot never reaches slot n.
#   2. Proposition 1's indicator partitions the 583 actives at the omitted RNA terminal into 323 the
#      law says are rescued, ALL of which are, and 260 it says are at risk. Zero counterexamples.
#
# 2026-10-22, ROUND 49, A DEFECT OF THIS INSTRUMENT, FIXED. These were printed in the same block as
# numbers measured from the input, with nothing marking the difference, and the only arithmetic
# enforced on them was 323 + 260 = 583, which holds BY CONSTRUCTION. A peer seat consequently told
# two other seats that this gate independently confirms the 583. It does not, and did not. Two of
# these are now MEASURED from the input file, the rest are labelled as deposit-derived with the
# producer that re-earns them, and the printout separates the two.
DEPOSIT_DERIVED = dict(dna_bulge_actives=6425, dna_at_slot_n=0,
                       rna_bulge_actives=17267, rna_at_pos_0=0, rna_at_terminal=583,
                       law_lossless=323)
# Re-earned by: make_sites129.py over the CHANGE-seq deposit, which prints all of these from
# 4,936,278 rows and reproduces the input file set-identically at record-set SHA-256
# dcc09aa424db55c216e356619f89a29120c737fc6f9fcd9b2bf1fb3aa9b0e479. Nothing below that line is
# available from the 260 records this script holds, so this script does not claim it.

# Added 2026-10-21 (round 19). The same comparison read as a design-time screen. A guide can lose a
# site to this omission only if g_n != g_(n-1), so the exposure is excludable before any scan. The
# claim is a NECESSARY condition, not a sufficient one, and the gate holds both halves: every guide
# that carries an exclusively-missed site must be in the exposed set, and the protected set must be
# empty of them. If a protected guide ever carries one, the proposition is falsified.
SCREEN = dict(panel=110, exposed=73, protected=37, carrying=40,
              carrying_in_exposed=40, carrying_in_protected=0)

# Added 2026-10-21 (round 20). How much the 129 are worth, which BOUNDS the claim rather than
# supporting it. All 129 sit at exactly four mismatches, which exclusivity forces, and that is the
# weakest activity stratum. The gate holds the bound so the paper cannot drift into implying these
# are typical sites: if the missed median ever rises above its stratum median, the limitation text
# is wrong and must be rewritten.
ACTIVITY = dict(missed_all_at_4mm=True, stratum_4mm_n=7994, stratum_4mm_median=24,
                missed_median=14, missed_pct_in_stratum=31, missed_at_or_above_stratum_median=37,
                missed_pct_overall=25, missed_in_top_decile=1, missed_max_reads=724)
# Round 27: a single-base bulge inside a homopolymer run has several equally valid placements, so
# the omitted-position classification is only as good as the aligner's gap placement. These hold the
# audit of that. The 129 are forced, not tie-broken, and the reason is exclusivity itself.
PLACEMENT = dict(reconciled=260, unique_terminal_129=129, tied_interior_129=0,
                 tied_interior_131=31, non_ngg_129=0, non_ngg_260=23)

SITES = os.environ.get("SITES129",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "sites129.json"))

COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
def rc(s): return "".join(COMP.get(c, "N") for c in reversed(s))
def mm(a, b): return sum(1 for x, y in zip(a, b) if x != y)

def bulgefree_best(sp, window):
    n, best = len(sp), 99
    for seq in (window, rc(window)):
        for i in range(0, max(0, len(seq) - (n + 3)) + 1):
            w = seq[i:i + n + 3]
            if len(w) < n + 3 or "N" in w: continue
            if not (w[n + 1] == "G" and w[n + 2] == "G"): continue
            best = min(best, mm(sp, w[:n]))
    return best

def main():
    if not os.path.exists(SITES):
        # 2026-10-21: this returned 0 with "SKIP", so the gate guarding the paper's headline passed
        # while checking nothing, and an exit-code sweep could not tell it from a real pass. Its
        # sibling verify_guide_screen.py already exits non-zero on the same absence. A gate whose
        # input is missing reports an INCOMPLETE RELEASE, never a pass.
        print(f"GATE FAILED: site records absent at {SITES}, so the 129-site claim cannot be\n"
              f"  checked here. This is an incomplete release, not a finding about the claim.\n"
              f"  Regenerate from the Lazzarotto 2020 tables, or set SITES129 to the records.")
        return 2
    sites = json.load(open(SITES))
    bad = []
    if len(sites) != EXPECT["at_omitted"]:
        bad.append(f"sites at an omitted position {len(sites)} != {EXPECT['at_omitted']}")
    excl = []
    for h in sites:
        sp, n = h["spacer"], len(h["spacer"])
        core = h["target"][:n - 1]
        if len(core) != n - 1: continue
        emitted = [sp[:i] + sp[i + 1:] for i in range(1, n - 1)]
        if min(mm(core, e) for e in emitted) > 4:
            excl.append(h)
    if len(excl) != EXPECT["exclusively_missed"]:
        bad.append(f"exclusively missed {len(excl)} != {EXPECT['exclusively_missed']}")
    gs = {h["spacer"] for h in excl}
    if len(gs) != EXPECT["guides"]:
        bad.append(f"distinct guides {len(gs)} != {EXPECT['guides']}")
    bf = [bulgefree_best(h["spacer"], h["window"]) for h in excl]
    w4 = sum(1 for b in bf if b <= 4); w6 = sum(1 for b in bf if b <= 6)
    if w4 != EXPECT["bulgefree_within_4"]:
        bad.append(f"bulge-free reachable at <=4mm: {w4} != {EXPECT['bulgefree_within_4']}; "
                   f"the exclusivity claim fails if this is nonzero")
    if w6 != EXPECT["bulgefree_within_6"]:
        bad.append(f"bulge-free reachable at <=6mm: {w6} != {EXPECT['bulgefree_within_6']}")
    rs = sorted(h["reads"] for h in excl)
    if int(rs[0]) != EXPECT["reads_min"] or int(rs[-1]) != EXPECT["reads_max"]:
        bad.append(f"read range {rs[0]:.0f}-{rs[-1]:.0f} != "
                   f"{EXPECT['reads_min']}-{EXPECT['reads_max']}")
    if bad:
        print("GATE FAILED:"); [print("  " + b) for b in bad]; return 2
    print("GATE PASSED: the missed-actives chain reproduces end to end\n")
    print(f"  at an omitted position                 {len(sites)}")
    print(f"  also reachable by an emitted pattern   {len(sites)-len(excl)}")
    print(f"  EXCLUSIVELY MISSED at 4 mismatches     {len(excl)}  over {len(gs)} guides")
    # 2026-10-21: the exclusivity mechanism, tested by count rather than argued. A site reached
    # below budget is reached by an emitted pattern too, so no sub-boundary site can be exclusive.
    below = [h for h in sites if int(float(h["mm"])) < 4]
    at_b  = [h for h in sites if int(float(h["mm"])) == 4]
    excl_ids = {id(h) for h in excl}
    below_excl = [h for h in below if id(h) in excl_ids]
    strat = {}
    for h in sites:
        k = int(float(h["mm"])); strat[k] = strat.get(k, 0) + 1
    print(f"  stratification by deposit mismatch count  "
          + ", ".join(f"{k}:{strat[k]}" for k in sorted(strat)))
    print(f"  below the budget {len(below)} of which exclusively missed {len(below_excl)} "
          f"(must be 0); at the boundary {len(at_b)}")
    if below_excl:
        # 2026-10-22: this appended to `bad`, which is read only by the `if bad:` above, so the
        # append happened after the gate had already printed GATE PASSED and nothing ever looked
        # at it again. Forced true on real data, the gate printed "must be 0", reported 1, and
        # exited 0. It returns here instead, because what it guards is the appendix's stated
        # mechanism for exclusivity arising only at the budget boundary.
        print(f"GATE FAILED: {len(below_excl)} sites below the budget are exclusively missed, so "
              f"exclusivity is NOT confined to the budget boundary and the stated mechanism is "
              f"false"); return 2
    print(f"  reachable bulge-free at <=4 mismatches {w4}   at <=6: {w6}")
    print(f"  nearest bulge-free alignment: median   {sorted(bf)[len(bf)//2]} mismatches")
    print(f"  reads {rs[0]:.0f} to {rs[-1]:.0f}, median {rs[len(rs)//2]:.0f}")
    # ---- placement audit: is the omitted-position classification forced or tie-broken? ----
    def ham(x, y): return sum(1 for u, v in zip(x, y) if u != v)
    rec = tie_int_129 = uniq_129 = tie_int_131 = 0
    nn129 = nn260 = 0
    exclset = {id(h) for h in excl}
    for h in sites:
        sp = h["spacer"]; nn = len(sp); core = h["target"][:nn - 1]; pam = h["target"][nn - 1:]
        costs = [ham(sp[:j] + sp[j + 1:], core) for j in range(nn)]
        best = min(costs); am = [j for j, c in enumerate(costs) if c == best]
        pam_mm = sum(1 for c in pam[1:] if c != "G")
        if best + pam_mm == int(float(h["mm"])): rec += 1
        if not (pam[1] == "G" and pam[2] == "G"): nn260 += 1
        interior = bool(set(am) - {0, nn - 1})
        if id(h) in exclset:
            if interior: tie_int_129 += 1
            elif len(am) == 1: uniq_129 += 1
            if not (pam[1] == "G" and pam[2] == "G"): nn129 += 1
        elif interior:
            tie_int_131 += 1
    P = PLACEMENT
    if rec != P["reconciled"]:
        print(f"GATE FAILED: the deposit's recorded mismatch count reconciles for {rec} of "
              f"{P['reconciled']} rows; the appendix says all of them, so the claim that the "
              f"deposit's gap placement is the mismatch-minimising one no longer holds"); return 2
    if tie_int_129 != P["tied_interior_129"]:
        print(f"GATE FAILED: {tie_int_129} of the exclusively-missed sites now have an interior "
              f"placement of equal cost; the appendix claims 0 and says the terminal placement is "
              f"FORCED by exclusivity. If this is nonzero the robustness argument is broken, not "
              f"merely imprecise"); return 2
    if uniq_129 != P["unique_terminal_129"]:
        print(f"GATE FAILED: unique terminal minimiser holds for {uniq_129} of "
              f"{P['unique_terminal_129']}"); return 2
    if tie_int_131 != P["tied_interior_131"]:
        print(f"GATE FAILED: tied-interior count among the non-exclusive rows is {tie_int_131}, "
              f"appendix says {P['tied_interior_131']}"); return 2
    if nn129 != P["non_ngg_129"] or nn260 != P["non_ngg_260"]:
        print(f"GATE FAILED: non-NGG PAMs: {nn129} in the 129 (appendix says {P['non_ngg_129']}), "
              f"{nn260} in the 260 (appendix says {P['non_ngg_260']})"); return 2
    print(f"  placement: deposit count reconciles on {rec}/{len(sites)} rows; the {len(excl)} "
          f"exclusively missed all have a UNIQUE terminal minimiser")
    print(f"             {tie_int_131} of the {len(sites)-len(excl)} non-exclusive rows do have an "
          f"interior tie; non-NGG PAMs {nn129} of {len(excl)} and {nn260} of {len(sites)}")

    a = DEPOSIT_DERIVED

    # MEASURED, here, from the input file. The at-risk count and the counterexample count are
    # properties of the 260 records this script holds, so it computes them instead of asserting them.
    at_risk = [h for h in sites if h["spacer"][-1] != h["spacer"][-2]]
    rescued_in_file = [h for h in sites if h["spacer"][-1] == h["spacer"][-2]]
    if len(at_risk) != len(sites):
        print(f"GATE FAILED: {len(rescued_in_file)} of {len(sites)} held records have equal last two "
              "spacer bases, so the law says they are RESCUED and they are in the at-risk file. "
              "That is a counterexample to the proposition, not a bookkeeping slip"); return 2
    if a["law_lossless"] + len(at_risk) != a["rna_at_terminal"]:
        print("GATE FAILED: the deposit-derived rescued count plus the MEASURED at-risk count does "
              f"not reach the deposit-derived terminal count: {a['law_lossless']} + {len(at_risk)} "
              f"!= {a['rna_at_terminal']}. One of the three is stale"); return 2
    if a["dna_at_slot_n"] != 0:
        print("GATE FAILED: a measured active now sits at the DNA slot n; the asymmetry claim moves"); return 2
    print()
    print(f"  MEASURED here, from {SITES}: all {len(at_risk)} records are at risk under the "
          f"indicator (g_n != g_(n-1)), over {len({h['spacer'] for h in at_risk})} distinct guides, "
          f"{len(rescued_in_file)} counterexamples")
    print("  ASSERTED from the deposit, NOT measured by this script "
          "(re-earn with make_sites129.py; see DEPOSIT_DERIVED):")
    print(f"    {a['dna_bulge_actives']:,} DNA-bulge actives, "
          f"{a['dna_at_slot_n']} at slot n -> the DNA omission costs nothing measured")
    print(f"    {a['rna_bulge_actives']:,} RNA-bulge actives, "
          f"{a['rna_at_pos_0']} at position 0, {a['rna_at_terminal']} at the terminal")
    print(f"    {a['law_lossless']} rescued by the law, which with the {len(at_risk)} measured above "
          f"accounts for all {a['rna_at_terminal']}")
    c = SCREEN
    # 2026-10-22: these three conditions compared SCREEN's constants to each other and never to the
    # input, so they could fire only if someone edited this file, never because the data moved. The
    # falsifiable half is measurable from the records this script already holds: a guide is
    # protected iff g_n == g_(n-1), and no guide carrying an exclusively-missed site may be one.
    carrying_protected = sorted(g for g in gs if g[-1] == g[-2])
    if carrying_protected:
        print(f"GATE FAILED: {len(carrying_protected)} guide(s) the law protects carry an "
              f"exclusively-missed site, e.g. {carrying_protected[0]}; the proposition is "
              f"falsified, not merely imprecise"); return 2
    if len(gs) != c["carrying"]:
        print(f"GATE FAILED: {len(gs)} guides carry a missed site, SCREEN says "
              f"{c['carrying']}"); return 2
    # panel, exposed and protected need all 110 spacers, which this script does not hold. They are
    # measured by verify_guide_screen.py from the DNA census index; asserted here, not checked.
    if c["exposed"] + c["protected"] != c["panel"]:
        print("GATE FAILED: the screen's two classes do not partition the panel"); return 2
    print(f"  design-time screen: {c['exposed']} exposed / {c['protected']} protected of "
          f"{c['panel']}; all {c['carrying']} guides carrying a missed site are exposed, "
          f"{c['carrying_in_protected']} protected")
    print(f"             necessary not sufficient: {c['exposed']-c['carrying']} exposed guides "
          f"carry none, so it excludes {c['protected']} guides and predicts no others")
    v = ACTIVITY
    # 2026-10-22: line 141 above prints the median of the missed sites' reads computed from the
    # input, and the block below printed ACTIVITY["missed_median"]. Both read 14 and nothing held
    # them together, so a data change would have put two different medians in one report.
    measured_median = int(rs[len(rs) // 2])
    if measured_median != v["missed_median"]:
        print(f"GATE FAILED: the missed sites' median read count measured here is "
              f"{measured_median}, ACTIVITY says {v['missed_median']}; the paper's activity "
              f"bound is quoted from the stale one"); return 2
    # measurable here too: exclusivity forces every missed site to the budget boundary.
    if any(int(float(h["mm"])) != 4 for h in excl):
        print("GATE FAILED: a missed site sits off the mismatch boundary, so the stated mechanism "
              "for its low activity no longer holds"); return 2
    if not v["missed_all_at_4mm"]:
        print("GATE FAILED: the missed set is no longer confined to the mismatch boundary, so the "
              "stated mechanism for its low activity no longer holds"); return 2
    if v["missed_median"] >= v["stratum_4mm_median"]:
        print(f"GATE FAILED: missed median {v['missed_median']} is no longer below its stratum's "
              f"{v['stratum_4mm_median']}; the limitation text overstates the weakness"); return 2
    print(f"  activity bound: all missed sites at 4mm, the weakest stratum (n={v['stratum_4mm_n']:,}, "
          f"median {v['stratum_4mm_median']} reads)")
    print(f"             missed median {v['missed_median']} reads, percentile {v['missed_pct_in_stratum']} "
          f"within that stratum and {v['missed_pct_overall']} overall; "
          f"{v['missed_in_top_decile']} of 129 in the top decile, max {v['missed_max_reads']}")
    print(f"             -> real off-targets, predominantly weak ones; the paper says so")
    return 0

if __name__ == "__main__":
    sys.exit(main())
