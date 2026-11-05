#!/usr/bin/env python3 -I
"""Gate for the minimality theorem: the released single-RNA-bulge library is minimal at budget 4.

Verifies the deposited certificate rather than re-searching for it, the same split the CV+ flag
certificate uses: searching is expensive and unverifiable, checking is cheap and total. Every
witness is checked against EVERY other pattern in the library, not a neighbourhood of it.

Also holds the two scoping facts a reviewer will recompute and the direction the design-time screen
errs in. Exits non-zero on any disagreement, and reports an incomplete release if an input is gone.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
# The RELEASED copy first, the working-tree sibling only as a fallback. This read used to point at
# HERE/../producers_2026-10-20, which exists beside code_release in our tree and NOWHERE in the
# release, so the gate passed for us and failed for anyone who cloned the repository, reporting
# "incomplete release" when the release was complete and the path was wrong. Found 2026-10-22 by
# running the released gates from a clean copy of code_release plus certificates and nothing else,
# which is the only way this class of defect shows up. The two copies are byte-identical
# (sha256 8f8c03c4...), so preferring the released one changes no result.
_IDX_RELEASED = os.path.join(HERE, "data", "producers_2026-10-20", "bulge_census_full_idx_rna.json")
_IDX_SIBLING  = os.path.join(HERE, "..", "producers_2026-10-20", "bulge_census_full_idx_rna.json")
# Overridable for the same reason MINCERT and SITES129 below are. Without it, a referee with the
# verifier, the certificate and the index in ONE directory gets GATE FAILED naming a path that does
# not exist, because both fallbacks are project-relative. That is a packaging defect presenting as a
# failed claim, and it is what "our build rechecks it" costs when a reader tries to do it themselves.
IDX  = os.environ.get("RNA1_IDX",
                      _IDX_RELEASED if os.path.exists(_IDX_RELEASED) else _IDX_SIBLING)
CERT = os.environ.get("MINCERT", os.path.join(HERE, "data", "minimality_cert.json"))
SITES = os.environ.get("SITES129", os.path.join(HERE, "data", "sites129.json"))

M = 4
EXPECT = dict(library=1629, guides=110, emitted_current=1560,
              exposed_per_guide=73, exposed_panel_wide=69, shift_rescued=4,
              excl_per_guide=129, excl_per_guide_guides=40,
              excl_panel_wide=96, excl_panel_wide_guides=37)

def dd(a, b): return sum(1 for x, y in zip(a, b) if x != y)

def main():
    for p, what in ((IDX, "the deposited rna1 pattern index"), (CERT, "the minimality certificate"),
                    (SITES, "the measured site records")):
        if not os.path.exists(p):
            print(f"GATE FAILED: {what} is absent at {p}, so the minimality claim cannot be\n"
                  f"  checked here. This is an incomplete release, not a finding about the claim.")
            return 2

    idx = json.load(open(IDX))
    sp = {}
    for pat, owners in idx.items():
        for gi, s, _ in owners: sp[gi] = s
    spacers = [sp[i] for i in sorted(sp)]
    n = len(spacers[0])
    bad = []

    # 1. the library, rebuilt from the spacers rather than taken on trust
    L = sorted({g[:i] + g[i+1:] for g in spacers for i in range(n)})
    if L != sorted({p[:n-1] for p in idx}):
        bad.append("the complete single-deletion library does not reproduce the deposited index")
    if len(L) != EXPECT["library"] or len(spacers) != EXPECT["guides"]:
        bad.append(f"library {len(L)} over {len(spacers)} guides != published "
                   f"{EXPECT['library']}/{EXPECT['guides']}")
    N = len(L)

    # 2. Fact 1: every pattern has a neighbour at distance exactly 1. Lemma 1 rests on this, so a
    #    single exception would void the boundary-shell claim rather than weaken it.
    nn = [min(dd(L[i], L[j]) for j in range(N) if j != i) for i in range(N)]
    if set(nn) != {1}:
        bad.append(f"nearest-neighbour distances are {sorted(set(nn))}, not uniformly 1; Lemma 1 "
                   f"and therefore the boundary-shell claim do not hold as stated")

    # 3. the certificate: every pattern owns a private point, at radius exactly M
    cert = json.load(open(CERT))
    wit = cert["minimality_witnesses"]
    if sorted(wit) != L:
        bad.append("the certificate does not cover exactly the library")
    else:
        radii, margins = set(), set()
        for i, P in enumerate(L):
            x = wit[P]
            r = dd(x, P)
            others = min(dd(x, L[j]) for j in range(N) if j != i)
            radii.add(r); margins.add(others)
            if r > M:
                bad.append(f"witness for {P} is at radius {r} > budget {M}"); break
            if others <= M:
                bad.append(f"witness for {P} is covered by another pattern at distance {others}; "
                           f"it is not private and the theorem fails"); break
        else:
            if radii != {M}:
                bad.append(f"witness radii {sorted(radii)} != {{{M}}}; Lemma 1 forbids anything "
                           f"below the budget, so this is a contradiction, not a loosening")
            # tightness control: a witness at margin > M+1 would mean the search was not finding
            # the boundary and the shell claim is untested rather than confirmed
            if margins != {M+1}:
                bad.append(f"witness margins {sorted(margins)} != {{{M+1}}}: the private points are "
                           f"not on the boundary shell")

    # 4. the constructive blind spot, against the whole library the implementation emits
    W = sorted({g[:i] + g[i+1:] for g in spacers for i in range(0, n-1)})
    if len(W) != EXPECT["emitted_current"]:
        bad.append(f"emitted library {len(W)} != published {EXPECT['emitted_current']}")
    exposed = [g for g in spacers if g[-1] != g[-2]]
    missing = [g for g in exposed if g[:n-1] not in set(W)]
    if len(exposed) != EXPECT["exposed_per_guide"]:
        bad.append(f"exposed per guide {len(exposed)} != {EXPECT['exposed_per_guide']}")
    if len(missing) != EXPECT["exposed_panel_wide"]:
        bad.append(f"exposed panel-wide {len(missing)} != {EXPECT['exposed_panel_wide']}")
    if len(exposed) - len(missing) != EXPECT["shift_rescued"]:
        bad.append(f"shift-rescued guides {len(exposed)-len(missing)} != {EXPECT['shift_rescued']}")
    blind = cert["blind_spot_witnesses"]
    if sorted(blind) != sorted(missing):
        bad.append("the blind-spot certificate does not cover exactly the panel-wide exposed set")
    else:
        for g in missing:
            x = blind[g]; P = g[:n-1]
            if dd(x, P) > M:
                bad.append(f"blind-spot witness for {g} is outside the budget"); break
            near = min(dd(x, q) for q in W)
            if near <= M:
                bad.append(f"blind-spot witness for {g} IS nominated by an emitted pattern at "
                           f"distance {near}; the claim fails for this guide"); break

    # 5. the scoping of the headline, both readings, from the released records
    sites = json.load(open(SITES))
    def excl(comparison_for):
        out = [h for h in sites
               if min(dd(h["target"][:len(h["spacer"])-1], e)
                      for e in comparison_for(h["spacer"])) > M]
        return len(out), len({h["spacer"] for h in out})
    own = excl(lambda s: [s[:i] + s[i+1:] for i in range(1, len(s)-1)])
    pan = excl(lambda s: W)
    if own != (EXPECT["excl_per_guide"], EXPECT["excl_per_guide_guides"]):
        bad.append(f"per-guide exclusivity {own} != published "
                   f"({EXPECT['excl_per_guide']}, {EXPECT['excl_per_guide_guides']})")
    if pan != (EXPECT["excl_panel_wide"], EXPECT["excl_panel_wide_guides"]):
        bad.append(f"panel-wide exclusivity {pan} != published "
                   f"({EXPECT['excl_panel_wide']}, {EXPECT['excl_panel_wide_guides']})")
    if pan[0] >= own[0]:
        bad.append("panel-wide exclusivity is not strictly below per-guide, so the two readings do "
                   "not differ and the appendix table is vacuous")

    if bad:
        print("GATE FAILED:"); [print("  " + b) for b in bad]; return 2

    print("GATE PASSED: minimality of the released library reproduces end to end\n")
    print(f"  library (complete single deletions, 110 guides)   {N}")
    print(f"  every pattern has a neighbour at distance          1   (all {N}, no exceptions)")
    print(f"  => Lemma 1: no private point below radius {M}")
    print(f"  private witness at radius exactly {M}               {len(wit)} of {N}, each verified")
    print(f"     against all {N-1} others; margin exactly {M+1} throughout")
    print(f"  MINIMAL at budget {M}: no pattern is redundant\n")
    print(f"  emitted by the current implementation              {len(W)}")
    print(f"  exposed guides: per guide {len(exposed)}, panel-wide {len(missing)} "
          f"({len(exposed)-len(missing)} covered by a one-base-shift partner)")
    print(f"  explicit sequences it cannot nominate at {M} mismatches   {len(blind)}")
    print(f"  measured actives exclusively missed: per guide {own[0]} over {own[1]} guides,")
    print(f"     panel-wide {pan[0]} over {pan[1]} guides (the floor under the headline)")
    return 0

if __name__ == "__main__":
    sys.exit(main())
