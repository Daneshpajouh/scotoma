#!/usr/bin/env python3 -I
"""Producer for the enumeration-scheme lemmas: class invariant, run identity, pattern counts.

Recomputes every published figure in the scheme section from the released pattern indexes and
EXITS NON-ZERO without printing if any disagrees. Run from the RECOMB_2027 directory.
"""
import json, sys
from collections import defaultdict, Counter

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
FILES = {'dna1': 'bulge_census_full_idx_dna.json', 'rna1': 'bulge_census_full_idx_rna.json',
         'dna2': 'bulge2_census_idx_dna2.json',    'rna2': 'bulge2_census_idx_rna2.json',
         'mixed': 'bulgemix_census_idx_mixed.json'}
SIZES = {'dna1': (1, 0), 'rna1': (0, 1), 'dna2': (2, 0), 'rna2': (0, 2), 'mixed': (1, 1)}
EXPECT_LEN = {'dna1': 24, 'rna1': 22, 'dna2': 25, 'rna2': 21, 'mixed': 23}
EXPECT_N   = {'dna1': 2200, 'rna1': 1629, 'dna2': 2200, 'rna2': 1571, 'mixed': 30951}
SHARED     = {'dna1': 0, 'rna1': 5, 'dna2': 0, 'rna2': 0, 'mixed': 95}  # patterns from >1 guide
SPACER_LEN, PAM_LEN = 20, 3

def runs(s):
    return 1 + sum(1 for a, b in zip(s, s[1:]) if a != b)

# The collision law (Proposition: which guides collide, and how many patterns they share).
# Two guides' pattern sets meet iff their spacers are related by a one-base shift, and then the
# intersection is exactly 1 in the single-deletion class and exactly n-1 = 19 in the mixed class.
# Verified exhaustively over all C(110,2) = 5,995 pairs, both directions, 2026-10-21.
SHIFT_PAIRS = 5
ISECT = {'dna1': None, 'dna2': None, 'rna2': None, 'rna1': 1, 'mixed': 19}


def shifted(a, b):
    """One spacer is a single base followed by the other's first n-1 bases."""
    return a[1:] == b[:-1] or b[1:] == a[:-1]


def collision_law(idx, spacer_of):
    """Returns (problems, n_shift_pairs). Checks the law in BOTH directions per class."""
    import itertools
    gs = sorted(spacer_of)
    shift = {tuple(sorted(q)) for q in itertools.combinations(gs, 2)
             if shifted(spacer_of[q[0]], spacer_of[q[1]])}
    probs = []
    if len(shift) != SHIFT_PAIRS:
        probs.append(f"one-base-shift pairs: {len(shift)} != published {SHIFT_PAIRS}")
    for k, want in ISECT.items():
        coll = {}
        for pat, owners in idx[k].items():
            os_ = sorted({o[0] for o in owners})
            for q in itertools.combinations(os_, 2):
                coll[q] = coll.get(q, 0) + 1
        got = set(coll)
        if want is None:
            if got:
                probs.append(f"{k}: expected no collisions, found {len(got)} pair(s)")
            continue
        if got - shift:
            probs.append(f"{k}: colliding pairs that are NOT one-base shifts: {sorted(got - shift)}")
        if shift - got:
            probs.append(f"{k}: one-base-shift pairs that do NOT collide: {sorted(shift - got)}")
        sizes = sorted(set(coll.values()))
        if sizes != [want]:
            probs.append(f"{k}: intersection sizes {sizes} != published [{want}]")
    return probs, len(shift)

# ---------------------------------------------------------------------------------------------
# Added 2026-10-21 (round 7). The signature check above is stated for the CANONICAL UNPADDED form.
# The 2016 wrapper does NOT use that form: Cas-OFFinder needs one fixed pattern length, so the
# wrapper pads every pattern with leading wildcards. In that encoding (length, #wildcards) does
# NOT recover (d, r) -- the no-bulge and DNA-bulge classes share a signature. The paper's remark
# offers itself as the justification published multi-class counts take silently, so the encoding
# those counts were actually computed in has to be checked, not just ours. What summation needs is
# weaker and does hold there: the pattern SETS are pairwise disjoint, separated by wildcard
# placement rather than count. Both facts are gated below; the appendix quotes them.
WRAPPER_SIG_COLLISIONS = {1: {('dna1', 'nobulge')}, 2: {('dna1', 'dna2', 'nobulge')}}

def wrapper_patterns(target, len_pam, bulge_dna, bulge_rna):
    """Pattern generation of hyugel/cas-offinder-bulge, non-reversed branch, transcribed from the
    published source (file `cas-offinder-bulge`, repo read 2026-10-21)."""
    out = defaultdict(set)
    out['nobulge'].add('N' * bulge_dna + target + 'N' * len_pam)
    for b in range(1, bulge_dna + 1):
        for i in range(1, len(target)):
            out[f'dna{b}'].add('N' * (bulge_dna - b) + target[:i] + 'N' * b
                               + target[i:] + 'N' * len_pam)
    for b in range(1, bulge_rna + 1):
        for i in range(1, len(target) - b):
            out[f'rna{b}'].add('N' * (bulge_dna + b) + target[:i] + target[i + b:] + 'N' * len_pam)
    return out


def wrapper_encoding_checks(spacers):
    """The signature must COLLAPSE in the wrapper encoding, and the sets must stay disjoint."""
    probs = []
    for budget, expected in WRAPPER_SIG_COLLISIONS.items():
        for sp in spacers:
            pats = wrapper_patterns(sp, PAM_LEN, budget, budget)
            by_sig = defaultdict(set)
            for cls, S in pats.items():
                shapes = {(len(x), x.count('N')) for x in S}
                if len(shapes) != 1:
                    probs.append(f"wrapper b={budget} {sp[:6]}..: {cls} not one shape"); continue
                by_sig[next(iter(shapes))].add(cls)
            got = {tuple(sorted(c)) for c in by_sig.values() if len(c) > 1}
            want = {tuple(sorted(c)) for c in expected}
            if got != want:
                probs.append(f"wrapper b={budget} {sp[:6]}..: signature collisions {sorted(got)} "
                             f"!= expected {sorted(want)}")
            names = sorted(pats)
            for a in range(len(names)):
                for bb in range(a + 1, len(names)):
                    inter = pats[names[a]] & pats[names[bb]]
                    if inter:
                        probs.append(f"wrapper b={budget} {sp[:6]}..: {names[a]} and {names[bb]} "
                                     f"share {len(inter)} pattern(s); sets must be disjoint")
    return probs


# Added 2026-10-21 (round 8). Closed form for the 2016 wrapper's single-bulge pattern deficit.
# The paper previously asserted BOTH RNA terminals are lossy omissions. That is false: a terminal
# deletion coincides with an interior one exactly when the terminal base lies in a run with its
# neighbour, so the omission is lossless for 40/110 guides at the start and 37/110 at the end.
# Gated here because the claim is now a law rather than a one-panel measurement.
DEFICIT_LOSSLESS_START, DEFICIT_LOSSLESS_END, DEFICIT_LOSSLESS_BOTH = 40, 37, 12
DEFICIT_PANEL_FULL, DEFICIT_PANEL_WRAPPER = 3834, 3581
# third-party validation: issue #2 of hyugel/cas-offinder-bulge, guide outside our panel
ISSUE2_GUIDE, ISSUE2_DNA, ISSUE2_RNA = "GGCCGACCTGTCGCTGACGC", 19, 16


def deficit_law(spacers):
    """wrapper DNA = n-1; wrapper RNA = r(g) - 1[g1!=g2] - 1[gn!=g(n-1)]. Checked per guide."""
    probs = []
    ls = le = both = 0
    full = wrap = 0
    for g in list(spacers) + [ISSUE2_GUIDE]:
        pats = wrapper_patterns(g, PAM_LEN, 1, 1)
        ds = 1 if g[0] != g[1] else 0
        de = 1 if g[-1] != g[-2] else 0
        pred_rna, pred_dna = runs(g) - ds - de, len(g) - 1
        if len(pats['rna1']) != pred_rna:
            probs.append(f"deficit law RNA: {g} actual {len(pats['rna1'])} != {pred_rna}")
        if len(pats['dna1']) != pred_dna:
            probs.append(f"deficit law DNA: {g} actual {len(pats['dna1'])} != {pred_dna}")
        if g == ISSUE2_GUIDE:
            if len(pats['dna1']) != ISSUE2_DNA or len(pats['rna1']) != ISSUE2_RNA:
                probs.append(f"issue #2 guide: ({len(pats['dna1'])},{len(pats['rna1'])}) != "
                             f"published ({ISSUE2_DNA},{ISSUE2_RNA})")
            continue
        ls += (ds == 0); le += (de == 0); both += (ds == 0 and de == 0)
        full += len(g) + runs(g); wrap += pred_dna + pred_rna
    for name, got, want in (("lossless start", ls, DEFICIT_LOSSLESS_START),
                            ("lossless end", le, DEFICIT_LOSSLESS_END),
                            ("lossless both", both, DEFICIT_LOSSLESS_BOTH),
                            ("panel full patterns", full, DEFICIT_PANEL_FULL),
                            ("panel wrapper patterns", wrap, DEFICIT_PANEL_WRAPPER)):
        if got != want:
            probs.append(f"deficit law {name}: {got} != published {want}")
    # non-vacuity control: the indicator terms must actually vary across the panel
    vals = {1 + (0 if g[0] == g[1] else 1) + (0 if g[-1] == g[-2] else 1) for g in spacers}
    if len(vals) < 2:
        probs.append(f"deficit law VACUOUS: only one deficit value on the panel ({vals})")
    return probs



# Added 2026-10-21 (round 10). The deficit law generalises to every bulge size, with a proof that
# uses no run structure: every wrapper RNA pattern begins with g_1 and ends with g_n, so the omitted
# leading word (which begins with g_{b+1}) and trailing word (which ends with g_{n-b}) are present
# exactly when those characters repeat. Gated over b = 1..10 with the converse checked per case.
DEFICIT_B_RANGE = range(1, 11)

def deficit_law_general_b(spacers):
    probs = []
    for b in DEFICIT_B_RANGE:
        vals = set()
        for g in spacers:
            n = len(g)
            wr = {g[:i] + g[i+b:] for i in range(1, n - b)}
            fr = {g[:i] + g[i+b:] for i in range(0, n - b + 1)}
            wd = {g[:i] + 'N'*b + g[i:] for i in range(1, n)}
            fd = {g[:i] + 'N'*b + g[i:] for i in range(1, n + 1)}
            ds = 0 if g[0] == g[b] else 1
            de = 0 if g[n-1] == g[n-1-b] else 1
            if len(wr) != len(fr) - ds - de:
                probs.append(f"deficit b={b} {g}: RNA {len(wr)} != {len(fr)}-{ds}-{de}")
            if len(wd) != n - 1 or len(fd) != n:
                probs.append(f"deficit b={b} {g}: DNA {len(wd)}/{len(fd)} != {n-1}/{n}")
            # converse control: omitted terminal words genuinely absent iff indicator is 1
            if ((g[b:] not in wr) + (g[:n-b] not in wr)) != ds + de:
                probs.append(f"deficit b={b} {g}: converse fails")
            # the two structural claims the proof rests on
            if not all(c[0] == g[0] and c[-1] == g[n-1] for c in wr):
                probs.append(f"deficit b={b} {g}: a wrapper pattern breaks the flank claim")
            vals.add(1 + ds + de)
        if len(vals) < 2:
            probs.append(f"deficit b={b}: VACUOUS, one deficit value only ({vals})")
    return probs


# Added 2026-10-21 (round 14). Proposition 1 carries no PAM-orientation clause, but every earlier
# check used only the wrapper's forward-PAM branch. The reversed branch (leading PAM, Cas12a family)
# pads on the other side, and there the DEGENERATE SLOT FLIPS: slot 0 for a trailing PAM, slot n for
# a leading one. The count is unchanged, so the proposition stands, but the appendix proof had to be
# generalised. Gated here over both orientations and three spacer lengths so neither can regress.
ORIENT_LENGTHS = (20, 23, 24)     # 23/24 cover Cas12a-family spacers
ORIENT_B = (1, 2, 3)

def _patterns(g, len_pam, bd, b, reversed_pam):
    if reversed_pam:
        nob = "N"*len_pam + g + "N"*bd
        slot = lambda i: "N"*len_pam + g[:i] + "N"*b + g[i:] + "N"*(bd-b)
        rna  = lambda i: "N"*len_pam + g[:i] + g[i+b:] + "N"*(bd+b)
    else:
        nob = "N"*bd + g + "N"*len_pam
        slot = lambda i: "N"*(bd-b) + g[:i] + "N"*b + g[i:] + "N"*len_pam
        rna  = lambda i: "N"*(bd+b) + g[:i] + g[i+b:] + "N"*len_pam
    n = len(g)
    return (nob,
            {slot(i) for i in range(1, n)},                  # wrapper DNA
            {slot(i) for i in range(0, n+1)} - {nob},        # full DNA, degenerate slot removed
            {rna(i) for i in range(1, n-b)},                 # wrapper RNA
            {rna(i) for i in range(0, n-b+1)})               # full RNA

def orientation_checks(spacers):
    probs = []
    for n in ORIENT_LENGTHS:
        panel = [(g*2)[:n] for g in spacers]
        for reversed_pam in (False, True):
            tag = "reversed" if reversed_pam else "forward"
            for b in ORIENT_B:
                degen_slots = set()
                for g in panel:
                    nob, dw, df, rw, rf = _patterns(g, PAM_LEN, b, b, reversed_pam)
                    ds = 0 if g[0] == g[b] else 1
                    de = 0 if g[n-1] == g[n-1-b] else 1
                    if len(df) != n:
                        probs.append(f"orient {tag} n={n} b={b}: full DNA {len(df)} != {n}")
                    if len(dw) != n - 1:
                        probs.append(f"orient {tag} n={n} b={b}: wrapper DNA {len(dw)} != {n-1}")
                    if len(rw) != len(rf) - ds - de:
                        probs.append(f"orient {tag} n={n} b={b}: RNA {len(rw)} != {len(rf)}-{ds}-{de}")
                    # which slot is degenerate, recorded so the flip itself is gated
                    s0 = _patterns(g, PAM_LEN, b, b, reversed_pam)
                    degen_slots.add(0 if nob == (("N"*PAM_LEN + g[:0] + "N"*b + g[0:] + "N"*0)
                                                if reversed_pam else
                                                ("N"*0 + g[:0] + "N"*b + g[0:] + "N"*PAM_LEN)) else n)
                want = {n} if reversed_pam else {0}
                if degen_slots != want:
                    probs.append(f"orient {tag} n={n} b={b}: degenerate slot {degen_slots} != {want}")
    return probs


def main():
    bad, idx = [], {}
    for k, f in FILES.items():
        idx[k] = json.load(open(f"{D}/{f}"))
        shapes = Counter((len(p), p.count('N')) for p in idx[k])
        if len(shapes) != 1:
            bad.append(f"{k}: patterns are not a single shape: {dict(shapes)}"); continue
        (L, N) = next(iter(shapes))
        d, r = SIZES[k]
        if L != SPACER_LEN + PAM_LEN + d - r: bad.append(f"{k}: length {L} != 23+{d}-{r}")
        if N != PAM_LEN + d:                  bad.append(f"{k}: wildcards {N} != 3+{d}")
        if L != EXPECT_LEN[k]:                bad.append(f"{k}: length {L} != published {EXPECT_LEN[k]}")
        if len(idx[k]) != EXPECT_N[k]:        bad.append(f"{k}: {len(idx[k])} patterns != published {EXPECT_N[k]}")
    # Guide identifiability: cross-guide pattern sharing, per class
    for k in FILES:
        sh = sum(1 for rows in idx[k].values() if len({r[1] for r in rows}) > 1)
        if sh != SHARED[k]:
            bad.append(f"{k}: {sh} patterns shared between guides != published {SHARED[k]}")
    # Lemma 1: (length, wildcards) separates all five classes AND the no-bulge pattern (23, 3)
    sigs = {k: (EXPECT_LEN[k], PAM_LEN + SIZES[k][0]) for k in FILES}
    sigs['nobulge'] = (SPACER_LEN + PAM_LEN, PAM_LEN)
    if len(set(sigs.values())) != len(sigs):
        bad.append(f"class signatures are NOT distinct: {sigs}")
    # Lemma 2: distinct single-RNA-deletion patterns per guide == homopolymer runs in its spacer
    spacers = sorted({row[1] for rows in idx['dna1'].values() for row in rows})
    per_guide = defaultdict(set)
    for pat, rows in idx['rna1'].items():
        for row in rows: per_guide[row[1]].add(pat)
    held = sum(1 for sp in spacers if runs(sp) == len(per_guide.get(sp, ())))
    if held != len(spacers):
        bad.append(f"run identity holds for only {held} of {len(spacers)} guides")
    if len(spacers) != 110:
        bad.append(f"expected 110 guides, got {len(spacers)}")
    shared = sum(len(v) for v in per_guide.values()) - len(idx['rna1'])
    if shared != 5:
        bad.append(f"patterns shared between guides = {shared}, published 5")
    spacer_of = {row[0]: row[1] for rows in idx['dna1'].values() for row in rows}
    law_problems, n_shift = collision_law(idx, spacer_of)
    bad.extend(law_problems)
    # Lemma 1b: in the WRAPPER's padded encoding the signature collapses, sets stay disjoint
    bad.extend(wrapper_encoding_checks(spacers))
    # Lemma 4: closed form for the wrapper's single-bulge deficit, + third-party validation
    bad.extend(deficit_law(spacers))
    # Lemma 5: the same law for every bulge size, with the converse and flank claims
    bad.extend(deficit_law_general_b(spacers))
    # Lemma 6: both PAM orientations and Cas12a spacer lengths, incl. the degenerate-slot flip
    bad.extend(orientation_checks(spacers))

    if bad:
        print("GATE FAILED:"); [print("  " + b) for b in bad]; return 2
    print("GATE PASSED: all scheme lemmas and pattern counts reproduced\n")
    print(f"{'class':<9}{'patterns':>9}{'len':>5}{'#N':>4}{'constrained':>12}")
    for k in ('dna2', 'dna1', 'mixed', 'rna1', 'rna2'):
        L, N = sigs[k]; print(f"{k:<9}{len(idx[k]):>9,}{L:>5}{N:>4}{L-N:>12}")
    L, N = sigs['nobulge']; print(f"{'nobulge':<9}{'-':>9}{L:>5}{N:>4}{L-N:>12}")
    print(f"\nLemma 1  (length, wildcards) distinct across all six classes: True (UNPADDED form)")
    print(f"Lemma 1b in the 2016 wrapper's padded encoding the pair does NOT separate: "
          f"nobulge/dna1 collide at a one-base budget and nobulge/dna1/dna2 at two,")
    print(f"         checked on all {len(spacers)} spacers; pattern SETS pairwise disjoint in "
          f"both budgets, so summing per-class counts stays licensed")
    print(f"Lemma 2  distinct RNA-deletion patterns == homopolymer runs: {held}/{len(spacers)} guides")
    print(f"         mean runs per 20-mer {sum(runs(s) for s in spacers)/len(spacers):.2f}"
          f"  vs random-sequence 1+19(3/4) = {1+19*0.75:.2f}")
    print(f"         patterns shared between two guides: {shared}")
    print(f"Lemma 3  colliding guide pairs == one-base-shift pairs, both directions: "
          f"{n_shift} pairs of {len(spacer_of)*(len(spacer_of)-1)//2}")
    print(f"         intersection sizes exactly 1 (rna1) and {ISECT['mixed']} (mixed); "
          f"zero in every DNA class")
    print(f"Lemma 4  wrapper deficit = 1 + 1[g1!=g2] + 1[gn!=g(n-1)] patterns, exact on "
          f"{len(spacers)}/{len(spacers)} guides")
    print(f"         terminal omission lossless: start {DEFICIT_LOSSLESS_START}, "
          f"end {DEFICIT_LOSSLESS_END}, both {DEFICIT_LOSSLESS_BOTH} of {len(spacers)} guides")
    print(f"         panel-wide {DEFICIT_PANEL_WRAPPER:,} of {DEFICIT_PANEL_FULL:,} patterns, "
          f"deficit {DEFICIT_PANEL_FULL-DEFICIT_PANEL_WRAPPER} = "
          f"{100*(DEFICIT_PANEL_FULL-DEFICIT_PANEL_WRAPPER)/DEFICIT_PANEL_FULL:.2f}% of patterns")
    print(f"         predicts third-party output in wrapper issue #2: "
          f"DNA {ISSUE2_DNA}, RNA {ISSUE2_RNA} on a guide outside the panel")
    print(f"Lemma 5  the same deficit law at b = {DEFICIT_B_RANGE.start}..{DEFICIT_B_RANGE.stop-1}: "
          f"{len(spacers)*len(DEFICIT_B_RANGE):,} cases, zero disagreements")
    print(f"         converse holds per case, and every wrapper RNA pattern begins with g_1 and "
          f"ends with g_n, which is what the proof uses")
    print(f"Lemma 6  the deficit holds for BOTH PAM orientations at n = "
          f"{', '.join(map(str, ORIENT_LENGTHS))} and b = {', '.join(map(str, ORIENT_B))}")
    print(f"         the degenerate slot FLIPS with orientation: slot 0 for a trailing PAM, "
          f"slot n for a leading one; the count does not change")
    return 0

if __name__ == '__main__':
    sys.exit(main())
