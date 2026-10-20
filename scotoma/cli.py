"""Command line interface.

Four subcommands:

``universe``  build a reference's target-window universe once, reuse it everywhere
``scan``      the matched-arm blind-spot measurement for one haplotype and one panel
``rescue``    how much of a haplotype's blind sequence a second reference recovers
``budget``    the validation budget implied by a candidate count

Every subcommand writes machine-readable JSON when asked, and prints a short human summary
either way. Gate failures are reported and set a non-zero exit status rather than being folded
into the numbers.
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from .blindspot import (DEFAULT_KEEP_PRESENT, blindspot, gc_matched_arms, interpret_gates,
                        matched_arms, partition_haplotype, replicate_ratios, stat_matched_arms)
from .budget import as_fraction, budget_rows
from .coding import (DEFAULT_PAM, PAM_SCOPES, iter_fasta, oriented_sites, shuffle_panel,
                     spacer_gc, window_complexity)
from .coding import encode_spacer
from .panel import load_panel, panel_from_windows
from .universe import (build_universe, contains, load_universe,
                       reference_candidate_counts, save_universe)

THRESHOLDS = (0, 1, 2, 3, 4, 5, 6)


def _write_json(path, payload):
    if not path:
        return

    def default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        return str(o)

    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=default)
    print(f"\nwrote {path}")


def cmd_universe(a) -> int:
    codes, mult, meta = build_universe(a.fasta, label=a.label,
                                      primary_only=a.primary_only, pam=a.pam)
    save_universe(a.out_npz, codes, mult, meta)
    print(json.dumps(meta, indent=2))
    print(f"\nwrote {a.out_npz}")
    _write_json(a.out_json, meta)
    return 0


def _check_pam_scope(universe_meta, requested, allow_assume) -> int:
    """Refuse to compare a haplotype extracted under one PAM scope to a universe built under another.

    NGG alone admits one window in sixteen per strand; NGG or NAG admits one in eight. So a site
    count under the wider scope is about twice the narrower one, and partitioning an NGG+NAG
    haplotype against an NGG universe would label every NAG site "reference-absent" when the
    reference simply was not asked about them. That is a denominator error of the exact kind that
    has already produced a three-order-of-magnitude overstatement in this project, so it is a
    hard gate rather than a warning.
    """
    stored = universe_meta.get("pam")
    if stored is None:
        if not allow_assume:
            print("GATE FAIL: this universe records no PAM scope, so it cannot be checked "
                  f"against the requested {requested}. Rebuild it with `scotoma universe`, or pass "
                  "--assume-universe-pam to proceed on your own assertion that the scopes match.",
                  file=sys.stderr)
            return 2
        print(f"WARNING: universe records no PAM scope. Proceeding on the asserted scope "
              f"{requested}; this is unverified.")
        return 0
    if stored != requested:
        print(f"GATE FAIL: universe was built under PAM scope {stored!r} but extraction was "
              f"requested under {requested!r}. These counts are not comparable; rebuild the "
              "universe or change --pam.", file=sys.stderr)
        return 2
    return 0


def cmd_scan(a) -> int:
    universe_codes, universe_meta = load_universe(a.universe)
    print(f"universe {universe_meta.get('label', a.universe)}: "
          f"{universe_codes.size:,} distinct canonical windows, "
          f"PAM scope {universe_meta.get('pam', '<unrecorded>')}")
    rc = _check_pam_scope(universe_meta, a.pam, a.assume_universe_pam)
    if rc:
        return rc

    guides, panel_report = load_panel(a.panel, changeseq=a.changeseq)
    print(f"panel: {panel_report['distinct_guides']} distinct guides from "
          f"{panel_report['rows_considered']} rows"
          + (f", {panel_report['rejected']} rejected" if panel_report["rejected"] else ""))
    if guides.size == 0:
        print("GATE FAIL: the panel is empty", file=sys.stderr)
        return 2
    if a.shuffle_guides:
        from .coding import spacer_gc as _gc
        before = _gc(guides)
        guides = shuffle_panel(guides, seed=a.seed)
        after = _gc(guides)
        if not (before.sum() == after.sum()):
            print("GATE FAIL: shuffling changed the panel's base composition", file=sys.stderr)
            return 2
        print(f"NEGATIVE CONTROL: panel shuffled, composition preserved "
              f"(total spacer GC {int(before.sum())} before and after)")

    absent, present_pool, stats = partition_haplotype(
        a.assembly, universe_codes, seed=a.seed, keep_present=a.keep_present, pam=a.pam)
    print(f"windows: {stats['windows_total']:,}; reference-present "
          f"{stats['windows_reference_present']:,} ({stats['present_fraction']:.4f}); "
          f"blind {stats['absent_arm_size']:,} ({stats['absent_fraction']:.4f})")

    if a.panel_from_blind:
        # The decisive discriminator. A panel drawn from reference-PRESENT windows is both genuine
        # genomic sequence and reference-provenanced, so it cannot separate the two explanations.
        # A panel drawn from the BLIND arm is genuine genomic sequence with the provenance
        # reversed. If the ratio inverts above 1, the effect is a selection relationship between
        # guides and the arm they came from. If it stays below 1, reference-absent sequence is
        # candidate-poor for real genomic guides whichever arm supplied them.
        drawn, draw_report = panel_from_windows(absent, int(guides.size),
                                                seed=a.seed, gc_match_to=guides)
        keep = ~np.isin(absent, drawn)
        removed = int(absent.size - int(keep.sum()))
        absent = absent[keep]
        guides = drawn
        print(f"CONTROL: panel drawn from {draw_report['pool']:,} reference-ABSENT windows, "
              f"{draw_report['n']} guides, GC-matched (max per-bin difference "
              f"{draw_report['max_per_bin_difference']}, mean spacer GC "
              f"{draw_report['mean_spacer_gc_target']:.4f} target against "
              f"{draw_report['mean_spacer_gc_drawn']:.4f} drawn)")
        print(f"         {removed:,} drawn windows removed from the blind arm to avoid "
              f"self-matching; blind arm now {absent.size:,}")
        if draw_report.get("gc_shortfall_by_bin"):
            print(f"         GC shortfall by bin: {draw_report['gc_shortfall_by_bin']}")

    if a.panel_from_reference:
        # Discriminating control: hold "genuine genomic sequence" fixed and vary only which guides.
        # The drawn windows themselves are removed from the visible pool, because a guide is at
        # distance zero from its own window and that self-match is an artifact of the draw. Their
        # near-cognate neighbourhoods are deliberately left in place: that a reference-derived
        # guide has its near-cognates in reference-present sequence is the mechanism under test.
        drawn, draw_report = panel_from_windows(present_pool, int(guides.size),
                                                seed=a.seed, gc_match_to=guides)
        keep = ~np.isin(present_pool, drawn)
        removed = int(present_pool.size - int(keep.sum()))
        present_pool = present_pool[keep]
        guides = drawn
        print(f"CONTROL: panel drawn from {draw_report['pool']:,} reference-present windows, "
              f"{draw_report['n']} guides, GC-matched (max per-bin difference "
              f"{draw_report['max_per_bin_difference']}, mean spacer GC "
              f"{draw_report['mean_spacer_gc_target']:.4f} target against "
              f"{draw_report['mean_spacer_gc_drawn']:.4f} drawn)")
        print(f"         {removed:,} drawn windows removed from the visible pool to avoid "
              f"self-matching; pool now {present_pool.size:,}")
        if draw_report.get("gc_shortfall_by_bin"):
            print(f"         GC shortfall by bin: {draw_report['gc_shortfall_by_bin']}")

    gc_report = None
    try:
        if a.match_on == "gc" or a.gc_matched:
            blind_arm, visible_arm, gc_report = stat_matched_arms(
                absent, present_pool, spacer_gc, seed=a.seed, label="spacer GC")
        elif a.match_on == "complexity":
            blind_arm, visible_arm, gc_report = stat_matched_arms(
                absent, present_pool, window_complexity, seed=a.seed,
                label="distinct 3-mers per window")
        else:
            blind_arm, visible_arm = matched_arms(absent, present_pool, seed=a.seed)
    except ValueError as exc:
        print(f"GATE FAIL: {exc}", file=sys.stderr)
        return 2
    print(f"matched arms: {blind_arm.size:,} windows each")
    if gc_report is not None:
        ident = gc_report.get("distributions_identical", gc_report.get("histograms_identical"))
        print(f"matched on {gc_report.get('statistic', 'spacer GC')}: "
              f"scale {gc_report['scale']:.4f} (limiting bin "
              f"{gc_report.get('limiting_bin', gc_report.get('limiting_gc_bin'))}), "
              f"max per-bin difference {gc_report['max_per_bin_difference']}, mean "
              f"{gc_report.get('mean_blind', gc_report.get('mean_spacer_gc_blind')):.6f} blind "
              f"against {gc_report.get('mean_visible', gc_report.get('mean_spacer_gc_visible')):.6f} "
              f"visible")
        if not ident:
            print("GATE FAIL: the matched distributions are not identical, so the control does "
                  "not hold", file=sys.stderr)
            return 2

    result = blindspot(blind_arm, visible_arm, guides, thresholds=THRESHOLDS,
                       partition_stats=stats,
                       control_pool=None if a.weak_control else present_pool)
    if a.shuffle_guides or a.panel_from_reference or a.panel_from_blind:
        # A synthetic panel has no on-target anywhere in the genome, so the on-target recovery
        # gate is meaningless here and its expectation is inverted: finding a zero-mismatch hit
        # would mean the shuffle had failed. Rewrite the gate rather than fail a run for doing
        # exactly what it was asked to do.
        hits = result["gates"]["positive_control_zero_mismatch_hits"]
        result["gates"]["positive_control_passed"] = bool(hits == 0)
        result["gates"]["positive_control_source"] = (
            "INVERTED for a synthetic panel: it must find no on-target, since a "
            "permuted spacer has none and a drawn window was removed from its own "
            "arm; " + result["gates"]["positive_control_source"])
        result["gates"]["synthetic_panel"] = True
    result["universe"] = universe_meta
    result["panel_report"] = panel_report
    result["assembly"] = a.assembly
    result["gc_matching"] = gc_report

    print(f"\npositive control: {result['gates']['positive_control_zero_mismatch_hits']} "
          f"zero-mismatch hits in the {result['gates']['positive_control_source']} "
          f"({result['gates']['positive_control_pool_size']:,} windows)")
    print(f"\n{'mm':>4} {'blind':>12} {'visible':>12} {'ratio':>10} {'per guide':>12}")
    for t in THRESHOLDS:
        print(f"{t:>4} {result['counts_blind'][t]:>12,} {result['counts_visible'][t]:>12,} "
              f"{result['ratio'][t]:>10.3f} {result['burden_per_guide'][t]:>12.2f}")

    if a.replicates > 1:
        reps = replicate_ratios(blind_arm, present_pool, guides,
                                replicates=a.replicates, seed=a.seed, thresholds=THRESHOLDS)
        result["replicates"] = reps
        print(f"\nvisible-arm redraws, {a.replicates} replicates "
              f"(separates draw variance from donor variance)")
        print(f"{'mm':>4} {'blind':>10} {'ratio mean':>11} {'sd':>8} {'min':>8} {'max':>8}")
        for t in THRESHOLDS:
            srow = reps["ratio_summary"][t]
            print(f"{t:>4} {reps['counts_blind'][t]:>10,} {srow['mean']:>11.4f} "
                  f"{srow['sd']:>8.4f} {srow['min']:>8.4f} {srow['max']:>8.4f}")

    problems = interpret_gates(result)
    _write_json(a.out_json, result)
    if problems:
        print("\nGATE FAILURES", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 2
    print("\nall gates passed")
    return 0


def cmd_rescue(a) -> int:
    first_codes, first_meta = load_universe(a.universe)
    second_codes, second_meta = load_universe(a.second_universe)
    print(f"first  reference {first_meta.get('label', a.universe)}: {first_codes.size:,} windows")
    print(f"second reference {second_meta.get('label', a.second_universe)}: "
          f"{second_codes.size:,} windows")
    for meta, name in ((first_meta, a.universe), (second_meta, a.second_universe)):
        rc = _check_pam_scope(meta, a.pam, a.assume_universe_pam)
        if rc:
            print(f"  (the offending universe is {name})", file=sys.stderr)
            return rc

    n_windows = 0
    n_blind_first = 0
    n_rescued = 0
    for _name, bases in iter_fasta(a.assembly):
        sites = oriented_sites(bases, pam=a.pam)
        if sites.size == 0:
            continue
        n_windows += int(sites.size)
        blind = sites[~contains(first_codes, sites)]
        n_blind_first += int(blind.size)
        if blind.size:
            n_rescued += int(contains(second_codes, blind).sum())

    blind_fraction = (n_blind_first / n_windows) if n_windows else float("nan")
    rescue_rate = (n_rescued / n_blind_first) if n_blind_first else float("nan")
    residual = n_blind_first - n_rescued
    residual_fraction = (residual / n_windows) if n_windows else float("nan")

    payload = {
        "assembly": a.assembly,
        "first_reference": first_meta,
        "second_reference": second_meta,
        "windows_total": n_windows,
        "blind_to_first": n_blind_first,
        "blind_fraction_first": blind_fraction,
        "rescued_by_second": n_rescued,
        "rescue_rate": rescue_rate,
        "residual_blind": residual,
        "residual_blind_fraction": residual_fraction,
    }
    print(f"\nwindows                      {n_windows:,}")
    print(f"blind to first reference     {n_blind_first:,}  ({blind_fraction:.6f})")
    print(f"rescued by second reference  {n_rescued:,}  ({rescue_rate:.4f} of blind)")
    print(f"residual blind               {residual:,}  ({residual_fraction:.6f} of windows)")
    _write_json(a.out_json, payload)
    return 0


def cmd_budget(a) -> int:
    rows = budget_rows(a.n, k_values=a.k, delta1=a.delta1)
    confidence = 1 - as_fraction(a.delta1)
    print(f"N = {a.n:,} candidate sites, delta1 = {a.delta1}\n")
    print(f"{'k':>5} {'m*':>12} {'m*/N':>10}")
    for r in rows:
        print(f"{r['k']:>5} {r['m_star']:>12,} {r['fraction_float']:>10.4f}")
    print("\nA row certifies 'fewer than k cleaved sites among the enumerated candidates' at")
    print(f"confidence {float(confidence):.4f}. It licenses 'none' only if an unclean guide")
    print("cleaves at least k enumerated sites, which is a claim about the enzyme, not the")
    print("sampling. Candidates your reference cannot represent are outside N and outside the")
    print("guarantee; that is what `scotoma scan` measures.")
    _write_json(a.out_json, {
        "n_sites": a.n,
        "delta1": str(a.delta1),
        "rows": [{"k": r["k"], "m_star": r["m_star"], "fraction": str(r["fraction"]),
                  "fraction_float": r["fraction_float"],
                  "miss_probability": str(r["miss_probability"])} for r in rows],
    })
    return 0


def cmd_reference_count(a) -> int:
    """What a conventional reference-anchored nominator would report for one guide."""
    universe_codes, universe_mult, meta = load_universe(a.universe, with_mult=True)
    rc = _check_pam_scope(meta, a.pam, a.assume_universe_pam)
    if rc:
        return rc
    guide = encode_spacer(a.spacer.strip().upper())
    if guide is None:
        print(f"GATE FAIL: {a.spacer!r} is not a 20-base ACGT spacer", file=sys.stderr)
        return 2
    print(f"universe {meta.get('label', a.universe)}: {universe_codes.size:,} canonical windows, "
          f"PAM scope {meta.get('pam', '<unrecorded>')}")
    print(f"spacer {a.spacer.strip().upper()}\n")

    out = reference_candidate_counts(universe_codes, universe_mult, guide,
                                    thresholds=tuple(a.thresholds))
    print(f"{'mm':>4} {'distinct windows':>18} {'genomic sites':>15}")
    for t in a.thresholds:
        print(f"{t:>4} {out['distinct'][t]:>18,} {out['sites'][t]:>15,}")
    print("\n'genomic sites' is the count a validation campaign is sized against: it sums")
    print("both-strand multiplicities, so a window occurring twice counts twice.")
    _write_json(a.out_json, {
        "universe": meta, "spacer": a.spacer.strip().upper(), "pam": a.pam,
        "distinct": out["distinct"], "sites": out["sites"],
    })
    return 0


def build_parser():
    p = argparse.ArgumentParser(
        prog="scotoma",
        description="Off-target candidate enumeration in sequence a reference cannot represent.")
    sub = p.add_subparsers(dest="command", required=True)

    u = sub.add_parser("universe", help="build a reference's target-window universe")
    u.add_argument("--fasta", nargs="+", required=True, help="reference FASTA files")
    u.add_argument("--label", required=True, help='name recorded in the output, e.g. "GRCh38"')
    u.add_argument("--out-npz", required=True)
    u.add_argument("--out-json")
    u.add_argument("--pam", choices=PAM_SCOPES, default=DEFAULT_PAM,
                   help="PAM scope to extract. NGG alone is one window in sixteen per "
                        "strand; NGG+NAG is one in eight, so counts under the two are "
                        "not comparable")
    u.add_argument("--primary-only", action="store_true",
                   help="skip records whose name contains '_' or '.'; a naming convention, not a "
                        "guarantee, and it removes sequence from the universe")
    u.set_defaults(func=cmd_universe)

    s = sub.add_parser("scan", help="matched-arm blind-spot measurement")
    s.add_argument("--assembly", required=True, help="haplotype-resolved assembly FASTA")
    s.add_argument("--universe", required=True, help="reference universe .npz")
    s.add_argument("--panel", required=True, help="guide panel file")
    s.add_argument("--pam", choices=PAM_SCOPES, default=DEFAULT_PAM,
                   help="PAM scope to extract. NGG alone is one window in sixteen per "
                        "strand; NGG+NAG is one in eight, so counts under the two are "
                        "not comparable")
    s.add_argument("--changeseq", action="store_true",
                   help="read the panel as a CHANGE-seq style table rather than one guide per line")
    s.add_argument("--keep-present", type=float, default=DEFAULT_KEEP_PRESENT,
                   help="retention probability for reference-present windows while streaming")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--panel-from-blind", action="store_true",
                   help="the decisive discriminator: draw the panel from the "
                        "reference-ABSENT arm, GC-matched. Real genomic sequence with "
                        "the provenance reversed. A ratio above 1 means the effect is "
                        "a guide-to-source selection relationship; below 1 means "
                        "reference-absent sequence is candidate-poor regardless")
    s.add_argument("--panel-from-reference", action="store_true",
                   help="discriminating control: replace the panel with real "
                        "reference-present windows drawn at random and GC-matched to "
                        "the supplied panel. Separates 'the clinical panel was built "
                        "against this reference' from 'any real genomic guide "
                        "resembles typical genomic sequence'. Use with "
                        "--keep-present 1.0")
    s.add_argument("--shuffle-guides", action="store_true",
                   help="negative control: permute each spacer's bases, preserving "
                        "composition exactly. If the density deficit is a property of "
                        "the arms it must survive; if it vanishes, the deficit is "
                        "specific to the real guide sequences, which matters because "
                        "published panels were designed against the reference")
    s.add_argument("--match-on", choices=["none", "gc", "complexity"], default="none",
                   help="stratify the visible-arm draw on a per-window statistic. "
                        "'gc' matches spacer GC content; 'complexity' matches the "
                        "number of distinct 3-mers per window, which GC matching "
                        "cannot control and which separates addressability from "
                        "repeat composition")
    s.add_argument("--gc-matched", action="store_true",
                   help="draw the visible arm to match the blind arm's spacer-GC "
                        "histogram bin by bin. This is the control that separates "
                        "addressability from composition. Use with --keep-present 1.0, "
                        "since matching every bin needs the unthinned pool")
    s.add_argument("--replicates", type=int, default=1,
                   help="redraw the visible arm this many times and report the ratio's "
                        "mean, sd and range. The blind arm is deterministic and is "
                        "scored once, so replicates are cheap. Use this rather than "
                        "quoting a cross-donor range, which conflates draw variance "
                        "with donor variance")
    s.add_argument("--weak-control", action="store_true",
                   help="run the positive control on the subsampled visible arm instead of the "
                        "full reference-present pool. Cheaper by one pass, but unreliable when "
                        "the blind arm is small, so a correct run can fail it by luck")
    s.add_argument("--assume-universe-pam", action="store_true",
                   help="proceed when a universe records no PAM scope, asserting it "
                        "matches --pam. Unverified; prefer rebuilding the universe")
    s.add_argument("--out-json")
    s.set_defaults(func=cmd_scan)

    r = sub.add_parser("rescue", help="how much blind sequence a second reference recovers")
    r.add_argument("--assembly", required=True)
    r.add_argument("--universe", required=True, help="first reference universe .npz")
    r.add_argument("--second-universe", required=True, help="second reference universe .npz")
    r.add_argument("--pam", choices=PAM_SCOPES, default=DEFAULT_PAM,
                   help="PAM scope to extract. NGG alone is one window in sixteen per "
                        "strand; NGG+NAG is one in eight, so counts under the two are "
                        "not comparable")

    r.add_argument("--assume-universe-pam", action="store_true",
                   help="proceed when a universe records no PAM scope, asserting it "
                        "matches --pam. Unverified; prefer rebuilding the universe")
    r.add_argument("--out-json")
    r.set_defaults(func=cmd_rescue)

    rcp = sub.add_parser("reference-count",
                         help="candidates for one guide as a reference-based nominator sees them")
    rcp.add_argument("--universe", required=True)
    rcp.add_argument("--spacer", required=True, help="20-base spacer")
    rcp.add_argument("--thresholds", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5, 6])
    rcp.add_argument("--pam", choices=PAM_SCOPES, default=DEFAULT_PAM)
    rcp.add_argument("--assume-universe-pam", action="store_true")
    rcp.add_argument("--out-json")
    rcp.set_defaults(func=cmd_reference_count)

    b = sub.add_parser("budget", help="validation budget from a candidate count")
    b.add_argument("--n", type=int, required=True, help="number of enumerated candidate sites")
    b.add_argument("--k", type=int, nargs="+", default=[1, 5, 10, 25])
    b.add_argument("--delta1", default="0.05",
                   help='exact confidence level as a string, e.g. "0.05" or "1/20"')
    b.add_argument("--out-json")
    b.set_defaults(func=cmd_budget)

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
