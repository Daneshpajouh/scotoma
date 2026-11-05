"""The measurement: how many off-target candidates live in target windows absent from a reference.

The design is a matched-arm comparison within a single genome, which is what makes the result
interpretable:

- **arm A**, the blind arm: every SpCas9 target window in a haplotype assembly whose canonical
  code is absent from the reference universe.
- **arm B**, the visible arm: an equal-sized uniform random sample of that same haplotype's
  windows that *are* present in the reference.

Both arms are scored against the same guide panel, so guide composition cancels exactly. Both
come from the same individual, so population structure cancels. The arms are equal in size, so
the ratio of candidate counts is a ratio of densities with no normalisation step to get wrong.

**Why the arms must be matched on size rather than normalised after the fact.** An earlier
version of this comparison varied the *guides* instead of the targets, pitting a clinical panel
against random genomic protospacers. Clinical guides turned out to be about 15 GC points richer
than random protospacers, so every such comparison measured composition rather than
addressability. Varying the target set while holding the panel fixed is the only axis on which
the question is answerable, and this module only implements that axis.

Nothing here is a statement about cleavage. Arms are counted in *candidate sites by sequence
similarity*. Whether a candidate is cut depends on chromatin, expression, delivery and enzyme
kinetics, none of which are inputs.
"""
from __future__ import annotations

import numpy as np

from .coding import (DEFAULT_PAM, SPACER_LEN, canonical, iter_fasta, oriented_sites,
                     spacer_gc, window_complexity)
from .distance import density_per_million, min_mismatches, threshold_counts
from .universe import contains

DEFAULT_THRESHOLDS = (0, 1, 2, 3, 4, 5, 6)
#: Fraction of reference-present windows retained while streaming. The blind arm is a few per
#: cent of windows, so this leaves the visible pool several times larger than it needs to be
#: while keeping peak memory bounded.
DEFAULT_KEEP_PRESENT = 0.12
#: A haplotype assembly whose windows are not overwhelmingly reference-present is not a human
#: haplotype against a human reference, and the comparison would be meaningless.
MIN_PRESENT_FRACTION = 0.90
#: Cap on the positive control's pool. The control asks only whether the panel's own on-targets
#: are recoverable from reference-present sequence, so a random subsample large enough to expect
#: several hits is as good as the whole pool. With ~110 on-targets among ~296M windows, 40M
#: expects about 15. Set to None to score the entire pool.
CONTROL_MAX_WINDOWS = 40_000_000


def partition_haplotype(fasta_path, universe_codes, seed: int = 0,
                        keep_present: float = DEFAULT_KEEP_PRESENT,
                        pam: str = DEFAULT_PAM):
    """Split a haplotype's target windows into reference-absent and a thinned reference-present pool.

    Streams the assembly one record at a time, so memory is bounded by the kept arrays rather
    than by the genome.

    Parameters
    ----------
    fasta_path : str
        Haplotype assembly, gzipped or not.
    universe_codes : numpy.ndarray
        Sorted unique canonical codes from :func:`scotoma.universe.load_universe`.
    seed : int
        Seed for the thinning draw. Recorded in the returned stats so a run is reproducible.
    keep_present : float
        Retention probability for reference-present windows, in ``(0, 1]``.

    Returns
    -------
    tuple
        ``(absent, present_pool, stats)``. Both arrays hold **oriented** codes with
        multiplicity, as a site that recurs is genuinely several sites.
    """
    if not 0 < keep_present <= 1:
        raise ValueError("keep_present must be in (0, 1]")
    rng = np.random.default_rng(seed)
    absent_parts, present_parts = [], []
    n_windows = 0
    n_present = 0
    n_records = 0
    for _name, bases in iter_fasta(fasta_path):
        sites = oriented_sites(bases, pam=pam)
        if sites.size == 0:
            continue
        n_records += 1
        present_mask = contains(universe_codes, sites)
        n_windows += int(sites.size)
        n_present += int(present_mask.sum())
        absent_parts.append(sites[~present_mask])
        present = sites[present_mask]
        if present.size:
            if keep_present >= 1.0:
                present_parts.append(present)
            else:
                present_parts.append(present[rng.random(present.size) < keep_present])

    absent = np.concatenate(absent_parts) if absent_parts else np.empty(0, np.uint64)
    present_pool = np.concatenate(present_parts) if present_parts else np.empty(0, np.uint64)
    present_fraction = (n_present / n_windows) if n_windows else float("nan")
    stats = {
        "records_with_sites": n_records,
        "windows_total": n_windows,
        "windows_reference_present": n_present,
        "present_fraction": present_fraction,
        "absent_fraction": (1.0 - present_fraction) if n_windows else float("nan"),
        "absent_arm_size": int(absent.size),
        "present_pool_size": int(present_pool.size),
        "keep_present": keep_present,
        "pam": pam,
        "seed": seed,
    }
    return absent, present_pool, stats


def matched_arms(absent, present_pool, seed: int = 0):
    """Draw a visible arm the same size as the blind arm, without replacement.

    Raises
    ------
    ValueError
        If the thinned visible pool is smaller than the blind arm, which means ``keep_present``
        was set too low for this genome. The message says what to raise it to.
    """
    absent = np.asarray(absent, dtype=np.uint64)
    present_pool = np.asarray(present_pool, dtype=np.uint64)
    if present_pool.size < absent.size:
        raise ValueError(
            f"visible pool ({present_pool.size:,}) is smaller than the blind arm "
            f"({absent.size:,}); raise keep_present and rerun"
        )
    rng = np.random.default_rng(seed + 1)
    if present_pool.size == absent.size:
        present_arm = present_pool
    else:
        present_arm = present_pool[rng.choice(present_pool.size, absent.size, replace=False)]
    assert present_arm.size == absent.size
    return absent, present_arm


def stat_matched_arms(absent, present_pool, stat_fn, seed: int = 0, label: str = "statistic",
                      n_bins=None):
    """Draw a visible arm whose distribution of ``stat_fn`` matches the blind arm's, bin by bin.

    The general form of :func:`gc_matched_arms`. ``stat_fn`` maps an array of codes to an integer
    statistic, so the same stratified draw controls GC content, linguistic complexity, or any
    other per-window integer property. Where a bin is short every bin is scaled by the single
    most limiting ratio, so the distributions still match exactly and only the common size falls.

    Returns ``(blind_arm, visible_arm, report)``; the report's ``distributions_identical`` is the
    gate and must be true.
    """
    absent = np.asarray(absent, dtype=np.uint64)
    present_pool = np.asarray(present_pool, dtype=np.uint64)
    rng = np.random.default_rng(seed + 2)

    s_blind = np.asarray(stat_fn(absent))
    s_pool = np.asarray(stat_fn(present_pool))
    if n_bins is None:
        lo = int(min(s_blind.min(initial=0), s_pool.min(initial=0)))
        hi = int(max(s_blind.max(initial=0), s_pool.max(initial=0)))
        bins = np.arange(lo, hi + 1)
    else:
        bins = np.arange(int(n_bins))

    want = np.array([(s_blind == b).sum() for b in bins], dtype=np.int64)
    have = np.array([(s_pool == b).sum() for b in bins], dtype=np.int64)
    needed = want > 0
    if not needed.any():
        raise ValueError("the blind arm is empty, so there is nothing to match")
    with np.errstate(divide="ignore", invalid="ignore"):
        headroom = np.where(needed, have / np.maximum(want, 1), np.inf)
    scale = float(min(1.0, headroom[needed].min()))
    limiting_bin = int(bins[needed][int(np.argmin(headroom[needed]))])
    if scale <= 0:
        raise ValueError(
            f"{label} bin {limiting_bin} is required by the blind arm but absent from the "
            "visible pool, so no matched comparison is possible on this input"
        )

    take = np.floor(want * scale).astype(np.int64)
    blind_idx, visible_idx = [], []
    for b, n_take in zip(bins, take):
        if n_take <= 0:
            continue
        b_pos = np.flatnonzero(s_blind == b)
        p_pos = np.flatnonzero(s_pool == b)
        blind_idx.append(b_pos if n_take == b_pos.size
                         else rng.choice(b_pos, n_take, replace=False))
        visible_idx.append(p_pos if n_take == p_pos.size
                           else rng.choice(p_pos, n_take, replace=False))
    blind_arm = absent[np.concatenate(blind_idx)]
    visible_arm = present_pool[np.concatenate(visible_idx)]

    sb, sv = np.asarray(stat_fn(blind_arm)), np.asarray(stat_fn(visible_arm))
    per_bin_diff = max(abs(int((sb == b).sum()) - int((sv == b).sum())) for b in bins)
    report = {
        "statistic": label,
        "scale": scale,
        "limiting_bin": limiting_bin,
        "arm_size": int(blind_arm.size),
        "mean_blind": float(sb.mean()) if sb.size else float("nan"),
        "mean_visible": float(sv.mean()) if sv.size else float("nan"),
        "max_per_bin_difference": int(per_bin_diff),
        "distributions_identical": bool(per_bin_diff == 0),
    }
    return blind_arm, visible_arm, report


def gc_matched_arms(absent, present_pool, seed: int = 0):
    """Draw a visible arm whose spacer-GC histogram matches the blind arm's, bin by bin.

    This is the control that separates addressability from composition. A plain size-matched draw
    leaves the two arms free to differ in GC content, and GC drives close-match density on its
    own, so a deficit could in principle be a composition artefact. Matching the histogram
    exactly removes that explanation.

    For each GC value in ``0 .. SPACER_LEN`` the blind arm's count is the target. If every bin can
    be served from the pool the arms are identical in GC by construction and the blind arm is used
    whole. If some bin is short, **every** bin is scaled down by the single most limiting ratio,
    so the histograms still match exactly and only the common size falls. Scaling one bin alone
    would preserve the size at the cost of the matching, which would defeat the control.

    Returns
    -------
    tuple
        ``(blind_arm, visible_arm, report)``. The report carries the scale actually applied, the
        limiting GC bin, the per-arm mean spacer GC, and the maximum per-bin discrepancy, which is
        the gate: it must be 0.
    """
    absent = np.asarray(absent, dtype=np.uint64)
    present_pool = np.asarray(present_pool, dtype=np.uint64)
    rng = np.random.default_rng(seed + 2)

    gc_blind = spacer_gc(absent)
    gc_pool = spacer_gc(present_pool)
    bins = np.arange(SPACER_LEN + 1)
    want = np.array([(gc_blind == b).sum() for b in bins], dtype=np.int64)
    have = np.array([(gc_pool == b).sum() for b in bins], dtype=np.int64)

    # The limiting bin is the one with the least headroom among bins the blind arm needs.
    needed = want > 0
    if not needed.any():
        raise ValueError("the blind arm is empty, so there is nothing to match")
    with np.errstate(divide="ignore", invalid="ignore"):
        headroom = np.where(needed, have / np.maximum(want, 1), np.inf)
    scale = float(min(1.0, headroom[needed].min()))
    limiting_bin = int(bins[needed][int(np.argmin(headroom[needed]))])
    if scale <= 0:
        raise ValueError(
            f"GC bin {limiting_bin} is required by the blind arm but absent from the visible "
            "pool, so no GC-matched comparison is possible on this input"
        )

    take = np.floor(want * scale).astype(np.int64)
    blind_idx, visible_idx = [], []
    for b, n_take in zip(bins, take):
        if n_take <= 0:
            continue
        b_pos = np.flatnonzero(gc_blind == b)
        p_pos = np.flatnonzero(gc_pool == b)
        blind_idx.append(b_pos if n_take == b_pos.size
                         else rng.choice(b_pos, n_take, replace=False))
        visible_idx.append(p_pos if n_take == p_pos.size
                           else rng.choice(p_pos, n_take, replace=False))
    blind_arm = absent[np.concatenate(blind_idx)]
    visible_arm = present_pool[np.concatenate(visible_idx)]

    gb, gv = spacer_gc(blind_arm), spacer_gc(visible_arm)
    per_bin_diff = max(abs(int((gb == b).sum()) - int((gv == b).sum())) for b in bins)
    report = {
        "scale": scale,
        "limiting_gc_bin": limiting_bin,
        "arm_size": int(blind_arm.size),
        "mean_spacer_gc_blind": float(gb.mean()) / SPACER_LEN if gb.size else float("nan"),
        "mean_spacer_gc_visible": float(gv.mean()) / SPACER_LEN if gv.size else float("nan"),
        "max_per_bin_difference": int(per_bin_diff),
        "histograms_identical": bool(per_bin_diff == 0),
    }
    return blind_arm, visible_arm, report


def blindspot(absent_arm, present_arm, guide_codes, thresholds=DEFAULT_THRESHOLDS,
              partition_stats=None, control_pool=None):
    """Score both arms against the panel and return the comparison.

    ``control_pool`` is the unsampled reference-present pool, and passing it is strongly
    preferred. The positive control asks whether the pipeline recovers the panel's own
    on-targets from reference-present sequence. Asked of the *visible arm*, that question is
    unreliable: the arm is a uniform subsample matched to the blind arm, so when the blind arm is
    a few per cent of the genome an on-target survives the draw only a few per cent of the time,
    and a correct pipeline fails the control by luck. Asked of the full pool it is a real test.
    When ``control_pool`` is omitted the control falls back to the visible arm and says so in
    ``gates["positive_control_source"]``, so a weak control is never silently reported as a
    strong one.

    Returns
    -------
    dict
        ``counts_blind`` and ``counts_visible`` by threshold, their densities per million
        sites, the ``ratio`` of blind to visible density at each threshold, the per-guide
        ``burden``, and a ``gates`` block.

    The ratio is the headline. A ratio above 1 means reference-absent sequence is enriched for
    guide-proximal candidates relative to reference sequence from the same genome. At or below 1
    means it is not, and the blind spot is a matter of quantity rather than of elevated density.
    """
    absent_arm = np.asarray(absent_arm, dtype=np.uint64)
    present_arm = np.asarray(present_arm, dtype=np.uint64)
    guide_codes = np.asarray(guide_codes, dtype=np.uint64)
    n_guides = int(guide_codes.size)

    best_blind = min_mismatches(absent_arm, guide_codes)
    best_visible = min_mismatches(present_arm, guide_codes)
    counts_blind = threshold_counts(best_blind, thresholds)
    counts_visible = threshold_counts(best_visible, thresholds)

    ratio = {}
    for t in thresholds:
        t = int(t)
        if counts_visible[t] == 0:
            ratio[t] = float("nan")
        else:
            # Arms are equal in size, so a count ratio is already a density ratio. Divide by the
            # size ratio anyway, so the number stays correct if a caller passes unmatched arms.
            size_ratio = (absent_arm.size / present_arm.size) if present_arm.size else float("nan")
            ratio[t] = (counts_blind[t] / counts_visible[t]) / size_ratio if size_ratio else float("nan")

    burden = {}
    for t in thresholds:
        t = int(t)
        burden[t] = (counts_blind[t] / n_guides) if n_guides else float("nan")

    if control_pool is not None:
        control_pool = np.asarray(control_pool, dtype=np.uint64)
        # The control only needs enough pool to expect several on-target hits. Scoring an
        # unthinned pool against the whole panel costs several times the rest of the job: at
        # keep_present=1.0 the pool is ~296M windows against ~12M per arm, so the control
        # dominated the runtime for no extra assurance. Cap it, and say what was used.
        if CONTROL_MAX_WINDOWS and control_pool.size > CONTROL_MAX_WINDOWS:
            rng = np.random.default_rng(12345)
            idx = rng.choice(control_pool.size, CONTROL_MAX_WINDOWS, replace=False)
            control_subset = control_pool[idx]
            control_source = (f"reference-present pool, capped random subsample of "
                              f"{CONTROL_MAX_WINDOWS:,} of {control_pool.size:,}")
        else:
            control_subset = control_pool
            control_source = "full reference-present pool"
        control_hits = int((min_mismatches(control_subset, guide_codes) == 0).sum())
        control_size = int(control_subset.size)
    else:
        control_hits = int(counts_visible.get(0, 0))
        control_source = "visible arm only (weak: see blindspot docstring)"
        control_size = int(present_arm.size)

    gates = {
        "arms_equal_size": bool(absent_arm.size == present_arm.size),
        "panel_non_empty": n_guides > 0,
        "positive_control_source": control_source,
        "positive_control_pool_size": control_size,
        "positive_control_zero_mismatch_hits": control_hits,
        "positive_control_passed": bool(control_hits > 0),
        "blind_arm_has_no_zero_mismatch_hit": bool(counts_blind.get(0, 0) == 0),
    }
    if partition_stats is not None:
        pf = partition_stats.get("present_fraction", float("nan"))
        gates["present_fraction"] = pf
        gates["present_fraction_passed"] = bool(pf == pf and pf > MIN_PRESENT_FRACTION)

    return {
        "n_guides": n_guides,
        "arm_size": int(absent_arm.size),
        "counts_blind": counts_blind,
        "counts_visible": counts_visible,
        "density_blind": density_per_million(counts_blind, absent_arm.size),
        "density_visible": density_per_million(counts_visible, present_arm.size),
        "ratio": ratio,
        "burden_per_guide": burden,
        "gates": gates,
        "partition": partition_stats or {},
    }


def replicate_ratios(absent_arm, present_pool, guide_codes, replicates: int = 5, seed: int = 0,
                     thresholds=DEFAULT_THRESHOLDS):
    """Redraw the visible arm several times to separate draw variance from donor variance.

    A single run reports one ratio per donor, so a range quoted across donors conflates two
    different sources of spread: how much people differ, and how much one person's number moves
    when the visible arm is redrawn. Quoting the first when you measured their sum overstates the
    biological signal.

    Only the visible arm is sampled. The blind arm is every reference-absent window, with no
    randomness in it, so it is scored **once** and reused across replicates. That makes R
    replicates cost about one extra visible-arm scoring each rather than R full runs.

    Returns
    -------
    dict
        Per threshold: the blind count, each replicate's visible count and ratio, and the mean,
        sample standard deviation and range of the ratio over replicates.
    """
    absent_arm = np.asarray(absent_arm, dtype=np.uint64)
    present_pool = np.asarray(present_pool, dtype=np.uint64)
    guide_codes = np.asarray(guide_codes, dtype=np.uint64)
    if replicates < 1:
        raise ValueError("replicates must be at least 1")
    if present_pool.size < absent_arm.size:
        raise ValueError(
            f"visible pool ({present_pool.size:,}) is smaller than the blind arm "
            f"({absent_arm.size:,}); raise keep_present and rerun"
        )

    counts_blind = threshold_counts(min_mismatches(absent_arm, guide_codes), thresholds)

    per_replicate = []
    for r in range(replicates):
        rng = np.random.default_rng(seed + 1000 + r)
        if present_pool.size == absent_arm.size:
            visible = present_pool
        else:
            visible = present_pool[rng.choice(present_pool.size, absent_arm.size, replace=False)]
        counts_visible = threshold_counts(min_mismatches(visible, guide_codes), thresholds)
        ratios = {}
        for t in thresholds:
            t = int(t)
            ratios[t] = (counts_blind[t] / counts_visible[t]) if counts_visible[t] else float("nan")
        per_replicate.append({"seed": seed + 1000 + r, "counts_visible": counts_visible,
                              "ratio": ratios})

    summary = {}
    for t in thresholds:
        t = int(t)
        vals = [rep["ratio"][t] for rep in per_replicate]
        finite = [v for v in vals if v == v]
        if finite:
            mean = sum(finite) / len(finite)
            if len(finite) > 1:
                var = sum((v - mean) ** 2 for v in finite) / (len(finite) - 1)
                sd = var ** 0.5
            else:
                sd = float("nan")
            summary[t] = {"mean": mean, "sd": sd, "min": min(finite), "max": max(finite),
                          "n_finite": len(finite)}
        else:
            summary[t] = {"mean": float("nan"), "sd": float("nan"), "min": float("nan"),
                          "max": float("nan"), "n_finite": 0}

    return {
        "replicates": replicates,
        "arm_size": int(absent_arm.size),
        "pool_size": int(present_pool.size),
        "counts_blind": counts_blind,
        "per_replicate": per_replicate,
        "ratio_summary": summary,
    }


def interpret_gates(result) -> list:
    """Human-readable failures from a :func:`blindspot` result; empty means every gate passed.

    The positive control deserves a word. The visible arm is drawn from reference-present
    sequence, so it should contain each guide's own on-target at distance zero, subject to the
    thinning. Zero such hits means the distance kernel is not finding sites it must find, and
    every other number in the result is then untrustworthy. The blind arm, conversely, must
    contain *no* distance-zero hit: an on-target that is absent from the reference would mean
    the partition is wrong.
    """
    g = result.get("gates", {})
    problems = []
    if not g.get("arms_equal_size", True):
        problems.append("arms are not equal in size; the ratio is not a density ratio")
    if not g.get("panel_non_empty", True):
        problems.append("the guide panel is empty")
    if not g.get("positive_control_passed", True):
        problems.append(
            "positive control failed: no distance-zero hit in "
            f"{g.get('positive_control_source', 'the control set')} "
            f"({g.get('positive_control_pool_size', 0):,} windows), so the pipeline is not "
            "recovering the panel's own on-targets from reference-present sequence. Check that "
            "the panel and the universe come from the same reference"
        )
    if not g.get("blind_arm_has_no_zero_mismatch_hit", True):
        problems.append(
            "a distance-zero hit appears in the blind arm, which should be impossible: "
            "an on-target cannot be absent from the reference the panel was designed on"
        )
    if "present_fraction_passed" in g and not g["present_fraction_passed"]:
        problems.append(
            f"only {g.get('present_fraction', float('nan')):.4f} of windows are "
            f"reference-present, below the {MIN_PRESENT_FRACTION} floor; this does not look "
            "like a human haplotype against a human reference"
        )
    return problems
