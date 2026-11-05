"""End-to-end tests of the matched-arm measurement on a synthetic genome we control.

The construction throughout: a "reference" of random sequence, and a "haplotype" that is the
reference plus a novel segment. Windows in the novel segment are reference-absent by
construction, so the partition has a known right answer, and a panel drawn from reference
windows gives the positive control a known right answer too.
"""
import os

import numpy as np
import pytest

from scotoma.blindspot import (blindspot, interpret_gates, matched_arms, partition_haplotype)
from scotoma.coding import decode_window, encode_spacer, oriented_sites
from scotoma.universe import build_universe, contains, load_universe, save_universe


@pytest.fixture
def genome(tmp_path, write_fasta, random_sequence):
    """A reference, a haplotype extending it, and the reference's universe on disk."""
    ref_seq = random_sequence(60_000, seed=21)
    novel_seq = random_sequence(6_000, seed=99)
    ref_path = write_fasta("ref.fa", [("chr1", ref_seq)])
    hap_path = write_fasta("hap.fa", [("chr1", ref_seq), ("novel", novel_seq)])

    codes, mult, meta = build_universe([ref_path], label="REF")
    universe_path = str(tmp_path / "REF.npz")
    save_universe(universe_path, codes, mult, meta)
    return {
        "ref_seq": ref_seq,
        "novel_seq": novel_seq,
        "ref_path": ref_path,
        "hap_path": hap_path,
        "universe_path": universe_path,
    }


def _panel_from_reference(ref_seq, n=12):
    """Guides whose on-targets are reference sequence, so the positive control must fire."""
    sites = oriented_sites(np.frombuffer(ref_seq.encode(), dtype=np.uint8))
    chosen = sites[:: max(1, sites.size // n)][:n]
    return np.array([encode_spacer(decode_window(int(c))[:20]) for c in chosen], dtype=np.uint64)


def test_partition_finds_exactly_the_novel_windows(genome):
    universe, _ = load_universe(genome["universe_path"])
    absent, pool, stats = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)

    novel_sites = oriented_sites(np.frombuffer(genome["novel_seq"].encode(), dtype=np.uint8))
    expected_absent = int((~contains(universe, novel_sites)).sum())

    assert stats["absent_arm_size"] == absent.size
    assert absent.size == expected_absent
    assert stats["windows_total"] == stats["windows_reference_present"] + absent.size


def test_every_blind_window_really_is_absent_and_every_pooled_one_present(genome):
    universe, _ = load_universe(genome["universe_path"])
    absent, pool, _ = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    assert not contains(universe, absent).any()
    assert contains(universe, pool).all()


def test_a_haplotype_identical_to_the_reference_has_an_empty_blind_arm(genome):
    universe, _ = load_universe(genome["universe_path"])
    absent, pool, stats = partition_haplotype(genome["ref_path"], universe, keep_present=1.0)
    assert absent.size == 0
    assert stats["present_fraction"] == 1.0
    assert pool.size == stats["windows_total"]


def test_matched_arms_are_equal_in_size(genome):
    universe, _ = load_universe(genome["universe_path"])
    absent, pool, _ = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    blind, visible = matched_arms(absent, pool, seed=0)
    assert blind.size == visible.size == absent.size


def test_matched_arms_refuses_when_the_visible_pool_is_too_small(genome):
    """Thinning too aggressively must fail loudly, not quietly produce a mismatched ratio."""
    universe, _ = load_universe(genome["universe_path"])
    absent, _, _ = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    tiny_pool = np.array([], dtype=np.uint64)
    with pytest.raises(ValueError, match="smaller than the blind arm"):
        matched_arms(absent, tiny_pool)


def test_visible_arm_is_drawn_without_replacement(genome):
    universe, _ = load_universe(genome["universe_path"])
    absent, pool, _ = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    _, visible = matched_arms(absent, pool, seed=3)
    # Sampling without replacement cannot invent duplicates beyond those already in the pool.
    _, pool_counts = np.unique(pool, return_counts=True)
    _, vis_counts = np.unique(visible, return_counts=True)
    assert vis_counts.max() <= pool_counts.max()


def test_positive_control_fires_and_the_blind_arm_has_no_on_target(genome):
    """The core end-to-end check.

    A panel built from reference windows must find its own on-targets in the visible arm, and
    must find none at distance zero in the blind arm: an on-target cannot be absent from the
    very reference the guide was read off.
    """
    universe, _ = load_universe(genome["universe_path"])
    guides = _panel_from_reference(genome["ref_seq"])
    absent, pool, stats = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    blind, visible = matched_arms(absent, pool, seed=0)
    result = blindspot(blind, visible, guides, partition_stats=stats, control_pool=pool)

    assert result["gates"]["positive_control_zero_mismatch_hits"] > 0
    assert result["gates"]["positive_control_passed"] is True
    assert result["gates"]["positive_control_source"].startswith("full")
    assert result["counts_blind"][0] == 0
    assert result["gates"]["blind_arm_has_no_zero_mismatch_hit"] is True
    assert interpret_gates(result) == []


def test_the_control_on_the_full_pool_finds_every_on_target(genome):
    """Each of the panel's on-targets is reference-present, so all of them must be recovered."""
    universe, _ = load_universe(genome["universe_path"])
    guides = _panel_from_reference(genome["ref_seq"], n=8)
    absent, pool, stats = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    blind, visible = matched_arms(absent, pool, seed=0)
    result = blindspot(blind, visible, guides, partition_stats=stats, control_pool=pool)
    # Every guide was read off a distinct reference window, so each contributes at least one
    # zero-mismatch hit to the unsampled pool.
    assert result["gates"]["positive_control_zero_mismatch_hits"] >= guides.size


def test_the_weak_control_is_weak_and_is_labelled_as_such(genome):
    """Regression test for a defect found by running the CLI, not the library.

    The control was originally computed on the subsampled visible arm. The blind arm here is
    about 9% of the genome, so the visible arm is a 9% draw from the present pool and the
    panel's on-targets usually miss it: a correct pipeline failed the control by luck. The weak
    form is retained behind a flag, and must announce itself so a weak control is never read as
    a strong one.
    """
    universe, _ = load_universe(genome["universe_path"])
    guides = _panel_from_reference(genome["ref_seq"], n=8)
    absent, pool, stats = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    blind, visible = matched_arms(absent, pool, seed=0)

    strong = blindspot(blind, visible, guides, partition_stats=stats, control_pool=pool)
    weak = blindspot(blind, visible, guides, partition_stats=stats)

    assert "weak" in weak["gates"]["positive_control_source"]
    assert weak["gates"]["positive_control_pool_size"] == visible.size
    assert strong["gates"]["positive_control_pool_size"] == pool.size
    assert pool.size > visible.size
    # The strong control sees at least as many on-targets as the subsample can.
    assert (strong["gates"]["positive_control_zero_mismatch_hits"]
            >= weak["gates"]["positive_control_zero_mismatch_hits"])


def test_gates_report_a_failed_positive_control(genome):
    """A panel of guides absent from the reference must trip the control, not pass silently."""
    universe, _ = load_universe(genome["universe_path"])
    alien = np.array([encode_spacer("A" * 20), encode_spacer("C" * 20)], dtype=np.uint64)
    absent, pool, stats = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    blind, visible = matched_arms(absent, pool, seed=0)
    result = blindspot(blind, visible, alien, partition_stats=stats, control_pool=pool)
    assert result["gates"]["positive_control_passed"] is False
    problems = interpret_gates(result)
    assert any("positive control" in p for p in problems)


def test_counts_are_cumulative_in_both_arms(genome):
    universe, _ = load_universe(genome["universe_path"])
    guides = _panel_from_reference(genome["ref_seq"])
    absent, pool, stats = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    blind, visible = matched_arms(absent, pool, seed=0)
    result = blindspot(blind, visible, guides, partition_stats=stats)
    for arm in ("counts_blind", "counts_visible"):
        values = [result[arm][t] for t in (0, 1, 2, 3, 4, 5, 6)]
        assert values == sorted(values)


def test_ratio_is_the_count_ratio_when_arms_are_matched(genome):
    """With equal arms the density ratio must reduce to the plain count ratio."""
    universe, _ = load_universe(genome["universe_path"])
    guides = _panel_from_reference(genome["ref_seq"])
    absent, pool, stats = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    blind, visible = matched_arms(absent, pool, seed=0)
    result = blindspot(blind, visible, guides, partition_stats=stats)
    for t in (4, 5, 6):
        if result["counts_visible"][t]:
            expected = result["counts_blind"][t] / result["counts_visible"][t]
            assert abs(result["ratio"][t] - expected) < 1e-9


def test_burden_is_candidates_per_guide(genome):
    universe, _ = load_universe(genome["universe_path"])
    guides = _panel_from_reference(genome["ref_seq"])
    absent, pool, stats = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    blind, visible = matched_arms(absent, pool, seed=0)
    result = blindspot(blind, visible, guides, partition_stats=stats)
    n = result["n_guides"]
    assert n == guides.size
    for t in (4, 6):
        assert abs(result["burden_per_guide"][t] - result["counts_blind"][t] / n) < 1e-9


def test_seed_makes_the_run_reproducible(genome):
    universe, _ = load_universe(genome["universe_path"])
    guides = _panel_from_reference(genome["ref_seq"])

    def run(seed):
        absent, pool, stats = partition_haplotype(
            genome["hap_path"], universe, seed=seed, keep_present=0.5)
        blind, visible = matched_arms(absent, pool, seed=seed)
        return blindspot(blind, visible, guides, partition_stats=stats)["counts_visible"]

    assert run(7) == run(7)


def test_keep_present_is_validated(genome):
    universe, _ = load_universe(genome["universe_path"])
    for bad in (0.0, -0.1, 1.5):
        with pytest.raises(ValueError, match="keep_present"):
            partition_haplotype(genome["hap_path"], universe, keep_present=bad)


def test_present_fraction_gate_trips_on_an_unrelated_genome(write_fasta, random_sequence,
                                                            tmp_path):
    """A haplotype that is not a variant of the reference must not be silently compared to it."""
    ref = write_fasta("r.fa", [("chr1", random_sequence(20_000, seed=31))])
    other = write_fasta("o.fa", [("chr1", random_sequence(20_000, seed=77))])
    codes, mult, meta = build_universe([ref], label="R")
    upath = str(tmp_path / "R.npz")
    save_universe(upath, codes, mult, meta)
    universe, _ = load_universe(upath)

    absent, pool, stats = partition_haplotype(other, universe, keep_present=1.0)
    assert stats["present_fraction"] < 0.9
    blind, visible = matched_arms(absent, pool[:absent.size] if pool.size >= absent.size else pool,
                                  seed=0) if pool.size >= absent.size else (None, None)
    if blind is not None:
        result = blindspot(blind, visible, np.array([encode_spacer("A" * 20)], dtype=np.uint64),
                           partition_stats=stats)
        assert result["gates"]["present_fraction_passed"] is False
        assert any("floor" in p for p in interpret_gates(result))


def test_control_pool_is_capped_and_says_so(genome, monkeypatch):
    """A capped control must announce the cap, so a reader never mistakes it for the full pool."""
    import importlib
    bs = importlib.import_module("scotoma.blindspot")
    universe, _ = load_universe(genome["universe_path"])
    guides = _panel_from_reference(genome["ref_seq"], n=8)
    absent, pool, stats = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    blind, visible = matched_arms(absent, pool, seed=0)

    monkeypatch.setattr(bs, "CONTROL_MAX_WINDOWS", max(1, pool.size // 2))
    capped = bs.blindspot(blind, visible, guides, partition_stats=stats, control_pool=pool)
    assert "capped random subsample" in capped["gates"]["positive_control_source"]
    assert capped["gates"]["positive_control_pool_size"] == pool.size // 2

    monkeypatch.setattr(bs, "CONTROL_MAX_WINDOWS", None)
    whole = bs.blindspot(blind, visible, guides, partition_stats=stats, control_pool=pool)
    assert whole["gates"]["positive_control_source"] == "full reference-present pool"
    assert whole["gates"]["positive_control_pool_size"] == pool.size
    # Capping may only ever reduce the hit count, never raise it.
    assert (capped["gates"]["positive_control_zero_mismatch_hits"]
            <= whole["gates"]["positive_control_zero_mismatch_hits"])


def test_capping_does_not_change_the_measured_ratio(genome, monkeypatch):
    """The cap touches only the control, never the arms, so no reported ratio may move."""
    import importlib
    bs = importlib.import_module("scotoma.blindspot")
    universe, _ = load_universe(genome["universe_path"])
    guides = _panel_from_reference(genome["ref_seq"], n=8)
    absent, pool, stats = partition_haplotype(genome["hap_path"], universe, keep_present=1.0)
    blind, visible = matched_arms(absent, pool, seed=0)
    monkeypatch.setattr(bs, "CONTROL_MAX_WINDOWS", max(1, pool.size // 3))
    a = bs.blindspot(blind, visible, guides, partition_stats=stats, control_pool=pool)
    monkeypatch.setattr(bs, "CONTROL_MAX_WINDOWS", None)
    b = bs.blindspot(blind, visible, guides, partition_stats=stats, control_pool=pool)
    assert a["counts_blind"] == b["counts_blind"]
    assert a["counts_visible"] == b["counts_visible"]
    assert a["ratio"] == b["ratio"]
