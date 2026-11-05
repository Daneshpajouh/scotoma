"""Visible-arm replicates, which separate draw variance from donor variance.

A single run gives one ratio per donor. Quoting a range across donors then reports the *sum* of
two spreads, how much people differ plus how much one person's number moves on a redraw, as if it
were the first alone. These tests pin that only the visible arm is resampled and that the blind
arm, which has no randomness in it, is identical across replicates.
"""
import numpy as np
import pytest

from scotoma.blindspot import matched_arms, partition_haplotype, replicate_ratios
from scotoma.coding import decode_window, encode_spacer, oriented_sites
from scotoma.universe import build_universe, load_universe, save_universe


@pytest.fixture
def genome(tmp_path, write_fasta, random_sequence):
    ref_seq = random_sequence(80_000, seed=61)
    ref = write_fasta("ref.fa", [("chr1", ref_seq)])
    hap = write_fasta("hap.fa", [("chr1", ref_seq), ("novel", random_sequence(8_000, seed=62))])
    codes, mult, meta = build_universe([ref], label="R")
    upath = str(tmp_path / "R.npz")
    save_universe(upath, codes, mult, meta)
    universe, _ = load_universe(upath)
    sites = oriented_sites(np.frombuffer(ref_seq.encode(), dtype=np.uint8))
    guides = np.array([encode_spacer(decode_window(int(c))[:20])
                       for c in sites[:: max(1, sites.size // 10)][:10]], dtype=np.uint64)
    absent, pool, stats = partition_haplotype(hap, universe, keep_present=1.0)
    return {"absent": absent, "pool": pool, "guides": guides, "stats": stats}


def test_blind_counts_are_identical_across_replicates(genome):
    """The blind arm is not sampled, so nothing about it may move with the seed."""
    a = replicate_ratios(genome["absent"], genome["pool"], genome["guides"],
                         replicates=4, seed=0)
    b = replicate_ratios(genome["absent"], genome["pool"], genome["guides"],
                         replicates=4, seed=999)
    assert a["counts_blind"] == b["counts_blind"]


def test_each_replicate_uses_a_distinct_seed(genome):
    out = replicate_ratios(genome["absent"], genome["pool"], genome["guides"],
                           replicates=5, seed=7)
    seeds = [rep["seed"] for rep in out["per_replicate"]]
    assert len(set(seeds)) == 5


def test_replicate_count_is_respected(genome):
    out = replicate_ratios(genome["absent"], genome["pool"], genome["guides"],
                           replicates=3, seed=0)
    assert out["replicates"] == 3
    assert len(out["per_replicate"]) == 3
    for t in (4, 6):
        assert out["ratio_summary"][t]["n_finite"] <= 3


def test_summary_statistics_bracket_the_replicates(genome):
    out = replicate_ratios(genome["absent"], genome["pool"], genome["guides"],
                           replicates=6, seed=3)
    for t in (5, 6):
        vals = [rep["ratio"][t] for rep in out["per_replicate"] if rep["ratio"][t] == rep["ratio"][t]]
        if len(vals) < 2:
            continue
        s = out["ratio_summary"][t]
        assert s["min"] <= s["mean"] <= s["max"]
        assert abs(s["min"] - min(vals)) < 1e-12
        assert abs(s["max"] - max(vals)) < 1e-12
        assert s["sd"] >= 0


def test_identical_pool_and_arm_sizes_give_zero_spread(genome):
    """With no room to sample, every replicate must be the same draw and the sd must be 0."""
    absent = genome["absent"]
    pool = absent.copy()          # pool exactly the size of the blind arm
    out = replicate_ratios(absent, pool, genome["guides"], replicates=4, seed=0)
    for t in (0, 1, 2, 3, 4, 5, 6):
        s = out["ratio_summary"][t]
        if s["n_finite"] > 1:
            assert s["sd"] == pytest.approx(0.0, abs=1e-12)
            assert s["min"] == pytest.approx(s["max"])


def _same(x, y):
    """Structural equality that treats NaN as equal to NaN.

    Plain ``==`` cannot be used here: a threshold where the visible arm has no candidate yields a
    NaN ratio, and NaN != NaN, so two bit-identical results compare unequal. The first version of
    this test failed for that reason and not because anything was nondeterministic.
    """
    if isinstance(x, dict):
        return isinstance(y, dict) and x.keys() == y.keys() and all(_same(x[k], y[k]) for k in x)
    if isinstance(x, (list, tuple)):
        return len(x) == len(y) and all(_same(u, v) for u, v in zip(x, y))
    if isinstance(x, float) and x != x:
        return isinstance(y, float) and y != y
    return x == y


def test_reproducible_under_the_same_seed(genome):
    a = replicate_ratios(genome["absent"], genome["pool"], genome["guides"], replicates=3, seed=5)
    b = replicate_ratios(genome["absent"], genome["pool"], genome["guides"], replicates=3, seed=5)
    assert _same(a, b)


def test_the_nan_is_where_the_visible_arm_is_empty(genome):
    """Pin the reason the previous test needed NaN-aware comparison, so it is not mistaken later."""
    out = replicate_ratios(genome["absent"], genome["pool"], genome["guides"],
                           replicates=2, seed=5)
    for rep in out["per_replicate"]:
        for t, ratio in rep["ratio"].items():
            if ratio != ratio:
                assert rep["counts_visible"][t] == 0, (
                    f"threshold {t} gave NaN with a non-zero visible count, which would be a bug"
                )


def test_a_too_small_pool_is_refused(genome):
    with pytest.raises(ValueError, match="smaller than the blind arm"):
        replicate_ratios(genome["absent"], np.array([], dtype=np.uint64), genome["guides"])


def test_zero_replicates_is_refused(genome):
    with pytest.raises(ValueError, match="at least 1"):
        replicate_ratios(genome["absent"], genome["pool"], genome["guides"], replicates=0)


def test_single_replicate_reports_no_standard_deviation(genome):
    out = replicate_ratios(genome["absent"], genome["pool"], genome["guides"], replicates=1)
    for t in (4, 6):
        sd = out["ratio_summary"][t]["sd"]
        assert sd != sd or out["ratio_summary"][t]["n_finite"] == 0


def test_the_matched_arms_draw_is_one_of_the_replicate_draws_in_distribution(genome):
    """Sanity: a single matched_arms ratio should sit inside the replicate range, loosely.

    Not an exact identity, since matched_arms uses its own seed offset. The check is only that
    both procedures sample the same pool to the same size.
    """
    blind, visible = matched_arms(genome["absent"], genome["pool"], seed=0)
    assert blind.size == visible.size == genome["absent"].size
    out = replicate_ratios(genome["absent"], genome["pool"], genome["guides"],
                           replicates=4, seed=0)
    assert out["arm_size"] == genome["absent"].size
