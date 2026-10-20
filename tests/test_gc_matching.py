"""Spacer-GC computation and the GC-matched control.

The control is what separates addressability from composition: GC drives close-match density on
its own, so a deficit between two arms that differ in GC could be a composition artefact. Matching
the histogram exactly removes that explanation, and the gate is that the match is exact.
"""
import numpy as np
import pytest

from scotoma.blindspot import gc_matched_arms, partition_haplotype
from scotoma.coding import (SPACER_LEN, decode_window, encode_spacer, encode_window, spacer_gc)
from scotoma.universe import build_universe, load_universe, save_universe


def _gc_by_decoding(code):
    """Brute-force reference implementation: decode, take the spacer, count G and C."""
    spacer = decode_window(int(code))[:SPACER_LEN]
    return sum(1 for b in spacer if b in "GC")


def test_gc_of_hand_built_spacers():
    assert int(spacer_gc(np.array([encode_spacer("A" * 20)], dtype=np.uint64))[0]) == 0
    assert int(spacer_gc(np.array([encode_spacer("T" * 20)], dtype=np.uint64))[0]) == 0
    assert int(spacer_gc(np.array([encode_spacer("G" * 20)], dtype=np.uint64))[0]) == 20
    assert int(spacer_gc(np.array([encode_spacer("C" * 20)], dtype=np.uint64))[0]) == 20
    assert int(spacer_gc(np.array([encode_spacer("ACGT" * 5)], dtype=np.uint64))[0]) == 10


def test_gc_matches_a_brute_force_decoder(random_sequence):
    """The campaign validated its GC counter against a decoder on 20,000 windows; same idea here."""
    rng = np.random.default_rng(7)
    spacers = ["".join(rng.choice(list("ACGT"), size=20)) for _ in range(3000)]
    codes = np.array([encode_window(s + "AGG") for s in spacers], dtype=np.uint64)
    fast = spacer_gc(codes)
    slow = np.array([_gc_by_decoding(c) for c in codes.tolist()], dtype=np.int8)
    assert np.array_equal(fast, slow)


def test_gc_ignores_the_pam():
    """A GC-rich PAM must not inflate the spacer statistic."""
    spacer = "A" * 20
    a = encode_window(spacer + "AGG")
    b = encode_window(spacer + "CGG")
    codes = np.array([a, b], dtype=np.uint64)
    assert list(int(x) for x in spacer_gc(codes)) == [0, 0]


def test_gc_is_in_range_on_random_windows(random_sequence):
    from scotoma.coding import oriented_sites
    seq = random_sequence(20_000, seed=9)
    gc = spacer_gc(oriented_sites(np.frombuffer(seq.encode(), dtype=np.uint8)))
    assert gc.size > 0
    assert gc.min() >= 0 and gc.max() <= SPACER_LEN


# --- the matched draw ---

def _arms(blind_spacers, pool_spacers):
    blind = np.array([encode_window(s + "AGG") for s in blind_spacers], dtype=np.uint64)
    pool = np.array([encode_window(s + "AGG") for s in pool_spacers], dtype=np.uint64)
    return blind, pool


def test_histograms_are_identical_after_matching():
    # Blind arm: two windows at GC 0 and one at GC 20.
    blind_spacers = ["A" * 20, "T" * 20, "G" * 20]
    # Pool: plenty of both.
    pool_spacers = (["A" * 20] * 5 + ["T" * 20] * 5 + ["G" * 20] * 5 + ["C" * 20] * 5
                    + ["ACGT" * 5] * 5)
    blind, pool = _arms(blind_spacers, pool_spacers)
    b, v, report = gc_matched_arms(blind, pool, seed=0)
    assert report["histograms_identical"] is True
    assert report["max_per_bin_difference"] == 0
    assert report["scale"] == 1.0
    assert b.size == v.size == 3


def test_mean_spacer_gc_is_equal_in_both_arms():
    rng = np.random.default_rng(3)
    blind_spacers = ["".join(rng.choice(list("ACGT"), size=20)) for _ in range(400)]
    pool_spacers = ["".join(rng.choice(list("ACGT"), size=20)) for _ in range(8000)]
    blind, pool = _arms(blind_spacers, pool_spacers)
    b, v, report = gc_matched_arms(blind, pool, seed=1)
    assert report["histograms_identical"] is True
    assert abs(report["mean_spacer_gc_blind"] - report["mean_spacer_gc_visible"]) < 1e-12


def test_a_short_bin_scales_every_bin_so_the_match_stays_exact():
    """Scaling only the short bin would keep the size and break the control, so all bins scale."""
    blind_spacers = ["A" * 20] * 10 + ["G" * 20] * 10
    pool_spacers = ["A" * 20] * 100 + ["G" * 20] * 5      # GC-20 bin is the limit, 5 of 10
    blind, pool = _arms(blind_spacers, pool_spacers)
    b, v, report = gc_matched_arms(blind, pool, seed=0)
    assert report["limiting_gc_bin"] == 20
    assert report["scale"] == pytest.approx(0.5)
    assert report["histograms_identical"] is True
    # Both bins were halved, so 5 + 5, not 10 + 5.
    assert b.size == v.size == 10
    assert int((spacer_gc(b) == 0).sum()) == 5
    assert int((spacer_gc(b) == 20).sum()) == 5


def test_a_missing_bin_is_refused_rather_than_silently_dropped():
    blind_spacers = ["G" * 20] * 4
    pool_spacers = ["A" * 20] * 50          # no GC-20 window at all
    blind, pool = _arms(blind_spacers, pool_spacers)
    with pytest.raises(ValueError, match="no GC-matched comparison is possible"):
        gc_matched_arms(blind, pool, seed=0)


def test_an_empty_blind_arm_is_refused():
    _, pool = _arms([], ["A" * 20] * 10)
    with pytest.raises(ValueError, match="blind arm is empty"):
        gc_matched_arms(np.array([], dtype=np.uint64), pool, seed=0)


def test_matching_is_reproducible_under_a_seed():
    rng = np.random.default_rng(4)
    blind_spacers = ["".join(rng.choice(list("ACGT"), size=20)) for _ in range(200)]
    pool_spacers = ["".join(rng.choice(list("ACGT"), size=20)) for _ in range(4000)]
    blind, pool = _arms(blind_spacers, pool_spacers)
    a1, b1, r1 = gc_matched_arms(blind, pool, seed=5)
    a2, b2, r2 = gc_matched_arms(blind, pool, seed=5)
    assert np.array_equal(a1, a2) and np.array_equal(b1, b2)
    assert r1 == r2


def test_on_a_synthetic_genome_the_control_runs_end_to_end(write_fasta, random_sequence, tmp_path):
    """The GC-matched run needs the unthinned pool, which is why keep_present must be able to be 1."""
    ref_seq = random_sequence(60_000, seed=51)
    ref = write_fasta("ref.fa", [("chr1", ref_seq)])
    hap = write_fasta("hap.fa", [("chr1", ref_seq), ("novel", random_sequence(5_000, seed=52))])
    codes, mult, meta = build_universe([ref], label="R")
    upath = str(tmp_path / "R.npz")
    save_universe(upath, codes, mult, meta)
    universe, _ = load_universe(upath)

    absent, pool, _ = partition_haplotype(hap, universe, keep_present=1.0)
    b, v, report = gc_matched_arms(absent, pool, seed=0)
    assert report["histograms_identical"] is True
    assert b.size == v.size > 0
    assert abs(report["mean_spacer_gc_blind"] - report["mean_spacer_gc_visible"]) < 1e-12
