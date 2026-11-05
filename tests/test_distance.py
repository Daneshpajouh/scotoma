"""Distance tests. The kernel must count spacer mismatches and must ignore the PAM."""
import numpy as np

from scotoma.coding import encode_spacer, encode_window
from scotoma.distance import (NO_MATCH, density_per_million, min_mismatches, spacer_mismatches,
                            threshold_counts)


def _hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def test_identical_window_is_distance_zero():
    code = encode_window("ACGT" * 5 + "TGG")
    assert int(spacer_mismatches(np.array([code], dtype=np.uint64), code)[0]) == 0


def test_pam_is_excluded_from_the_distance():
    """Two windows with the same spacer and different PAMs are the same target to a guide."""
    spacer = "ACGTACGTACGTACGTACGT"
    a = encode_window(spacer + "AGG")
    b = encode_window(spacer + "TGG")
    c = encode_window(spacer + "CGG")
    codes = np.array([a, b, c], dtype=np.uint64)
    assert list(spacer_mismatches(codes, a)) == [0, 0, 0]


def test_distance_equals_hamming_on_the_spacer():
    guide_spacer = "A" * 20
    guide = encode_spacer(guide_spacer)
    cases = [
        "A" * 20,
        "C" + "A" * 19,
        "CC" + "A" * 18,
        "ACGTACGTACGTACGTACGT",
        "T" * 20,
    ]
    codes = np.array([encode_window(s + "AGG") for s in cases], dtype=np.uint64)
    observed = list(int(x) for x in spacer_mismatches(codes, guide))
    expected = [_hamming(s, guide_spacer) for s in cases]
    assert observed == expected


def test_distance_is_symmetric_in_its_arguments():
    x = encode_spacer("ACGTACGTACGTACGTACGT")
    y = encode_spacer("ACGTACGTACGTACGTACGA")
    assert int(spacer_mismatches(np.array([y], dtype=np.uint64), x)[0]) == \
           int(spacer_mismatches(np.array([x], dtype=np.uint64), y)[0])


def test_maximum_distance_is_the_spacer_length():
    a = encode_spacer("A" * 20)
    t = encode_spacer("T" * 20)
    assert int(spacer_mismatches(np.array([t], dtype=np.uint64), a)[0]) == 20


def test_min_over_panel_takes_the_closest_guide():
    target = encode_window("A" * 20 + "AGG")
    far = encode_spacer("T" * 20)
    near = encode_spacer("C" + "A" * 19)      # one mismatch from the target
    codes = np.array([target], dtype=np.uint64)
    assert int(min_mismatches(codes, np.array([far], dtype=np.uint64))[0]) == 20
    assert int(min_mismatches(codes, np.array([far, near], dtype=np.uint64))[0]) == 1


def test_empty_panel_gives_no_match():
    codes = np.array([encode_window("A" * 20 + "AGG")], dtype=np.uint64)
    assert int(min_mismatches(codes, np.array([], dtype=np.uint64))[0]) == int(NO_MATCH)


def test_empty_target_set_returns_empty():
    guides = np.array([encode_spacer("A" * 20)], dtype=np.uint64)
    assert min_mismatches(np.array([], dtype=np.uint64), guides).size == 0


def test_chunking_does_not_change_the_answer():
    rng = np.random.default_rng(5)
    spacers = ["".join(rng.choice(list("ACGT"), size=20)) for _ in range(400)]
    codes = np.array([encode_window(s + "AGG") for s in spacers], dtype=np.uint64)
    guides = np.array([encode_spacer(s) for s in spacers[:7]], dtype=np.uint64)
    whole = min_mismatches(codes, guides, chunk=10 ** 9)
    split = min_mismatches(codes, guides, chunk=13)
    assert np.array_equal(whole, split)


def test_threshold_counts_are_cumulative():
    best = np.array([0, 1, 1, 3, 6, 20], dtype=np.int8)
    counts = threshold_counts(best, thresholds=(0, 1, 2, 3, 4, 5, 6))
    assert counts[0] == 1
    assert counts[1] == 3
    assert counts[2] == 3
    assert counts[3] == 4
    assert counts[6] == 5
    for a, b in zip(range(0, 6), range(1, 7)):
        assert counts[a] <= counts[b]


def test_density_scales_by_arm_size():
    counts = {4: 50}
    assert density_per_million(counts, 1_000_000)[4] == 50.0
    assert density_per_million(counts, 500_000)[4] == 100.0


def test_density_of_an_empty_arm_is_nan_not_an_error():
    value = density_per_million({4: 0}, 0)[4]
    assert value != value
