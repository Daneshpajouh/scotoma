"""Reference-side candidate counting, the incumbent quantity.

This is what a conventional reference-anchored nominator reports and therefore the `N` a
validation campaign is sized against. The tests pin the orientation recovery, which is the only
subtle part: a universe stores canonical codes and has lost the strand, so the admissible oriented
forms must be recovered before any distance is taken.
"""
import numpy as np
import pytest

from scotoma.coding import decode_window, encode_spacer, encode_window, oriented_sites, revcomp64
from scotoma.universe import build_universe, reference_candidate_counts


def _bases(text):
    return np.frombuffer(text.encode(), dtype=np.uint8)


def _universe_from(seq, pam="NGG"):
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".fa")
    with os.fdopen(fd, "w") as fh:
        fh.write(">chr1\n" + seq + "\n")
    return build_universe([path], label="t", pam=pam)


def test_a_guides_own_on_target_is_found_at_distance_zero():
    spacer = "ACGTACGTACGTACGTACGT"
    seq = spacer + "AGG"
    codes, mult, _ = _universe_from(seq)
    out = reference_candidate_counts(codes, mult, encode_spacer(spacer), thresholds=(0, 1, 4))
    assert out["distinct"][0] == 1
    assert out["sites"][0] == 1


def test_orientation_is_recovered_for_a_reverse_strand_site():
    """The decisive case: the stored canonical code is not the oriented one.

    ``T*20 + TGG`` is a forward site whose reverse complement sorts lower, so the universe stores
    the reverse complement. Counting distances on the stored code would miss the guide entirely.
    """
    spacer = "T" * 20
    seq = spacer + "TGG"
    codes, mult, _ = _universe_from(seq)
    # Confirm the premise: what is stored is NOT the oriented form.
    oriented = oriented_sites(_bases(seq))[0]
    assert int(codes[0]) != int(oriented)
    assert int(codes[0]) == int(revcomp64(oriented))

    out = reference_candidate_counts(codes, mult, encode_spacer(spacer), thresholds=(0,))
    assert out["distinct"][0] == 1, "orientation was not recovered"


def test_a_window_with_no_valid_pam_in_either_orientation_is_never_counted():
    """Such a window cannot be in a universe, but the guard must hold if one is passed in."""
    bogus = encode_window("A" * 20 + "AAA")       # no GG in either orientation
    codes = np.array(sorted([bogus]), dtype=np.uint64)
    mult = np.array([1], dtype=np.int32)
    out = reference_candidate_counts(codes, mult, encode_spacer("A" * 20), thresholds=(0, 6, 20))
    assert out["distinct"][20] == 0


def test_counts_are_cumulative_and_monotone():
    rng = np.random.default_rng(2)
    seq = "".join(rng.choice(list("ACGT"), size=30_000))
    codes, mult, _ = _universe_from(seq)
    guide = encode_spacer("".join(rng.choice(list("ACGT"), size=20)))
    out = reference_candidate_counts(codes, mult, guide)
    dvals = [out["distinct"][t] for t in (0, 1, 2, 3, 4, 5, 6)]
    svals = [out["sites"][t] for t in (0, 1, 2, 3, 4, 5, 6)]
    assert dvals == sorted(dvals)
    assert svals == sorted(svals)
    # Sites counts multiplicity, so it can never be below the distinct count.
    for d, s in zip(dvals, svals):
        assert s >= d


def test_sites_counts_multiplicity_and_distinct_does_not():
    spacer = "ACGTACGTACGTACGTACGT"
    window = spacer + "AGG"
    seq = window + "TTTTTTTTTTTT" + window        # the same window twice
    codes, mult, meta = _universe_from(seq)
    assert codes.size == 1 and int(mult[0]) == 2
    out = reference_candidate_counts(codes, mult, encode_spacer(spacer), thresholds=(0,))
    assert out["distinct"][0] == 1
    assert out["sites"][0] == 2


def test_chunking_does_not_change_the_answer():
    rng = np.random.default_rng(6)
    seq = "".join(rng.choice(list("ACGT"), size=60_000))
    codes, mult, _ = _universe_from(seq)
    guide = encode_spacer("".join(rng.choice(list("ACGT"), size=20)))
    whole = reference_candidate_counts(codes, mult, guide, chunk=10 ** 9)
    split = reference_candidate_counts(codes, mult, guide, chunk=777)
    assert whole == split


def test_agrees_with_a_brute_force_scan_of_the_same_sequence():
    """Independent check: enumerate oriented sites directly and score them, bypassing the universe.

    This is the cross-validation that matters, because it does not share the orientation-recovery
    code path with the function under test.
    """
    from scotoma.distance import min_mismatches
    rng = np.random.default_rng(14)
    seq = "".join(rng.choice(list("ACGT"), size=120_000))
    codes, mult, _ = _universe_from(seq)
    guide = encode_spacer("".join(rng.choice(list("ACGT"), size=20)))

    via_universe = reference_candidate_counts(codes, mult, guide)
    oriented = oriented_sites(_bases(seq))
    best = min_mismatches(oriented, np.array([guide], dtype=np.uint64))
    for t in (0, 1, 2, 3, 4, 5, 6):
        brute_sites = int((best <= t).sum())
        assert via_universe["sites"][t] == brute_sites, (
            f"threshold {t}: universe path {via_universe['sites'][t]} against brute force "
            f"{brute_sites}"
        )


def test_mult_may_be_omitted():
    spacer = "ACGTACGTACGTACGTACGT"
    codes, _, _ = _universe_from(spacer + "AGG")
    out = reference_candidate_counts(codes, None, encode_spacer(spacer), thresholds=(0,))
    assert out["distinct"][0] == out["sites"][0] == 1


# --- panel-wide path, and the metadata robustness fix ---

def test_panel_counts_match_the_single_guide_function_guide_for_guide():
    """The fast panel path must agree exactly with the one-guide path it optimises."""
    from scotoma.universe import panel_reference_counts
    rng = np.random.default_rng(21)
    seq = "".join(rng.choice(list("ACGT"), size=80_000))
    codes, mult, _ = _universe_from(seq)
    guides = np.array([encode_spacer("".join(rng.choice(list("ACGT"), size=20)))
                       for _ in range(6)], dtype=np.uint64)

    panel = panel_reference_counts(codes, mult, guides, chunk=3_000)
    for gi, g in enumerate(guides.tolist()):
        one = reference_candidate_counts(codes, mult, g)
        assert panel[gi]["distinct"] == one["distinct"], f"guide {gi} distinct"
        assert panel[gi]["sites"] == one["sites"], f"guide {gi} sites"


def test_panel_counts_are_chunk_invariant():
    from scotoma.universe import panel_reference_counts
    rng = np.random.default_rng(22)
    seq = "".join(rng.choice(list("ACGT"), size=40_000))
    codes, mult, _ = _universe_from(seq)
    guides = np.array([encode_spacer("".join(rng.choice(list("ACGT"), size=20)))
                       for _ in range(3)], dtype=np.uint64)
    assert (panel_reference_counts(codes, mult, guides, chunk=10 ** 9)
            == panel_reference_counts(codes, mult, guides, chunk=511))


def test_panel_progress_callback_is_called_once_per_chunk():
    from scotoma.universe import panel_reference_counts
    rng = np.random.default_rng(23)
    seq = "".join(rng.choice(list("ACGT"), size=20_000))
    codes, mult, _ = _universe_from(seq)
    guides = np.array([encode_spacer("A" * 20)], dtype=np.uint64)
    seen = []
    panel_reference_counts(codes, mult, guides, chunk=1_000,
                           progress=lambda i, n: seen.append((i, n)))
    assert seen, "progress was never called"
    assert seen[-1][0] == seen[-1][1], "final call should report the last chunk"
    assert [i for i, _ in seen] == list(range(1, len(seen) + 1))


def test_bytes_metadata_is_recovered_not_fatal(tmp_path):
    """Regression: a real GS_CHM13 file stores meta as bytes and killed a cluster job.

    str(b'{...}') is "b'{...}'", which fails JSON parsing at character 1.
    """
    import json as _json
    from scotoma.universe import load_universe
    out = str(tmp_path / "bytes_meta.npz")
    np.savez_compressed(
        out,
        codes=np.array([1, 2, 3], dtype=np.uint64),
        mult=np.array([1, 1, 1], dtype=np.int32),
        meta=np.array([_json.dumps({"label": "B", "pam": "NGG"}).encode("utf-8")], dtype=object),
    )
    codes, meta = load_universe(out)
    assert codes.size == 3
    assert meta["label"] == "B"
    assert meta["pam"] == "NGG"


def test_unreadable_metadata_degrades_instead_of_raising(tmp_path):
    """Metadata is descriptive; a corrupt block must not destroy a loaded universe."""
    from scotoma.universe import load_universe
    out = str(tmp_path / "bad_meta.npz")
    np.savez_compressed(
        out,
        codes=np.array([1, 2, 3], dtype=np.uint64),
        mult=np.array([1, 1, 1], dtype=np.int32),
        meta=np.array([b"\xff\xfe not json at all"], dtype=object),
    )
    codes, meta = load_universe(out)
    assert codes.size == 3
    assert "meta_unreadable" in meta
