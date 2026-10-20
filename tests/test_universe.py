"""Universe tests. Membership must be strand-free; distance must not be."""
import numpy as np
import pytest

from scotoma.coding import canonical, encode_window, oriented_sites, revcomp64
from scotoma.universe import (build_universe, contains, load_universe, save_universe, site_count)


def test_build_counts_a_known_tiny_genome(write_fasta):
    # One forward site only.
    path = write_fasta("tiny.fa", [("chr1", "A" * 20 + "AGG")])
    codes, mult, meta = build_universe([path], label="tiny")
    assert codes.size == 1
    assert int(mult[0]) == 1
    assert meta["label"] == "tiny"
    assert meta["distinct_canonical"] == 1
    assert meta["target_sites_total"] == 1
    assert meta["records"] == 1


def test_build_rejects_a_genome_with_no_usable_sequence(write_fasta):
    path = write_fasta("empty.fa", [("chr1", "ACGT")])
    with pytest.raises(ValueError):
        build_universe([path], label="empty")


def test_repeated_window_gets_multiplicity(write_fasta):
    window = "A" * 20 + "AGG"
    spacer_gap = "TTTTTTTTTTTT"
    path = write_fasta("rep.fa", [("chr1", window + spacer_gap + window)])
    codes, mult, meta = build_universe([path], label="rep")
    assert codes.size == 1
    assert int(mult[0]) == 2
    assert meta["max_sites"] == 2
    assert meta["single_site"] == 0


def test_membership_is_strand_free(write_fasta, random_sequence):
    """A site the reference holds on one strand must be found when queried from the other."""
    seq = random_sequence(3000, seed=1)
    path = write_fasta("ref.fa", [("chr1", seq)])
    codes, _ = load_universe(_save(write_fasta, path, "u.npz"))
    sites = oriented_sites(np.frombuffer(seq.encode(), dtype=np.uint8))
    assert sites.size > 0
    assert contains(codes, sites).all()
    # Query the same sites as they would appear from the opposite strand.
    assert contains(codes, revcomp64(sites)).all()


def test_absent_window_is_reported_absent(write_fasta, random_sequence):
    seq = random_sequence(3000, seed=2)
    path = write_fasta("ref.fa", [("chr1", seq)])
    codes, _ = load_universe(_save(write_fasta, path, "u.npz"))
    # A homopolymer window is not going to occur in 3 kb of random sequence.
    alien = np.array([encode_window("A" * 20 + "AGG")], dtype=np.uint64)
    assert not contains(codes, alien).any()


def test_save_and_load_roundtrip(tmp_path, write_fasta, random_sequence):
    seq = random_sequence(2000, seed=3)
    path = write_fasta("ref.fa", [("chr1", seq)])
    codes, mult, meta = build_universe([path], label="rt")
    out = str(tmp_path / "rt.npz")
    save_universe(out, codes, mult, meta)
    loaded, loaded_mult, loaded_meta = load_universe(out, with_mult=True)
    assert np.array_equal(loaded, codes)
    assert np.array_equal(loaded_mult, mult)
    assert loaded_meta == meta


def test_stored_codes_are_sorted_and_unique(write_fasta, random_sequence):
    seq = random_sequence(5000, seed=4)
    path = write_fasta("ref.fa", [("chr1", seq)])
    codes, _, _ = build_universe([path], label="s")
    assert (codes[1:] > codes[:-1]).all()


def test_load_rejects_an_unsorted_file(tmp_path):
    """Every membership test binary-searches, so an unsorted file must fail loudly."""
    import json
    out = str(tmp_path / "bad.npz")
    np.savez_compressed(out,
                        codes=np.array([5, 3, 9], dtype=np.uint64),
                        mult=np.array([1, 1, 1], dtype=np.int32),
                        meta=np.array([json.dumps({"label": "bad"})], dtype=object))
    with pytest.raises(ValueError, match="not sorted"):
        load_universe(out)


def test_site_count_returns_multiplicity_and_zero_for_absent(write_fasta):
    window = "A" * 20 + "AGG"
    path = write_fasta("rep.fa", [("chr1", window + "TTTTTTTTTTTT" + window)])
    codes, mult, meta = build_universe([path], label="rep")
    present = np.array([encode_window(window)], dtype=np.uint64)
    absent = np.array([encode_window("C" * 20 + "AGG")], dtype=np.uint64)
    assert int(site_count(codes, mult, present)[0]) == 2
    assert int(site_count(codes, mult, absent)[0]) == 0


def test_contains_on_an_empty_universe_is_all_false():
    empty = np.array([], dtype=np.uint64)
    q = np.array([encode_window("A" * 20 + "AGG")], dtype=np.uint64)
    assert contains(empty, q).tolist() == [False]


def test_two_fasta_files_are_pooled(write_fasta):
    a = write_fasta("a.fa", [("chr1", "A" * 20 + "AGG")])
    b = write_fasta("b.fa", [("chr2", "C" * 20 + "AGG")])
    codes, _, meta = build_universe([a, b], label="pooled")
    assert meta["records"] == 2
    assert codes.size == 2


def test_gzipped_input_reads_the_same(write_fasta, random_sequence):
    seq = random_sequence(1500, seed=6)
    plain = write_fasta("p.fa", [("chr1", seq)])
    gz = write_fasta("p.fa.gz", [("chr1", seq)])
    codes_plain, _, _ = build_universe([plain], label="p")
    codes_gz, _, _ = build_universe([gz], label="p")
    assert np.array_equal(codes_plain, codes_gz)


def test_primary_only_skips_decorated_record_names(write_fasta):
    path = write_fasta("alt.fa", [("chr1", "A" * 20 + "AGG"),
                                  ("chr1_KI270706v1_random", "C" * 20 + "AGG")])
    kept, _, meta_kept = build_universe([path], label="k", primary_only=True)
    allrec, _, meta_all = build_universe([path], label="a", primary_only=False)
    assert meta_kept["skipped_records"] == 1
    assert kept.size == 1
    assert allrec.size == 2


def _save(write_fasta, fasta_path, name):
    """Build a universe from a FASTA and persist it next to that file."""
    import os
    codes, mult, meta = build_universe([fasta_path], label="u")
    out = os.path.join(os.path.dirname(fasta_path), name)
    save_universe(out, codes, mult, meta)
    return out
