"""Codec tests, written from what the encoding must mean rather than from how it is implemented."""
import numpy as np
import pytest

from scotoma.coding import (CODE_TABLE, K, PAM_LEN, SPACER_LEN, canonical, decode_window,
                          encode_spacer, encode_window, oriented_sites, revcomp64, spacer_of)


def _bases(text):
    return np.frombuffer(text.encode(), dtype=np.uint8)


def test_constants_are_consistent():
    assert K == 23
    assert PAM_LEN == 3
    assert SPACER_LEN == 20
    assert SPACER_LEN + PAM_LEN == K


def test_roundtrip_encode_decode():
    for window in ("A" * 20 + "AGG", "ACGT" * 5 + "TGG", "G" * 20 + "CGG"):
        code = encode_window(window)
        assert code is not None
        assert decode_window(code) == window


def test_encode_rejects_bad_input():
    assert encode_window("ACGT") is None                      # too short
    assert encode_window("A" * 24) is None                     # too long
    assert encode_window("N" + "A" * 19 + "AGG"[1:]) is None   # contains N
    assert encode_spacer("A" * 19) is None
    assert encode_spacer("N" * 20) is None


def test_encode_spacer_puts_the_spacer_where_distance_reads_it():
    spacer = "ACGTACGTACGTACGTACGT"
    code = encode_spacer(spacer)
    assert code is not None
    assert spacer_of(code) == spacer
    # The PAM bits are left zero, which decodes to AAA.
    assert decode_window(code) == spacer + "AAA"


def test_revcomp_is_an_involution():
    codes = np.array([encode_window("ACGT" * 5 + "TGG"),
                      encode_window("A" * 20 + "AGG"),
                      encode_window("TTGCAACGTTGCAACGTTGC" + "AGG")], dtype=np.uint64)
    assert np.array_equal(revcomp64(revcomp64(codes)), codes)


def test_revcomp_matches_a_hand_computed_case():
    # revcomp("CCT" + "T"*20) reverses to "T"*20 + "TCC" and complements to "A"*20 + "AGG".
    code = encode_window("CCT" + "T" * 20)
    assert decode_window(int(revcomp64(code))) == "A" * 20 + "AGG"


def test_canonical_is_strand_free():
    code = encode_window("ACGT" * 5 + "TGG")
    rc = int(revcomp64(code))
    assert int(canonical(code)) == int(canonical(rc)) == min(code, rc)


def test_forward_site_is_a_window_ending_in_gg():
    sites = oriented_sites(_bases("A" * 20 + "AGG"))
    assert sites.size == 1
    assert decode_window(int(sites[0])) == "A" * 20 + "AGG"


def test_reverse_site_is_a_window_starting_cc_and_comes_back_oriented():
    """A reverse-strand site appears as CC... in the forward sequence.

    It must be returned reverse-complemented, so that what Scotoma hands downstream always reads
    [spacer]NGG. This is the distinction the whole package turns on.
    """
    sites = oriented_sites(_bases("CCT" + "T" * 20))
    assert sites.size == 1
    oriented = decode_window(int(sites[0]))
    assert oriented == "A" * 20 + "AGG"
    assert oriented[-2:] == "GG"


def test_forward_site_is_returned_oriented_even_when_canonical_differs():
    """A case where the oriented and canonical forms are *not* equal.

    The obvious reverse-strand example above cannot catch a canonicalisation bug, because for
    that window the oriented form happens to be the smaller of the two and so equals the
    canonical form. Here ``T*20 + TGG`` has reverse complement ``CCA + A*20``, which is
    numerically smaller, so canonical != oriented and returning the wrong one shows.
    """
    window = "T" * 20 + "TGG"
    sites = oriented_sites(_bases(window))
    assert sites.size == 1
    oriented = int(sites[0])
    assert oriented == encode_window(window)
    assert oriented != int(canonical(oriented))
    assert decode_window(oriented)[-2:] == "GG"


def test_reverse_site_is_returned_oriented_even_when_canonical_differs():
    """The reverse-strand counterpart, chosen so the canonical form is the genomic window.

    ``CC + G*20 + A`` is a reverse-strand site whose oriented form is ``T + C*20 + GG``. The
    genomic window sorts lower, so canonicalising here would hand back a code that does not read
    [spacer]NGG at all.
    """
    window = "CC" + "G" * 20 + "A"
    assert len(window) == K
    sites = oriented_sites(_bases(window))
    assert sites.size == 1
    oriented = decode_window(int(sites[0]))
    assert oriented == "T" + "C" * 20 + "GG"
    assert int(sites[0]) != int(canonical(int(sites[0])))


def test_a_window_on_both_strands_counts_twice():
    """Ending GG and beginning CC are two different binding sites, not one."""
    window = "CC" + "A" * 19 + "GG"
    assert len(window) == K
    sites = oriented_sites(_bases(window))
    assert sites.size == 2
    forward = encode_window(window)
    assert sorted(int(x) for x in sites) == sorted([forward, int(revcomp64(forward))])


def test_windows_containing_n_are_dropped():
    # The only GG-ending window here also contains an N, so nothing is extracted.
    assert oriented_sites(_bases("A" * 10 + "N" + "A" * 9 + "AGG")).size == 0


def test_short_sequence_yields_nothing():
    assert oriented_sites(_bases("ACGT")).size == 0
    assert oriented_sites(_bases("")).size == 0


def test_every_extracted_site_has_a_valid_pam():
    rng = np.random.default_rng(11)
    seq = "".join(rng.choice(list("ACGT"), size=4000))
    for code in oriented_sites(_bases(seq)).tolist():
        assert decode_window(int(code))[-2:] == "GG"


def test_code_table_maps_lowercase_like_uppercase():
    assert CODE_TABLE[ord("a")] == CODE_TABLE[ord("A")] == 0
    assert CODE_TABLE[ord("t")] == CODE_TABLE[ord("T")] == 3
    assert CODE_TABLE[ord("N")] == 255


def test_lowercase_sequence_is_extracted():
    """Soft-masked assembly sequence is real sequence and must not be silently skipped."""
    upper = oriented_sites(_bases("A" * 20 + "AGG"))
    lower = oriented_sites(_bases(("A" * 20 + "AGG").lower()))
    assert np.array_equal(upper, lower)
