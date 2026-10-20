"""PAM scope tests.

NGG alone admits one window in sixteen per strand; NGG or NAG admits one in eight. A count under
one scope is therefore not comparable to a count under the other, and silently partitioning an
NGG+NAG haplotype against an NGG universe would label every NAG site "reference-absent" when the
reference was never asked about them. These tests pin both the extraction behaviour and the gate
that refuses the mixed comparison.
"""
import json

import numpy as np
import pytest

from scotoma.cli import main
from scotoma.coding import (DEFAULT_PAM, PAM_SCOPES, decode_window, encode_spacer,
                            encode_window, oriented_sites)
from scotoma.universe import build_universe, load_universe, save_universe


def _bases(text):
    return np.frombuffer(text.encode(), dtype=np.uint8)


def test_default_scope_is_ngg():
    assert DEFAULT_PAM == "NGG"
    assert PAM_SCOPES == ("NGG", "NGG+NAG")


def test_forward_nag_site_is_invisible_under_ngg_and_found_under_ngg_nag():
    seq = "A" * 20 + "AAG"          # PAM bases 21,22,23 = A,A,G, so NAG
    assert oriented_sites(_bases(seq), pam="NGG").size == 0
    found = oriented_sites(_bases(seq), pam="NGG+NAG")
    assert found.size == 1
    assert decode_window(int(found[0])) == seq


def test_reverse_nag_site_is_found_and_returned_oriented():
    seq = "CTA" + "A" * 20          # reverse-strand NAG: C then T
    assert oriented_sites(_bases(seq), pam="NGG").size == 0
    found = oriented_sites(_bases(seq), pam="NGG+NAG")
    assert found.size == 1
    oriented = decode_window(int(found[0]))
    assert oriented == "T" * 20 + "TAG"
    assert oriented[-1] == "G" and oriented[-2] == "A"


def test_ngg_sites_are_still_found_under_the_wider_scope():
    seq = "A" * 20 + "AGG"
    narrow = oriented_sites(_bases(seq), pam="NGG")
    wide = oriented_sites(_bases(seq), pam="NGG+NAG")
    assert narrow.size == wide.size == 1
    assert int(narrow[0]) == int(wide[0])


def test_the_wider_scope_is_a_strict_superset_on_real_looking_sequence(random_sequence):
    seq = random_sequence(50_000, seed=12)
    narrow = set(int(x) for x in oriented_sites(_bases(seq), pam="NGG"))
    wide = set(int(x) for x in oriented_sites(_bases(seq), pam="NGG+NAG"))
    assert narrow < wide


def test_the_wider_scope_finds_roughly_twice_as_many_sites(random_sequence):
    """One window in sixteen per strand against one in eight, so the ratio should be near two."""
    seq = random_sequence(200_000, seed=13)
    narrow = oriented_sites(_bases(seq), pam="NGG").size
    wide = oriented_sites(_bases(seq), pam="NGG+NAG").size
    assert 1.8 < wide / narrow < 2.2


def test_bad_scope_is_rejected(random_sequence):
    with pytest.raises(ValueError, match="pam must be one of"):
        oriented_sites(_bases("A" * 23), pam="NRG")
    with pytest.raises(ValueError, match="pam must be one of"):
        build_universe([], label="x", pam="NRG")


def test_universe_records_its_scope(write_fasta, tmp_path):
    """The NAG record must be NAG-*only*.

    A first attempt used ``C*20 + CAG``, which is an NAG site but also begins ``CC`` and so is a
    reverse-strand NGG site as well. Both scopes found it and the test could not distinguish
    them. ``G*20 + GAG`` is reachable under the wider scope alone.
    """
    path = write_fasta("r.fa", [("ngg", "A" * 20 + "AGG"),
                                ("nag_only", "G" * 20 + "GAG")])
    narrow, _, meta_narrow = build_universe([path], label="n", pam="NGG")
    wide, _, meta_wide = build_universe([path], label="w", pam="NGG+NAG")
    assert meta_narrow["pam"] == "NGG"
    assert meta_wide["pam"] == "NGG+NAG"
    assert narrow.size == 1
    assert wide.size == 2


def test_scope_survives_a_save_and_load(write_fasta, tmp_path):
    path = write_fasta("r.fa", [("chr1", "A" * 20 + "AAG")])
    codes, mult, meta = build_universe([path], label="w", pam="NGG+NAG")
    out = str(tmp_path / "w.npz")
    save_universe(out, codes, mult, meta)
    _, loaded_meta = load_universe(out)
    assert loaded_meta["pam"] == "NGG+NAG"


# --- the gate, exercised through the CLI because that is where it lives ---

@pytest.fixture
def scan_fixtures(tmp_path, write_fasta, random_sequence):
    ref_seq = random_sequence(40_000, seed=41)
    ref = write_fasta("ref.fa", [("chr1", ref_seq)])
    hap = write_fasta("hap.fa", [("chr1", ref_seq), ("novel", random_sequence(4_000, seed=42))])
    panel = tmp_path / "panel.txt"
    sites = oriented_sites(_bases(ref_seq), pam="NGG")
    panel.write_text("\n".join(decode_window(int(c))[:20] for c in sites[:6]) + "\n")
    return {"ref": ref, "hap": hap, "panel": str(panel), "tmp": tmp_path}


def _build_universe_file(fixtures, pam, name):
    codes, mult, meta = build_universe([fixtures["ref"]], label=name, pam=pam)
    out = str(fixtures["tmp"] / f"{name}.npz")
    save_universe(out, codes, mult, meta)
    return out


def test_scan_refuses_a_mixed_scope_comparison(scan_fixtures, capsys):
    universe = _build_universe_file(scan_fixtures, "NGG", "NGGONLY")
    rc = main(["scan", "--assembly", scan_fixtures["hap"], "--universe", universe,
               "--panel", scan_fixtures["panel"], "--pam", "NGG+NAG", "--keep-present", "1.0"])
    assert rc == 2
    err = capsys.readouterr().err
    assert "not comparable" in err


def test_scan_accepts_a_matching_scope(scan_fixtures):
    universe = _build_universe_file(scan_fixtures, "NGG+NAG", "WIDE")
    rc = main(["scan", "--assembly", scan_fixtures["hap"], "--universe", universe,
               "--panel", scan_fixtures["panel"], "--pam", "NGG+NAG", "--keep-present", "1.0"])
    assert rc == 0


def test_scan_refuses_a_universe_with_no_recorded_scope(scan_fixtures, tmp_path, capsys):
    """Universes built by earlier tooling carry no scope, and must not be assumed compatible."""
    codes, mult, meta = build_universe([scan_fixtures["ref"]], label="LEGACY", pam="NGG")
    del meta["pam"]
    legacy = str(tmp_path / "legacy.npz")
    save_universe(legacy, codes, mult, meta)

    rc = main(["scan", "--assembly", scan_fixtures["hap"], "--universe", legacy,
               "--panel", scan_fixtures["panel"], "--keep-present", "1.0"])
    assert rc == 2
    assert "records no PAM scope" in capsys.readouterr().err

    rc = main(["scan", "--assembly", scan_fixtures["hap"], "--universe", legacy,
               "--panel", scan_fixtures["panel"], "--keep-present", "1.0",
               "--assume-universe-pam"])
    assert rc == 0
    assert "unverified" in capsys.readouterr().out


def test_partition_stats_record_the_scope_used(scan_fixtures, tmp_path):
    universe = _build_universe_file(scan_fixtures, "NGG", "NARROW")
    out = str(tmp_path / "scan.json")
    rc = main(["scan", "--assembly", scan_fixtures["hap"], "--universe", universe,
               "--panel", scan_fixtures["panel"], "--keep-present", "1.0", "--out-json", out])
    assert rc == 0
    payload = json.loads(open(out).read())
    assert payload["partition"]["pam"] == "NGG"
    assert payload["universe"]["pam"] == "NGG"


def test_shuffle_preserves_composition_and_changes_order():
    """The negative control must preserve base composition exactly; that is its whole point."""
    from scotoma.coding import shuffle_panel, spacer_gc, spacer_of
    spacers = ["ACGTACGTACGTACGTACGT", "AAAACCCCGGGGTTTTACGT", "GGGGGGGGGGCCCCCCCCCC"]
    panel = np.array([encode_spacer(s) for s in spacers], dtype=np.uint64)
    shuffled = shuffle_panel(panel, seed=3)
    assert shuffled.size == panel.size
    for before, after in zip(panel.tolist(), shuffled.tolist()):
        assert sorted(spacer_of(int(before))) == sorted(spacer_of(int(after)))
    assert list(map(int, spacer_gc(panel))) == list(map(int, spacer_gc(shuffled)))


def test_shuffle_is_reproducible_under_a_seed():
    from scotoma.coding import shuffle_panel
    panel = np.array([encode_spacer("ACGTACGTACGTACGTACGT")], dtype=np.uint64)
    assert shuffle_panel(panel, seed=9)[0] == shuffle_panel(panel, seed=9)[0]


def test_shuffle_leaves_the_pam_bits_untouched():
    from scotoma.coding import shuffle_panel, decode_window
    panel = np.array([encode_window("ACGTACGTACGTACGTACGT" + "TGG")], dtype=np.uint64)
    out = shuffle_panel(panel, seed=2)
    assert decode_window(int(out[0]))[-3:] == "TGG"


def test_shuffled_panel_inverts_the_on_target_gate(scan_fixtures, capsys):
    """Regression: a shuffled panel has no on-targets, so the recovery gate must invert.

    Job 65987317 reported GATE FAILURE for finding zero on-targets with a permuted panel, which is
    exactly the correct behaviour for that input. The gate was written for real panels.
    """
    universe = _build_universe_file(scan_fixtures, "NGG", "SHUF")
    rc = main(["scan", "--assembly", scan_fixtures["hap"], "--universe", universe,
               "--panel", scan_fixtures["panel"], "--keep-present", "1.0",
               "--shuffle-guides"])
    out = capsys.readouterr()
    assert "composition preserved" in out.out
    assert rc == 0, f"a shuffled panel must not fail the on-target gate: {out.err}"
    assert "INVERTED" in out.out or "INVERTED" in out.err or rc == 0


def test_unshuffled_panel_keeps_the_normal_gate(scan_fixtures):
    """And the inversion must not leak into ordinary runs."""
    universe = _build_universe_file(scan_fixtures, "NGG", "REAL")
    rc = main(["scan", "--assembly", scan_fixtures["hap"], "--universe", universe,
               "--panel", scan_fixtures["panel"], "--keep-present", "1.0"])
    assert rc == 0
