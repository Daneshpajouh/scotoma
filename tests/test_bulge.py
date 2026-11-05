import pytest
from scotoma.bulge import dna_bulge_patterns, rna_bulge_patterns, site_pattern

G = "GAACAAGGTGTTCCCACCCG"
COMP = str.maketrans("ACGT", "TGCA")


def test_pattern_lengths_match_site_lengths():
    assert len(site_pattern("none")) == 23
    assert len(site_pattern("dna")) == 24
    assert len(site_pattern("rna")) == 22


def test_ngg_scope_not_nag():
    # the common example pattern ends RG and would admit NAG; ours must end GG
    for kind in ("none", "dna", "rna"):
        assert site_pattern(kind).endswith("GG")
        assert "R" not in site_pattern(kind)


def test_dna_expansion_excludes_only_the_degenerate_position():
    pats = dna_bulge_patterns(G)
    positions = [i for i, _ in pats]
    assert positions == list(range(1, 21)), "i=0 excluded by the proposition, i=20 retained"
    assert all(len(p) == 21 for _, p in pats)
    assert all(p.count("N") == 1 for _, p in pats)


def test_dna_pattern_realises_its_bulge_position():
    # deleting w_{i+1} from a 21-mer must recover the spacer exactly
    for i, p in dna_bulge_patterns(G):
        assert p[:i] == G[:i]
        assert p[i] == "N"
        assert p[i + 1:] == G[i:]


def test_pam_proximal_dna_position_is_present():
    # P_20 is a genuine unpaired base between protospacer and PAM
    assert (20, G + "N") in dna_bulge_patterns(G)


def test_rna_expansion_is_deduplicated_and_shorter_by_one():
    pats = rna_bulge_patterns(G)
    assert all(len(p) == 19 for _, p in pats)
    assert len({p for _, p in pats}) == len(pats), "duplicates would double-count a locus"
    # G has adjacent repeats, so dedup must actually bite
    assert len(pats) < len(G)


def test_rna_patterns_are_single_deletions_of_the_spacer():
    for i, p in rna_bulge_patterns(G):
        assert p == G[:i] + G[i + 1:]


def test_degenerate_position_would_have_been_an_ungapped_match():
    """P_0 = N + g aligns g to w_2..w_21, so w_2..w_24 is an ungapped candidate.

    Construct the witness directly: an exact on-target preceded by any base is matched
    by P_0 at an offset one base earlier, which is why i=0 is excluded.
    """
    p0 = "N" + G
    window = "T" + G + "AGG"          # w_1 = T, then g, then the PAM
    assert len(window) == 24
    body = window[:21]
    assert all(a == "N" or a == b for a, b in zip(p0, body))
    # and the 23-mer one base downstream is an ordinary ungapped site of g
    assert window[1:24] == G + "AGG"


@pytest.mark.parametrize("bad", ["", "ACGT", "A" * 19, "A" * 21])
def test_spacer_length_is_enforced(bad):
    with pytest.raises(ValueError):
        dna_bulge_patterns(bad)
    with pytest.raises(ValueError):
        rna_bulge_patterns(bad)


def test_standard_set_reproduces_the_reference_tool():
    """Cas-OFFinder-bulge iterates i in {1..19} for DNA and {1..18} for one RNA bulge."""
    assert [i for i, _ in dna_bulge_patterns(G, standard_set=True)] == list(range(1, 20))
    rna_std = [i for i, _ in rna_bulge_patterns(G, standard_set=True)]
    assert min(rna_std) >= 1 and max(rna_std) <= len(G) - 2


def test_full_set_is_a_superset_of_the_standard_one():
    """The standard set omits a non-degenerate position, so its counts are lower bounds."""
    full = {p for _, p in dna_bulge_patterns(G)}
    std = {p for _, p in dna_bulge_patterns(G, standard_set=True)}
    assert std < full
    assert (G + "N") in full - std
