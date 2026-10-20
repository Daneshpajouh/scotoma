"""Single-base bulge enumeration by wildcard expansion.

A bulge-inclusive off-target census is usually taken to need a gapped aligner, or
an index over 22-mers and 24-mers alongside the 23-mer one. Neither is necessary.
A single-base bulge can be reduced to exact matching against a pattern carrying one
wildcard, because a wildcard is precisely an unconstrained position and so models a
bulged base, and because a wildcard cannot contribute a mismatch the budget keeps its
usual meaning over the real spacer bases.

Definitions. A spacer is g in {A,C,G,T}^20 and an SpCas9 site is a 23-mer ending NGG.

  ungapped       w_1..w_23 with w_22 = w_23 = G and d(w_1..w_20, g) <= k.

  DNA bulge      w_1..w_24 with w_23 = w_24 = G and some j in {1..21} such that
                 deleting w_j from w_1..w_21 leaves a 20-mer within k of g. The
                 bulged base is w_j, an extra base in the DNA.

  RNA bulge      w_1..w_22 with w_21 = w_22 = G and some i in {1..20} such that
                 deleting g_i from g leaves a 19-mer within k of w_1..w_19.

The expansion. For i in {0..20} let P_i = g_1..g_i N g_{i+1}..g_20, of length 21.
Matching P_i against w_1..w_21 leaves w_{i+1} free and aligns the rest of g in order,
which is exactly "delete w_{i+1}, then compare to g". So P_i realises bulge position
j = i+1, and P_0..P_20 realise all 21 positions.

PROPOSITION (the PAM-distal position is degenerate). Every match of P_0 coincides with
an ungapped candidate one base downstream.

  Proof. A P_0 match on w_1..w_21 with w_23 = w_24 = G leaves w_1 free and gives
  d(w_2..w_21, g) <= k. The 23-mer w_2..w_24 then has first twenty bases w_2..w_21,
  within k of g, and last two w_23 = w_24 = G. So w_2..w_24 is an ungapped candidate
  at budget k. QED

P_0 therefore carries no bulge information, and since its reported start differs by one
base a locus-level deduplication does not remove it: it would be scored as a new
candidate and inflate the census. i = 0 is excluded for that reason, which is necessary
rather than cosmetic.

By contrast P_20 = g_1..g_20 N requires w_23 = w_24 = G while the ungapped reading of
w_1..w_23 would require w_22 = w_23 = G, which P_20 does not imply. P_20 is NOT
degenerate: it is a genuine unpaired DNA base between protospacer and PAM, and it is
retained. The RNA terminals are likewise not provably degenerate and are retained.

Cost. 20 patterns per guide for one DNA bulge and at most 20 for one RNA bulge, against
a 40-fold scan cost for building and querying 22-mer and 24-mer universes.
"""

SPACER_LEN = 20


def dna_bulge_patterns(spacer):
    """Expanded 21-mers for one DNA bulge, as (bulge_position, pattern) pairs.

    Position i=0 is excluded by the proposition in this module's docstring: it is
    degenerate with the ungapped candidate one base downstream. i=20 is retained
    because it is not degenerate.
    """
    if len(spacer) != SPACER_LEN:
        raise ValueError(f"spacer must be {SPACER_LEN} bases, got {len(spacer)}")
    return [(i, spacer[:i] + "N" + spacer[i:]) for i in range(1, len(spacer) + 1)]


def rna_bulge_patterns(spacer):
    """Expanded 19-mers for one RNA bulge, as (deleted_index, pattern) pairs.

    Deduplicated: deleting either of two identical adjacent bases yields the same
    pattern, and a duplicate would double-count the locus.
    """
    if len(spacer) != SPACER_LEN:
        raise ValueError(f"spacer must be {SPACER_LEN} bases, got {len(spacer)}")
    seen, out = set(), []
    for i in range(len(spacer)):
        p = spacer[:i] + spacer[i + 1:]
        if p not in seen:
            seen.add(p)
            out.append((i, p))
    return out


def site_pattern(kind):
    """The fixed-position pattern a matcher scans for, NGG scope.

    'none' is the ordinary 23-mer, 'dna' the 24-mer admitting one DNA bulge, 'rna'
    the 22-mer admitting one RNA bulge. Note that the widely copied example pattern
    N*21+RG admits NAG as well as NGG and so silently changes the PAM scope.
    """
    if kind == "none":
        return "N" * 21 + "GG"
    if kind == "dna":
        return "N" * 22 + "GG"
    if kind == "rna":
        return "N" * 20 + "GG"
    raise ValueError(f"kind must be none, dna or rna, got {kind!r}")
