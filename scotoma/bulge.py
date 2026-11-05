"""Single-base bulge enumeration by wildcard expansion.

The construction here is **not ours**. It is Cas-OFFinder-bulge (Jeongbin Park and
Sangsu Bae, v1.2 beta, 2016-07-06, https://github.com/hyugel/cas-offinder-bulge),
which enumerates bulge-bearing off-target candidates by inserting wildcards into the
spacer for DNA bulges and deleting spacer bases for RNA bulges, then handing the
expanded patterns to an exact matcher. A wildcard is an unconstrained position, so it
models a bulged base, and it cannot contribute a mismatch, so the budget keeps its
usual meaning over the real spacer bases. No gapped aligner is required. The original
handles bulges of arbitrary size; this module covers the single-base case.

What this module adds is one observation about the position set.

Definitions. A spacer is g in {A,C,G,T}^20 and an SpCas9 site is a 23-mer ending NGG.
Write P_i = g_1..g_i N g_{i+1}..g_20 for the DNA expansion at position i in {0..20}.

The reference implementation iterates i in {1..19}, excluding both ends.

  i = 0 must be excluded. A P_0 match on w_1..w_21 with w_23 = w_24 = G leaves w_1
  free and puts w_2..w_21 within the budget of g, so w_2..w_24 is an ordinary ungapped
  candidate: P_0 carries no bulge information. Because its reported start differs by
  one base, a locus-level deduplication does not remove it, and it would be counted as
  a new candidate.

  i = 20 need not be. P_20 requires w_23 = w_24 = G, whereas the ungapped reading of
  w_1..w_23 would require w_22 = w_23 = G, which P_20 does not imply. It describes a
  genuine unpaired DNA base between protospacer and PAM. The RNA terminals are likewise
  not degenerate.

So bulge counts obtained with the standard position set are LOWER BOUNDS. These
functions default to the full set with only i = 0 removed, and standard_set=True
reproduces the reference tool's set for comparability.
"""

SPACER_LEN = 20


def dna_bulge_patterns(spacer, standard_set=False):
    """Expanded 21-mers for one DNA bulge, as (bulge_position, pattern) pairs.

    i=0 is always excluded: it is degenerate with the ungapped candidate one base
    downstream. i=20 is retained by default because it is not degenerate; pass
    standard_set=True to reproduce Cas-OFFinder-bulge's set, which omits it.
    """
    if len(spacer) != SPACER_LEN:
        raise ValueError(f"spacer must be {SPACER_LEN} bases, got {len(spacer)}")
    last = len(spacer) if standard_set else len(spacer) + 1
    return [(i, spacer[:i] + "N" + spacer[i:]) for i in range(1, last)]


def rna_bulge_patterns(spacer, standard_set=False):
    """Expanded 19-mers for one RNA bulge, as (deleted_index, pattern) pairs.

    Deduplicated: deleting either of two identical adjacent bases yields the same
    pattern, and a duplicate would double-count the locus. standard_set=True
    reproduces Cas-OFFinder-bulge's set, which omits both terminals.
    """
    if len(spacer) != SPACER_LEN:
        raise ValueError(f"spacer must be {SPACER_LEN} bases, got {len(spacer)}")
    seen, out = set(), []
    rng = range(1, len(spacer) - 1) if standard_set else range(len(spacer))
    for i in rng:
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
