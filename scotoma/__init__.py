"""Scotoma: off-target candidate enumeration in target windows absent from a reference genome.

Every variant-aware off-target screen in current use is anchored to a linear reference. It can
substitute bases, and recent tools can represent indels and population variants, but a candidate
site in sequence that is *absent* from the reference is not expressible as a variant of it, so
such a site is invisible by construction rather than by oversight.

Scotoma enumerates there. Given a haplotype-resolved assembly, a reference, and a guide panel, it
partitions the assembly's SpCas9 target windows into those the reference contains and those it
does not, scores both against the same panel, and reports the candidate density of each. The
comparison is within one genome and the arms are matched on size, so guide composition and
population structure cancel.

Typical use::

    from scotoma import load_universe, load_panel, partition_haplotype, matched_arms, blindspot

    universe, _ = load_universe("GS_GRCh38.npz")
    guides, _ = load_panel("panel.txt")
    absent, pool, stats = partition_haplotype("HG00235.hap1.fa.gz", universe)
    blind, visible = matched_arms(absent, pool)
    result = blindspot(blind, visible, guides, partition_stats=stats)

The one rule to carry into any use of the lower-level pieces: mismatch distance is computed on
**oriented** codes and reference membership on **canonical** ones. :mod:`scotoma.coding` explains
why mixing them silently biases one arm of a comparison.
"""
from .blindspot import (DEFAULT_KEEP_PRESENT, DEFAULT_THRESHOLDS, blindspot, gc_matched_arms,
                        interpret_gates, matched_arms, partition_haplotype, replicate_ratios,
                        stat_matched_arms)
from .budget import budget_rows, census_sensitivity, m_star, miss_probability
from .coding import (DEFAULT_PAM, K, PAM_LEN, PAM_SCOPES, SPACER_LEN, canonical, decode_window,
                     encode_spacer, encode_window, iter_fasta, oriented_sites, revcomp64,
                     shuffle_panel, shuffle_spacer, spacer_gc, spacer_of,
                     window_complexity)
from .distance import (density_per_million, min_mismatches, spacer_mismatches, threshold_counts)
from .panel import load_panel, load_panel_changeseq, load_panel_text
from .universe import build_universe, contains, load_universe, save_universe, site_count

__version__ = "0.1.0"

__all__ = [
    "K", "PAM_LEN", "SPACER_LEN",
    "DEFAULT_PAM", "PAM_SCOPES",
    "canonical", "decode_window", "encode_spacer", "encode_window", "iter_fasta",
    "oriented_sites", "revcomp64", "shuffle_panel", "shuffle_spacer", "spacer_gc",
    "spacer_of", "window_complexity",
    "spacer_mismatches", "min_mismatches", "threshold_counts", "density_per_million",
    "build_universe", "save_universe", "load_universe", "contains", "site_count",
    "load_panel", "load_panel_text", "load_panel_changeseq",
    "partition_haplotype", "matched_arms", "gc_matched_arms", "blindspot", "interpret_gates",
    "replicate_ratios", "stat_matched_arms",
    "DEFAULT_THRESHOLDS", "DEFAULT_KEEP_PRESENT",
    "m_star", "miss_probability", "budget_rows", "census_sensitivity",
    "__version__",
]
