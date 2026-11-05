# The exactness dividend on an independent published benchmark

One command, one input file, no genome and no executable:

    python3 -I verify_exactness_dividend.py

Fourteen checks, two of them controls that must come out at the extremes (an alternating
spacer has no redundancy; a homopolymer spacer is nineteen of twenty redundant).

## What it establishes

A per-position emitter enumerates one deletion placement per spacer position, n per guide.
The number of distinct deletion words is the run count r, and r is also the exact minimum
number of equal-radius balls that cover the deletion union, so no family of any kind is
smaller. The difference n - r is work no correct emitter needs.

On the 61 spacers of the published tool comparison in `published_benchmark_guides.txt`:

    placements enumerated                 1,220
    distinct words required (exact)          891
    provably redundant                       329   = 27.0%
    worst guide  TTTTTTTGTTCCTCTTTTGG    12 of 20 redundant

The identity r = number of distinct deletion words is re-derived per guide by set
construction, not assumed.

Separately, 47 of the 61 guides (77.0%) have unequal last two spacer bases and are therefore
at risk under the terminal law; 14 are rescued by it.

## Provenance, which is second-hand

The preprint that reports this benchmark names a repository, `chopoff-benchmark`, that is not
public under its authors' account: a direct lookup, a code search for the filename it cites,
a repository-name search and both of the relevant laboratory organisations return nothing.
The guide list here is the copy carried by a public derivative repository, which still holds
the original project file. The only check available on the list itself is that the preprint
states 61 sgRNAs and the file holds exactly 61, all 20 nt and over the four letters, which
the verifier asserts. The companion file in that repository prints the counts 3,058 and
10,918, whose difference is the 7,860 the preprint's text states.

This limitation is stated in the paper as well as here. The arithmetic above depends on the
guide list alone and can be repeated by anyone from the same public file.

## What this directory does NOT claim

It does not claim that the terminal omission accounts for the published 7,860 missing
off-targets. The full list is not deposited, only examples, and the first printed example is
an alignment carrying a gap in both strings at once: measured against the reference the same
file prints, no one-sided bulge reaches the stated distance there, so that site is a
different limitation.

## The exact locator, so the recovered list can be re-fetched

A reviewer noted that the provenance was disclosed but the release named no locator, so the
calculation authenticates combinatorial properties rather than membership in the cited panel.
Replacing sixty of the sixty-one spacers while preserving each guide's equality pattern still
passed every assertion, which is correct behaviour for a combinatorial check and is exactly why
the locator belongs here:

    repository   github.com/rickbeeloo/sassy-crispr-bench
    commit       f498a4d96a313dbbd832f1c46fa94939743adff6
    path         data/curated_guides_wo_PAM.txt
    git blob     d6d13c0b6e4be6f13353ce385a950ac7e990cfce   (1,280 bytes)
    sha256       1e6c9b7d1a52bb20e3bdfacf886c4546e123df3a53703731b85c7aac7a76ec94

That repository is a public derivative; the repository the preprint names is not public under
its authors' account. The arithmetic here is therefore explicitly conditional on this recovered
list, not on a uniquely authenticated original, and the verifier asserts only the two facts the
list can support against the preprint's text: exactly 61 spacers, each 20 nt.
