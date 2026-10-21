"""Bulge-inclusive enumeration by pattern expansion, needing no bulge-aware aligner.

An SpCas9 target is 20nt spacer + NGG. For a single-base bulge:

  DNA bulge  the DNA target carries one base the guide skips, so the protospacer is 21nt
             and the site is 24nt. Expanded guide = the 20mer with an 'N' inserted at one
             position (21nt) + 'NNN'; pattern = 'N'*22 + 'GG'. The inserted N is a wildcard,
             which is exactly "this target base is bulged out and need not match", so it
             cannot contribute a mismatch.

  RNA bulge  the guide carries one base the DNA lacks, so the protospacer is 19nt and the
             site is 22nt. Expanded guide = the 20mer with one base deleted (19nt) + 'NNN';
             pattern = 'N'*20 + 'GG'.

Mismatches are therefore counted over the 20 real spacer bases (DNA case) or the 19 surviving
ones (RNA case), with the PAM excluded in both because the guide's last three characters are N.
"""
import sys

def dna_bulge_guides(spacer):
    """21nt guides: an N inserted at each INTERNAL position of the protospacer.

    Position 0 is excluded. Inserting N before the first spacer base leaves the 20 real bases
    matching protospacer positions 2..21 with the PAM at 22..24, which is exactly an ungapped
    site shifted one base, not a bulge. The toy validation found precisely that artefact: the
    exact on-target planted at offset 40 was re-found at offset 39 by the position-0 pattern,
    and because it lands at a different coordinate a locus-level dedup would have scored it as
    a new candidate and inflated the census.
    """
    out = []
    for i in range(1, len(spacer) + 1):
        g = spacer[:i] + "N" + spacer[i:]
        out.append((i, g))
    return out

def rna_bulge_guides(spacer):
    """19nt guides: one spacer base deleted, at any of the 20 positions. Deduplicated, since
    deleting either of two identical adjacent bases gives the same guide.

    All 20 deletion positions are used, including both terminals. That is deliberate and it is
    what the paper's undercount result rests on: the reference implementation of
    Cas-OFFinder-bulge restricts RNA deletions to internal positions, and restricting them is
    precisely the omission we measure, at 7.44% of single-RNA-bulge loci. The DNA case does
    exclude position 0, because inserting a wildcard before the first spacer base genuinely
    re-finds an ungapped site at a shifted coordinate; a terminal deletion does not have that
    property, since it shortens rather than shifts the match.

    An earlier version of this docstring claimed the terminals were excluded here too. That
    described the reference position set, not this code, and it was corrected on 2026-10-20.
    """
    seen, out = set(), []
    for i in range(0, len(spacer)):
        g = spacer[:i] + spacer[i+1:]
        if g not in seen:
            seen.add(g)
            out.append((i, g))
    return out

def write_input(path, refdir, kind, spacers, mm):
    if kind == "dna":
        pattern, gen, plen = "N" * 22 + "GG", dna_bulge_guides, 24
    elif kind == "rna":
        pattern, gen, plen = "N" * 20 + "GG", rna_bulge_guides, 22
    else:
        pattern, gen, plen = "N" * 21 + "GG", (lambda s: [(0, s)]), 23
    n = 0
    index = []
    with open(path, "w") as f:
        f.write(refdir + "\n")
        f.write(pattern + "\n")
        for gi, sp in enumerate(spacers):
            for pos, g in gen(sp):
                line = g + "NNN"
                assert len(line) == plen, (kind, len(line), plen)
                f.write(f"{line} {mm}\n")
                index.append((line, gi, sp, pos))
                n += 1
    return pattern, n, index

if __name__ == "__main__":
    import json
    refdir, kind, panel, out, mmS, idxout = sys.argv[1:7]
    mm = int(mmS)
    if panel.endswith(".json"):
        spacers = [r["spacer"] for r in json.load(open(panel))["rows"]]
    else:
        spacers = [l.strip() for l in open(panel) if l.strip()]
    pat, n, index = write_input(out, refdir, kind, spacers, mm)
    # the expanded-guide string is not unique across guides, so keep every owner
    owners = {}
    for line, gi, sp, pos in index:
        owners.setdefault(line, []).append([gi, sp, pos])
    json.dump(owners, open(idxout, "w"))
    print(f"{kind}: pattern {pat} ({len(pat)}), {n} expanded guides from {len(spacers)} spacers, "
          f"{len(owners)} distinct expanded strings, mm<={mm}")
