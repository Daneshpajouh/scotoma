#!/usr/bin/env python3
"""The minimum number of concrete patterns a complete single-deletion search must use.

    python3 -I min_patterns.py GUIDES.txt [--budget M]          certificate per guide
    python3 -I min_patterns.py GUIDES.txt --audit EMITTED.tsv   verdict on a tool's output

A one-base RNA bulge makes a guide searchable only through its distinct single-deletion words.
For a guide g of length n the run count

    r(g) = 1 + |{ i : g_i != g_{i+1} }|

is computable in one pass and equals the number of those distinct words. The companion's
Theorem 5 (universal low-budget floor) says that for every n >= 20 and every budget m <= 9
no family of fewer than r equal-radius Hamming balls centred in the four-letter alphabet
covers the union those words must search. So on that range r is not a heuristic target: it is
the exact minimum, emitting fewer is an enumeration defect, and emitting more costs output
without covering anything new.

Outside that range this program reports the run count and says the minimum is NOT PROVED. It
does not guess. Centres are four-letter words only; wildcards are not admitted (companion
notation: X is an unused concrete letter, never a wildcard).
"""
import csv, sys

PROVED_MIN_LENGTH, PROVED_MAX_BUDGET = 20, 9
LETTERS = set("ACGT")


def runs(g):
    return 1 + sum(1 for i in range(len(g) - 1) if g[i] != g[i + 1])


def deletion_words(g):
    """The distinct words a single deletion can produce. Its size is r(g)."""
    return {g[:i] + g[i + 1:] for i in range(len(g))}


def proved(n, m):
    return n >= PROVED_MIN_LENGTH and m <= PROVED_MAX_BUDGET


def certificate(g, m):
    g = g.upper()
    if not g or any(c not in LETTERS for c in g):
        raise ValueError("guide must be a nonempty word over ACGT")
    n, r = len(g), runs(g)
    ok = proved(n, m)
    return dict(guide=g, n=n, r=r,
                minimum=r if ok else None,
                status="PROVED" if ok else "NOT PROVED",
                # a position-enumerating emitter spends one pattern per deletion position;
                # the distinct-word count is what a complete search actually needs
                positions=n,
                redundancy=n - r)


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    path, m, audit = argv[0], 4, None
    i = 1
    while i < len(argv):
        if argv[i] == "--budget":
            m = int(argv[i + 1]); i += 2
        elif argv[i] == "--audit":
            audit = argv[i + 1]; i += 2
        else:
            print("unknown argument %s" % argv[i]); return 2
    if not isinstance(m, int) or m < 0:
        print("REJECT: budget must be a nonnegative integer, got %r" % (m,))
        return 2
    guides = [l.strip().upper() for l in open(path) if l.strip()
              and not l.startswith("#")]
    rows = [certificate(g, m) for g in guides]

    if audit is None:
        print("# minimum concrete patterns for a complete single-deletion search, budget %d" % m)
        print("\t".join(["guide", "n", "r", "minimum", "status", "positions", "redundancy"]))
        for c in rows:
            print("\t".join(str(c[k]) for k in
                            ("guide", "n", "r", "minimum", "status", "positions", "redundancy")))
        tot = sum(c["r"] for c in rows)
        print("# %d guides, %d patterns required in total, %d of %d proved"
              % (len(rows), tot, sum(1 for c in rows if c["minimum"] is not None), len(rows)))
        return 0

    emitted = {}
    with open(audit, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            emitted.setdefault(r["guide"].upper(), []).append(r["pattern"].upper())
    print("\t".join(["guide", "r", "emitted_rows", "distinct_emitted", "verdict"]))
    bad = 0
    for c in rows:
        g = c["guide"]
        got = emitted.get(g, [])
        words = deletion_words(g)
        distinct = {p for p in got if p in words}
        foreign = [p for p in got if p not in words]
        if foreign:
            v = "INVALID: %d patterns are not deletion words of this guide" % len(foreign)
        elif len(distinct) < c["r"]:
            # Word-set shortfall. "Loses sites" is reserved for budgets where the
            # companion proves K=r, so omitting a classical word is a cover defect.
            if c["minimum"] is not None:
                v = "INCOMPLETE: %d of %d words, loses sites" % (len(distinct), c["r"])
            else:
                v = "INCOMPLETE_LIBRARY: %d of %d classical words" % (len(distinct), c["r"])
        elif len(got) > c["r"]:
            v = "REDUNDANT: %d rows for %d words" % (len(got), c["r"])
        else:
            v = "EXACT"
        bad += v != "EXACT"
        print("\t".join([g, str(c["r"]), str(len(got)), str(len(distinct)), v]))
    print("# %d of %d guides are not served exactly" % (bad, len(rows)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
