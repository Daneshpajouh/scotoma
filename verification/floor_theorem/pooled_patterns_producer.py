#!/usr/bin/env python3
"""Where the 1,629 pooled patterns come from: rebuild data/pooled_patterns_1629.txt from the guides.

WHY THIS EXISTS. pooled_volume_bound.py reads data/pooled_patterns_1629.txt, and until 2026-10-23
nothing in the release said where that list came from or rebuilt it. A reader could check the
bound against the file, but not that the file is the paper's panel. A released list should not be
its own authority, so this script regenerates it from the guides and prints its digest.

THE SOURCE. The 110 guides of the CHANGE-seq measurement table of Lazzarotto et al. (2020), file
CHANGEseq_CR_Lazzarotto_2020_dataset.csv, the file name the cluster scripts read with
--changeseq. The copy used here is datasets/CHANGEseq/include_on_targets/ in files/datasets.zip of
the OrensteinLab/CRISPR-Bulge repository, fetched through Git LFS (the Zenodo record
10.5281/zenodo.10902445 holds only the LFS pointer). That table is 878,760,836 bytes with SHA-256
68f3e1d7eb9e9e18e3ddc859752be8d156378a3e681b0fff329fe679366adc11. It is a third-party deposit and
is not redistributed here. Any copy with the same sgRNA column gives the same output; the script
prints the input's digest so a run log shows which copy was read.

WHAT IT DOES. Each distinct value of the sgRNA column is a 23-letter window ending in NGG. Its
first 20 letters are the guide. The patterns are the distinct single-base deletions of each guide
(length 19), pooled over all guides, one per line in sorted order. The script checks:

  1. 110 distinct guides over ACGT, each from a 23-letter window ending in NGG.
  2. 1,634 patterns counted guide by guide, 1,629 distinct: the libraries are NOT pairwise
     disjoint. Exactly 5 patterns are shared, each by two guides offset by one base, and they are
     listed with the guides that share them.
  3. The rebuilt list's SHA-256 equals the digest below, and its bytes equal the released file.
  4. The control: dropping one guide that has a pattern of its own must change the rebuilt
     digest, or check 3 proves nothing.

Usage:  pooled_patterns_producer.py DEPOSIT.csv [--patterns RELEASED.txt] [--write OUT.txt]
Default released file: data/pooled_patterns_1629.txt beside this script. --write saves the rebuilt
list for comparison; it refuses to write over the released file.
Exit 0 only when every check above holds; 1 otherwise; 2 if it could not run.
"""
import argparse
import csv
import hashlib
import os
import sys
from collections import defaultdict

N_GUIDES, LEN, N_WITH_MULTIPLICITY, N_PATTERNS, N_SHARED = 110, 20, 1634, 1629, 5
EXPECTED_SHA256 = "b749f0869f5d4a97d1e6c6c48f750694dbec54460846ad4730b5cba76136fa02"
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "data", "pooled_patterns_1629.txt")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def patterns(g):
    """The distinct single-deletion reads of g: one per maximal run."""
    return {g[:i] + g[i + 1:] for i in range(len(g))}


def render(guides):
    """The pattern file for these guides, exactly as released: sorted, one per line."""
    union = set().union(*(patterns(g) for g in guides))
    return "".join(p + "\n" for p in sorted(union)).encode("ascii")


def offset_by_one(a, b):
    return a[1:] == b[:-1] or b[1:] == a[:-1]


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("deposit")
    parser.add_argument("--patterns", default=DEFAULT)
    parser.add_argument("--write")
    args = parser.parse_args(argv[1:])      # a usage error exits 2, as documented

    try:
        windows = set()
        with open(args.deposit, newline="") as fh:
            for row in csv.DictReader(fh):
                windows.add(row["sgRNA"].strip().upper())
        released_bytes = open(args.patterns, "rb").read()
    except OSError as exc:
        print(f"NOT RUN: {exc}")
        return 2
    except KeyError:
        print(f"NOT RUN: {args.deposit} has no sgRNA column")
        return 2
    if not windows:
        print(f"NOT RUN: {args.deposit} has no rows")
        return 2

    print(f"deposit  {args.deposit}")
    print(f"         SHA-256 {sha256_file(args.deposit)}")
    guides = sorted({w[:LEN] for w in windows})
    owners = defaultdict(list)
    for g in guides:
        for p in patterns(g):
            owners[p].append(g)
    shared = {p: gs for p, gs in owners.items() if len(gs) > 1}
    total = sum(len(gs) for gs in owners.values())
    built = render(guides)
    built_sha = hashlib.sha256(built).hexdigest()

    # Control: the same build without one guide must not reproduce the digest. The dropped guide
    # must have a pattern no other guide has, or dropping it could not change anything.
    own = [g for g in guides if any(len(owners[p]) == 1 for p in patterns(g))]
    dropped = own[0] if own else guides[0]
    rest = [g for g in guides if g != dropped]
    control_sha = hashlib.sha256(render(rest)).hexdigest()
    lost = len(patterns(dropped) - set().union(*(patterns(g) for g in rest)))

    checks = [
        (f"{N_GUIDES} distinct guides over ACGT",
         len(guides) == len(windows) == N_GUIDES and all(set(g) <= set("ACGT") for g in guides),
         f"{len(windows)} windows, {len(guides)} guides, letters {''.join(sorted(set(''.join(guides))))}"),
        ("every window 23 letters ending NGG",
         all(len(w) == LEN + 3 and w[LEN:] == "NGG" for w in windows),
         f"lengths {sorted({len(w) for w in windows})}, tails {sorted({w[LEN:] for w in windows})}"),
        (f"{N_WITH_MULTIPLICITY:,} patterns guide by guide", total == N_WITH_MULTIPLICITY,
         f"{total:,}"),
        (f"{N_PATTERNS:,} distinct patterns", len(owners) == N_PATTERNS, f"{len(owners):,}"),
        (f"{N_SHARED} shared, each by two guides offset by one base",
         len(shared) == N_SHARED
         and all(len(gs) == 2 and offset_by_one(*gs) for gs in shared.values()),
         f"{len(shared)} shared"),
        ("digest of the rebuilt list", built_sha == EXPECTED_SHA256, built_sha),
        ("rebuilt bytes equal the released file", built == released_bytes,
         f"released {hashlib.sha256(released_bytes).hexdigest()[:8]}, {len(released_bytes):,} bytes"),
        ("control: dropping one guide changes the digest", control_sha != built_sha,
         f"drop {dropped}: {lost} patterns lost, digest {control_sha[:8]}"),
    ]

    width = max(len(name) for name, _, _ in checks)
    for name, ok, detail in checks:
        print(f"  {name:<{width}}  {'ok ' if ok else 'NO '}  {detail}")
    print()
    print("shared patterns, with the guides whose libraries both contain them:")
    for p in sorted(shared):
        print(f"  {p}  {'  '.join(shared[p])}")

    if args.write:
        if os.path.exists(args.write) and any(
                os.path.exists(r) and os.path.samefile(args.write, r) for r in (args.patterns, DEFAULT)):
            print(f"\nNOT WRITTEN: {args.write} is a released file")
            return 2
        with open(args.write, "wb") as fh:
            fh.write(built)
        print(f"\nwrote {args.write}")

    good = all(ok for _, ok, _ in checks)
    print()
    print("VERDICT:", "PASS" if good else "FAIL")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
