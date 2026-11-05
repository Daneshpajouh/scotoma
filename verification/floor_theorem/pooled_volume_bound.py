#!/usr/bin/env python3
"""The volume lower bound on covering the pooled panel, computed from the released patterns.

WHY THIS EXISTS. Until 2026-10-24 the body printed "no fewer than 839 centres of any provenance
here at m = 4, 51.5% of the library", and that figure had no producer on record. Reconstructing it
gives

    ceil(286,613,932 / 341,716) = 839 = 51.5% of 1,629,

so it rests on a LOWER BOUND of 286,613,932 for the size of the union of the patterns' radius-4
balls. The union's exact size, computed here from the released pattern list, is 381,985,082, which
gives 1,118 and 68.6%. Both are valid lower bounds on the covering number, since a cover needs at
least ceil(|union| / |ball|) balls; 839 is simply the weaker of the two. Once this script had
actually been run the body was changed to print the stronger 1,118 (68.6%), so the figure the paper
quotes is now the exact count below. Neither reaches the 1,629 the certificate proves, which is the
point the appendix makes: no volume argument gets to the answer.

WHAT IT COSTS, SAID UP FRONT. The union is held as a set of packed 38-bit integers, one per word.
Measured on this machine 2026-10-24: 187 s, and a peak resident set of 27.4 GiB with a peak
memory footprint of 33.7 GiB, for the 1,629 patterns. An earlier note in this file said about
11 GB; that was wrong by 2.7x and would defeat a reader on a 16 GB machine, so the figure here
is now the one /usr/bin/time -l reported. That is
what an exact count costs here; there is no cheaper exact route in this script, and a sampled
estimate would not be a bound.

THE BALL SIZE IS COMPUTED, NOT QUOTED. |B_4| at length 19 over four letters is the binomial sum
sum_{r<=4} C(19,r) * 3^r = 341,716.

Usage:  pooled_volume_bound.py PATTERNS.txt [--expect-union N] [--expect-bound K]
Prints the input's SHA-256, so a run log shows which list was checked.
Exit 0 when the computation completes and matches any --expect values given; 1 on a mismatch;
2 if it could not run.
"""
import hashlib
import itertools
import math
import sys
import time

CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
RADIUS = 4


def ball_size(n, radius, q=4):
    return sum(math.comb(n, r) * (q - 1) ** r for r in range(radius + 1))


def exact_union(patterns, radius):
    """Every word within `radius` of some pattern, counted once, as packed integers."""
    n = len(patterns[0])
    union = set()
    for word in patterns:
        packed = sum(CODE[c] << (2 * i) for i, c in enumerate(word))
        for r in range(radius + 1):
            for pos in itertools.combinations(range(n), r):
                cleared = packed
                for p in pos:
                    cleared &= ~(3 << (2 * p))
                variants = [cleared]
                for p in pos:
                    original = (packed >> (2 * p)) & 3
                    variants = [v | (c << (2 * p)) for v in variants
                                for c in range(4) if c != original]
                union.update(variants)
    return union


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    if not args:
        print(__doc__)
        return 2
    try:
        raw = open(args[0], "rb").read()
    except OSError as exc:
        print(f"NOT RUN: {exc}")
        return 2
    # Read as bytes so the digest is of the file itself; split lines as text-mode open() would.
    text = raw.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    patterns = [l.strip() for l in text.split("\n") if l.strip()]
    print(f"input  {args[0]}")
    print(f"       SHA-256 {hashlib.sha256(raw).hexdigest()}")
    if not patterns or len({len(p) for p in patterns}) != 1:
        print("NOT RUN: patterns must be non-empty and all the same length")
        return 2

    n = len(patterns[0])
    b = ball_size(n, RADIUS)
    print(f"{len(patterns)} patterns of length {n}; |B_{RADIUS}| = {b:,} (computed, not quoted)")
    print("counting the exact union; about 190 s and a peak resident set near 27 GiB here",
          flush=True)

    t0 = time.time()
    union = exact_union(patterns, RADIUS)
    bound = -(-len(union) // b)
    print(f"exact union: {len(union):,} words in {time.time() - t0:.0f}s")
    print(f"volume bound ceil(union / ball) = {bound:,}"
          f"   = {bound / len(patterns):.1%} of the {len(patterns):,} patterns")

    bad = 0
    for flag, got, label in (("--expect-union", len(union), "union"),
                             ("--expect-bound", bound, "bound")):
        if flag in argv:
            want = int(argv[argv.index(flag) + 1].replace(",", ""))
            ok = got == want
            print(f"  {label}: expected {want:,}, got {got:,} -> {'match' if ok else 'MISMATCH'}")
            bad += not ok
    print()
    print("VERDICT:", "PASS" if not bad else "FAIL")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
