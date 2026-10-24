#!/usr/bin/env python3
"""The panel behind the body's 1,628, rebuilt and checked end to end from its guide list.

WHAT THE BODY CLAIMS. Right after the pooled theorem, which proves the real 110-guide panel needs
1,629 radius-4 balls, the text adds that the value is a property of that panel and not of the
shape: "there are 110 guides of length 20 with pairwise disjoint libraries, run counts 9 to 20 and
1,629 patterns, whose union of radius-4 balls is covered by 1,628."

WHY THIS FILE EXISTS. The original producer could not run as shipped: it opened a bare
"synth_guides.txt" from the current directory, while the list is named space_synthetic_guides.txt,
and it printed its findings and exited 0 whatever they said. It was never released.

CORRECTION, 2026-10-23. An earlier version of this file said the list "could not be found" and that
two of the producer's checks were false FOR ANY INPUT. BOTH STATEMENTS WERE WRONG. The list exists
under the other name; a search for it used a broad pattern on a narrow path and a narrow pattern on
the broad one, so it was missed. And the panel is written over {0,1,2,3}, the relabelled alphabet
the proof works in, so the producer's membership checks were consistent with its real input all
along: p and all sixteen witnesses ARE in that library. The checks only failed against a rebuilt
panel written over {A,C,G,T}. The original list is what is released here.

WHAT THIS CHECKS, ALL OF IT, IN ONE RUN:

  1. 110 guides, each of length 20.
  2. Libraries pairwise disjoint: the total pattern count equals the union's size.
  3. 1,629 distinct patterns.
  4. Run counts spanning 9 to 20.
  5. The redundant pattern p and its 16 witnesses are IN the library. This is the step the original
     could not perform, and it is what makes the rest a statement about THIS panel.
  6. B_4(p) lies inside the union of the 16 witnesses' balls, exhaustively over all 341,716 words.
  7. The control: dropping one witness must leave part of B_4(p) uncovered, or 6 proves nothing.

Together 5, 6 and 7 give the claim: p's ball is already covered by other patterns of the library,
so the remaining 1,628 patterns' balls still cover the whole union.

Usage:  pooled_counterexample.py [GUIDES.txt]
Default input: data/synthetic_panel_110.txt beside this script. Prints the input's SHA-256, so a
run log shows which list was checked.
Exit 0 only when every check above holds; 1 otherwise; 2 if it could not run.
"""
import hashlib
import itertools
import os
import sys
from collections import Counter

LEN, N_GUIDES, TARGET_PATTERNS, RMIN, RMAX, RADIUS = 20, 110, 1629, 9, 20, 4
ALPHABET = "0123"   # the relabelled alphabet the panel and the proof use
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "data", "synthetic_panel_110.txt")
P = "0101010101010101010"                     # the redundant pattern, in those letters
SWAP = {"0": "1", "1": "0"}
WITNESSES = [P[:j] + SWAP[P[j]] + P[j + 1:] for j in range(16)]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def run_count(g):
    return 1 + sum(g[i] != g[i - 1] for i in range(1, len(g)))


def patterns(g):
    """The distinct single-deletion reads of g: one per maximal run."""
    out, i = set(), 0
    while i < len(g):
        j = i
        while j + 1 < len(g) and g[j + 1] == g[i]:
            j += 1
        out.add(g[:i] + g[i + 1:])
        i = j + 1
    return out


def ball(word, radius):
    out = set()
    for r in range(radius + 1):
        for pos in itertools.combinations(range(len(word)), r):
            choices = [[c for c in ALPHABET if c != word[p]] for p in pos]
            for subs in itertools.product(*choices):
                x = list(word)
                for p, c in zip(pos, subs):
                    x[p] = c
                out.add("".join(x))
    return out


def main(argv):
    path = argv[1] if len(argv) > 1 else DEFAULT
    try:
        raw = open(path, "rb").read()
    except OSError as exc:
        print(f"NOT RUN: {exc}")
        return 2
    # Read as bytes so the digest is of the file itself; split lines as text-mode open() would.
    text = raw.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    guides = [l.strip() for l in text.split("\n") if l.strip()]
    print(f"input  {path}")
    print(f"       SHA-256 {hashlib.sha256(raw).hexdigest()}")

    library, total = set(), 0
    for g in guides:
        p = patterns(g)
        library |= p
        total += len(p)
    counts = Counter(run_count(g) for g in guides)

    checks = [
        (f"{N_GUIDES} guides", len(guides) == N_GUIDES, f"{len(guides)}"),
        (f"every guide length {LEN}", all(len(g) == LEN for g in guides),
         f"{sorted({len(g) for g in guides})}"),
        ("libraries pairwise disjoint", len(library) == total, f"{total - len(library)} shared"),
        (f"{TARGET_PATTERNS} distinct patterns", len(library) == TARGET_PATTERNS, f"{len(library)}"),
        (f"run counts span {RMIN} to {RMAX}", min(counts) == RMIN and max(counts) == RMAX,
         f"{min(counts)}..{max(counts)}"),
        ("p is in the library", P in library, f"{P in library}"),
        ("all 16 witnesses in the library", all(w in library for w in WITNESSES),
         f"{sum(w in library for w in WITNESSES)} of 16"),
    ]

    b = ball(P, RADIUS)
    missed = [x for x in b if min(hamming(x, w) for w in WITNESSES) > RADIUS]
    control = [x for x in b if min(hamming(x, w) for w in WITNESSES[1:]) > RADIUS]
    checks.append((f"B_{RADIUS}(p) inside the witnesses' union", not missed,
                   f"{len(missed)} uncovered of {len(b):,}"))
    checks.append(("control: dropping one witness uncovers part of it", bool(control),
                   f"{len(control)} uncovered"))

    width = max(len(name) for name, _, _ in checks)
    for name, ok, detail in checks:
        print(f"  {name:<{width}}  {'ok ' if ok else 'NO '}  {detail}")

    good = all(ok for _, ok, _ in checks)
    print()
    if good:
        print(f"So p's ball is covered by other patterns of this library: the remaining "
              f"{TARGET_PATTERNS - 1:,} balls still cover the union.")
    print("VERDICT:", "PASS" if good else "FAIL")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
