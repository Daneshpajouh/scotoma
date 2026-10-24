#!/usr/bin/env python3
"""Check the redundancy that lets a pooled library be covered by fewer balls than it has patterns.

WHAT THE PAPER CLAIMS, AND WHICH HALF THIS CHECKS. The body says there are 110 guides of length 20
with pairwise disjoint libraries, run counts 9 to 20 and 1,629 patterns, whose union of radius-4
balls is covered by 1,628. That claim has two parts:

  1. A MECHANISM: some pattern of such a library is redundant, meaning its own radius-4 ball lies
     entirely inside the union of other patterns' radius-4 balls, so dropping it loses nothing and
     1,628 balls suffice. THIS SCRIPT CHECKS THAT, exhaustively, over the whole ball.
  2. A CONSTRUCTION: that a 110-guide panel with pairwise disjoint libraries and exactly 1,629
     patterns exists containing such a configuration. THIS SCRIPT DOES NOT CHECK THAT; see
     pooled_counterexample.py, which does, on the released guide list.

     CORRECTION, 2026-10-23: an earlier version of this paragraph said that guide list "could not
     be found". It could. It is space_synthetic_guides.txt, and the search that missed it used a
     broad pattern on a narrow path and a narrow pattern on the broad one. It is released now.

Saying so is the point. The earlier producer, controls/space_pooled_checks.py, could not run as
shipped: it opened a bare "synth_guides.txt" from the current directory with no path and no
argument, it printed its findings and exited 0 whatever they were, and neither it nor its input was
released. A reader following the paper's own availability clause could check none of this.

THE MECHANISM, CONCRETELY. Take the length-19 pattern p = 0101...0 and the 16 words w_1..w_16 that
differ from it in exactly one of its first 16 coordinates. Every word within Hamming distance 4 of
p is also within distance 4 of some w_j, so B_4(p) is contained in the union of the B_4(w_j) and p
is redundant in any library containing all of them.

THE CONTROL IS PART OF THE RUN, NOT A SEPARATE STEP. Dropping w_1 must leave part of B_4(p)
uncovered; if it does not, the containment test is incapable of failing and its pass means nothing.
Both assertions are checked and both must hold for a PASS.

Usage: pooled_redundancy.py
Exit 0 when the containment holds AND the control fails as it must; 1 otherwise.
"""
import itertools
import sys

ALPHABET = "0123"
RADIUS = 4
PATTERN = "0101010101010101010"      # length 19, the relabelled chain's shape
WITNESS_COORDS = 16                  # w_j flips coordinate j, for j = 0..15


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def ball(word, radius, alphabet):
    """Every word within `radius` substitutions of `word`, as a set."""
    out = set()
    for r in range(radius + 1):
        for positions in itertools.combinations(range(len(word)), r):
            choices = [[c for c in alphabet if c != word[p]] for p in positions]
            for subs in itertools.product(*choices):
                x = list(word)
                for p, c in zip(positions, subs):
                    x[p] = c
                out.add("".join(x))
    return out


def uncovered(target_ball, centres, radius):
    """Words of the ball that no centre reaches within radius."""
    return [x for x in target_ball if min(hamming(x, c) for c in centres) > radius]


def main():
    flip = {"0": "1", "1": "0"}
    witnesses = [PATTERN[:j] + flip[PATTERN[j]] + PATTERN[j + 1:]
                 for j in range(WITNESS_COORDS)]
    b = ball(PATTERN, RADIUS, ALPHABET)
    print(f"pattern length {len(PATTERN)}, alphabet {len(ALPHABET)} letters, radius {RADIUS}")
    print(f"|B_{RADIUS}(p)| = {len(b)}")

    missed = uncovered(b, witnesses, RADIUS)
    print(f"uncovered by w_1..w_{WITNESS_COORDS}: {len(missed)}"
          f"   (0 means p is redundant, so its library needs one fewer ball)")

    # The control. It must fail, or the test above cannot fail either.
    missed_control = uncovered(b, witnesses[1:], RADIUS)
    print(f"CONTROL, drop w_1:          {len(missed_control)}"
          f"   (must be > 0, else the containment test proves nothing)")

    contained = not missed
    control_fires = bool(missed_control)
    good = contained and control_fires
    print()
    print("VERDICT:", "PASS" if good else "FAIL")
    if not good:
        print(f"   containment {contained}, control fires {control_fires}")
    print("Checks the mechanism only. The 110-guide, 1,629-pattern construction is checked by "
          "pooled_counterexample.py, on the released list data/synthetic_panel_110.txt.")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
