#!/usr/bin/env python3 -I
"""Gate for the scattered-transversal theorem: the construction holds at every length, not only at 20.

The theorem is stated for every spacer length n >= 20, while the construction is written with a
cycle taken mod 19, a modulus that looks tied to n = 20. The proof that carries it to larger n is a
monotonicity argument, not a finite check: every support lies inside a fixed 19-coordinate window,
coordinates outside that window are overwritten by no witness, and an untouched coordinate can only
ADD to a distance. So the value computed on the cycle is a lower bound at every length, attained at
n = 20 where the chain fills the cycle.

A reader is entitled to disbelieve that sentence, so this recomputes the two quantities the proof
needs, directly from the construction as the paper specifies it, at six lengths. It takes no input
files and depends on nothing in the tree, so it runs from a clean copy.

A gate that cannot fail proves nothing, so this ends with a negative control: one overwritten letter
is reverted and the same checks must then report a violation. If the control passes, the gate says
so and exits non-zero, because a checker that accepts a broken construction has not checked anything.

Construction, from the body's method section:
  tau_k(s) = 1 + ((k + s - 1) mod 19), offsets 4, 6, 9, 12, 16 carrying letters 2, 3, 2, 3, 2;
  budget m uses the first m offsets, so m = 4 omits s = 16;
  pattern D_k = 1^k 0^(N-k) with N = n - 1, and W_k is D_k overwritten at each tau_k(s).

Exits 0 when every length meets the requirement and the control fails as it must, 1 otherwise.
"""
import sys

OFFSETS = [(4, 2), (6, 3), (9, 2), (12, 3), (16, 2)]
LENGTHS = (20, 21, 25, 30, 40, 60)
CYCLE = 19


def build(n, r, m, revert=None):
    """Patterns and witnesses at length n for run count r and budget m."""
    width = n - 1
    pats, wits = [], []
    for k in range(r):
        pattern = [1] * k + [0] * (width - k)
        witness = list(pattern)
        for i, (offset, letter) in enumerate(OFFSETS[:m]):
            position = (k + offset - 1) % CYCLE
            # revert names an offset to leave unwritten, which is the negative control.
            if i != revert:
                witness[position] = letter
        pats.append(pattern)
        wits.append(witness)
    return pats, wits


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def measure(n, r, m, revert=None):
    """(radii seen, minimum pairwise distance, minimum private distance)."""
    pats, wits = build(n, r, m, revert)
    radii = {hamming(wits[k], pats[k]) for k in range(r)}
    pairwise = min(hamming(wits[k], wits[j])
                   for k in range(r) for j in range(k + 1, r))
    private = min(hamming(wits[k], pats[j])
                  for k in range(r) for j in range(r) if j != k)
    return radii, pairwise, private


def main():
    failures = 0
    print("A 2m-scattered transversal needs radius exactly m, pairwise > 2m and private > m.")
    print(f"{'m':>2} {'n':>4} {'r':>4} {'radius':>8} {'pairwise':>9} {'private':>8}  result")
    for m in (4, 5):
        need_pair, need_priv = 2 * m + 1, m + 1
        for n in LENGTHS:
            # The largest run count the length admits, which is the hardest case.
            radii, pairwise, private = measure(n, n, m)
            ok = radii == {m} and pairwise >= need_pair and private >= need_priv
            failures += not ok
            print(f"{m:>2} {n:>4} {n:>4} {str(sorted(radii)):>8} {pairwise:>9} {private:>8}  "
                  f"{'holds' if ok else f'FAILS, needs {need_pair} and {need_priv}'}")

    # The control runs at m=5 and prints numbers that would PASS at m=4, directly beneath a table
    # whose last rows are m=5. A reader checking the control against the wrong budget calls it
    # broken, which happened on first reading, so every number here carries its budget and its
    # threshold on the same line.
    cm, cn = 5, 20
    need_pair, need_priv = 2 * cm + 1, cm + 1
    print(f"\nNegative control at m={cm}, n={cn}: leave the first offset unwritten, "
          f"which must break the separation.")
    radii, pairwise, private = measure(cn, cn, cm, revert=0)
    broke = not (radii == {cm} and pairwise >= need_pair and private >= need_priv)
    print(f"  m={cm}: radius {sorted(radii)} (needs [{cm}]), pairwise {pairwise} "
          f"(needs >= {need_pair}), private {private} (needs >= {need_priv})")
    if broke:
        print("  the control fails as it must, so these checks can tell a good construction "
              "from a broken one")
    else:
        print("  THE CONTROL PASSED. This gate accepts a broken construction and proves nothing.")
        failures += 1

    print(f"\n{failures} failure(s).")
    if failures == 0:
        print("The construction meets the requirement at every length checked, so the theorem's "
              "reach beyond n = 20 rests on the argument above and not on a check at one length.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
