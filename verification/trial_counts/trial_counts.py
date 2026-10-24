"""Recompute the trial counts of the appendix's block-deletion paragraphs (TRIAL-COUNTS-PRODUCER).

Independent of the paper's own checker, which is not deposited. Definitions, 1-indexed g_1..g_n:

  B_b(g)   block-deletion set: g with one contiguous block of b bases deleted, every start 0..n-b
  W_b(g)   the wrapper's RNA set: the same, starts 1..n-b-1 only (both flanks non-empty);
           this is T-34 expansion.py rna_patterns(scheme="wrapper"), slots range(1, n - b)
  ball     the b-deletion ball: every subsequence of length n - b

Checks:
  A  54,000 = 18 x 3,000: at n = 20, for b = 1..18, every word of B_b(g) lies within Hamming
     distance 1 of a word of W_b(g); 3,000 seeded random spacers per b.
  B  20,000: |W_b(g)| = 1 + |{i in [2, n-b-1] : g_i != g_{i+b}}|, with the classical
     |B_b(g)| = 1 + |{i in [1, n-b] : g_i != g_{i+b}}| checked on the same cases.
  C  1,100 = 110 guides x b = 1..10, forward and reversed reading: the deficit |B_b| - |W_b|
     equals 1[g_1 != g_{b+1}] + 1[g_n != g_{n-b}]; every W word begins with g_1 and ends with g_n;
     the omitted terminal words absent from W number exactly the indicator sum.
Controls: the ball version of A must fail for b >= 2; a checker with threshold 0 must fail;
|B_2| = 6 and a ball of 21 for TAAAGACAAT.
Exit 1 on any failed check or control.
"""

import itertools
import random
import sys

if len(sys.argv) < 2 or sys.argv[1].startswith("-"):
    sys.exit("usage: trial_counts.py sgrna_panel110.txt   "
             "(the 110 guide sequences, one per line, each ending NGG)")
PANEL = sys.argv[1]
SEED = 20261024


def block_set(g, b, lo, hi):
    return {g[:i] + g[i + b:] for i in range(lo, hi)}


def B(g, b):
    return block_set(g, b, 0, len(g) - b + 1)


def W(g, b):
    return block_set(g, b, 1, len(g) - b)


def ball(g, b):
    n = len(g)
    return {"".join(g[j] for j in keep) for keep in itertools.combinations(range(n), n - b)}


def hd(x, y):
    return sum(1 for p, q in zip(x, y) if p != q)


def covered(words, emitted, limit):
    return all(x in emitted or min(hd(x, w) for w in emitted) <= limit for x in words)


def diffs(g, b, lo, hi):
    # count of 1-indexed i in [lo, hi] with g_i != g_{i+b}
    return sum(1 for i in range(lo, hi + 1) if g[i - 1] != g[i - 1 + b])


def main():
    rng = random.Random(SEED)
    bad = []

    # A
    n = 20
    total_a = held_a = 0
    for b in range(1, 19):
        for _ in range(3000):
            g = "".join(rng.choice("ACGT") for _ in range(n))
            total_a += 1
            held_a += covered(B(g, b), W(g, b), 1)
    print(f"A  block-deletion proposition, n = 20, b = 1..18, {total_a // 18:,} spacers per b: "
          f"holds in {held_a:,} of {total_a:,}")
    if held_a != total_a:
        bad.append("A")

    ball_fail = {}
    for b in (2, 3):
        fails = 0
        for _ in range(200):
            g = "".join(rng.choice("ACGT") for _ in range(n))
            fails += not covered(ball(g, b), W(g, b), 1)
        ball_fail[b] = fails
    print(f"   control, same statement over the ball, 200 spacers each: fails b=2 {ball_fail[2]}, b=3 {ball_fail[3]}")
    if not all(ball_fail.values()):
        bad.append("A-ball-control")
    zero_fail = sum(1 for _ in range(200)
                    if not covered(B(g := "".join(rng.choice("ACGT") for _ in range(n)), 1), W(g, 1), 0))
    print(f"   control, threshold 0 instead of 1, b = 1, 200 spacers: fails {zero_fail}")
    if zero_fail == 0:
        bad.append("A-zero-control")
    t = "TAAAGACAAT"
    print(f"   control, TAAAGACAAT b = 2: |B_2| = {len(B(t, 2))}, ball = {len(ball(t, 2))}")
    if (len(B(t, 2)), len(ball(t, 2))) != (6, 21):
        bad.append("A-example")

    # B
    total_b = held_b = held_classic = 0
    for k in range(20000):
        alpha = "ACGT" if k % 2 == 0 else "AC"
        n = rng.randint(3, 30)
        b = rng.randint(1, n - 2)
        g = "".join(rng.choice(alpha) for _ in range(n))
        total_b += 1
        held_b += len(W(g, b)) == 1 + diffs(g, b, 2, n - b - 1)
        held_classic += len(B(g, b)) == 1 + diffs(g, b, 1, n - b)
    print(f"B  closed form of |W_b|, n = 3..30, b = 1..n-2, half over ACGT and half over AC: "
          f"0 failures needed, {total_b - held_b} failures in {total_b:,}; classical |B_b| "
          f"{total_b - held_classic} failures")
    if held_b != total_b or held_classic != total_b:
        bad.append("B")

    # C
    guides = [line.strip() for line in open(PANEL) if line.strip()]
    spacers = [s[:-3] for s in guides]
    assert len(spacers) == 110 and all(len(s) == 20 and q.endswith("NGG") for s, q in zip(spacers, guides))
    for arm, read in (("forward", lambda s: s), ("reversed", lambda s: s[::-1])):
        cases = agree = ends = converse = 0
        law_b1 = 0
        values = {b: set() for b in range(1, 11)}
        for s in spacers:
            g = read(s)
            n = len(g)
            for b in range(1, 11):
                cases += 1
                bw, ww = B(g, b), W(g, b)
                ind = (g[0] != g[b]) + (g[-1] != g[n - 1 - b])
                deficit = len(bw) - len(ww)
                agree += deficit == ind
                ends += all(w[0] == g[0] and w[-1] == g[-1] for w in ww)
                terminal = {g[b:], g[:n - b]}
                converse += sum(1 for x in terminal if x not in ww) == ind
                values[b].add(deficit)
                if b == 1:
                    law_b1 += deficit == ind
        spans = sum(1 for b in values if values[b] == {0, 1, 2})
        print(f"C  {arm}: {cases:,} cases, deficit = indicator sum in {agree:,}, structural ends in {ends:,}, "
              f"converse count in {converse:,}; b = 1 law {law_b1} of 110; b values spanning deficits "
              f"0, 1, 2: {spans} of 10")
        if not (cases == agree == ends == converse == 1100 and law_b1 == 110):
            bad.append("C-" + arm)

    print("RESULT " + ("PASS" if not bad else "FAIL " + " ".join(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
