# The four trial counts, recomputed rather than quoted

This folder reproduces four figures the paper prints as the sizes of checks it ran: `54,000`,
`3,000`, `20,000` and `1,100`. They are counts of cases tested, so the only honest way to support
them is to run the cases again and count.

One file, standard library only, tested with Python 3.9.6. About a second.

## What it checks

| figure | what it counts |
|---|---|
| `54,000` | the b-block deletion proposition at `n = 20`, `b = 1` to `18`, 3,000 spacers for each: every word of `B_b(g)` within Hamming distance 1 of a wrapper word |
| `3,000` | the spacers drawn per `b`, so that `54,000 = 18 x 3,000` |
| `20,000` | the closed form of the wrapper's RNA total, `1 + #{i in [2, n-b-1] : g_i != g_{i+b}}`, over `n = 3` to `30` and `b = 1` to `n-2`, half the draws over four letters and half over two |
| `1,100` | the 110 panel guides at `b = 1` to `10`, checking that the deficit equals the indicator sum, in both PAM orientations |

## What to run

The panel is the 110 guide sequences, one per line, each ending `NGG`. It is not redistributed
here; it is the guide list of the study the paper analyses.

```
python3 -I trial_counts.py sgrna_panel110.txt
```

A pass prints one line per part and ends `RESULT PASS`, exit 0. Without an argument it exits 1 with
a usage line rather than guessing a path.

## The controls are the point

A count of 54,000 successes out of 54,000 attempts is worthless on its own: a check that cannot
fail will report that every time. Three controls run beside it and each is meant to fail.

- The same statement taken over the ball instead of the sphere **fails at `b = 2` and `b = 3`**,
  200 of 200 each.
- A threshold of 0 instead of 1 **fails 184** of 200 at `b = 1`.
- A named spacer is shown with `|B_2| = 6` against a ball of 21, so the two quantities are visibly
  different objects rather than two names for one.

If those three stop failing, the main count has stopped meaning anything.

## What this does not establish

It counts cases, not correctness of the propositions those cases test. The proofs are in the paper;
this says the checks the paper describes were run at the sizes it states, and that they pass at
those sizes.

The draws are seeded, so the counts are the same on every run. That makes the figures reproducible
and it also means a reader is checking the same sample we did, not an independent one.
