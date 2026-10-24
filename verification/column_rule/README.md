# The column rule, and the counting bound that settles m = 6

This folder checks the appendix paragraph "A counting bound that settles m = 6 without any search",
which rules out every scattered transversal, cyclic or otherwise, at `N = 4m - 1` for `m >= 6`.

The whole argument turns on one inequality about a single coordinate, so that inequality is checked
here by enumeration rather than taken on the argument's word.

One file, standard library only, tested with Python 3.9.6. It reads nothing and finds nothing
outside itself, so you can run it from any working directory.

## What it checks

| File | What it computes | Paper claim it backs |
|---|---|---|
| `check_column_rule.py` | Four parts. (A) The column rule over every column at `w = 2` to `6`: every split of ones and zeros, every choice of which entries are overwritten, and every letter each can take over the four-letter alphabet. (B) A control: the same rule with its slope weakened from `w-2` to `w-3`. (C) The structure of five consecutive staircase words. (D) The arithmetic from those two facts to the bound. | "The column rule it rests on was verified here over every column at w = 2 to 6 ... 36,400 in all, with no violation; weakening its slope from w-2 to w-3 produces 2,615 violations", and the derivation `10(2m+1) <= 20 + 15m + (N-4)`, that is `N >= 5m-6`. |

## What to run, and what a pass looks like

```
python3 -I check_column_rule.py
```

About half a minute. A pass prints four `PASS` lines and `ALL PARTS PASS`, and exits 0:

- **[A]** `36,400 columns at w = 2 to 6, 0 violations, 55 (w,j,t) cases attain the bound`. The
  second number is the one that matters; the third says the bound is attained often enough not to
  be loose.
- **[B]** `slope w-3 instead of w-2 gives 2615 violations`. This part **fails if that count is
  zero**, because a rule nobody can violate proves nothing. It is the control, and it is meant to
  find violations.
- **[C]** `five consecutive staircase words always split in exactly four columns contributing
  4+6+6+4 = 20, with N-4 constant`, checked at several lengths and every starting index.
- **[D]** the arithmetic, and that `5m-6` exceeds `4m-1` exactly from `m = 6`.

## What this does not establish

The bound is a lower bound on length. It says no transversal exists below `5m-6`; it does not say
one exists at or above it, and the witnesses that show existence are checked elsewhere, in
`verification/uniqueness/`.

It also does not give the classification at `m = 5`, which is a separate argument: this count rules
lengths out, and says nothing about which overwrite sets achieve the ones that remain.
