# The witness family is forced: the enumeration behind it

This folder checks the appendix paragraph "The witness family, enumerated", which says that the
overwrite set the floor theorem's construction prints is the only one there is, up to reflection
and exchange of the two labels.

One file, standard library only, tested with Python 3.9.6. It finds nothing outside its own folder,
so you can run it from any working directory.

## What it checks

| File | What it computes | Paper claim it backs |
|---|---|---|
| `check_uniqueness.py` | For a length `N` and budget `m`, every pair of overwrite set and labelling, counting those whose witness family is `2m`-scattered. At `N = 19`, `m = 5` it also prints each valid family and its gap multiset. | "At `N = 19`, `m = 5` the overwrite set is unique up to reflection and label exchange, four of 372,096", and the surrounding counts at the obstruction boundary. |

## What to run, and what a pass looks like

```
python3 -I check_uniqueness.py
```

Searches all 372,096 pairs at `N = 19`, `m = 5` and prints the four valid families, then checks the
boundary counts at `m = 2`, `3`, `4` and the slack case at `N = 20`. A pass ends with
`VERDICT: PASS` and exits 0. Expect about a minute.

The four families it returns are the support `{4, 6, 9, 12, 16}` the construction prints, its
reflection `{4, 8, 11, 14, 16}`, and both labellings of each. All four share the gap multiset
`{2, 3, 3, 4, 7}`.

```
python3 -I check_uniqueness.py --wide
```

Adds `m = 6` at its own boundary `N = 23`, which returns none. This is the slow one: it searches
6.5 million pairs and takes several minutes.

```
python3 -I check_uniqueness.py --self-test
```

The controls. It plants two defects, a budget one too low and a length one short, and each must
move the count away from four. The length-one-short case returns zero, which is the bound
`N >= 4m - 1` appearing on its own rather than being assumed.

## What this does not establish

The enumeration covers the **cyclic** overwrite rule with two labels, which is the construction the
paper uses. It says nothing about whether some other, non-cyclic transversal exists at a budget
where the cyclic rule fails. The paper states the claim with that scope and so does this checker.
