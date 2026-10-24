# The CFD gap penalty is asymmetric: the four-arm measurement

This folder holds the script behind the appendix sentence that begins "What it does with them is
asymmetric, and we measured it", and the output of the run that sentence reports.

The claim is that CFD's gap penalty helps on RNA bulges and hurts on DNA bulges. That is a
measurement, not an argument, so the thing a reader needs is the program that made it.

One file, standard library only, tested with Python 3.9.6.

## What it checks

| File | What it computes | Paper claim it backs |
|---|---|---|
| `cfd_gap_arms.py` | Four scoring arms over every single-bulge row of the panel, stratified by kind and total edits, with paired DeLong intervals and Holm across the strata. | "On 1,749,342 single-bulge rows, holding total edits fixed and comparing against the same score with the gap simply dropped, the RNA-bulge penalty adds +0.2749, +0.1395, +0.0810 and +0.0927 AUROC at two to five edits, every interval excluding zero after Holm. The DNA-bulge penalty subtracts 0.069 to 0.139 at three and four edits under either attribution convention." |
| `expected_output.txt` | The output of the run on 2026-10-23 that those figures come from. | The numbers above, line for line. |

The four arms are the point. Arm B drops the gap column and scores the mismatches at their true
guide positions, which is what a tool that ignores bulges computes. Arm C adds one factor for the
gap. An RNA bulge leaves the unpaired base on the guide, where it has a definite position, so there
is one answer. A DNA bulge leaves it on the target, with no guide position at all, so the script
reports both attribution conventions, C after and C before, and their difference is the measurable
cost of choosing one. C minus B is then the gap penalty's own contribution, which is the quantity
in the sentence.

Two further arms are controls and are meant to look wrong. R reversed reads the gap table from the
other end, and PI permuted sends each position to a seeded random other position. Neither should
beat arm B. On the deposited panel neither ever does: across the seven strata that clear the floor,
R reversed subtracts between 0.169 and 0.449, and PI permuted between 0.013 and 0.249. Every one of
the fourteen control measurements has the sign a scrambled position table should have, and the four
gains in the paper are the only arm that goes the other way.

## What you need first

Unlike the other folders here, this one does not carry its input. It reads the CHANGE-seq panel of
Lazzarotto et al. 2020 (Nature Biotechnology 38, 1317 to 1327, doi 10.1038/s41587-020-0555-7),
which is that study's deposit and is not ours to redistribute. Download their dataset CSV, then
pass its path:

```
python3 -I cfd_gap_arms.py CHANGEseq_CR_Lazzarotto_2020_dataset.csv
```

The columns it reads are `Align.sgRNA`, `Align.off-target`, `Align.#Bulges`, `Align.#Mismatches`
and `reads`. Rows with anything other than exactly one bulge are skipped, which is the class the
claim is about.

That study's data is distributed in more than one arrangement, under names that differ by little,
so check you have the same one before concluding anything from a difference. The check is free and
needs no hash from us: the second line of output reports how many single-bulge rows were scored,
and on the table these figures come from it reads exactly `1749342`. A different count means a
different arrangement of the deposit, not a failure to reproduce.

If no row at all is scored the script stops with `GATE FAILED` and exits 1, because scoring nothing
is an invocation error and should not be reported as a clean run.

## What a pass looks like

The first line prints the permutation seed, so the control arm is reproducible rather than merely
random. The second reports how many rows were scored. On the deposited panel that is
`single-bulge rows scored: 1749342   unscorable: 0   zeroed by an IUPAC base: 9063`, and anything
else means a different table.

Then one line per stratum, and the Holm table. Compare it against `expected_output.txt`. The four
figures in the paper are the `C after` rows of the RNA block at two, three, four and five edits.

A stratum with fewer than twenty actives or twenty inactives prints a note instead of an interval
and takes no test. That floor was fixed before the run, and the script says so where it applies
rather than dropping those strata silently.

Two of the twenty four tests do not clear Holm: DNA at five edits under C after, and RNA at two
edits under the permutation control. Both are in the output and neither is in the paper's claim,
which is confined to the strata that do clear.

## What this does not establish

It measures ranking, not biology. A higher AUROC means the score orders this panel's measured
actives above its inactives better; it does not mean a site is a genuine off-target, and no row
here was assayed by us.

It also says nothing about whether the gap penalty's magnitude is right. The comparison is against
the same score with the gap dropped, so it answers whether applying the penalty is better than
ignoring the bulge, for each kind, and nothing further.

The DNA result depends on the attribution convention only in size, not in sign: both conventions
subtract at three and four edits. The script reports both so that a reader can see that, rather
than taking our choice of one.
