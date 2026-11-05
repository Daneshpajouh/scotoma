# A printed private witness at budget six, on the clinical guide sg1617

Paper B proves that for every guide of length at least twenty and every budget at most nine the
exact minimum number of equal-radius patterns is the run count r, with the competing centres
unrestricted. That is a statement about all centre placements, proved by weighting and by integer
obstructions, and it is not proved by exhibiting a witness.

This directory answers a different and smaller question, asked by an external reviewer: make the
consequence tangible. Show one concrete set of r-1 centres, at a budget inside the range only
this paper certifies, with a printed word that those r-1 balls do not cover, on a guide whose
off-target list a reader already has a reason to care about.

## The witness

    guide     CTAACAGTTGCTTTTATCAC   sg1617, the BCL11A +58 erythroid enhancer guide
                                     used in the exa-cel trials, and the guide Cancellieri
                                     et al. searched at up to six mismatches with two bulges
    runs      15                     C T AA C A G TT G C TTTT A T C A C
    r         15                     distinct single-deletion words, one per run
    omitted   CTAACAGTTGCTTTTATCA    the terminal deletion, which the released build never emits
    witness   ACCGGCGTTGCTTTTATCA

    distance from the witness to the omitted word        6
    distance from the witness to each of the other 14    7 or more

So the witness lies in U_6, because it is within the budget of a library word, and no radius-6
ball of the other fourteen library words contains it. Dropping one deletion word at budget six
therefore loses coverage on this guide, and the loss is a nineteen-character string a reader can
check by hand in a few minutes.

## What this is not

It is not the theorem. The theorem is about arbitrary centres, and a search that replaced the
library with fourteen cleverly placed centres elsewhere in sequence space is exactly the case a
witness argument cannot settle. It is also not a claim about any assembly: the witness is a word
of sequence space, and whether it occurs in a genome is a separate question this directory does
not address. The companion paper reports the genomic figure, 129 measured-active sites on 40
guides at budget four.

## Replay

    python3 -I verify_budget_six_witness.py

15 claims, all reproduced. Everything is recomputed from the twenty-character guide string: the
checker rebuilds the deletion library and the run count rather than reading any stored list.

Controls, each of which had to be able to fail:

- single-coordinate forgeries of the witness lose the property, and 45 of them do
- restoring the omitted word covers the witness, so the property is about the omission
- the achieved separation of 7 attains the triangle bound exactly, since the nearest competing
  word sits at distance 1 from the omitted one, so no witness within 6 of it can do better
- an exhaustive sweep over every six-coordinate flip construction confirms the printed witness is
  optimal rather than a lucky find

The first construction tried reached only 5 and would have supported no claim at all. The
exhaustive sweep is in the checker so that the difference is auditable.
