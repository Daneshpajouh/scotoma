# Cas-OFFinder version-3 tags: single-base RNA-bulge positions emitted (lead, 2026-10-31 13:5x PDT)

Source: git clone of github.com/snugel/cas-offinder (into the scratchpad, never executed).
Tags (commit, commit date): 3.0.0b 8ffdc8b4 (2021-01-24), 3.0.0b2 a5ed8d1b (2021-01-25),
3.0.0b3 1b7d8031 (2021-07-29), develop 0a9ac005 (2022-08-31). cas-offinder.cpp is byte-identical
between 3.0.0b and 3.0.0b2.

RNA-bulge loop, verbatim in rna_loop_3.0.0b.txt and rna_loop_3.0.0b3.txt:
- 3.0.0b, 3.0.0b2: `for (j = 1; j < sline[0].find_last_not_of('N') + 1 - i; j++)`
- 3.0.0b3, develop: `for (j = sline[0].find_first_not_of('N'); j < sline[0].find_last_not_of('N') + 1 - i; j++)`
For a 20-nt spacer followed by NNN and i = 1 (exact emulation of the loop bounds):
- 3.0.0b and 3.0.0b2 emit deletions at 0-based positions 1..18 (18 patterns): positions 0 and 19,
  the first and the last spacer base, are never deleted.
- 3.0.0b3 and develop emit 0..18 (19 patterns): position 19, the last (PAM-proximal) base, is never
  deleted. The 3.0.0b3 change restored the first position and left the last one out.
Consequence: every published version-3 tag omits the terminal RNA position; the two earlier betas
also omit the first. Evidence state: H1/H2 (exact static reading of the loop plus emulation of its
bounds). Executing the parser of 3.0.0b and 3.0.0b2 (as the paper already did for b3/develop) would
raise it; not yet done. Whether an omitted position costs a pattern depends on the guide's run
structure (deleting either base of a run gives the same word), as the paper's closed form states.
