# Scotoma

Off-target candidate enumeration directly in haplotype-resolved assemblies, including in target
windows absent from a reference.

Every variant-aware off-target screen in current use is anchored to a linear reference. The good
ones substitute bases, handle indels, and incorporate population variants from large panels. But a
candidate site lying in sequence that is **absent** from the reference is not expressible as a
variant of that reference, so it is invisible by construction rather than by oversight. FDA's
Casgevy postmarketing requirement, for instance, frames in-silico off-target reanalysis over
variants above an allele-frequency threshold in continental population groups, and the operative
criterion is defined on match *identity*, by reducing mismatches or lowering gaps. Sequence with no
reference coordinate does not enter that definition.

Scotoma enumerates there. Given a haplotype-resolved assembly, a reference, and a guide panel, it
partitions the assembly's SpCas9 target windows into those the reference contains and those it does
not, scores both against the same panel, and reports the candidate density of each.

**It counts candidate sites by sequence similarity. It does not predict cleavage.** Whether a
candidate is cut depends on chromatin, expression, delivery and enzyme kinetics, none of which are
inputs here. Nothing in this repository is a safety claim about any guide or therapy.

## Install

Python 3.9 or newer and numpy. Nothing else.

```
pip install -e .
```

Or run it in place, with no install at all:

```
python3 -m scotoma.cli --help
```

## The four commands

### 1. Build a reference universe, once

A universe is a sorted array of the distinct canonical target windows a reference contains, plus
each one's both-strand multiplicity. Building it once turns every later membership question into a
binary search.

```
scotoma universe --fasta GRCh38.fa.gz --label GRCh38 --out-npz GS_GRCh38.npz
scotoma universe --fasta chm13v2.0.fa.gz --label CHM13 --out-npz GS_CHM13.npz
```

The output is a plain `.npz` holding `codes`, `mult` and a JSON `meta`, readable by anything that
can read numpy. `--primary-only` skips records whose names contain `_` or `.`; it is off by default,
because that is a naming convention rather than a guarantee and silently dropping sequence from a
reference would inflate every absence count downstream.

### 2. Measure the blind spot

```
scotoma scan --assembly HG00235.hap1.fa.gz \
           --universe GS_GRCh38.npz \
           --panel panel.txt \
           --out-json HG00235.hap1.blindspot.json
```

The panel is one guide per line, either 20-mer spacers or full 23-mer windows. `--changeseq` reads
a CHANGE-seq style measurement table instead and takes the panel from its zero-mismatch,
zero-bulge on-target rows.

Output is a table of candidate counts in each arm at each mismatch threshold and the per-guide
burden:

```
  mm        blind      visible    per guide
   4       21,603       26,311         8.72
```

The tool also prints the ratio of the two arms. **Do not interpret that ratio as a property of the
sequence.** We measured it under four guide panels with the arms, arm sizes, universe, PAM scope and
spacer-GC distribution held identical, and it moves monotonically from 0.26 to 2.64 with nothing
changed but where the panel's guides came from. A matched-arm candidate-density comparison is
dominated by the panel's provenance. The per-guide burden and the out-of-reach share do not involve
a comparison between arms and are unaffected.

### 3. Ask what a second reference recovers

```
scotoma rescue --assembly HG00235.hap1.fa.gz \
             --universe GS_GRCh38.npz \
             --second-universe GS_CHM13.npz
```

Reports the share of the haplotype's reference-blind windows that a second reference contains, and
the residual blind fraction that survives the migration.

### 3b. What a reference-based nominator would report

```
scotoma reference-count --universe GS_GRCh38.npz --spacer CTAACAGTTGCTTTTATCAC
```

Counts candidates for one guide **in the reference**, at every mismatch threshold, reporting both
distinct canonical windows and multiplicity-summed genomic sites. This is the incumbent quantity:
the `N` a validation campaign is sized against.

### 4. Turn a candidate count into a validation budget

The count below is an NGG-scope figure. Counts under NGG and under NGG-or-NAG are not
interchangeable, so pass a count from the same PAM scope as the rest of your analysis.

```
scotoma budget --n 11096 --k 1 5 10 25 --delta1 0.05
```

A screen that assays a subset of a guide's candidate sites and finds nothing has bounded, not
eliminated, the risk. `m*(N, k, delta1)` is the smallest number of sites that must be assayed before
missing all of `k` genuinely cleaved sites has probability at most `delta1`. Every sample size is an
exact integer search in rational arithmetic; a float `--delta1` is refused outright, because 0.05 is
not 1/20 in binary and a sample size off by one at a published confidence level is a defect rather
than a rounding detail.

## The design, and why it is the only axis that answers the question

The measurement is a matched-arm comparison inside a single genome.

- **arm A, blind**: every target window in the assembly whose canonical code is absent from the
  reference universe.
- **arm B, visible**: an equal-sized uniform random sample of that same assembly's windows that
  *are* present in the reference.

Both arms are scored against the same panel, so guide composition cancels exactly. Both come from
the same individual, so population structure cancels. The arms are equal in size, so the ratio of
counts is a ratio of densities with no normalisation step to get wrong.

This matters because the obvious alternative is wrong. An earlier version of this comparison varied
the *guides* instead of the targets, pitting a clinical panel against random genomic protospacers.
Clinical guides turned out to be about 15 GC points richer than random protospacers, so every such
comparison measured composition rather than addressability, and four separate null models built that
way had to be discarded. Varying the target set while holding the panel fixed is the only axis on
which the question is answerable, and Scotoma implements only that axis.

### Oriented versus canonical, the one rule to carry

Two operations need two different representations of a site, and conflating them silently corrupts
a comparison:

- **Membership** in a reference is strand-free, so it is tested on the canonical code,
  `min(code, revcomp(code))`.
- **Mismatch distance** to a spacer is strand-specific, so it must be computed on the oriented
  code, the one that actually reads `[spacer]NGG`.

About 58% of canonical codes are the reverse complement of their oriented form, so a canonical code
begins with `CC` roughly half the time. That is a compositional constraint no real spacer has. Using
canonical codes for distance imposes it on whichever arm was canonicalised and not on the other.
`scotoma.coding.oriented_sites` returns oriented codes and canonicalisation happens only inside the
membership test. Two tests in `tests/test_coding.py` are built specifically on windows where the
oriented and canonical forms differ, because the obvious test cases are ones where they coincide and
so cannot catch the bug.

### Controls available on `scan`

| flag | what it tests |
|---|---|
| `--match-on gc` | draws the visible arm to match the blind arm's spacer-GC histogram bin by bin |
| `--match-on complexity` | matches instead on distinct 3-mers per window, which GC matching cannot control and which separates addressability from repeat composition |
| `--shuffle-guides` | permutes each spacer, preserving composition exactly: tests whether a result needs the real guide sequences |
| `--panel-from-reference` | replaces the panel with real reference-present windows, GC-matched and drawn from distinct low-multiplicity codes |
| `--panel-from-blind` | the same, drawn from the reference-absent arm, so provenance is reversed |
| `--replicates N` | redraws the visible arm N times, separating draw variance from donor variance |

Under the three synthetic-panel flags the on-target recovery gate is **inverted**, since a permuted
spacer has no on-target anywhere and a drawn window is removed from its own arm. Drawn panels are
taken from distinct codes rather than window instances, because sampling instances selects sequences
in proportion to how often they occur and yields a panel of high-copy repeats.

### Gates

`scotoma scan` exits non-zero on gate failure rather than folding a broken run into the numbers.

| gate | what it catches |
|---|---|
| arms equal in size | a ratio that is not a density ratio |
| positive control: zero-mismatch hits in the reference-present pool | the panel and the universe are not from the same reference, or the distance kernel is not recovering on-targets |
| no zero-mismatch hit in the blind arm | a broken partition; an on-target cannot be absent from the reference it was read off |
| at least 90% of windows reference-present | this is not a human haplotype against a human reference |

The positive control runs on the **full** reference-present pool by default. Running it on the
subsampled visible arm, which `--weak-control` still allows, is unreliable: the arm is matched to the
blind arm, so when the blind arm is a few per cent of the genome an on-target survives the draw only
a few per cent of the time and a correct pipeline fails by luck. That defect was found by running
the CLI rather than the library, and `tests/test_blindspot.py` pins it.

## One import note

`scotoma/__init__.py` exports a function named `blindspot`, which shadows the same-named submodule
on the package, so `import scotoma.blindspot as bs` binds the **function**. This is conventional
Python, the same shape as `datetime.datetime`. To reach the module, for example to override
`CONTROL_MAX_WINDOWS`, use `importlib.import_module("scotoma.blindspot")`.

## Tests

```
python3 -m pytest tests/ -q
```

142 tests, no network and no genome on disk: every fixture is synthetic, with the reference and the
haplotype constructed so the partition has a known right answer. One of the 142 is a cross-check
against the optional `m0_search` companion and skips when that module is not importable, so a clean
clone reports 141 passed and 1 skipped.

The suite was checked by mutation rather than by its own green run. Three deliberate bugs were
introduced and all three were caught: including the PAM in the distance (4 failures), returning
canonical instead of oriented codes (4 failures), and dropping canonicalisation from the membership
test (11 failures). Several tests exist because a defect reached a real run first and is now pinned:
the positive control evaluated on a subsampled arm, a shuffled panel failing an on-target gate it
cannot satisfy, bytes-encoded universe metadata, and a panel drawn over window instances rather
than distinct codes. The budget arithmetic is additionally cross-checked, when the companion
`m0_search` module is importable, against an independently written exact implementation; the two
agree on every case tested.

## Limitations, stated rather than buried

- Candidate sites are counted by **sequence similarity only**. No cleavage, no chromatin, no
  activity prediction.
- An NGG PAM is assumed. Other PAMs and other nucleases are not implemented.
- The `scan` pipeline counts mismatches only, and a site reachable only through a bulge is not
  enumerated by it. Single-base bulges are modelled separately, by `scotoma.bulge`: it expands a
  spacer into bulge-bearing patterns using the Cas-OFFinder-bulge construction of Jeongbin Park and
  Sangsu Bae and hands those to an exact matcher, which is the route the bulge figures come from.
  Bulges of more than one base, and indels elsewhere in the protospacer, are not enumerated at all.
- The comparison is **per haplotype**. Haplotypes share most of their sequence, so per-haplotype
  figures do not add to a per-person figure, and the union across an individual's two haplotypes is
  not computed.
- Statistical power falls away at tight thresholds. At three mismatches and below, counts are small
  enough that a single genome supports no useful inference.
- `m*` certifies "fewer than k cleaved sites among the enumerated candidates". It licenses "none"
  only under the further assumption that an unclean guide cleaves at least k enumerated sites, which
  is a claim about the enzyme that this code cannot check. Candidates in windows absent from the reference
  are outside `N` and outside the guarantee, which is the whole point of `scotoma scan`.

## License

MIT. See `LICENSE`.

## Minimum covering number: verification scripts

`verification/floor_theorem/` holds self-contained scripts that check every
quantitative claim in the minimality result, including the counterexamples that
bound its parameter range. They need only the standard library and print what
they checked. See that directory's README for the table of what each one
establishes.
