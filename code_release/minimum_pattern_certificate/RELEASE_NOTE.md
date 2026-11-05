# Guide-level minimum pattern certificate

Local staging only; not pushed to the public deposit.

## What a tool author runs

```sh
python3 -I min_patterns.py GUIDES.txt --budget 4
python3 -I min_patterns.py GUIDES.txt --budget 4 --audit EMITTED.tsv
python3 -I verify_minimum_pattern_certificate.py
python3 -I mutations.py
python3 -I audit_mutations.py
```

On every guide of length at least 20 and budget at most 9, the run count
`r = 1 + |{i : g_i != g_{i+1}}|` is reported as the exact minimum number of
equal-radius DNA-centred patterns for a complete single-deletion search
(companion Theorem 5). Outside that range the tool reports the run count and
`NOT PROVED`; it does not guess.

`--audit` classifies an emitted pattern list per guide as `EXACT`, `INCOMPLETE`
(loses sites), `REDUNDANT`, or `INVALID`.

## Producer and checkers

| File | Role |
|---|---|
| `min_patterns.py` | Producer |
| `verify_minimum_pattern_certificate.py` | Checker (51 claims on deposited panels + DNA-only exhaustive mechanism tests) |
| `mutations.py` | Six planted input corruptions; all rejected |
| `audit_mutations.py` / `audit_fixtures/` | Four emitted-set corruptions flipping each verdict |

## Limit recorded, not papered over

At settleable small lengths, every sound DNA-centre set found had size exactly `r`.
Extra centres hit the searchable union but are unsound. A larger symmetry-reduced
search remains open. The field value of this package is the O(n) certificate and
audit, not a re-proof of Theorem 5.

## Private-witness measurement (academic_lead_r6, 1 Nov ~18:20 Pacific)

At budget 4, every omission of one deletion word on all 171 deposited guides
leaves an explicit DNA private witness in U_4. This shows the classical library
is minimal among its own subsets; it is not a re-proof of K over unrestricted
centres.

```sh
python3 -I private_witnesses.py --budget 4
python3 -I verify_private_witnesses.py private_witnesses_m4.json
python3 -I mutations_private_witnesses.py
```

Recorded negative (unchanged): at settleable lengths the sound DNA-centre set
still has size exactly r; |sound|>r remains open for symmetry-reduced search.

## |sound| > r settleable witness (academic_lead_r6, 1 Nov ~19:10 Pacific)

**Settled:** there exist settleable DNA instances with |sound| > r, Um a proper
nonempty subset of the space, and m > n-3 (outside the companion unique-sound
range). Primary witness: guide `AACCG`, n=5, m=3, r=3, |sound|=4 (extra `ACCC`),
Um=220/256, explicit DNA 2-cover `{AAAC, ACCG}`.

```sh
python3 -I sound_gt_r_witness.py
python3 -I verify_sound_gt_r_witness.py
python3 -I mutations_sound_gt_r.py
```

**Negative (precise):** no settleable instance found with |sound|>r *and* no
r-1 DNA cover of Um. Outside the theorem range, every checked |sound|>r
witness admitted an explicit (r-1)-cover. The certificate's "no r-1 cover"
exhaustive teeth therefore stay on |sound|=r instances with m <= n-3.
