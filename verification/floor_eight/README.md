# The floor at eight: the weight certificate for twenty runs

These files check the body's theorem "The floor at eight, for every guide of length twenty": at
mismatch budget m = 8, the complete single-deletion library of a spacer of length 20 with r runs
cannot be covered by fewer than r radius-8 balls, centred anywhere.

After the relabelling of the general floor theorem (checked by
`verification/floor_theorem/canonical_relabelling_lemma.py`), the library with r runs is the
staircase 1^k 0^(19-k), k < r. Write U_r for the union of its radius-8 balls. The certificate is a
weight on 46,100 words of U_20, of mass 30348/1573, that gives every radius-8 ball, centred
anywhere, at most one unit. So U_20 needs 20 balls. For r < 20 the body gives the one-line step: a
cover of U_r by r - 1 balls, together with the balls at the other 20 - r staircase words, would
cover U_20 with 19. The scripts here check the certificate; that step is not a computation.

## Files

| File | What it is | What it checks |
|---|---|---|
| `payload.bin` | The certificate: 564 templates over {0, 1, *}, 5 bytes each (a base-3 index and an integer coefficient over 1573). A `*` stands for either letter outside {0, 1}, with the weight split evenly. SHA-256 2ab76785c7602e5c9bddfbcb135e1b1f1d3326637de2bd1b7836c08a96ceff4e. | Nothing; it is the data. |
| `dual20.js` | The exhaustive check, in Node.js. It carries its own copy of the payload and refuses to run unless that copy's SHA-256 equals the one above. | Every one of the 3^19 = 1,162,261,467 centre classes, by an exact transform in modular arithmetic: the largest ball mass is exactly 1 (integer score 402,688 = 256 x 1573), at centre class 1001010000001000000, and the mass is 30348/1573, so 1573 x 19 falls short of it by 461. |
| `independent_check.py` | A second implementation in Python that shares nothing with `dual20.js` but the payload. Standard library only. | By direct expansion and Hamming distance, with exact fractions: every support word lies in U_20, within 8 of the staircase (it prints the distances it finds; `dual20.js` asserts they are exactly 8); the mass is 30348/1573; the ball at the attaining centre holds exactly 1; and a **sample** of 140 other centres (the staircase and 120 random ones) holds at most 1. **The sample is not a sweep and this script does not establish the maximum; `dual20.js` does.** Control A raises one coefficient that touches the attaining ball and must push it above 1; control B raises one that does not and must leave it at 1. |
| `r7_certificate.py` | A smaller certificate for the seven-run staircase at N = 19, m = 8, written into the script. | Mass exactly 7, and every centre class (64 binary prefixes by 14 counts, 896 in all) at most 1, so seven balls are needed. |
| `impossibility.py` | Two volume-type bounds, written into the script. | The binary-support bound at N = 19, m = 8 first falls to r - 1 at r = 13 (so that support cannot certify 13 runs or more); and the unrestricted bound 4^19 / V_12 = 68719476736/12017379241, below 6, the "5.72 at m = 12" ceiling of the body. |
| `SHA256SUMS` | Digests of the five files above. | `shasum -a 256 -c SHA256SUMS` |

## Running them

    shasum -a 256 -c SHA256SUMS
    python3 independent_check.py
    python3 r7_certificate.py
    python3 impossibility.py
    node dual20.js

The Python scripts need Python 3.9 or later and read files only from their own folder, so they can
be run from anywhere. All three refuse to run under `python3 -O`, which would switch their
assertions off, and exit 2. `independent_check.py` takes about 10 s and also exits 2 if
`payload.bin` is missing.

`dual20.js` needs Node.js 18.14 or later (worker threads, a shared array and `os.availableParallelism`). **It holds one 32-bit
value per centre class, about 5 GB of memory at peak**, and uses up to 24 threads, taking about 15 s on 24. It reads no file. It
exits 1 if any check fails.
