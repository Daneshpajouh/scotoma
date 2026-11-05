#!/usr/bin/env python3
"""The sound-pattern theorem on the real panel: steps (iii) and (iv) of its proof, for every guide.

THE STATEMENT. For a guide with complete single-deletion library L (words of length N = 19), r = |L|,
and U the union of the radius-m balls around L with m = 4 <= N - 2, the only sound centres are the r
members of L. sound_patterns.py checks the theorem exhaustively at small lengths with no input. This
script checks the deployed panel for the centres that use only library letters, which are what steps
(iii) and (iv) of the proof handle. It does NOT enumerate the panel's wildcard patterns or the
centres carrying a letter no library word uses (steps (i) and (ii)); those are one-line arguments,
checked exhaustively at small lengths by sound_patterns.py, and their panel count is 5^19.

WHAT IT CHECKS, FOR EACH OF THE 110 GUIDES.
  * Every centre that uses only library letters at each coordinate: at a run-boundary coordinate
    either of its two letters, elsewhere the one letter every library word carries. That is
    2^(r-1) centres per guide, 3,697,920 over the panel, of which 1,634 are library words.
  * Each of the other 3,696,286 is shown unsound by the witness of the proof's step (iii): with u
    the first boundary coordinate carrying the right-hand letter and v the last carrying the
    left-hand one (u < v, since the centre is not a library word), change the first m coordinates
    other than u and v to a letter no library word carries there. The result is checked by direct
    distance: within m of the centre, at least m + 1 from every member of L, so outside U.
  * Each of the 1,634 library words owns a private point, built as the proof's step (iv) builds it
    and checked the same way.
  * The hypothesis: m <= N - 2 holds for every guide (4 <= 17). The weaker form m <= N - r + 1 fails
    for r >= N - m + 2 = 17 at N = 19, and the script prints how many guides that is.

THE SOURCE. The guides come from the CHANGE-seq table read the way pooled_patterns_producer.py reads
it (the sgRNA column, first 20 letters, 23-letter windows ending NGG). That table is a third-party
deposit and is not redistributed; its source and SHA-256 are in pooled_patterns_producer.py.

CONTROL, RUN EVERY TIME. The same check with one coordinate fewer changed must report failures, both
for the witnesses and for the private points.

Usage:  sound_panel.py DEPOSIT.csv
Python 3.9 or later, standard library only, about 90 s. Prints the table's SHA-256.
Exit 0 when every check holds and the control fails as it should; 1 otherwise; 2 if it could not run.
"""
import csv
import hashlib
import sys
import time

M, LEN = 4, 20
N = LEN - 1
ALPHA = "ACGT"
CODE = {c: i for i, c in enumerate(ALPHA)}
LOW_BITS = int("01" * N, 2)
EXPECTED = {"guides": 110, "centres": 3697920, "library": 1634, "others": 3696286, "weaker": 21}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def pack(word):
    return sum(CODE[c] << (2 * i) for i, c in enumerate(word))


def dist(x, y):
    v = x ^ y
    return bin((v | (v >> 1)) & LOW_BITS).count("1")


def check_guide(g, nchange):
    """Returns (centres, library words among them, others, witness failures, private failures)."""
    lib = sorted({g[:i] + g[i + 1:] for i in range(LEN)})
    lib_packed = [pack(w) for w in lib]
    lib_set = set(lib_packed)
    bnd = [p for p in range(N) if g[p] != g[p + 1]]
    k = len(bnd)
    free = [next(x for x in ALPHA if x not in (g[p], g[p + 1])) for p in range(N)]
    free_packed = [CODE[free[p]] << (2 * p) for p in range(N)]
    clear = [3 << (2 * p) for p in range(N)]
    left = pack(g[:N])                                  # every boundary coordinate at its left letter
    delta = [(CODE[g[p]] ^ CODE[g[p + 1]]) << (2 * p) for p in bnd]

    def witness(c, keep_out):
        pos = [p for p in range(N) if p not in keep_out][:nchange]
        w = c
        for p in pos:
            w = (w & ~clear[p]) | free_packed[p]
        return w

    centres = in_lib = wfail = 0
    for b in range(1 << k):
        c = left
        for t in range(k):
            if (b >> t) & 1:
                c ^= delta[t]
        centres += 1
        if c in lib_set:
            in_lib += 1
            continue
        ones = [t for t in range(k) if (b >> t) & 1]
        zeros = [t for t in range(k) if not (b >> t) & 1]
        u, v = min(ones), max(zeros)
        if u > v:                                       # sorted but not in L: the structure is wrong
            wfail += 1
            continue
        w = witness(c, {bnd[u], bnd[v]})
        if dist(w, c) > M or min(dist(w, l) for l in lib_packed) <= M:
            wfail += 1

    pfail = 0
    for j, D in enumerate(lib_packed):
        uv = set()
        for i, other in enumerate(lib_packed):
            if i != j and dist(D, other) == 1:
                diff = D ^ other
                uv.add(next(p for p in range(N) if (diff >> (2 * p)) & 3))
        w = witness(D, uv)
        if dist(w, D) > M or any(dist(w, l) <= M for i, l in enumerate(lib_packed) if i != j):
            pfail += 1
    return centres, in_lib, centres - in_lib, wfail, pfail, len(lib)


def main(argv):
    if len(argv) < 2:
        print("NOT RUN: give the CHANGE-seq table, as for pooled_patterns_producer.py")
        return 2
    try:
        windows = set()
        with open(argv[1], newline="") as fh:
            for row in csv.DictReader(fh):
                windows.add(row["sgRNA"].strip().upper())
    except OSError as exc:
        print(f"NOT RUN: {exc}")
        return 2
    except KeyError:
        print(f"NOT RUN: {argv[1]} has no sgRNA column")
        return 2
    print(f"{argv[1]}  SHA-256 {sha256_file(argv[1])}")
    guides = sorted({w[:LEN] for w in windows})
    problems = []
    if not (len(guides) == len(windows) == EXPECTED["guides"]
            and all(len(w) == LEN + 3 and w[LEN:] == "NGG" for w in windows)
            and all(set(g) <= set(ALPHA) for g in guides)):
        print("FAIL: the table does not give 110 guides of 20 letters over ACGT in 23-letter NGG windows")
        return 1

    t0 = time.time()
    totals = [0, 0, 0, 0, 0, 0]
    for g in guides:
        row = check_guide(g, M)
        for i, x in enumerate(row):
            totals[i] += x
        if row[0] != 2 ** (row[5] - 1) or row[1] != row[5]:
            problems.append(f"{g}: {row[0]} centres and {row[1]} library words, against 2^(r-1) and r = {row[5]}")
    centres, in_lib, others, wfail, pfail, rsum = totals
    print(f"guides {len(guides)}  m {M}  N {N}  hypothesis m <= N-2: {M <= N - 2}")
    print(f"library-letter centres {centres:,} (sum 2^(r-1))  library words among them {in_lib:,} (sum r {rsum:,})")
    print(f"other centres {others:,}  shown unsound by the step (iii) witness: {others - wfail:,}  failures {wfail}")
    print(f"private points {in_lib:,}  failures {pfail}")
    r_of = [1 + sum(g[p] != g[p + 1] for p in range(N)) for g in guides]
    print(f"weaker form m <= N-r+1 fails at r >= {N - M + 2}: {sum(r >= N - M + 2 for r in r_of)} guides"
          f" | r >= 18: {sum(r >= 18 for r in r_of)} | largest r {max(r_of)}")
    for key, got in (("centres", centres), ("library", in_lib), ("others", others)):
        if got != EXPECTED[key]:
            problems.append(f"{key}: {got:,}, expected {EXPECTED[key]:,}")
    weaker = sum(r >= N - M + 2 for r in r_of)
    if weaker != EXPECTED["weaker"]:
        problems.append(f"weaker form excludes {weaker} guides, expected {EXPECTED['weaker']}")
    if wfail or pfail:
        problems.append(f"{wfail} witness failures and {pfail} private-point failures")

    # Control: one coordinate fewer changed must break both constructions somewhere.
    c_w = c_p = 0
    for g in guides:
        row = check_guide(g, M - 1)
        c_w += row[3]
        c_p += row[4]
    print(f"control, {M - 1} coordinates changed: witness failures {c_w:,}  private-point failures {c_p:,}")
    if not c_w or not c_p:
        problems.append("control: one change fewer did not break the constructions")

    print(f"seconds {time.time() - t0:.0f}")
    for p in problems:
        print(f"FAIL: {p}")
    print("VERDICT:", "FAIL" if problems else "PASS")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
