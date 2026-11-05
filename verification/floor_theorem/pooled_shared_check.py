#!/usr/bin/env python3
"""The five shared patterns: the pooled certificate counts each once, never guide by guide.

WHY THIS EXISTS. The 110 guides' single-deletion libraries are not pairwise disjoint: 1,634 patterns
counted guide by guide, 1,629 distinct, 5 shared, each by two guides offset by one base. The pooled
theorem is a bound over the 1,629 distinct patterns. A proof that slipped between the two counts
would be wrong by five, so this script checks the certificate on exactly those five, with the
guides that share them. The whole certificate (every word in the union, no radius-4 ball holding
three) is checked by pooled_floor_certificate.py; this script does not repeat that.

THE SOURCE. The guides come from the CHANGE-seq table read the way pooled_patterns_producer.py reads
it (the sgRNA column, first 20 letters). That table is a third-party deposit and is not
redistributed here; its source and SHA-256 are in pooled_patterns_producer.py. The pattern list and
the certificate words are the released data/pooled_patterns_1629.txt and data/W_3258.txt, and both
must match their released digests.

WHAT IT CHECKS.
  1. 110 guides; 1,634 patterns guide by guide; 1,629 distinct, byte for byte the released list.
  2. 5 shared patterns, each in exactly two guides offset by one base.
  3. Each shared pattern owns exactly 2 certificate words, not 4: each at distance exactly 4 from it
     and at least 5 from every other distinct pattern, including the rest of both sharing guides.
     The ball centred at each shared pattern carries mass exactly 1.
  4. In the guide-by-guide view (1,634 entries) exactly the 10 words of the shared patterns fall
     inside two entries' balls: a count of 1,634 would fail the certificate's privacy there.
  5. 3,258 distinct words, so the mass is exactly 1,629, a Fraction, one unit per distinct pattern.
It then prints, as reasoning rather than as a check, why 1,634 cannot be the bound: the 1,629
patterns cover the union themselves, and a guide-by-guide certificate would need two more words at
distance 4 from each shared pattern, putting mass 2 in the ball centred there.

REMOVAL CONTROLS, RUN EVERY TIME. Each reruns checks 1 to 5 on altered input and must see a named
check fail: one sharing guide dropped (check 2, the count); one sharing guide replaced by a guide
that shares one pattern but is not offset by one base, five shared kept (check 2, the offset);
two extra words planted at one shared pattern (checks 3 and 4); the two words of one shared
pattern dropped (checks 3 and 5); one of its words replaced by a word at distance 3 (check 3, with
check 5 still holding); one replaced by a word at distance 4 from it and from a second pattern
(check 3, the privacy clause, with check 5 still holding). Every pattern here has a neighbour at
distance 1, so no word below distance 4 is private, and the distance-3 control cannot separate the
"exactly 4" clause from the privacy clause.

Usage:  pooled_shared_check.py DEPOSIT.csv [DATA_DIR]   (default: the data folder beside this script)
Python 3.9 or later, standard library only, about 30 s. Prints the SHA-256 of every input it reads.
Exit 0 only when every check and control holds; 1 otherwise; 2 if it could not run.
"""
import csv
import hashlib
import itertools
import os
import sys
from collections import defaultdict
from fractions import Fraction

RADIUS, LEN = 4, 20
N_GUIDES, N_WITH_MULTIPLICITY, N_PATTERNS, N_SHARED = 110, 1634, 1629, 5
PATTERNS_SHA256 = "b749f0869f5d4a97d1e6c6c48f750694dbec54460846ad4730b5cba76136fa02"
WORDS_SHA256 = "7213bb411954d2247e6df1e0f31533f6bf88fa87499048e9e7f14e36fe052a6e"
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
LOW_BITS = int("01" * (LEN - 1), 2)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def patterns(g):
    return {g[:i] + g[i + 1:] for i in range(len(g))}


def offset_by_one(a, b):
    return a[1:] == b[:-1] or b[1:] == a[:-1]


def pack(word):
    return sum(CODE[c] << (2 * i) for i, c in enumerate(word))


def dist(x, y):
    v = x ^ y
    return bin((v | (v >> 1)) & LOW_BITS).count("1")


def flips(word, k, skip=None):
    """Words that change `word` at k positions (never at `skip`), one fixed letter swap each."""
    swap = {"A": "C", "C": "A", "G": "T", "T": "G"}
    spots = [i for i in range(len(word)) if i != skip]
    for pos in itertools.combinations(spots, k):
        yield "".join(swap[c] if i in pos else c for i, c in enumerate(word))


def evaluate_shared(guides):
    """The shared patterns of these guides, and nothing else."""
    owners = defaultdict(list)
    for g in guides:
        for p in patterns(g):
            owners[p].append(g)
    return {p: gs for p, gs in owners.items() if len(gs) > 1}


def evaluate(guides, words, released_bytes, show=False, windows=None):
    """Checks 1 to 5 on these inputs; returns {check name: (ok, detail)} and the shared patterns."""
    owners = defaultdict(list)
    for g in guides:
        for p in patterns(g):
            owners[p].append(g)
    distinct = sorted(owners)
    shared = {p: gs for p, gs in owners.items() if len(gs) > 1}
    total = sum(len(gs) for gs in owners.values())
    rendered = "".join(p + "\n" for p in distinct).encode("ascii")
    packed = [pack(p) for p in distinct]
    index = {p: i for i, p in enumerate(distinct)}
    mult = [len(owners[p]) for p in distinct]

    rows_d, inside_entries = [], []
    for w in words:
        pw = pack(w)
        d = [dist(pw, q) for q in packed]
        rows_d.append(d)
        inside_entries.append(sum(mult[k] for k, x in enumerate(d) if x <= RADIUS))

    # Check 3 reads each shared pattern's ball: the words inside it, their distance to it, and their
    # distance to the nearest OTHER distinct pattern. The ball's mass is half the words inside it.
    shared_ok, rows, shared_words = True, [], set()
    for p in sorted(shared):
        i = index[p]
        mine = [k for k, d in enumerate(rows_d) if d[i] <= RADIUS]
        shared_words.update(mine)
        at = sorted({rows_d[k][i] for k in mine})
        nxt = [min(x for j, x in enumerate(rows_d[k]) if j != i) for k in mine]
        ball_mass = Fraction(1, 2) * len(mine)
        ok = len(mine) == 2 and at == [RADIUS] and all(x > RADIUS for x in nxt) and ball_mass == 1
        shared_ok = shared_ok and ok
        rows.append(f"{p}  {' '.join(shared[p])}  {len(mine):>5}  {at}  {min(nxt) if nxt else '-'}  {ball_mass}")
    in_two = {k for k, n in enumerate(inside_entries) if n == 2}
    mass = Fraction(1, 2) * len(set(words))

    checks = {
        "1 guides, counts, released list": (
            len(guides) == N_GUIDES and total == N_WITH_MULTIPLICITY and len(distinct) == N_PATTERNS
            and rendered == released_bytes
            and (windows is None or (len(windows) == len(guides)
                                     and all(len(w) == LEN + 3 and w[LEN:] == "NGG" for w in windows))),
            f"{len(guides)} guides, {total:,} guide by guide, {len(distinct):,} distinct"),
        "2 five shared, two guides each, offset by one": (
            len(shared) == N_SHARED and all(len(gs) == 2 and offset_by_one(*gs) for gs in shared.values()),
            f"{len(shared)} shared, at most {max(mult)} guides per pattern"),
        "3 each shared: 2 words, at 4, others at 5+, ball mass 1": (
            shared_ok and len(shared) > 0, f"{len(rows)} shared patterns read"),
        "4 inside two entries in the 1,634 view: exactly the shared words": (
            len(in_two) == 2 * N_SHARED and in_two == shared_words,
            f"{len(in_two)} words, expected the {2 * N_SHARED} of the shared patterns"),
        "5 distinct words, mass exactly 1,629": (
            len(set(words)) == len(words) and mass == N_PATTERNS, f"{len(set(words))} distinct, mass {mass}"),
    }
    if show:
        print("shared pattern        guides                                     words  own  next  ball mass")
        for r in rows:
            print(r)
    return checks, shared


def main(argv):
    if len(argv) < 2:
        print("NOT RUN: give the CHANGE-seq table, as for pooled_patterns_producer.py")
        return 2
    here = os.path.dirname(os.path.abspath(__file__))
    data = argv[2] if len(argv) > 2 else os.path.join(here, "data")
    released_path = os.path.join(data, "pooled_patterns_1629.txt")
    words_path = os.path.join(data, "W_3258.txt")
    try:
        windows = set()
        with open(argv[1], newline="") as fh:
            for row in csv.DictReader(fh):
                windows.add(row["sgRNA"].strip().upper())
        released_bytes = open(released_path, "rb").read()
        words = open(words_path, encoding="ascii").read().split()
    except (OSError, UnicodeDecodeError) as exc:
        print(f"NOT RUN: {exc}")
        return 2
    except KeyError:
        print(f"NOT RUN: {argv[1]} has no sgRNA column")
        return 2
    digests = {}
    for path in (argv[1], released_path, words_path):
        digests[path] = sha256_file(path)
        print(f"{os.path.basename(path)}  SHA-256 {digests[path]}")

    guides = sorted({w[:LEN] for w in windows})
    problems = []
    if digests[released_path] != PATTERNS_SHA256:
        problems.append("the pattern list is not the released one")
    if digests[words_path] != WORDS_SHA256:
        problems.append("the word list is not the released one")
    malformed = []
    if not all(len(g) == LEN and set(g) <= set(CODE) for g in guides):
        malformed.append(f"a guide is not {LEN} letters over ACGT")
    if not all(len(w) == LEN - 1 and set(w) <= set(CODE) for w in words):
        malformed.append(f"a word is not {LEN - 1} letters over ACGT")
    if malformed:
        for p in problems + malformed:
            print(f"FAIL: {p}")
        return 1

    print()
    checks, shared = evaluate(guides, words, released_bytes, show=True, windows=windows)

    print()
    print("why 1,634 cannot be the bound (reasoning from the checks above, not a further check):")
    print(f"  the {N_PATTERNS:,} patterns cover the union of their own balls, so no bound exceeds {N_PATTERNS:,};")
    print("  a guide-by-guide certificate needs two more words at distance 4 from each shared pattern,")
    print("  so the ball centred there, which check 3 finds carrying mass 1, would carry mass 2.")

    # Removal controls: alter one input, rerun checks 1 to 5, and require the named checks to fail.
    p0 = sorted(shared)[0] if shared else None
    controls = []
    if p0 is not None:
        drop_guide = [g for g in guides if g != shared[p0][0]]
        c, _ = evaluate(drop_guide, words, released_bytes)
        controls.append(("one sharing guide dropped: check 2 fails", not c["2 five shared, two guides each, offset by one"][0]))
        flip = {"A": "C", "C": "A", "G": "T", "T": "G"}
        others = [pack(q) for q in set().union(*(patterns(g) for g in guides)) if q != p0]
        extra = [w for w in ("".join(flip[ch] if s <= i < s + RADIUS else ch for i, ch in enumerate(p0))
                             for s in range(LEN - RADIUS))
                 if w not in words and all(dist(pack(w), q) > RADIUS for q in others)][:2]
        c, _ = evaluate(guides, words + extra, released_bytes)
        controls.append(("two words planted at a shared pattern: checks 3 and 4 fail",
                         len(extra) == 2 and not c["3 each shared: 2 words, at 4, others at 5+, ball mass 1"][0]
                         and not c["4 inside two entries in the 1,634 view: exactly the shared words"][0]))
        mine = {w for w in words if dist(pack(w), pack(p0)) <= RADIUS}
        c, _ = evaluate(guides, [w for w in words if w not in mine], released_bytes)
        controls.append(("a shared pattern's words dropped: checks 3 and 5 fail",
                         not c["3 each shared: 2 words, at 4, others at 5+, ball mass 1"][0]
                         and not c["5 distinct words, mass exactly 1,629"][0]))

        # Each clause on its own. A sharing pair that is NOT offset by one base, with the shared
        # count kept at five: only the offset clause of check 2 can fail.
        h = shared[p0][1]
        swapped = None
        for k, letter in itertools.product(range(LEN), "ACGT"):
            g2 = h[:k] + letter + h[k + 1:]
            if g2 == h or g2 in guides or offset_by_one(g2, h):
                continue
            trial = [g for g in guides if g != shared[p0][0]] + [g2]
            if len(trial) == N_GUIDES and len(evaluate_shared(trial)) == N_SHARED:
                swapped = trial
                break
        c, s2 = evaluate(swapped or guides, words, released_bytes)
        controls.append(("a sharing pair not offset by one, five shared kept: check 2 fails",
                         swapped is not None and len(s2) == N_SHARED
                         and not c["2 five shared, two guides each, offset by one"][0]))
        own0 = sorted(mine)
        # A word at distance 3 from the shared pattern replaces one of its two: count and mass stay.
        # Every pattern here has a neighbour at distance 1, so no word below distance 4 is private,
        # and this control exercises the "exactly 4" and "others at 5+" clauses together.
        near3 = next((w for w in flips(p0, RADIUS - 1) if w not in words), None)
        c, _ = evaluate(guides, [w for w in words if w != own0[0]] + [near3] if near3 else words,
                        released_bytes)
        controls.append(("a word at distance 3 replaces one: check 3 fails, check 5 holds",
                         near3 is not None and not c["3 each shared: 2 words, at 4, others at 5+, ball mass 1"][0]
                         and c["5 distinct words, mass exactly 1,629"][0]))
        # A word at distance 4 from the shared pattern AND from its nearest other pattern replaces
        # one of its two: only the "others at 5+" clause of check 3 can fail.
        q = min((q for q in set().union(*(patterns(g) for g in guides)) if q != p0),
                key=lambda q: dist(pack(q), pack(p0)))
        j = next(i for i in range(len(p0)) if p0[i] != q[i])
        third = next(L for L in "ACGT" if L not in (p0[j], q[j]))
        base = p0[:j] + third + p0[j + 1:]
        tie = next((w for w in flips(base, RADIUS - 1, skip=j) if w not in words
                    and dist(pack(w), pack(p0)) == RADIUS and dist(pack(w), pack(q)) == RADIUS), None)
        c, _ = evaluate(guides, [w for w in words if w != own0[0]] + [tie] if tie else words, released_bytes)
        controls.append(("a word at 4 from a second pattern replaces one: check 3 fails",
                         tie is not None and not c["3 each shared: 2 words, at 4, others at 5+, ball mass 1"][0]
                         and c["5 distinct words, mass exactly 1,629"][0]))
    else:
        controls.append(("controls need a shared pattern; none found", False))

    print()
    width = max(len(n) for n in list(checks) + [n for n, _ in controls])
    for name, (ok, detail) in checks.items():
        print(f"  {name:<{width}}  {'ok ' if ok else 'NO '}  {detail}")
    for name, ok in controls:
        print(f"  {name:<{width}}  {'ok ' if ok else 'NO '}")
    for p in problems:
        print(f"FAIL: {p}")
    good = not problems and all(ok for ok, _ in checks.values()) and all(ok for _, ok in controls)
    print()
    print("VERDICT:", "PASS" if good else "FAIL")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
