#!/usr/bin/env python3
"""The pooled floor at 1,629: checks the two-word certificate from the released files alone.

WHAT IT PROVES. The body's pooled theorem says the union of the radius-4 balls around the panel's
1,629 distinct patterns cannot be covered by fewer than 1,629 radius-4 balls, with centres anywhere.
The argument is the weighted covering bound: put weight 1/2 on each of 3,258 words. If every word lies
in the union and no radius-4 ball, centred anywhere, holds three of them, then no ball carries more
than one unit, so any cover needs at least 3,258 / 2 = 1,629 balls. The 1,629 patterns are themselves a
cover, so the covering number is exactly 1,629.

THE SAME WEIGHTS BOUND THE FRACTIONAL OPTIMUM. They are a feasible solution of the dual of the
fractional covering program, of mass 1,629, so the fractional optimum is at least 1,629; it is at most
the covering number, so it is exactly 1,629 as well. One file carries both bounds, and the per-pattern
mass check below is the check that fails if the file is edited.

NO INTEGRAL TRANSVERSAL INSIDE THESE WORDS. The half weights cannot be rounded to one word per
pattern with every two chosen words more than 8 apart. That choice is a 2-SAT problem over the 20,602
pairs of words from different patterns that lie within 8 of each other; it is unsatisfiable, and for
every one of the 1,629 patterns its two choices fall in one strongly connected component of the
implication graph. This is a fact about these words, not about every possible certificate.

WHAT IT READS. data/pooled_patterns_1629.txt (the distinct patterns) and data/W_3258.txt (the
certificate words). Nothing else: no guide list, no guide membership. The five patterns that two
guides share by a one-base offset are single entries of the pattern list, so the bound counts each
once; a count of 1,634 taken guide by guide would be wrong, and the patterns themselves show it is
impossible.

HOW THE BALL TEST WORKS. A ball of radius 4 holding three words puts them pairwise within 8, so only
triples pairwise within 8 need a test. For each one an exact dynamic program over the 19 coordinates
asks whether any centre over ACGT lies within 4 of all three. All arithmetic is on integers and
Fractions.

CONTROLS, RUN EVERY TIME. (1) A planted third word at distance 4 from a pattern that already owns two
words must be caught by the same enumeration and test. (2) At radius 5 the same program must accept
some triple, so it is not a test that always rejects. (3) The 2-SAT solver must return a verified
satisfying assignment for two planted satisfiable instances with the same number of variables, one
random and one made of the real conflicting pairs that a hidden choice of words satisfies, and must
reject x and not x; and the clauses, read back as word choices, must say what the words say.

Usage:  pooled_floor_certificate.py [DATA_DIR]      (default: the data folder beside this script)
Prints the SHA-256 of each list it reads. Python 3.9 or later, standard library only, about 12 s.
Exit 0 when the certificate holds and every control behaves; 1 when a check fails; 2 if it could not
run (a missing or unreadable input).
"""
import hashlib
import os
import random
import sys
from fractions import Fraction

RADIUS = 4
PATTERNS_FILE = "pooled_patterns_1629.txt"
WORDS_FILE = "W_3258.txt"
PATTERNS_SHA256 = "b749f0869f5d4a97d1e6c6c48f750694dbec54460846ad4730b5cba76136fa02"
WORDS_SHA256 = "7213bb411954d2247e6df1e0f31533f6bf88fa87499048e9e7f14e36fe052a6e"
EXPECTED_TRIPLES = 7404  # a property of these words, not of the panel; other valid certificates differ
EXPECTED_CONFLICTS = 20602  # pairs of words from different patterns within 8, for the released words
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def read_list(path):
    with open(path, "rb") as f:
        raw = f.read()
    seqs = [line.strip() for line in raw.decode("ascii").splitlines() if line.strip()]
    return seqs, hashlib.sha256(raw).hexdigest()


def pack(word):
    return sum(CODE[c] << (2 * i) for i, c in enumerate(word))


def hamming(x, y, low_bits):
    v = x ^ y
    return bin((v | (v >> 1)) & low_bits).count("1")


def common_centre(a, b, c, radius):
    """True when some word over ACGT lies within `radius` of all three of a, b, c."""
    states = {(0, 0, 0)}
    for x, y, z in zip(a, b, c):
        steps = {(int(L != x), int(L != y), int(L != z)) for L in "ACGT"}
        states = {(p + s, q + t, r + u) for p, q, r in states for s, t, u in steps
                  if p + s <= radius and q + t <= radius and r + u <= radius}
        if not states:
            return False
    return True


def close_pairs(words, cut, low_bits):
    """For each word, the later words within `cut` of it."""
    packed = [pack(w) for w in words]
    near = [set() for _ in words]
    for i in range(len(words)):
        pi = packed[i]
        for j in range(i + 1, len(words)):
            if hamming(pi, packed[j], low_bits) <= cut:
                near[i].add(j)
    return near


def triples_within(words, cut, low_bits):
    near = close_pairs(words, cut, low_bits)
    for i in range(len(words)):
        for j in sorted(near[i]):
            for k in sorted(near[i] & near[j]):
                yield i, j, k


def ball_test(words, radius, low_bits):
    triples = list(triples_within(words, 2 * radius, low_bits))
    held = [t for t in triples if common_centre(*(words[i] for i in t), radius)]
    return triples, held


def two_sat(n_vars, clauses):
    """Clauses are pairs (a, b) meaning a or b; literal 2v + s, its negation is the literal ^ 1.
    Returns (satisfiable, an assignment of the even literals or None, variables whose two literals
    share a strongly connected component of the implication graph)."""
    n = 2 * n_vars
    out = [[] for _ in range(n)]
    back = [[] for _ in range(n)]
    for a, b in clauses:
        out[a ^ 1].append(b)
        out[b ^ 1].append(a)
        back[b].append(a ^ 1)
        back[a].append(b ^ 1)
    order, seen = [], [False] * n
    for s in range(n):
        if seen[s]:
            continue
        seen[s] = True
        stack = [(s, 0)]
        while stack:
            v, k = stack[-1]
            if k < len(out[v]):
                stack[-1] = (v, k + 1)
                u = out[v][k]
                if not seen[u]:
                    seen[u] = True
                    stack.append((u, 0))
            else:
                stack.pop()
                order.append(v)
    comp, c = [-1] * n, 0
    for s in reversed(order):
        if comp[s] != -1:
            continue
        comp[s] = c
        stack = [s]
        while stack:
            v = stack.pop()
            for u in back[v]:
                if comp[u] == -1:
                    comp[u] = c
                    stack.append(u)
        c += 1
    tied = sum(1 for v in range(n_vars) if comp[2 * v] == comp[2 * v + 1])
    if tied:
        return False, None, tied
    return True, [comp[2 * v] > comp[2 * v + 1] for v in range(n_vars)], 0


def satisfies(assign, clauses):
    def value(lit):
        return assign[lit >> 1] == (lit % 2 == 0)
    return all(value(a) or value(b) for a, b in clauses)


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    data = argv[1] if len(argv) > 1 else os.path.join(here, "data")
    try:
        patterns, p_sha = read_list(os.path.join(data, PATTERNS_FILE))
        words, w_sha = read_list(os.path.join(data, WORDS_FILE))
    except (OSError, UnicodeDecodeError) as exc:
        print(f"NOT RUN: cannot read the inputs in {data}: {exc}")
        print(f"supply {PATTERNS_FILE} and {WORDS_FILE} in that folder")
        return 2
    print(f"patterns {PATTERNS_FILE} sha256 {p_sha}")
    print(f"words    {WORDS_FILE} sha256 {w_sha}")

    failures = []
    if p_sha != PATTERNS_SHA256:
        failures.append("pattern list is not the released one")
    if w_sha != WORDS_SHA256:
        failures.append("word list is not the released one")
    n = len(patterns[0]) if patterns else 0
    alphabet_ok = all(len(s) == n and set(s) <= set(CODE) for s in patterns + words)
    if not alphabet_ok or n == 0:
        for f in failures:
            print(f"FAIL: {f}")
        print("FAIL: every line must be a word of one length over ACGT")
        return 1
    low_bits = int("01" * n, 2)

    distinct_patterns = len(set(patterns))
    distinct_words = len(set(words))
    print(f"patterns {len(patterns)} distinct {distinct_patterns} | words {len(words)} distinct {distinct_words}")
    if distinct_patterns != len(patterns) or distinct_words != len(words):
        failures.append("repeated lines")

    packed_p = [pack(p) for p in patterns]
    owned = [0] * len(patterns)
    owner_of = []
    outside_union = not_exactly_radius = in_two_balls = 0
    for w in words:
        pw = pack(w)
        d = [hamming(pw, q, low_bits) for q in packed_p]
        best = min(d)
        inside = [i for i, x in enumerate(d) if x <= RADIUS]
        if best > RADIUS:
            outside_union += 1
            owner_of.append(-1)
            continue
        if best != RADIUS:
            not_exactly_radius += 1
        if len(inside) > 1:
            in_two_balls += 1
        owned[d.index(best)] += 1
        owner_of.append(d.index(best))
    owning_two = sum(1 for c in owned if c == 2)
    print(f"words outside the union {outside_union} | own distance not exactly {RADIUS} {not_exactly_radius}"
          f" | inside two patterns' balls {in_two_balls} | patterns owning exactly two {owning_two} of {len(patterns)}")
    if outside_union:
        failures.append("a word lies outside the union, so its weight proves nothing")
    if not_exactly_radius:
        failures.append(f"{not_exactly_radius} words are not at distance exactly {RADIUS} from their pattern")
    if in_two_balls:
        failures.append(f"{in_two_balls} words lie inside two patterns' balls, so they are not private")
    if owning_two != len(patterns):
        failures.append(f"only {owning_two} of {len(patterns)} patterns own exactly two words")

    triples, held = ball_test(words, RADIUS, low_bits)
    mass = Fraction(1, 2) * len(words)
    print(f"triples pairwise within {2 * RADIUS}: {len(triples)} | held by one radius-{RADIUS} ball: {len(held)}")
    if held:
        failures.append(f"{len(held)} triples share a radius-{RADIUS} ball, so a ball can carry 3/2")
    if mass != len(patterns):
        failures.append(f"mass {mass} is not the pattern count {len(patterns)}")
    if w_sha == WORDS_SHA256 and len(triples) != EXPECTED_TRIPLES:
        failures.append(f"{len(triples)} triples within {2 * RADIUS}, expected {EXPECTED_TRIPLES} for the released words")

    # Control 1: a third word at distance RADIUS from a pattern that owns two must be caught.
    if 2 not in owned:
        failures.append("control 1: no pattern owns two words, so the planted case cannot be built")
        for f in failures:
            print(f"FAIL: {f}")
        return 1
    target = owned.index(2)
    p = patterns[target]
    flip = {"A": "C", "C": "A", "G": "T", "T": "G"}
    planted = None
    for start in range(n - RADIUS + 1):
        cand = "".join(flip[c] if start <= i < start + RADIUS else c for i, c in enumerate(p))
        if cand not in words:
            planted = cand
            break
    _, held_planted = ball_test(words + [planted], RADIUS, low_bits)
    caught = any(len(words) in t for t in held_planted)
    print(f"control 1, planted third word near pattern {target}: caught {caught}")
    if not caught:
        failures.append("control 1: the planted triple was not caught")

    # Control 2: at radius + 1 the program must accept some triple.
    accepted = sum(1 for t in triples if common_centre(*(words[i] for i in t), RADIUS + 1))
    print(f"control 2, same triples at radius {RADIUS + 1}: accepted {accepted} of {len(triples)}")
    if not accepted:
        failures.append("control 2: the ball test rejected everything at radius + 1")

    # No integral transversal inside these words: choose one word per pattern, every two chosen
    # words more than 2 * RADIUS apart. Variable v is pattern v; literal 2v + s picks its s-th word.
    if not failures:
        slot, seen_pattern = [], {}
        for o in owner_of:
            slot.append(seen_pattern.get(o, 0))
            seen_pattern[o] = seen_pattern.get(o, 0) + 1
        near = close_pairs(words, 2 * RADIUS, low_bits)
        pairs = [(i, j) for i in range(len(words)) for j in near[i] if owner_of[i] != owner_of[j]]
        clauses = [((2 * owner_of[i] + slot[i]) ^ 1, (2 * owner_of[j] + slot[j]) ^ 1) for i, j in pairs]
        sat, _, tied = two_sat(len(patterns), clauses)
        print(f"one word per pattern, pairwise more than {2 * RADIUS} apart: {len(clauses)} conflicting pairs,"
              f" satisfiable {sat}, patterns whose two choices share a component {tied}")
        if sat:
            failures.append("an integral transversal exists inside the words, so the README is wrong")
        if w_sha == WORDS_SHA256 and (len(clauses) != EXPECTED_CONFLICTS or tied != len(patterns)):
            failures.append(f"{len(clauses)} conflicting pairs and {tied} tied patterns, expected"
                            f" {EXPECTED_CONFLICTS} and {len(patterns)} for the released words")
        # Controls for the solver: planted satisfiable instances, and x and not x.
        rng = random.Random(1629)
        hidden = [rng.random() < 0.5 for _ in range(len(patterns))]
        planted_clauses = []
        for _ in range(len(clauses)):
            a, b = rng.randrange(2 * len(patterns)), rng.randrange(2 * len(patterns))
            if not satisfies(hidden, [(a, b)]):
                a ^= 1
            planted_clauses.append((a, b))
        p_sat, p_assign, _ = two_sat(len(patterns), planted_clauses)
        p_ok = p_sat and satisfies(p_assign, planted_clauses)
        u_sat, _, _ = two_sat(1, [(0, 0), (1, 1)])
        # The real conflict pairs, keeping only those a hidden choice of words satisfies: same
        # structure, satisfiable by construction.
        kept_pairs = [pr for pr, c in zip(pairs, clauses) if satisfies(hidden, [c])]
        kept = [c for pr, c in zip(pairs, clauses) if satisfies(hidden, [c])]
        k_sat, k_assign, _ = two_sat(len(patterns), kept)
        k_ok = k_sat and satisfies(k_assign, kept)

        # The encoding, read back as words rather than literals. A choice of words: pattern v takes
        # its slot-0 word when its even literal is true.
        def chosen(assign):
            return [(slot[i] == 0) == assign[owner_of[i]] for i in range(len(words))]
        hid = chosen(hidden)
        both_hidden = {pr for pr in pairs if hid[pr[0]] and hid[pr[1]]}
        encoding_ok = (both_hidden == set(pairs) - set(kept_pairs)
                       and all(sum(1 for i in range(len(words)) if hid[i] and owner_of[i] == v) == 1
                               for v in set(owner_of)))
        if k_ok:
            got = chosen(k_assign)
            encoding_ok = encoding_ok and not any(got[i] and got[j] for i, j in kept_pairs)
        print(f"control 3, the clause encoding read back as word choices: {encoding_ok}")
        if not encoding_ok:
            failures.append("control 3: the clauses do not say what the word choices say")
        print(f"control 3, planted satisfiable instances: random {p_ok}, real pairs a hidden choice satisfies"
              f" ({len(kept)} of {len(clauses)}) {k_ok}, each SAT with a verified assignment"
              f" | x and not x: UNSAT {not u_sat}")
        if not p_ok or not k_ok or u_sat:
            failures.append("control 3: the 2-SAT solver misread a planted instance")

    if failures:
        for f in failures:
            print(f"FAIL: {f}")
        return 1
    print(f"mass {mass} | every radius-{RADIUS} ball carries at most 1")
    print(f"=> no cover by fewer than {mass} balls; the {len(patterns)} patterns are a cover,"
          f" so the covering number is exactly {len(patterns)}")
    print(f"=> the weights are a feasible dual of mass {mass}, so the fractional optimum is at least"
          f" {mass}; it is at most the covering number, so it is exactly {len(patterns)} too")
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
