#!/usr/bin/env python3
"""Check the second-common-word criterion by literal four-letter enumeration."""
import itertools
import json


def deletions(word):
    return {word[:i] + word[i + 1:] for i in range(len(word))}


def check_pair(g, x):
    w = g[:-1]
    h = x + w
    common = deletions(g) & deletions(h)
    assert w in common and g != h
    shorter = deletions(w)
    candidates = {x + a for a in shorter} & {b + g[-1] for b in shorter}
    observed = len(common) == 2
    original = bool(candidates)
    repaired = bool(candidates - {w})
    assert len(common) <= 2
    assert repaired == observed, (g, h, common, candidates)
    assert candidates - {w} == common - {w}
    return observed, original != observed, sorted(common), sorted(candidates)


def main():
    totals = {"pairs": 0, "one_common_word": 0, "two_common_words": 0,
              "original_mismatches": 0, "repaired_mismatches": 0,
              "inserted_base_equals_first": 0}
    by_length = []
    for n in range(2, 7):
        count = bad = positive = 0
        for letters in itertools.product("ACGT", repeat=n):
            g = "".join(letters)
            for x in "ACGT":
                if x + g[:-1] == g:
                    continue
                observed, mismatch, _, _ = check_pair(g, x)
                count += 1
                bad += mismatch
                positive += observed
                totals["inserted_base_equals_first"] += x == g[0]
        by_length.append({"n": n, "pairs": count, "original_mismatches": bad,
                          "two_common_words": positive})
        totals["pairs"] += count
        totals["original_mismatches"] += bad
        totals["two_common_words"] += positive
        totals["one_common_word"] += count - positive
    assert totals["pairs"] == 21804
    assert totals["original_mismatches"] == 1224
    assert totals["one_common_word"] > 0 and totals["two_common_words"] > 0
    assert totals["inserted_base_equals_first"] > 0
    examples = []
    for g, x in (("ACAC", "C"), ("C" + "A" * 19, "C")):
        observed, mismatch, common, candidates = check_pair(g, x)
        examples.append({"g": g, "h": x + g[:-1], "common": common,
                         "candidate_equalities": candidates, "has_second_word": observed,
                         "original_wrong": mismatch})
    assert examples[0]["common"] == ["ACA", "CAC"]
    assert examples[1]["original_wrong"] and not examples[1]["has_second_word"]
    print(json.dumps({"totals": totals, "by_length": by_length, "examples": examples}, indent=2))
    print("COLLISION_CRITERION_VERIFIED: 21804 pairs, 1224 original failures, zero repaired failures")


if __name__ == "__main__":
    main()
