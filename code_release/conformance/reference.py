#!/usr/local/bin/python3 -I
"""Single-deletion library, standard library only. PAMs are not part of a guide."""
import argparse
import json


def validate(guide):
    if not isinstance(guide, str) or not guide or set(guide) - set("ACGT"):
        raise ValueError("guide must be a nonempty uppercase ACGT string without PAM")
    return guide


def run_count(guide):
    """O(n) time, O(1) auxiliary space; no deletion words are constructed."""
    validate(guide)
    return 1 + sum(guide[i] != guide[i - 1] for i in range(1, len(guide)))


def required_words(guide):
    """Delete at each run start: one word per run; O(n*r) materialized output."""
    validate(guide)
    return {guide[:i] + guide[i + 1:] for i in range(len(guide))
            if i == 0 or guide[i] != guide[i - 1]}


def required(guide):
    """Return (required distinct count, required word set)."""
    return run_count(guide), required_words(guide)


def certificate(guide, budget=4):
    if type(budget) is not int or budget < 0:
        raise ValueError("budget must be a nonnegative integer")
    r, words = required(guide)
    return {"guide": guide, "r": r, "words": sorted(words), "budget": budget,
            "sound_uniqueness_applies": budget <= len(guide) - 3,
            "paper_a_unrestricted_floor_applies": len(guide) >= 20 and budget <= 5}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("guide")
    parser.add_argument("--budget", type=int, default=4)
    parser.add_argument("--count-only", action="store_true")
    args = parser.parse_args()
    try:
        if args.budget < 0:
            raise ValueError("budget must be nonnegative")
        result = run_count(args.guide) if args.count_only else certificate(args.guide, args.budget)
    except ValueError as exc:
        parser.error(str(exc))
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
