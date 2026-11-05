#!/usr/local/bin/python3 -I
"""Check a single-deletion word family. Exit: exact=0, nonexact=1, malformed=2."""
import argparse
import csv
import json
from pathlib import Path
import runpy
import sys
from collections import Counter

REF = runpy.run_path(str(Path(__file__).with_name("reference.py")))



def _no_duplicate_keys(pairs):
    """json.loads keeps the LAST value for a repeated key, silently.

    The README promises that duplicate guide IDs are rejected, and without this hook they are not:
    an emitted file carrying the same guide twice has its first list discarded, so an invalid row
    can be hidden behind a valid one and the suite reports "exact". Every JSON read here goes
    through this, so the contract the README states is the contract the code enforces.
    """
    seen = {}
    for key, value in pairs:
        if key in seen:
            raise ValueError("duplicate key in JSON input: " + str(key))
        seen[key] = value
    return seen


def _load_json(text):
    return json.loads(text, object_pairs_hook=_no_duplicate_keys)

def classify(guide, emitted):
    """Keep coverage, invalidity and repeated rows separate, even when they coexist."""
    if not isinstance(emitted, (list, tuple, set, frozenset)):
        raise ValueError("emitted patterns must be a collection of strings")
    if any(not isinstance(word, str) for word in emitted):
        raise ValueError("every pattern must be a string")
    r, expected = REF["required"](guide)
    counts = Counter(emitted)
    observed = set(counts)
    missing = sorted(expected - observed)
    invalid = sorted(observed - expected)
    repeats = sum(count - 1 for count in counts.values())
    status = ("invalid" if invalid else "incomplete" if missing else
              "redundant" if repeats else "exact")
    return {"guide": guide, "r": r, "emitted_rows": sum(counts.values()),
            "distinct_emitted": len(observed), "valid_distinct": len(observed & expected),
            "missing_words": missing, "invalid_words": invalid,
            "redundant_rows": repeats,
            "multiplicities": {w: c for w, c in sorted(counts.items()) if c > 1},
            "word_set_exact": not missing and not invalid, "status": status}


def panel(guides, emitted, budget=4):
    if type(budget) is not int or budget < 0:
        raise ValueError("budget must be a nonnegative integer")
    if not guides:
        raise ValueError("guide panel cannot be empty")
    if not isinstance(emitted, dict) or set(emitted) - set(guides):
        raise ValueError("emitted guide IDs must belong to the guide panel")
    rows = {gid: classify(g, emitted.get(gid, [])) for gid, g in guides.items()}
    statuses = Counter(row["status"] for row in rows.values())
    status = next((s for s in ("invalid", "incomplete", "redundant", "exact") if statuses[s]), None)
    return {"status": status, "guides": len(rows), "budget": budget,
            "required_words": sum(row["r"] for row in rows.values()),
            "emitted_rows": sum(row["emitted_rows"] for row in rows.values()),
            "valid_distinct_words": sum(row["valid_distinct"] for row in rows.values()),
            "missing_words": sum(len(row["missing_words"]) for row in rows.values()),
            "invalid_words": sum(len(row["invalid_words"]) for row in rows.values()),
            "redundant_rows": sum(row["redundant_rows"] for row in rows.values()),
            "guides_incomplete": sum(bool(row["missing_words"]) for row in rows.values()),
            "guide_status_counts": dict(sorted(statuses.items())),
            "sound_uniqueness_applies_to_all": all(budget <= len(g) - 3 for g in guides.values()),
            "per_guide": rows}


def read_guides(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["guide_id", "guide"]:
            raise ValueError("guide TSV header must be guide_id<TAB>guide")
        guides = {}
        for row in reader:
            if None in row or None in row.values() or not row["guide_id"]:
                raise ValueError("malformed guide row")
            if row["guide_id"] in guides:
                raise ValueError("duplicate guide ID")
            guides[row["guide_id"]] = REF["validate"](row["guide"])
        return guides


def read_emitted(path, single=False):
    path = Path(path)
    if path.suffix == ".json":
        data = _load_json(path.read_text())
        if single:
            if not isinstance(data, list):
                raise ValueError("single-guide JSON must be a list of strings")
            return {"guide": data}
        if not isinstance(data, dict) or any(not isinstance(v, list) for v in data.values()):
            raise ValueError("panel JSON must map guide IDs to lists of strings")
        return data
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["guide_id", "pattern"]:
            raise ValueError("emitted TSV header must be guide_id<TAB>pattern")
        result = {}
        for row in reader:
            if None in row or None in row.values() or not row["guide_id"]:
                raise ValueError("malformed emitted row")
            result.setdefault(row["guide_id"], []).append(row["pattern"])
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--guide", help="single uppercase ACGT spacer, no PAM")
    group.add_argument("--guides", type=Path, help="panel TSV: guide_id, guide")
    parser.add_argument("--emitted", type=Path, required=True, help="JSON or TSV word rows")
    parser.add_argument("--budget", type=int, default=4)
    args = parser.parse_args()
    try:
        guides = {"guide": args.guide} if args.guide is not None else read_guides(args.guides)
        result = panel(guides, read_emitted(args.emitted, args.guide is not None), args.budget)
    except (OSError, ValueError, TypeError) as exc:
        print(f"input error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "exact" else 1


if __name__ == "__main__":
    sys.exit(main())
