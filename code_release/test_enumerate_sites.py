#!/usr/bin/env python3
"""Tests for enumerate_sites.py. Run: python3 test_enumerate_sites.py (pytest also works).

The vectorised scan is compared with a naive window-by-window scan written independently
here, on random sequences with planted windows (both strands, NGG and NAG PAMs, 0 to 7
mismatches, bad PAMs, windows containing N, soft-masked bases)."""
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import enumerate_sites as E  # noqa: E402

SPACER = E.SPACER_EXACEL
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(COMP[c] for c in reversed(s))


def naive(seq, spacer, max_mm=6):
    """Window-by-window reference implementation, both strands, upper-cased. Coordinates are
    1-based on the forward strand: for '-' the window read 5' to 3' on the minus strand is the
    reverse complement of the forward bases start .. start + 22."""
    seq = seq.upper()
    out = []
    for strand, s in (("+", seq), ("-", rc(seq))):
        for i in range(len(s) - 22):
            w = s[i:i + 23]
            if set(w) - set("ACGT"):
                continue
            if w[21:] != "GG" and not (w[21] == "A" and w[22] == "G"):
                continue
            mm = sum(a != b for a, b in zip(w[:20], spacer))
            if mm <= max_mm:
                start = i + 1 if strand == "+" else len(seq) - i - 23 + 1
                out.append((strand, start, w, mm))
    return out


def mutate(spacer, n, rng):
    pos = rng.sample(range(20), n)
    s = list(spacer)
    for p in pos:
        s[p] = rng.choice([b for b in "ACGT" if b != s[p]])
    return "".join(s)


def random_genome(rng, length, plants):
    g = [rng.choice("ACGT") for _ in range(length)]
    for w in plants:
        i = rng.randrange(0, length - len(w))
        g[i:i + len(w)] = list(w)
    return "".join(g)


def test_matches_naive_scan():
    rng = random.Random(20261013)
    checked = 0
    for trial in range(6):
        plants = [SPACER + "AGG", mutate(SPACER, 2, rng) + "TGG", mutate(SPACER, 6, rng) + "CAG",
                  mutate(SPACER, 7, rng) + "AGG", SPACER + "ATG", rc(mutate(SPACER, 3, rng) + "GGG"),
                  rc(mutate(SPACER, 5, rng) + "AAG"), SPACER[:10] + "N" + SPACER[11:] + "TGG"]
        g = random_genome(rng, 4000, plants)
        g = g[:500] + g[500:700].lower() + g[700:]          # soft-masked stretch
        fast = E.scan_encoded(E.encode(g.encode()), E.encode(SPACER.encode()))
        slow = naive(g, SPACER)
        assert sorted(fast) == sorted(slow), (trial, len(fast), len(slow))
        checked += len(slow)
    return checked


def test_coordinates_roundtrip():
    """A planted window is found at the planted forward coordinates, on both strands."""
    rng = random.Random(7)
    plus = SPACER + "AGG"
    minus_fwd = rc(mutate(SPACER, 1, rng) + "TGG")      # forward-strand bases of a minus-strand site
    g = "".join(rng.choice("ACGT") for _ in range(300))
    g = g[:100] + plus + g[123:200] + minus_fwd + g[223:]
    recs = E.scan_encoded(E.encode(g.encode()), E.encode(SPACER.encode()), 1)
    got = {(s, st): (w, mm) for s, st, w, mm in recs}
    assert ("+", 101) in got and got[("+", 101)][0] == plus and got[("+", 101)][1] == 0
    assert ("-", 201) in got and got[("-", 201)][1] == 1 and g[200:223] == minus_fwd
    assert rc(got[("-", 201)][0]) == minus_fwd
    return 5


def test_summary_excludes_on_target():
    g = "T" * 50 + SPACER + "AGG" + "T" * 50 + mutate(SPACER, 2, random.Random(1)) + "TGG" + "T" * 50
    recs = E.scan_encoded(E.encode(g.encode()), E.encode(SPACER.encode()))
    s = E.summarise(recs, SPACER)
    assert s["on_target_windows"] == 1
    assert s["windows_total"] == s["off_target_windows"] + 1
    assert s["mismatch_histogram_windows"]["0"] == 1 and s["mismatch_histogram_windows"]["2"] == 1
    return 4


TESTS = [test_matches_naive_scan, test_coordinates_roundtrip, test_summary_excludes_on_target]

if __name__ == "__main__":
    failed = 0
    for t in TESTS:
        try:
            n = t()
            print(f"PASS  {t.__name__}  ({n} checks)")
        except Exception as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(TESTS) - failed} of {len(TESTS)} tests passed")
    sys.exit(1 if failed else 0)
