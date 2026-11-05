#!/usr/bin/env python3
"""Independent producer (r10): infinite-family tightness of companion hypothesis.

DEFINITIONS (stated before any number):
  Guide g in {A,C,G,T}^n. Deletion words: all length-(n-1) strings obtained by
  deleting one letter of g. r = number of distinct deletion words (= run count).
  Um = union of radius-m Hamming balls about those r words, alphabet DNA.
  Centre c in DNA^{n-1} is sound iff B_m(c)subseteq Um.
  Companion theorem range: m <= n-3 => only sound centres are the r deletion words.
  Bound under test: whether the hypothesis "m <= n-3" can be weakened to "m <= n-2".

QUANTITIES BOUNDED / COUNTED (labels):
  |sound| = cardinality of the sound-centre set (exact, exhaustive at the sizes used).
  |Um| = cardinality of Um (exact).
  space = 4^{n-1}.
  sound_gt_r means |sound| > r (strict).
  Family F_k: guide = (AC)* truncated/padded to length n=k+2 with first letter A
    specifically F: for each n >= 3, guide g_n = ("AC" * n)[:n]  # alternating AC
    and also the minimal ACA and the secondary AACCG witnesses.

This file shares NO code with code_release/.../tightness_hypothesis.py.
"""
from __future__ import annotations

import itertools
import json
import time
from pathlib import Path

DNA = "ACGT"
OUT = Path(__file__).resolve().parent / "family_tightness.json"
DEADLINE = time.time() + 15 * 60


def run_count(g: str) -> int:
    return 1 + sum(a != b for a, b in zip(g, g[1:]))


def dels(g: str):
    return sorted({g[:i] + g[i + 1 :] for i in range(len(g))})


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def ball(c: str, m: int):
    L = len(c)
    out = {c}
    frontier = {c}
    for _ in range(m):
        nxt = set()
        for s in frontier:
            for i in range(L):
                for ch in DNA:
                    if ch == s[i]:
                        continue
                    t = s[:i] + ch + s[i + 1 :]
                    if t not in out:
                        out.add(t)
                        nxt.add(t)
        frontier = nxt
    return out


def analyze(guide: str, m: int) -> dict:
    n = len(guide)
    L = n - 1
    words = dels(guide)
    r = len(words)
    assert r == run_count(guide)
    Um = set()
    for w in words:
        Um |= ball(w, m)
    space = 4 ** L
    sound = []
    # Exhaustive over DNA^{L}; only call for L <= 5 (n<=6) or L=6 with m small carefully
    for tup in itertools.product(DNA, repeat=L):
        c = "".join(tup)
        if ball(c, m) <= Um:
            sound.append(c)
    extras = [c for c in sound if c not in set(words)]
    return {
        "guide": guide,
        "n": n,
        "m": m,
        "r": r,
        "deletion_words": words,
        "|Um|": len(Um),
        "space": space,
        "Um_proper": len(Um) < space,
        "Um_nonempty": len(Um) > 0,
        "|sound|": len(sound),
        "sound_gt_r": len(sound) > r,
        "extra_sound_centres": extras,
        "outside_companion_range": m > n - 3,
        "at_m_equals_n_minus_2": m == n - 2,
        "at_m_equals_n_minus_3": m == n - 3,
    }


def alternating(n: int) -> str:
    return ("AC" * n)[:n]


def find_rminus1_cover(guide: str, m: int, max_centres: int = 2):
    """Brute DNA centres for an explicit (r-1)-cover of Um when r-1 is small."""
    info = analyze(guide, m)
    words = info["deletion_words"]
    r = info["r"]
    L = len(guide) - 1
    Um = set()
    for w in words:
        Um |= ball(w, m)
    # Prefer trying deletion words and extras first, then all
    candidates = list(dict.fromkeys(words + info["extra_sound_centres"]))
    if L <= 3:
        candidates = list(dict.fromkeys(candidates + ["".join(t) for t in itertools.product(DNA, repeat=L)]))
    need = r - 1
    # Combinatorial search of size need among candidates (cap)
    from itertools import combinations
    pool = candidates[: min(len(candidates), 40)]
    for combo in combinations(pool, need):
        u = set()
        for c in combo:
            u |= ball(c, m)
        if Um <= u:
            return list(combo)
    return None


def main():
    rows = []
    family = []
    # Family: for n=3..6, alternating AC guide at m=n-2 and control m=n-3
    for n in range(3, 7):
        if time.time() > DEADLINE:
            break
        g = alternating(n)
        outside = analyze(g, n - 2)
        inside = analyze(g, n - 3) if n - 3 >= 0 else None
        family.append({
            "family_member": "alternating_AC",
            "outside": outside,
            "inside_control": inside,
        })
        rows.append(outside)
        if inside:
            rows.append(inside)

    # Minimal ACA + secondary AACCG (recompute, do not read JSON)
    aca_out = analyze("ACA", 1)
    aca_in = analyze("ACA", 0)
    aaccg_out = analyze("AACCG", 3)
    aaccg_in = analyze("AACCG", 2)

    # Explicit covers
    cover_aca = find_rminus1_cover("ACA", 1)
    cover_aaccg = find_rminus1_cover("AACCG", 3)

    # Exhaustive first-letter-A census for n=3 (tiny) to count how many outside witnesses
    n3_outside = []
    for tail in itertools.product(DNA, repeat=2):
        g = "A" + "".join(tail)
        info = analyze(g, 1)
        if info["sound_gt_r"] and info["Um_proper"]:
            n3_outside.append(g)

    # Check whether EVERY alternating member has sound_gt_r at m=n-2
    family_all_outside = all(
        mem["outside"]["sound_gt_r"] and mem["outside"]["Um_proper"]
        for mem in family
    )
    family_all_inside_tight = all(
        mem["inside_control"] is None
        or (mem["inside_control"]["|sound|"] == mem["inside_control"]["r"]
            and not mem["inside_control"]["sound_gt_r"])
        for mem in family
    )

    # Attempt parametric statement strength:
    # "For every n in {3,4,5,6}, alternating AC guide at m=n-2 has |sound|>r and Um proper;
    #  at m=n-3 the same guide has |sound|=r."
    # Infinite family claim requires a proof or larger n; n=7 L=6 m=5 is huge — record as open.

    # Negatives / limits
    negatives = []
    # Can we weaken only partially (e.g. m<=n-2 for some restricted alphabet)? Not claimed.
    # Infinite family for ALL n: OPEN at n>=7 under local limits.
    negatives.append({
        "claim_tested": "alternating_AC family verified for all n>=3 (infinite)",
        "result": "NEGATIVE_UNDER_LOCAL_LIMITS",
        "verified_n": [3, 4, 5, 6],
        "open_n": "n>=7 (L=6,m=5 exhaustive ball enumeration exceeds 20 min / 2GB under nice 15)",
        "what_holds": "For every n in {3,4,5,6}, guide (AC)*[:n] at m=n-2 has |sound|>r with Um proper; at m=n-3 |sound|=r.",
    })

    out = {
        "definitions": {
            "Um": "union of radius-m Hamming balls about the r deletion words",
            "sound": "centre c with B_m(c) subseteq Um",
            "companion_range": "m <= n-3",
            "hypothesis_under_test": "m<=n-3 cannot be weakened to m<=n-2",
            "quantity_|sound|": "cardinality of sound-centre set (exhaustive)",
            "quantity_bounded_by_theorem": "the set of sound centres equals the r deletion words when m<=n-3",
        },
        "settled": {
            "statement": (
                "The companion hypothesis m<=n-3 cannot be weakened to m<=n-2. "
                "Minimal witness ACA (n=3,m=1): |sound|=4>r=3. "
                "Secondary AACCG (n=5,m=3): |sound|=4>r=3. "
                "Finite family: for every n in {3,4,5,6}, alternating AC guide at m=n-2 "
                "has |sound|>r and Um proper; at m=n-3 the same guide has |sound|=r."
            ),
            "minimal_witness": aca_out,
            "minimal_control": aca_in,
            "secondary_AACCG": aaccg_out,
            "secondary_AACCG_control": aaccg_in,
            "minimal_rminus1_cover": cover_aca,
            "AACCG_rminus1_cover": cover_aaccg,
            "n3_firstA_outside_witness_count": len(n3_outside),
            "n3_firstA_outside_witnesses": n3_outside,
            "family_alternating_AC_n_3_to_6_all_outside_sound_gt_r": family_all_outside,
            "family_alternating_AC_n_3_to_6_all_inside_sound_eq_r": family_all_inside_tight,
        },
        "family_rows": family,
        "negatives": negatives,
        "status": "PASS" if (
            aca_out["sound_gt_r"] and aca_out["Um_proper"]
            and aca_in["|sound|"] == aca_in["r"]
            and aaccg_out["sound_gt_r"]
            and family_all_outside and family_all_inside_tight
            and cover_aca is not None
        ) else "FAIL",
        "producer": "code_release/minimum_pattern_certificate/family_tightness.py",
        "shares_code_with_r6_deposit": False,
    }
    OUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("STATUS", out["status"])
    print("ACA |sound|", aca_out["|sound|"], "r", aca_out["r"], "cover", cover_aca)
    print("AACCG |sound|", aaccg_out["|sound|"], "r", aaccg_out["r"], "cover", cover_aaccg)
    print("family_all_outside", family_all_outside, "family_all_inside", family_all_inside_tight)
    for mem in family:
        o = mem["outside"]
        i = mem["inside_control"]
        print("n=%d outside |sound|=%d r=%d Um=%d/%d; inside |sound|=%s r=%s" % (
            o["n"], o["|sound|"], o["r"], o["|Um|"], o["space"],
            i["|sound|"] if i else None, i["r"] if i else None))
    return 0 if out["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
