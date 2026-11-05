#!/usr/bin/env python3
"""Producer: settleable |sound| > r witnesses outside the companion range m <= n-3.

Definitions (before numbers):
  Um = union of radius-m Hamming balls about the r deletion words of a guide.
  A DNA centre c is sound when B_m(c) subseteq Um.
  Companion theorem: for m <= n-3 the only sound centres are the r deletion words.
  This producer searches m > n-3 only, requires Um proper and nonempty, and
  records an explicit DNA (r-1)-cover when one exists.

Primary deposited witness (independently rechecked):
  guide AACCG, n=5, m=3 > n-3=2, r=3, |sound|=4, extra ACCC, Um=220/256,
  explicit 2-cover {AAAC, ACCG}.

    python3 -I sound_gt_r_witness.py
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LETTERS = "ACGT"

# Canonical witness the paper/certificate may cite
WITNESS = {
    "guide": "AACCG",
    "n": 5,
    "m": 3,
    "n_minus_3": 2,
    "r": 3,
    "deletion_words": ["AACC", "AACG", "ACCG"],
    "extra_sound_centres": ["ACCC"],
    "Um_size": 220,
    "space": 256,
    "explicit_rminus1_cover": ["AAAC", "ACCG"],
}


def runs(g: str) -> int:
    return 1 + sum(g[i] != g[i + 1] for i in range(len(g) - 1))


def deletion_words(g: str):
    return sorted({g[:i] + g[i + 1 :] for i in range(len(g))})


def ball(c: str, m: int) -> set[str]:
    out = {c}
    frontier = {c}
    for _ in range(m):
        nxt = set()
        for s in frontier:
            ch = list(s)
            for i in range(len(s)):
                o = ch[i]
                for let in LETTERS:
                    if let == o:
                        continue
                    ch[i] = let
                    t = "".join(ch)
                    ch[i] = o
                    if t not in out:
                        nxt.add(t)
                        out.add(t)
        frontier = nxt
    return out


def analyze(guide: str, m: int) -> dict:
    n = len(guide)
    L = n - 1
    words = deletion_words(guide)
    r = len(words)
    assert r == runs(guide)
    Um = set()
    for w in words:
        Um |= ball(w, m)
    space = 4**L
    sound = []
    for tup in itertools.product(LETTERS, repeat=L):
        c = "".join(tup)
        if ball(c, m) <= Um:
            sound.append(c)
    extras = [c for c in sound if c not in words]
    # search one explicit (r-1)-cover among DNA centres
    cover = None
    centres = ["".join(t) for t in itertools.product(LETTERS, repeat=L)]
    hitting = [c for c in centres if ball(c, m) & Um]
    for combo in itertools.combinations(hitting, r - 1):
        u = set()
        for c in combo:
            u |= ball(c, m)
            if Um <= u:
                cover = list(combo)
                break
        if cover:
            break
    return {
        "guide": guide,
        "n": n,
        "m": m,
        "n_minus_3": n - 3,
        "outside_companion_soundness_range": m > n - 3,
        "r": r,
        "deletion_words": words,
        "|Um|": len(Um),
        "space": space,
        "Um_proper": len(Um) < space,
        "Um_nonempty": len(Um) > 0,
        "|sound|": len(sound),
        "sound_gt_r": len(sound) > r,
        "extra_sound_centres": extras,
        "explicit_rminus1_cover": cover,
        "has_rminus1_DNA_cover": cover is not None,
    }


def main():
    # Verify canonical witness and a small panel of outside-range instances
    rows = []
    for g, m in [
        ("AACCG", 3),
        ("AACCGG", 4),
        ("AACCG", 2),  # inside theorem: expect |sound|==r
    ]:
        rows.append(analyze(g, m))
    primary = analyze(WITNESS["guide"], WITNESS["m"])
    ok = (
        primary["outside_companion_soundness_range"]
        and primary["Um_proper"]
        and primary["sound_gt_r"]
        and primary["|sound|"] == 4
        and primary["extra_sound_centres"] == ["ACCC"]
        and primary["has_rminus1_DNA_cover"]
        and set(primary["explicit_rminus1_cover"]) == set(WITNESS["explicit_rminus1_cover"])
    )
    inside = analyze("AACCG", 2)
    ok_inside = (not inside["outside_companion_soundness_range"]) and (
        inside["|sound|"] == inside["r"]
    )
    out = {
        "definition": (
            "sound iff B_m(c)subseteq Um; companion unique-sound theorem applies "
            "only for m<=n-3; this deposit settles |sound|>r outside that range"
        ),
        "primary_witness": primary,
        "control_inside_theorem_AACCG_m2": inside,
        "panel": rows,
        "status": "PASS" if ok and ok_inside else "FAIL",
        "settled_questions": {
            "exists_settleable_sound_gt_r_Um_proper_outside_theorem": True,
            "primary_also_has_explicit_DNA_rminus1_cover": True,
            "note": (
                "No settleable instance was found with |sound|>r AND no r-1 DNA "
                "cover; every outside-theorem |sound|>r witness checked so far "
                "admits an explicit (r-1)-cover. Certificate 'no r-1 cover' teeth "
                "therefore remain on |sound|=r instances inside m<=n-3."
            ),
        },
    }
    path = HERE / "sound_gt_r_witness.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in ("status", "settled_questions")}, indent=2))
    print("primary |sound|", primary["|sound|"], "r", primary["r"], "cover", primary["explicit_rminus1_cover"])
    print("wrote", path)
    return 0 if out["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
