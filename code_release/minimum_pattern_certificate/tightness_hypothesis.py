#!/usr/bin/env python3
"""Producer: companion unique-sound hypothesis tightness census.

Definitions (before numbers):
  Um = union of radius-m Hamming balls about the r deletion words of a guide of length n.
  Centre c in Sigma^{n-1} is sound when B_m(c) subseteq Um.
  Companion theorem: for m <= n-3 the only sound centres are the r deletion words.
  Hypothesis tightness: m <= n-3 cannot be weakened to m <= n-2.

Settled statement (this deposit):
  The hypothesis m <= n-3 cannot be weakened to m <= n-2.  Minimal first-letter-A
  witness under exhaustive enumeration for n <= 5: guide ACA at m = n-2 = 1 has
  |sound| = 4 > r = 3 with Um proper (12/16); the same guide at m = n-3 = 0 has
  |sound| = r = 3.  Previously deposited AACCG (n=5,m=3) remains a valid larger
  witness and is re-checked as a secondary.

Caps: nice 15, exhaustive for n<=5 (first-letter A), structured+random for n in {6,7}.
"""
from __future__ import annotations

import itertools
import json
import random
import time
from pathlib import Path

LETTERS = "ACGT"
HERE = Path(__file__).resolve().parent
OUT = HERE / "tightness_census.json"
DEADLINE = time.time() + 12 * 60


def runs(g: str) -> int:
    return 1 + sum(g[i] != g[i + 1] for i in range(len(g) - 1))


def deletion_words(g: str):
    return sorted({g[:i] + g[i + 1 :] for i in range(len(g))})


def encode(s: str) -> int:
    v = 0
    for ch in s:
        v = (v << 2) | LETTERS.index(ch)
    return v


def ball_indices(c_idx: int, L: int, m: int) -> set[int]:
    out = {c_idx}
    frontier = {c_idx}
    for _ in range(m):
        nxt = set()
        for idx in frontier:
            for pos in range(L):
                shift = 2 * (L - 1 - pos)
                cur = (idx >> shift) & 3
                for let in range(4):
                    if let == cur:
                        continue
                    nidx = idx ^ ((cur ^ let) << shift)
                    if nidx not in out:
                        nxt.add(nidx)
                        out.add(nidx)
        frontier = nxt
    return out


def decode(idx: int, L: int) -> str:
    chars = []
    for pos in range(L):
        shift = 2 * (L - 1 - pos)
        chars.append(LETTERS[(idx >> shift) & 3])
    return "".join(chars)


def analyze_idx(guide: str, m: int) -> dict:
    n = len(guide)
    L = n - 1
    words = deletion_words(guide)
    r = len(words)
    assert r == runs(guide)
    Um = set()
    for w in words:
        Um |= ball_indices(encode(w), L, m)
    space = 1 << (2 * L)
    sound = []
    for c_idx in range(space):
        if ball_indices(c_idx, L, m) <= Um:
            sound.append(c_idx)
    wordset = set(words)
    extras = [decode(i, L) for i in sound if decode(i, L) not in wordset]
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
    }


def guides_first_A(n: int):
    for tail in itertools.product(LETTERS, repeat=n - 1):
        yield "A" + "".join(tail)


def structured_guides(n: int, seed: int):
    rng = random.Random(seed)
    out = []
    for r in range(1, n + 1):
        base, rem = divmod(n, r)
        sizes = [base + (1 if i < rem else 0) for i in range(r)]
        s = "".join(LETTERS[i % 4] * sz for i, sz in enumerate(sizes))
        if s[0] != "A":
            mapping = {s[0]: "A", "A": s[0]}
            s = "".join(mapping.get(c, c) for c in s)
        out.append(s)
    for k in range(1, n):
        for other in "CGT":
            out.append("A" * k + other * (n - k))
    for _ in range(120):
        out.append("A" + "".join(rng.choice(LETTERS) for _ in range(n - 1)))
    # dedupe preserve order
    return list(dict.fromkeys(out))


def main():
    t0 = time.time()
    census = []
    witnesses_outside = []
    boundary_failures_inside = []
    timed_out = False

    for n in range(3, 7):  # n=7 deferred: L=6 m=5 too heavy under 20min local
        if time.time() > DEADLINE:
            timed_out = True
            break
        for m_mode, m in (("outside_n_minus_2", n - 2), ("inside_n_minus_3", n - 3)):
            if m < 0:
                continue
            if time.time() > DEADLINE:
                timed_out = True
                break
            if n <= 5:
                guide_iter = list(guides_first_A(n))
                mode = "exhaustive_first_A"
            else:
                guide_iter = structured_guides(n, 20261101 + n * 17 + m)
                mode = "structured_plus_120_random_first_A"
            n_guides = n_proper = n_gt = 0
            first_witness = None
            for g in guide_iter:
                if time.time() > DEADLINE:
                    timed_out = True
                    break
                n_guides += 1
                row = analyze_idx(g, m)
                if row["Um_proper"] and row["Um_nonempty"]:
                    n_proper += 1
                    if row["sound_gt_r"]:
                        n_gt += 1
                        rec = {
                            "guide": g,
                            "n": n,
                            "m": m,
                            "r": row["r"],
                            "|sound|": row["|sound|"],
                            "extra_sound_centres": row["extra_sound_centres"][:8],
                            "|Um|": row["|Um|"],
                            "space": row["space"],
                        }
                        if first_witness is None:
                            first_witness = rec
                        if m_mode.startswith("outside"):
                            witnesses_outside.append(rec)
                        else:
                            boundary_failures_inside.append(rec)
            cell = {
                "n": n,
                "m": m,
                "m_mode": m_mode,
                "enum_mode": mode,
                "guides_classified": n_guides,
                "Um_proper_nonempty": n_proper,
                "sound_gt_r_count": n_gt,
                "first_witness": first_witness,
                "wall_seconds_partial": round(time.time() - t0, 3),
            }
            census.append(cell)
            print("CELL", cell, flush=True)
            if timed_out:
                break

    primary = analyze_idx("ACA", 1)
    control = analyze_idx("ACA", 0)
    secondary = analyze_idx("AACCG", 3)
    secondary_ctrl = analyze_idx("AACCG", 2)

    uniq = {}
    for w in witnesses_outside:
        uniq[(w["guide"], w["m"])] = w
    witnesses_unique = sorted(uniq.values(), key=lambda w: (w["n"], w["m"], w["guide"]))
    minimal = witnesses_unique[0] if witnesses_unique else None
    smaller_than_aaccg = [w for w in witnesses_unique if w["n"] < 5]

    # Find explicit (r-1)-cover for primary if any (small: brute pairs)
    cover = None
    if primary["sound_gt_r"]:
        L = primary["n"] - 1
        words = primary["deletion_words"]
        Um = set()
        for w in words:
            Um |= ball_indices(encode(w), L, primary["m"])
        space = 1 << (2 * L)
        # try all pairs of DNA words
        found = False
        for i in range(space):
            if found:
                break
            bi = ball_indices(i, L, primary["m"])
            if not (bi & Um):
                continue
            for j in range(i + 1, space):
                if bi | ball_indices(j, L, primary["m"]) >= Um:
                    cover = [decode(i, L), decode(j, L)]
                    found = True
                    break

    result = {
        "definition": (
            "Companion unique-sound theorem: for m<=n-3 the only sound centres are the r "
            "deletion words. Hypothesis tightness: m<=n-3 cannot be weakened to m<=n-2."
        ),
        "settled_statement": (
            "The hypothesis m <= n-3 cannot be weakened to m <= n-2: the minimal "
            "first-letter-A witness under exhaustive enumeration at n<=5 is guide ACA "
            "at m = n-2 = 1, with |sound| = 4 > r = 3 and Um proper (12/16); at the "
            "theorem boundary m = n-3 = 0 the same guide has |sound| = r = 3. "
            "Secondary: AACCG at m=3=n-2 also has |sound|=4>r=3 (prior deposit)."
        ),
        "primary_minimal_witness": primary,
        "primary_explicit_rminus1_cover": cover,
        "control_at_theorem_boundary": control,
        "secondary_AACCG_outside": secondary,
        "secondary_AACCG_boundary": secondary_ctrl,
        "minimal_outside_witness_in_census": minimal,
        "outside_witnesses_with_n_lt_5": smaller_than_aaccg,
        "census": census,
        "outside_witnesses_unique_count": len(witnesses_unique),
        "outside_witnesses_unique_sample": witnesses_unique[:40],
        "inside_boundary_failures_count": len(boundary_failures_inside),
        "inside_boundary_failures": boundary_failures_inside[:20],
        "timed_out": timed_out,
        "wall_seconds": round(time.time() - t0, 3),
        "status": (
            "PASS"
            if (
                primary["sound_gt_r"]
                and primary["Um_proper"]
                and primary["outside_companion_soundness_range"]
                and not control["sound_gt_r"]
                and control["|sound|"] == control["r"]
                and len(boundary_failures_inside) == 0
            )
            else "FAIL"
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(
        "WROTE",
        OUT,
        "status",
        result["status"],
        "outside_unique",
        len(witnesses_unique),
        "inside_failures",
        len(boundary_failures_inside),
        "minimal",
        minimal,
        "cover",
        cover,
        "wall",
        result["wall_seconds"],
        flush=True,
    )
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
