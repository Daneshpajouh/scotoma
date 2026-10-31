#!/usr/bin/env python3
"""Independent exact checker for the three m=8 staircase certificates.

Definitions follow paper_v7/sections/03_method.tex: four-letter Hamming space,
staircase P_j = 1^j 0^{N-j}, U the union of radius-m balls about the first r
patterns, and a feasible packing a nonnegative weight on U whose mass in every
radius-m ball is at most 1. The certificate row format is the guarded shell
(g, j, t, f_left, f_right, v) with normalizer A. Exact integers only.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import combinations
from math import comb, lcm
from pathlib import Path

CERT_DIR = Path(__file__).resolve().parent
M_BUDGET = 8
MAX_WORKERS = 8


def pattern(j: int, N: int) -> list[int]:
    return [1 if i < j else 0 for i in range(N)]


def hamming(a: list[int], b: list[int]) -> int:
    return sum(x != y for x, y in zip(a, b))


def regions(j: int, g: int, N: int) -> tuple[list[int], list[int], list[int]]:
    left = list(range(max(0, j - g)))
    guard = list(range(max(0, j - g), min(N, j + g)))
    right = list(range(min(N, j + g), N))
    return left, guard, right


def shell_size(N: int, g: int, j: int, t: int, fl: int, fr: int, m: int) -> int:
    left, _, right = regions(j, g, N)
    L, R = len(left), len(right)
    s = m - t
    if t < 0 or s < 0 or fl < 0 or fr < 0:
        return 0
    if t > L or s > R or fl > t or fr > s or fl > g - 1 or fr > g - 1:
        return 0
    return comb(L, t) * comb(t, fl) * comb(R, s) * comb(s, fr) * (1 << (m - fl - fr))


def admissible_lcm(N: int, guards: list[int], m: int) -> tuple[int, int]:
    """LCM of shell sizes over every legal (g, j, t, f_left, f_right).

    j runs through the full staircase 0..N. Shells with j >= r are not given
    mass; they only enlarge the common denominator, exactly as a cleared
    denominator must.
    """
    H = 1
    count = 0
    for g in guards:
        for j in range(N + 1):
            left, _, right = regions(j, g, N)
            L, R = len(left), len(right)
            for t in range(m + 1):
                s = m - t
                if t > L or s > R:
                    continue
                for fl in range(g):
                    if fl > t:
                        continue
                    for fr in range(g):
                        if fr > s:
                            continue
                        S = comb(L, t) * comb(t, fl) * comb(R, s) * comb(s, fr) * (
                            1 << (m - fl - fr)
                        )
                        H = lcm(H, S)
                        count += 1
    return H, count


class Component:
    __slots__ = (
        "g", "j", "t", "fl", "fr", "v", "S", "W", "m", "N",
        "left", "guard", "right", "end", "table", "nDL", "nAR", "nDR", "nDG",
    )

    def __init__(self, g, j, t, fl, fr, v, S, W, m, N):
        self.g = g
        self.j = j
        self.t = t
        self.fl = fl
        self.fr = fr
        self.v = v
        self.S = S
        self.W = W
        self.m = m
        self.N = N
        self.left, self.guard, self.right = regions(j, g, N)
        self.end = len(self.left) + len(self.guard)
        self.table = None
        self.nDL = self.nAR = self.nDR = self.nDG = 0


class Certificate:
    def __init__(self, path: Path | None, data: dict):
        self.path = path
        self.N = int(data["N"])
        self.m = int(data["m"])
        self.r = int(data["r"])
        self.A = int(data["normalizer"])
        self.guards = [int(g) for g in data["guards"]]
        self.rows = [tuple(int(x) for x in row) for row in data["masses"]]
        self.H, self.n_admissible = admissible_lcm(self.N, self.guards, self.m)
        self.comps: list[Component] = []
        for g, j, t, fl, fr, v in self.rows:
            S = shell_size(self.N, g, j, t, fl, fr, self.m)
            if S <= 0 or g not in self.guards or not (0 <= j < self.r) or v <= 0:
                raise ValueError(f"illegal row {(g, j, t, fl, fr, v)} S={S}")
            if self.H % S != 0:
                raise ValueError("shell size does not divide admissible LCM")
            W = v * (self.H // S)
            self.comps.append(Component(g, j, t, fl, fr, v, S, W, self.m, self.N))
        self.mass_num = self.H * sum(c.v for c in self.comps)
        self.mass = Fraction(self.mass_num, self.A)
        self._arm_cache: dict = {}
        self._by_end: dict[int, list[int]] = defaultdict(list)
        for i, c in enumerate(self.comps):
            self._by_end[c.end].append(i)

    @classmethod
    def load(cls, path: Path) -> "Certificate":
        return cls(path, json.loads(path.read_text()))


def arm_ways(a: int, b: int, x: int, nedit: int, nflip: int, cache: dict) -> dict[int, int]:
    """Ways by Hamming distance for one arm.

    a positions agree with the pattern letter, b carry the opposite binary
    letter, x carry a letter in {2, 3}. Exactly nedit edits, nflip of them
    the opposite binary letter, the rest uniform in {2, 3}.
    """
    key = (a, b, x, nedit, nflip)
    hit = cache.get(key)
    if hit is not None:
        return hit
    u = nedit - nflip
    out: dict[int, int] = {}
    if u < 0 or nflip < 0 or a < 0 or b < 0 or x < 0:
        cache[key] = out
        return out
    for fA in range(nflip + 1):
        max_fB = nflip - fA
        for fB in range(max_fB + 1):
            fX = nflip - fA - fB
            if fX > x or fA > a or fB > b:
                continue
            for uA in range(u + 1):
                if fA + uA > a:
                    continue
                max_uB = u - uA
                for uB in range(max_uB + 1):
                    uX = u - uA - uB
                    if fB + uB > b or fX + uX > x:
                        continue
                    ways_pos = (
                        comb(a, fA)
                        * comb(a - fA, uA)
                        * comb(b, fB)
                        * comb(b - fB, uB)
                        * comb(x, fX)
                        * comb(x - fX, uX)
                    )
                    if ways_pos == 0:
                        continue
                    fixed = (
                        fA
                        + fX
                        + uA
                        + uB
                        + (b - fB - uB)
                        + (x - fX - uX)
                    )
                    base = ways_pos << (uA + uB)
                    for k in range(uX + 1):
                        d = fixed + k
                        out[d] = out.get(d, 0) + base * comb(uX, k)
    cache[key] = out
    return out


def gf_count(comp: Component, centre: list[int], cache: dict) -> int:
    P = pattern(comp.j, comp.N)
    gd = 0
    gz = go = gx = 0
    for i in comp.guard:
        ch = centre[i]
        gd += ch != P[i]
        if ch == 0:
            gz += 1
        elif ch == 1:
            go += 1
        elif ch == 2 or ch == 3:
            gx += 1
        else:
            raise ValueError(ch)
    if gd > comp.m:
        return 0
    la = lb = lx = 0
    for i in comp.left:
        ch = centre[i]
        if ch == 1:
            la += 1
        elif ch == 0:
            lb += 1
        else:
            lx += 1
    ra = rb = rx = 0
    for i in comp.right:
        ch = centre[i]
        if ch == 0:
            ra += 1
        elif ch == 1:
            rb += 1
        else:
            rx += 1
    if la + lb + lx != len(comp.left) or ra + rb + rx != len(comp.right):
        raise AssertionError("arm counts")
    left_w = arm_ways(la, lb, lx, comp.t, comp.fl, cache)
    right_w = arm_ways(ra, rb, rx, comp.m - comp.t, comp.fr, cache)
    limit = comp.m - gd
    total = 0
    for d1, w1 in left_w.items():
        if d1 > limit:
            continue
        for d2, w2 in right_w.items():
            if d1 + d2 <= limit:
                total += w1 * w2
    return total


def load_numerator(cert: Certificate, centre: list[int]) -> int:
    total = 0
    cache = cert._arm_cache
    for comp in cert.comps:
        total += gf_count(comp, centre, cache) * comp.W
    return total


def build_tables(cert: Certificate) -> None:
    cache = cert._arm_cache
    m = cert.m
    for comp in cert.comps:
        L = len(comp.left)
        R = len(comp.right)
        nAL, nDL = L + 1, L + 1
        nAR, nDR = R + 1, R + 1
        nDG = m + 1
        table = [0] * (nAL * nDL * nAR * nDR * nDG)
        comp.nDL, comp.nAR, comp.nDR, comp.nDG = nDL, nAR, nDR, nDG
        left_cache = {}
        right_cache = {}
        for al in range(nAL):
            for dl in range(nDL - al):
                xl = L - al - dl
                left_cache[(al, dl)] = arm_ways(al, dl, xl, comp.t, comp.fl, cache)
        for ar in range(nAR):
            for dr in range(nDR - ar):
                xr = R - ar - dr
                right_cache[(ar, dr)] = arm_ways(ar, dr, xr, m - comp.t, comp.fr, cache)
        for (al, dl), lw in left_cache.items():
            for (ar, dr), rw in right_cache.items():
                prefix = (((al * nDL + dl) * nAR + ar) * nDR + dr) * nDG
                for dg in range(nDG):
                    limit = m - dg
                    s = 0
                    if limit >= 0:
                        for d1, w1 in lw.items():
                            if d1 > limit:
                                continue
                            for d2, w2 in rw.items():
                                if d1 + d2 <= limit:
                                    s += w1 * w2
                    table[prefix + dg] = s * comp.W
        comp.table = table


def _lookup(comp: Component, al: int, dl: int, ar: int, dr: int, dg: int) -> int:
    L = len(comp.left)
    R = len(comp.right)
    if dg > comp.m or al < 0 or dl < 0 or ar < 0 or dr < 0:
        return 0
    if al + dl > L or ar + dr > R:
        return 0
    nDL, nAR, nDR, nDG = comp.nDL, comp.nAR, comp.nDR, comp.nDG
    return comp.table[((((al * nDL + dl) * nAR + ar) * nDR + dr) * nDG) + dg]


def dp_one_total(cert: Certificate, Z: int, O: int) -> tuple[int, list[int] | None]:
    """Maximum load numerator over centres with letter counts (Z, O, N-Z-O) in {0,1,2}."""
    N = cert.N
    X = N - Z - O
    if X < 0:
        return -1, None
    by_end = cert._by_end
    cur: dict[tuple[int, int, int], tuple[int, int]] = {(0, 0, 0): (0, -1)}
    # value, parent packed later; parent reconstructed by predecessor map per layer
    parents: list[dict] = []
    for p in range(N):
        nxt: dict[tuple[int, int, int], int] = {}
        pred: dict[tuple[int, int, int], tuple[int, int, int, int]] = {}
        active = by_end.get(p + 1, ())
        for (z, o, hist), (val, _) in cur.items():
            xused = p - z - o
            for letter, nz, no, nx_ok in (
                (0, z + 1, o, xused <= X),
                (1, z, o + 1, xused <= X),
                (2, z, o, xused + 1 <= X),
            ):
                if letter == 0 and z + 1 > Z:
                    continue
                if letter == 1 and o + 1 > O:
                    continue
                if letter == 2 and xused + 1 > X:
                    continue
                nh = ((hist << 2) | letter) & 0xFF
                key = (nz, no, nh)
                old = nxt.get(key)
                if old is None or val > old:
                    nxt[key] = val
                    pred[key] = (z, o, hist, letter)
        if active:
            for key, val in list(nxt.items()):
                z, o, hist = key
                xused = (p + 1) - z - o
                add = 0
                for idx in active:
                    comp = cert.comps[idx]
                    gl = len(comp.guard)
                    h = hist
                    letters = [0] * gl
                    for t in range(gl - 1, -1, -1):
                        letters[t] = h & 3
                        h >>= 2
                    gs = (p + 1) - gl
                    gz = go = gx = gd = 0
                    j = comp.j
                    for t, ch in enumerate(letters):
                        gd += ch != (1 if gs + t < j else 0)
                        if ch == 0:
                            gz += 1
                        elif ch == 1:
                            go += 1
                        else:
                            gx += 1
                    al = o - go
                    dl = z - gz
                    ar = Z - z
                    dr = O - o
                    add += _lookup(comp, al, dl, ar, dr, gd)
                nxt[key] = val + add
        parents.append(pred)
        cur = {k: (v, 0) for k, v in nxt.items()}
    if not cur:
        return -1, None
    (bz, bo, bh), (best, _) = max(cur.items(), key=lambda kv: kv[1][0])
    # reconstruct
    z, o, hist = bz, bo, bh
    letters_rev = []
    for p in range(N - 1, -1, -1):
        pz, po, ph, letter = parents[p][(z, o, hist)]
        letters_rev.append(letter)
        z, o, hist = pz, po, ph
    letters_rev.reverse()
    return best, letters_rev


_WORKER_CERT: Certificate | None = None


def _worker_init(path_str: str):
    global _WORKER_CERT
    cert = Certificate.load(Path(path_str))
    build_tables(cert)
    _WORKER_CERT = cert


def _worker_total(zo: tuple[int, int]) -> tuple[int, int, int, list[int] | None]:
    Z, O = zo
    assert _WORKER_CERT is not None
    val, word = dp_one_total(_WORKER_CERT, Z, O)
    return val, Z, O, word


def exhaustive_max(cert: Certificate, workers: int = MAX_WORKERS) -> dict:
    totals = [(Z, O) for Z in range(cert.N + 1) for O in range(cert.N - Z + 1)]
    best = -1
    best_word = None
    best_zo = None
    # Tables are rebuilt in each worker. The parent also builds them for checks.
    if cert.comps[0].table is None:
        build_tables(cert)
    ctx = ProcessPoolExecutor(
        max_workers=workers,
        initializer=_worker_init,
        initargs=(str(cert.path),),
    )
    try:
        for val, Z, O, word in ctx.map(_worker_total, totals, chunksize=1):
            if val > best:
                best = val
                best_word = word
                best_zo = (Z, O, cert.N - Z - O)
    finally:
        ctx.shutdown()
    return {
        "numerator": best,
        "A": cert.A,
        "load": str(Fraction(best, cert.A)),
        "within_one": best <= cert.A,
        "counts_012": best_zo,
        "word": best_word,
    }


def iter_skeletons(comp: Component):
    t, fl = comp.t, comp.fl
    s, fr = comp.m - comp.t, comp.fr
    for eL in combinations(comp.left, t):
        for fL in combinations(eL, fl):
            fLset = set(fL)
            for eR in combinations(comp.right, s):
                for fR in combinations(eR, fr):
                    yield eL, fLset, eR, set(fR)


def skeleton_word(comp: Component, eL, fLset, eR, fRset, neutral: int = 2) -> list[int]:
    word = pattern(comp.j, comp.N)
    for i in eL:
        word[i] = 0 if i in fLset else neutral
    for i in eR:
        word[i] = 1 if i in fRset else neutral
    return word


def support_audit(cert: Certificate, sample_neutral: int = 400) -> dict:
    """Direct Hamming check of every edit skeleton.

    Library patterns are binary. Both unused letters differ from 0 and from 1,
    so every neutral filling of a skeleton has the same distances to the
    library as the filling that writes 2. One representative per skeleton
    therefore decides membership of every support word.
    """
    checked = 0
    outside = 0
    bad_anchor = 0
    sample_mismatch = 0
    sampled = 0
    rng = random.Random(20261029)
    library = [pattern(j, cert.N) for j in range(cert.r)]
    for comp in cert.comps:
        part = comp.left + comp.guard + comp.right
        if part != list(range(cert.N)):
            raise AssertionError("regions do not partition the coordinates")
        anchor = library[comp.j]
        skeletons = 0
        for eL, fLset, eR, fRset in iter_skeletons(comp):
            skeletons += 1
            word = skeleton_word(comp, eL, fLset, eR, fRset, 2)
            d_anchor = hamming(word, anchor)
            if d_anchor != cert.m:
                bad_anchor += 1
            best = d_anchor
            if best > cert.m:
                # still compute the true minimum; anchor already exceeds
                best = min(hamming(word, P) for P in library)
            if best > cert.m:
                outside += 1
            checked += 1
            if sampled < sample_neutral and rng.randrange(2000) == 0:
                word3 = skeleton_word(comp, eL, fLset, eR, fRset, 3)
                d3 = [hamming(word3, P) for P in library]
                d2 = [hamming(word, P) for P in library]
                sampled += 1
                if d2 != d3:
                    sample_mismatch += 1
        expect = comb(len(comp.left), comp.t) * comb(comp.t, comp.fl)
        expect *= comb(len(comp.right), cert.m - comp.t) * comb(cert.m - comp.t, comp.fr)
        if skeletons != expect:
            raise AssertionError("skeleton count disagrees with binomial coefficient")
        if skeletons * (1 << (cert.m - comp.fl - comp.fr)) != comp.S:
            raise AssertionError("skeleton expansion disagrees with shell size")
    return {
        "skeletons": checked,
        "outside_U": outside,
        "anchor_distance_not_m": bad_anchor,
        "neutral_swap_samples": sampled,
        "neutral_swap_mismatches": sample_mismatch,
        "passed": outside == 0 and bad_anchor == 0 and sample_mismatch == 0,
    }


def literal_count(comp: Component, centre: list[int]) -> int:
    """Independent ball count: enumerate edit sets, DP only over neutral letters."""
    P = pattern(comp.j, comp.N)
    m = comp.m
    gd = sum(centre[i] != P[i] for i in comp.guard)
    if gd > m:
        return 0
    t, fl = comp.t, comp.fl
    s, fr = m - t, comp.fr
    left_set_base = comp.left
    right = comp.right
    cnt = 0
    for eL in combinations(left_set_base, t):
        eLset = set(eL)
        for fL in combinations(eL, fl):
            fLset = set(fL)
            for eR in combinations(right, s):
                eRset = set(eR)
                for fR in combinations(eR, fr):
                    fRset = set(fR)
                    dist = gd
                    for i in left_set_base:
                        if i not in eLset and centre[i] != P[i]:
                            dist += 1
                    for i in right:
                        if i not in eRset and centre[i] != P[i]:
                            dist += 1
                    for i in fLset:
                        dist += centre[i] != 0
                    for i in fRset:
                        dist += centre[i] != 1
                    if dist > m:
                        continue
                    budget = m - dist
                    dp = [1] + [0] * budget
                    neutrals = [i for i in eL if i not in fLset] + [
                        i for i in eR if i not in fRset
                    ]
                    for i in neutrals:
                        d2 = 0 if centre[i] == 2 else 1
                        d3 = 0 if centre[i] == 3 else 1
                        nd = [0] * (budget + 1)
                        for have, w in enumerate(dp):
                            if not w:
                                continue
                            if have + d2 <= budget:
                                nd[have + d2] += w
                            if have + d3 <= budget:
                                nd[have + d3] += w
                        dp = nd
                        if not any(dp):
                            break
                    cnt += sum(dp)
    return cnt


def expand_measure(comps: list[Component]) -> dict[tuple[int, ...], int]:
    weights: dict[tuple[int, ...], int] = defaultdict(int)
    for comp in comps:
        P = pattern(comp.j, comp.N)
        t, fl = comp.t, comp.fl
        s, fr = comp.m - t, comp.fr
        uL = t - fl
        uR = s - fr
        for eL in combinations(comp.left, t):
            for fL in combinations(eL, fl):
                fLset = set(fL)
                neutrals_L = [i for i in eL if i not in fLset]
                for eR in combinations(comp.right, s):
                    for fR in combinations(eR, fr):
                        fRset = set(fR)
                        neutrals_R = [i for i in eR if i not in fRset]
                        neutrals = neutrals_L + neutrals_R
                        nmask = len(neutrals)
                        for mask in range(1 << nmask):
                            word = P[:]
                            for i in fLset:
                                word[i] = 0
                            for i in fRset:
                                word[i] = 1
                            for b, i in enumerate(neutrals):
                                word[i] = 2 if ((mask >> b) & 1) == 0 else 3
                            weights[tuple(word)] += comp.W
    return weights


def ball_from_expansion(weights: dict, centre: list[int], m: int) -> int:
    total = 0
    for word, wt in weights.items():
        if hamming(list(word), centre) <= m:
            total += wt
    return total


def self_test() -> dict:
    rng = random.Random(7)
    failures = []

    def check(name, ok, detail=""):
        if not ok:
            failures.append(f"{name} {detail}")

    # Small full expansion versus generating functions, overlaps included.
    N, m = 5, 2
    raw_rows = []
    for g in (1, 2):
        for j in range(N + 1):
            left, _, right = regions(j, g, N)
            for t in range(m + 1):
                s = m - t
                if t > len(left) or s > len(right):
                    continue
                for fl in range(g):
                    if fl > t:
                        continue
                    for fr in range(g):
                        if fr > s:
                            continue
                        if rng.randrange(3) == 0:
                            raw_rows.append((g, j, t, fl, fr, rng.randint(1, 5)))
    # Force two overlapping-capable rows and one negative-control geometry later.
    if len(raw_rows) < 4:
        raw_rows.append((1, 2, 1, 0, 0, 3))
        raw_rows.append((1, 2, 0, 0, 0, 2))
    data = {
        "N": N,
        "m": m,
        "r": N + 1,
        "normalizer": 1,
        "guards": [1, 2],
        "masses": raw_rows,
    }
    # Normalizer 1 is fine for comparisons of numerators. H still scales W.
    cert = Certificate(None, data)
    # Expansion uses W, GF uses W. Compare on every centre of a subset plus all ternary.
    weights = expand_measure(cert.comps)
    expanded_S = sum(c.S for c in cert.comps)
    check("expansion word-mass", sum(weights.values()) == sum(c.S * c.W for c in cert.comps),
          f"got {sum(weights.values())}")
    # Every expanded word lies at distance m from its contributing anchors; union check.
    outside = 0
    library = [pattern(j, N) for j in range(N + 1)]
    for word in weights:
        if min(hamming(list(word), P) for P in library) > m:
            outside += 1
    check("small support", outside == 0, str(outside))
    cache = cert._arm_cache
    compared = 0
    for _ in range(30):
        centre = [rng.randrange(4) for _ in range(N)]
        gf = sum(gf_count(c, centre, cache) * c.W for c in cert.comps)
        brute = ball_from_expansion(weights, centre, m)
        compared += 1
        if gf != brute:
            failures.append(f"gf!=expansion centre {centre} {gf} {brute}")
            break
    # Literal enumerator versus GF on every component, several centres.
    for centre in (
        pattern(2, N),
        [2] * N,
        [3] * N,
        [0, 1, 2, 3, 1],
    ):
        for comp in cert.comps:
            g = gf_count(comp, centre, cache)
            lit = literal_count(comp, centre)
            if g != lit:
                failures.append(f"gf!=literal {centre} j={comp.j} g={comp.g} {g}!={lit}")
                break
        else:
            continue
        break
    # Symmetry: swapping 2 and 3 does not change the load.
    for _ in range(20):
        centre = [rng.randrange(4) for _ in range(N)]
        swapped = [2 if ch == 3 else ch for ch in centre]
        a = sum(gf_count(c, centre, cache) * c.W for c in cert.comps)
        b = sum(gf_count(c, swapped, cache) * c.W for c in cert.comps)
        if a != b:
            failures.append(f"symmetry {centre}")
            break
    # DP versus exhaustive ternary scan.
    build_tables(cert)
    brute_max = -1
    brute_at = None
    for mask in range(3 ** N):
        centre = []
        x = mask
        for _ in range(N):
            centre.append(x % 3)
            x //= 3
        val = 0
        # use tables through a direct gf to avoid depending on dp's lookup in the oracle
        val = sum(gf_count(c, centre, cache) * c.W for c in cert.comps)
        if val > brute_max:
            brute_max = val
            brute_at = centre
    dp_max = -1
    dp_at = None
    for Z in range(N + 1):
        for O in range(N - Z + 1):
            val, word = dp_one_total(cert, Z, O)
            if val > dp_max:
                dp_max = val
                dp_at = word
    check("dp equals ternary exhaustion", dp_max == brute_max, f"dp {dp_max} brute {brute_max}")
    if dp_at is not None:
        check(
            "dp word attains",
            sum(gf_count(c, dp_at, cache) * c.W for c in cert.comps) == dp_max,
        )
    # Four-letter max equals the ternary max: folding 3 to 2 preserves every load.
    four_max = -1
    for mask in range(4 ** N):
        centre = []
        x = mask
        for _ in range(N):
            centre.append(x % 4)
            x //= 4
        val = sum(gf_count(c, centre, cache) * c.W for c in cert.comps)
        if val > four_max:
            four_max = val
    check("four-letter exhaustion equals ternary", four_max == brute_max, f"{four_max}")
    # Pattern-centre identity: every shell word sits in its anchor ball.
    for comp in cert.comps:
        c0 = pattern(comp.j, N)
        check(
            f"anchor shell {comp.j},{comp.g},{comp.t}",
            gf_count(comp, c0, cache) == comp.S,
            str(gf_count(comp, c0, cache)),
        )
    return {
        "passed": not failures,
        "failures": failures,
        "small_N": N,
        "small_m": m,
        "rows": len(cert.comps),
        "expanded_occurrences": expanded_S,
        "gf_vs_expansion_centres": compared,
        "dp_max": dp_max,
        "brute_max": brute_max,
        "attaining": brute_at,
    }


def structured_and_search(cert: Certificate, workers: int = MAX_WORKERS) -> dict:
    rng = random.Random(11 + cert.N * 100 + cert.r)
    N = cert.N
    seeds: list[list[int]] = []
    for j in range(cert.r):
        seeds.append(pattern(j, N))
    seeds.append([0] * N)
    seeds.append([1] * N)
    seeds.append([2] * N)
    seeds.append([3] * N)
    for j in range(N):
        w = pattern(j, N)
        w[min(j, N - 1)] = 2
        seeds.append(w[:])
        if j + 1 < N:
            w2 = pattern(j, N)
            w2[j] = 3
            seeds.append(w2)
    for width in (1, 2, 3, 4):
        for start in range(0, N - width + 1, max(1, width)):
            for letter in (2, 3):
                w = pattern(min(cert.r - 1, start), N)
                for i in range(start, start + width):
                    w[i] = letter
                seeds.append(w)
    for comp in cert.comps:
        # one explicit support word per component, as a centre
        eL = tuple(comp.left[: comp.t])
        fL = set(eL[: comp.fl])
        eR = tuple(comp.right[: cert.m - comp.t])
        fR = set(eR[: comp.fr])
        seeds.append(skeleton_word(comp, eL, fL, eR, fR, 2))
        seeds.append(skeleton_word(comp, eL, fL, eR, fR, 3))
    for _ in range(48):
        seeds.append([rng.randrange(4) for _ in range(N)])
    # staircase midpoints: on the interval where P_a and P_b differ, take a prefix of the far pattern
    for a in range(cert.r):
        for b in (a + 1, a + 2, a + 4, a + 8):
            if b >= cert.r:
                continue
            w = pattern(a, N)
            span = range(min(a, b), max(a, b))
            for take, letter in ((len(span) // 2, None), (1, 2), (2, 3)):
                w2 = w[:]
                for i in list(span)[:take]:
                    w2[i] = pattern(b, N)[i] if letter is None else letter
                seeds.append(w2)

    def climb(seed: list[int]) -> tuple[int, list[int]]:
        cache = {}
        cur = seed[:]
        cur_val = load_numerator(cert, cur)
        for _step in range(12):
            improved = False
            best_local = cur_val
            best_word = cur
            order = list(range(N))
            rng.shuffle(order)
            for i in order:
                original = cur[i]
                for letter in range(4):
                    if letter == original:
                        continue
                    cur[i] = letter
                    val = load_numerator(cert, cur)
                    if val > best_local:
                        best_local = val
                        best_word = cur[:]
                        improved = True
                cur[i] = original
            # one coordinated pair
            i, k = rng.randrange(N), rng.randrange(N)
            saved = cur[i], cur[k]
            cur[i], cur[k] = rng.randrange(4), rng.randrange(4)
            val = load_numerator(cert, cur)
            if val > best_local:
                best_local = val
                best_word = cur[:]
                improved = True
            else:
                cur[i], cur[k] = saved
            if not improved:
                break
            cur, cur_val = best_word, best_local
        return cur_val, cur

    # De-duplicate seeds
    uniq = []
    seen = set()
    for s in seeds:
        t = tuple(s)
        if t not in seen:
            seen.add(t)
            uniq.append(s)

    best = -1
    best_word = None
    n_eval_seeds = len(uniq)
    # Climb serially but the load function is pure; split seeds across processes
    # by calling climb in the parent to reuse one arm cache. Cache grows and helps.
    over = []
    for seed in uniq:
        val, word = climb(seed)
        if val > best:
            best = val
            best_word = word
        if val > cert.A:
            over.append(word)
            break
    return {
        "seeds": n_eval_seeds,
        "best_numerator": best,
        "best_load": str(Fraction(best, cert.A)),
        "best_word": best_word,
        "refuted": best > cert.A,
        "overloads": len(over),
    }


def nstar(m: int) -> int:
    N = m + 2
    while Fraction(m, N - 2) + Fraction(m, N - m - 1) > 1:
        N += 1
        if N > 10 * m + 20:
            raise RuntimeError("N* search did not finish")
    return N


def closed_form_nstar(m: int) -> int:
    disc = 5 * m * m - 2 * m + 1
    # least k with 2k >= 3m+3 and (2k-3m-3)^2 >= disc
    k = (3 * m + 3 + math.isqrt(disc - 1) + 1) // 2
    # isqrt(disc-1)+1 = ceil(sqrt(disc))
    while True:
        gap = 2 * k - 3 * m - 3
        if gap >= 0 and gap * gap >= disc:
            prev = 2 * (k - 1) - 3 * m - 3
            if k == m + 2 or prev < 0 or prev * prev < disc:
                return k
        k += 1


def global_coverage() -> dict:
    m = 8
    k = nstar(m)
    k_form = closed_form_nstar(m)
    L22 = Fraction(m, 22 - 2) + Fraction(m, 22 - m - 1)
    L23 = Fraction(m, 23 - 2) + Fraction(m, 23 - m - 1)
    L21 = Fraction(m, 21 - 2) + Fraction(m, 21 - m - 1)
    # selector indicator at word length k-1 = 22, full-chain run count k = 23
    indicator = 1 if k * (L22 - 1) >= L22 else 0
    guide_threshold = k + indicator
    rescaled = Fraction(k, 1) / L22
    full_chain_target = k - 1  # need mass > (number of patterns) - 1 = 22
    boundary_ok = rescaled > full_chain_target
    # weighting test r(M-1) < M at the appendix's reported fractions
    M23 = Fraction(66, 65)
    M22 = Fraction(62, 57)
    weight_certifies_n23 = 23 * (M23 - 1) < M23
    weight_certifies_n22_full = 22 * (M22 - 1) < M22
    # R at n=22 from M=62/57: largest r with r(M-1) < M
    R22 = 0
    for r in range(1, 23):
        if r * (M22 - 1) < M22:
            R22 = r

    def reason(n: int, r: int) -> str:
        if not (1 <= r <= n):
            return "impossible"
        if n == 20:
            return "thm:eight"
        if r <= 20:
            return "padding of thm:eight"
        if n == 21 and r == 21:
            return "cell n=21 r=21"
        if n == 22 and r == 22:
            return "cell n=22 r=22"
        if n == 22 and r == 21:
            return "padding of cell (21,21); also step-down from cell (22,22)"
        if n == 23:
            return "selector rescaling at guide 23 (word length 22)"
        if n >= 24:
            return "unconditional floor at word length N=n-1>=23"
        return "UNCOVERED"

    uncovered = []
    for n in range(20, 41):
        for r in range(1, n + 1):
            if reason(n, r) == "UNCOVERED":
                uncovered.append((n, r))
    # general argument past 40: every larger n is in the n>=24 case or, for r<=20,
    # padding. n starts at 41 >= 24.
    return {
        "Nstar_word_length": k,
        "closed_form": k_form,
        "formula_matches_definition": k == k_form,
        "L_word22": str(L22),
        "L_word23": str(L23),
        "L_word21": str(L21),
        "L_word22_equals_66_65": L22 == Fraction(66, 65),
        "L_word21_equals_62_57": L21 == Fraction(62, 57),
        "L_word23_equals_20_21": L23 == Fraction(20, 21),
        "selector_indicator": indicator,
        "guide_threshold": guide_threshold,
        "rescaled_mass_guide23": str(rescaled),
        "rescaled_exceeds_22": boundary_ok,
        "weighting_test_guide23": weight_certifies_n23,
        "weighting_test_guide22_full_run": weight_certifies_n22_full,
        "R_from_M_at_guide22": R22,
        "uncovered_20_to_40": uncovered,
        "three_open_cells_before_certificates": [
            reason(21, 21),
            reason(22, 21),
            reason(22, 22),
        ],
    }


def literal_audit(cert: Certificate, centres: list[list[int]]) -> dict:
    cache = cert._arm_cache
    mismatches = []
    for centre in centres:
        for comp in cert.comps:
            g = gf_count(comp, centre, cache)
            lit = literal_count(comp, centre)
            if g != lit:
                mismatches.append({"j": comp.j, "g": comp.g, "t": comp.t, "gf": g, "literal": lit})
                if len(mismatches) >= 5:
                    return {"compared_ok": False, "mismatches": mismatches}
    return {"centres": len(centres), "components": len(cert.comps), "compared_ok": not mismatches, "mismatches": mismatches}


def make_overloaded(cert: Certificate) -> Certificate:
    data = {
        "N": cert.N,
        "m": cert.m,
        "r": cert.r,
        "normalizer": cert.A,
        "guards": cert.guards,
        "masses": [list(row) for row in cert.rows],
    }
    # One extra unit of mass on the first shell, before dividing by A.
    bump = cert.A // cert.H + 1
    data["masses"][0][5] += bump
    forged = Certificate(cert.path, data)
    return forged


def control_overload(cert: Certificate) -> dict:
    forged = make_overloaded(cert)
    anchor = pattern(forged.comps[0].j, forged.N)
    num = load_numerator(forged, anchor)
    component_mass = Fraction(forged.comps[0].v * forged.H, forged.A)
    return {
        "anchor_numerator": num,
        "A": forged.A,
        "anchor_load": str(Fraction(num, forged.A)),
        "component_mass": str(component_mass),
        "checker_rejects": num > forged.A,
        "mass_still_above": forged.mass > forged.r - 1,
    }


def control_outside_word(cert: Certificate) -> dict:
    word = [2] * cert.N
    library = [pattern(j, cert.N) for j in range(cert.r)]
    dists = [hamming(word, P) for P in library]
    return {
        "word": "2" * cert.N,
        "min_distance": min(dists),
        "budget": cert.m,
        "checker_rejects": min(dists) > cert.m,
    }


def certificate_paths() -> list[Path]:
    wanted = [
        "packing_N20_r21_m8.json",
        "packing_N21_r21_m8.json",
        "packing_N21_r22_m8.json",
    ]
    return [CERT_DIR / name for name in wanted]


def dump(obj) -> None:
    print(json.dumps(obj, indent=2, default=str))


def cmd_self_test(_args) -> int:
    result = self_test()
    dump(result)
    return 0 if result["passed"] else 1


def cmd_mass(_args) -> int:
    rows = []
    ok = True
    for path in certificate_paths():
        cert = Certificate.load(path)
        claimed = None
        raw = json.loads(path.read_text())
        if "total_mass" in raw:
            num_s, den_s = raw["total_mass"].split("/")
            claimed = Fraction(int(num_s), int(den_s))
        row = {
            "file": path.name,
            "N": cert.N,
            "n_guide": cert.N + 1,
            "r": cert.r,
            "components": len(cert.comps),
            "admissible_shells": cert.n_admissible,
            "H": str(cert.H),
            "mass": str(cert.mass),
            "exceeds_r_minus_1": cert.mass > cert.r - 1,
            "matches_printed_fraction": claimed == cert.mass,
            "sum_shell_sizes": sum(c.S for c in cert.comps),
        }
        ok = ok and row["exceeds_r_minus_1"]
        rows.append(row)
    dump({"passed": ok, "rows": rows})
    return 0 if ok else 1


def cmd_support(_args) -> int:
    ok = True
    rows = []
    for path in certificate_paths():
        cert = Certificate.load(path)
        audit = support_audit(cert)
        audit["file"] = path.name
        ok = ok and audit["passed"]
        rows.append(audit)
    dump({"passed": ok, "rows": rows})
    return 0 if ok else 1


def cmd_literal(_args) -> int:
    ok = True
    rows = []
    rng = random.Random(19)
    for path in certificate_paths():
        cert = Certificate.load(path)
        centres = [
            pattern(0, cert.N),
            pattern(cert.r // 2, cert.N),
            pattern(cert.r - 1, cert.N),
            [2] * cert.N,
            [3 if i % 2 == 0 else 2 for i in range(cert.N)],
            [rng.randrange(4) for _ in range(cert.N)],
        ]
        audit = literal_audit(cert, centres)
        audit["file"] = path.name
        # symmetry spot check
        mixed = [rng.choice((0, 1, 2, 3)) for _ in range(cert.N)]
        folded = [2 if ch == 3 else ch for ch in mixed]
        same = load_numerator(cert, mixed) == load_numerator(cert, folded)
        audit["symmetry_spot"] = same
        ok = ok and audit["compared_ok"] and same
        rows.append(audit)
    dump({"passed": ok, "rows": rows})
    return 0 if ok else 1


def cmd_search(_args) -> int:
    ok = True
    rows = []
    for path in certificate_paths():
        cert = Certificate.load(path)
        found = structured_and_search(cert, workers=1)
        found["file"] = path.name
        ok = ok and not found["refuted"]
        rows.append(found)
    dump({"refuted": not ok, "rows": rows})
    return 1 if not ok else 0


def cmd_dp(args) -> int:
    ok = True
    rows = []
    workers = min(MAX_WORKERS, args.workers)
    for path in certificate_paths():
        cert = Certificate.load(path)
        result = exhaustive_max(cert, workers=workers)
        # recompute winner by the direct sum
        if result["word"] is not None:
            direct = load_numerator(cert, result["word"])
            result["direct_numerator"] = direct
            result["direct_matches_dp"] = direct == result["numerator"]
            if not result["direct_matches_dp"]:
                ok = False
        result["file"] = path.name
        result["exceeds_one"] = result["numerator"] > cert.A
        if result["word"] is not None:
            audit = literal_audit(cert, [result["word"], pattern(cert.r // 2, cert.N)])
            result["literal_on_winner"] = audit["compared_ok"]
            ok = ok and audit["compared_ok"]
        ok = ok and not result["exceeds_one"]
        # do not dump huge ints without strings
        result["numerator"] = str(result["numerator"])
        result["A"] = str(result["A"])
        rows.append(result)
        dump({"partial": result["file"], "load": result["load"], "within_one": result["within_one"]})
    dump({"within_one": ok, "rows": rows})
    return 0 if ok else 1


def cmd_global(_args) -> int:
    result = global_coverage()
    result["passed"] = (
        result["formula_matches_definition"]
        and result["Nstar_word_length"] == 23
        and result["guide_threshold"] == 23
        and result["rescaled_exceeds_22"]
        and result["L_word22_equals_66_65"]
        and result["L_word21_equals_62_57"]
        and result["weighting_test_guide23"]
        and not result["weighting_test_guide22_full_run"]
        and result["uncovered_20_to_40"] == []
    )
    dump(result)
    return 0 if result["passed"] else 1


def cmd_control(_args) -> int:
    cert = Certificate.load(certificate_paths()[0])
    over = control_overload(cert)
    outside = control_outside_word(cert)
    passed = over["checker_rejects"] and outside["checker_rejects"]
    dump({"passed": passed, "overload": over, "outside_support": outside})
    return 0 if passed else 1


def cmd_verify(args) -> int:
    """Run the exact certificate checks and both positive and negative controls."""
    for name, check in (("self-test", cmd_self_test), ("mass", cmd_mass),
                        ("support", cmd_support), ("dp", cmd_dp),
                        ("control", cmd_control)):
        print("CHECK " + name, flush=True)
        status = check(args)
        if status:
            return status
    print("VERIFIED: all three certificates have support in U, mass > r-1, "
          "and maximum radius-8 load <= 1 over every four-letter centre.", flush=True)
    return 0


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)
    p_verify = sub.add_parser("verify")
    p_verify.add_argument("--workers", type=int, choices=range(1, MAX_WORKERS + 1), default=4)
    p_verify.set_defaults(func=cmd_verify)
    sub.add_parser("self-test").set_defaults(func=cmd_self_test)
    sub.add_parser("mass").set_defaults(func=cmd_mass)
    sub.add_parser("support").set_defaults(func=cmd_support)
    sub.add_parser("literal").set_defaults(func=cmd_literal)
    sub.add_parser("search").set_defaults(func=cmd_search)
    p_dp = sub.add_parser("dp")
    p_dp.add_argument("--workers", type=int, default=MAX_WORKERS)
    p_dp.set_defaults(func=cmd_dp)
    sub.add_parser("global").set_defaults(func=cmd_global)
    sub.add_parser("control").set_defaults(func=cmd_control)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
