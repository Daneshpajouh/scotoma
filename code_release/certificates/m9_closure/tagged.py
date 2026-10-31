"""Guarded shells with one tagged cell on each side of the guard.

The tag is the coordinate immediately outside the guard. It is fixed to the
pattern letter, flipped to the other binary letter, or replaced by a
non-binary letter. The far arm carries a uniform edit count. Every word in a
shell differs from its anchor pattern in exactly m coordinates, so the shell
sits inside U.
"""

from __future__ import annotations

from math import comb

from shells import build_poly, pattern


def shell_size(N: int, m: int, g: int, j: int, t: int, u: int, fl: int, fr: int, tagL: int, tagR: int) -> int:
    if g < 1 or j < 0 or j > N:
        return 0
    if tagL not in (-1, 0, 1, 2) or tagR not in (-1, 0, 1, 2):
        return 0
    hasL = j - g - 1 >= 0
    hasR = j + g < N
    if (tagL >= 0) != hasL or (tagR >= 0) != hasR:
        return 0
    extra = (tagL in (1, 2)) + (tagR in (1, 2))
    if t + u + extra != m or min(t, u, fl, fr) < 0 or fl > t or fr > u:
        return 0
    L = (j - g - 1) if hasL else max(0, j - g)
    R = (N - (j + g + 1)) if hasR else max(0, N - (j + g))
    if t > L or u > R:
        return 0
    ways = comb(L, t) * comb(t, fl) * comb(R, u) * comb(u, fr) * (1 << (t - fl + u - fr))
    if tagL == 2:
        ways *= 2
    if tagR == 2:
        ways *= 2
    return ways


def list_shells(N: int, m: int, guards: tuple[int, ...], flip_cap: int, r: int | None = None):
    if r is None:
        r = N + 1
    out = []
    for g in guards:
        for j in range(r):
            tagL_opts = (-1,) if j - g - 1 < 0 else (0, 1, 2)
            tagR_opts = (-1,) if j + g >= N else (0, 1, 2)
            for tagL in tagL_opts:
                for tagR in tagR_opts:
                    extra = (tagL in (1, 2)) + (tagR in (1, 2))
                    L = (j - g - 1) if tagL >= 0 else max(0, j - g)
                    R = (N - (j + g + 1)) if tagR >= 0 else max(0, N - (j + g))
                    for t in range(m + 1):
                        u = m - extra - t
                        if u < 0 or t > L or u > R:
                            continue
                        for fl in range(min(flip_cap, t) + 1):
                            for fr in range(min(flip_cap, u) + 1):
                                out.append((g, j, t, u, fl, fr, tagL, tagR))
    return out


def tag_multiplicity(tag: int, centre_letter: int, pattern_letter: int, m: int) -> list[int]:
    """Distance -> number of tag letters."""
    out = [0] * (m + 1)
    if tag < 0:
        out[0] = 1
        return out
    if tag == 0:
        out[0 if centre_letter == pattern_letter else 1] = 1
        return out
    if tag == 1:
        flip = 1 - pattern_letter
        out[0 if centre_letter == flip else 1] = 1
        return out
    if centre_letter >= 2:
        out[0] += 1
        out[1] += 1
    else:
        out[1] += 2
    return out


class Counter:
    def __init__(self, N: int, m: int, guards: tuple[int, ...], flip_cap: int, r: int | None = None):
        self.N = N
        self.m = m
        self.guards = tuple(guards)
        self.flip_cap = flip_cap
        self.r = N + 1 if r is None else r
        self.shells = list_shells(N, m, self.guards, flip_cap, self.r)
        self.sizes = [shell_size(N, m, *sh) for sh in self.shells]
        if any(s <= 0 for s in self.sizes):
            raise RuntimeError("shell list contains an empty shell")
        self.poly = build_poly(N, m, flip_cap)
        self._pref: dict[tuple, list[int]] = {}
        for key, row in self.poly.items():
            acc = [0]
            total = 0
            for v in row:
                total += v
                acc.append(total)
            self._pref[key] = acc

    def count(self, centre: tuple[int, ...] | list[int], index: int) -> int:
        g, j, t, u, fl, fr, tagL, tagR = self.shells[index]
        N = self.N
        m = self.m
        left_guard = max(0, j - g)
        right_guard = min(N, j + g)
        e = 0
        for i in range(left_guard, j):
            e += centre[i] != 1
        for i in range(j, right_guard):
            e += centre[i] != 0
        if e > m:
            return 0
        if tagL >= 0:
            far_hi = j - g - 1
            nearL = j - g - 1
        else:
            far_hi = left_guard
            nearL = None
        if tagR >= 0:
            r0 = j + g + 1
            nearR = j + g
        else:
            r0 = right_guard
            nearR = None
        ag = dg = cg = 0
        for i in range(far_hi):
            ch = centre[i]
            if ch == 1:
                ag += 1
            elif ch == 0:
                dg += 1
            else:
                cg += 1
        ar = dr = cr = 0
        for i in range(r0, N):
            ch = centre[i]
            if ch == 0:
                ar += 1
            elif ch == 1:
                dr += 1
            else:
                cr += 1
        left_tag = tag_multiplicity(tagL, centre[nearL] if nearL is not None else 0, 1, m)
        right_tag = tag_multiplicity(tagR, centre[nearR] if nearR is not None else 0, 0, m)
        left = self._pref.get((ag, dg, cg, t, fl))
        right = self._pref.get((ar, dr, cr, u, fr))
        if left is None or right is None:
            return 0
        total = 0
        for dl in range(m + 1):
            lv = left[dl + 1] - left[dl]
            if not lv:
                continue
            for dtl, ltv in enumerate(left_tag):
                if not ltv:
                    continue
                for dr in range(m - dl - dtl + 1):
                    rv = right[dr + 1] - right[dr]
                    if not rv:
                        continue
                    base = dl + dtl + dr
                    for dtr, rtv in enumerate(right_tag):
                        if rtv and base + dtr + e <= m:
                            total += lv * ltv * rv * rtv
        return total

    def numerator(self, centre: tuple[int, ...], coeffs: list[int]) -> int:
        total = 0
        for i, coeff in enumerate(coeffs):
            if coeff:
                total += coeff * self.count(centre, i)
        return total


def anchor_of(shell, N: int) -> tuple[int, ...]:
    return tuple(pattern(shell[1], N))
