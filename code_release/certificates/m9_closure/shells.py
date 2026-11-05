"""Guarded private shells and exact ball intersection counts.

A shell (g, j, t, f_left, f_right) fixes the guard coordinates
j-g, ..., j+g-1, clipped to the word, to the pattern P_j = 1^j 0^{N-j}.
It places t edits on the left arm and m-t edits on the right arm.
Exactly f_left and f_right of those edits are the opposite binary letter.
The rest are letter 2 or letter 3. Both flip counts are at most g-1.
"""

from __future__ import annotations

from math import comb, lcm


def pattern(j: int, N: int) -> list[int]:
    return [1 if i < j else 0 for i in range(N)]


def regions(j: int, g: int, N: int) -> tuple[int, int, int]:
    left = max(0, j - g)
    guard_end = min(N, j + g)
    return left, guard_end, N - guard_end


def shell_size(N: int, m: int, g: int, j: int, t: int, fl: int, fr: int) -> int:
    L, guard_end, R = regions(j, g, N)
    # regions returns (left arm length, guard end index, right arm length)
    s = m - t
    if min(g, t, s, fl, fr, L, R) < 0:
        return 0
    if t > L or s > R or fl > t or fr > s:
        return 0
    return comb(L, t) * comb(t, fl) * comb(R, s) * comb(s, fr) * (1 << (m - fl - fr))


def list_shells(N: int, m: int, guards: tuple[int, ...], r: int | None = None):
    if r is None:
        r = N + 1
    out = []
    for g in guards:
        for j in range(r):
            L = max(0, j - g)
            R = max(0, N - j - g)
            for t in range(m + 1):
                s = m - t
                if t > L or s > R:
                    continue
                for fl in range(min(g - 1, t) + 1):
                    for fr in range(min(g - 1, s) + 1):
                        out.append((g, j, t, fl, fr))
    return out


def admissible_lcm(N: int, m: int, guards: tuple[int, ...]) -> tuple[int, int]:
    H = 1
    count = 0
    for shell in list_shells(N, m, guards, r=N + 1):
        S = shell_size(N, m, *shell)
        H = lcm(H, S)
        count += 1
    return H, count


def build_poly(N: int, m: int, max_flip: int) -> dict:
    """Single-parent product of the one-position factors.

    Key (agree, disagree, nonbinary, edits, flips) -> distance multiplicities.
    """
    dist: dict[tuple[int, int, int, int, int], list[int]] = {}

    def put(a: int, b: int, c: int, s: int, f: int, d: int, v: int) -> None:
        if v <= 0 or s > m or f > max_flip or d > m:
            return
        key = (a, b, c, s, f)
        row = dist.get(key)
        if row is None:
            row = [0] * (m + 1)
            dist[key] = row
        row[d] += v

    put(0, 0, 0, 0, 0, 0, 1)
    layers: list[dict[tuple[int, int, int], list[tuple[int, int, int]]]] = [
        {(0, 0, 0): [(0, 0, 0)]}
    ]
    # Store sparse source triples per length, rebuilt from dist keys of that length.
    have = {(0, 0, 0)}
    for length in range(1, N + 1):
        new_have = []
        for a in range(length + 1):
            for b in range(length - a + 1):
                c = length - a - b
                if c > 0:
                    pa, pb, pc, kind = a, b, c - 1, 2
                elif b > 0:
                    pa, pb, pc, kind = a, b - 1, c, 1
                else:
                    pa, pb, pc, kind = a - 1, b, c, 0
                if (pa, pb, pc) not in have and length > 1:
                    raise RuntimeError("missing parent")
                sources = []
                # Gather every stored (s,f) row of the parent.
                # Parents exist because length grows by one and (0,0,0) starts the chain.
                prefix = (pa, pb, pc)
                # Scan possible s,f. Parent rows live in dist.
                for s in range(m + 1):
                    for f in range(max_flip + 1):
                        row = dist.get(prefix + (s, f))
                        if row:
                            sources.append((s, f, row))
                for s, f, row in sources:
                    for d, v in enumerate(row):
                        if not v:
                            continue
                        if kind == 0:
                            put(a, b, c, s, f, d, v)
                            put(a, b, c, s + 1, f, d + 1, 2 * v)
                            put(a, b, c, s + 1, f + 1, d + 1, v)
                        elif kind == 1:
                            put(a, b, c, s, f, d + 1, v)
                            put(a, b, c, s + 1, f, d + 1, 2 * v)
                            put(a, b, c, s + 1, f + 1, d, v)
                        else:
                            put(a, b, c, s, f, d + 1, v)
                            put(a, b, c, s + 1, f, d, v)
                            put(a, b, c, s + 1, f, d + 1, v)
                            put(a, b, c, s + 1, f + 1, d + 1, v)
                new_have.append((a, b, c))
        have.update(new_have)
    return dist


class Counter:
    def __init__(self, N: int, m: int, guards: tuple[int, ...], r: int | None = None):
        self.N = N
        self.m = m
        self.guards = tuple(guards)
        self.r = N + 1 if r is None else r
        self.shells = list_shells(N, m, self.guards, self.r)
        self.sizes = [shell_size(N, m, *sh) for sh in self.shells]
        self.poly = build_poly(N, m, max(self.guards) - 1)
        self._pref: dict[tuple, list[int]] = {}
        for key, row in self.poly.items():
            acc = [0]
            for v in row:
                acc.append(acc[-1] + v)
            self._pref[key] = acc

    def count(self, centre: tuple[int, ...] | list[int], index: int) -> int:
        g, j, t, fl, fr = self.shells[index]
        N = self.N
        left_end = max(0, j - g)
        right_start = min(N, j + g)
        ag = dg = cg = 0
        for i in range(left_end):
            ch = centre[i]
            if ch == 1:
                ag += 1
            elif ch == 0:
                dg += 1
            else:
                cg += 1
        ar = dr = cr = 0
        for i in range(right_start, N):
            ch = centre[i]
            if ch == 0:
                ar += 1
            elif ch == 1:
                dr += 1
            else:
                cr += 1
        e = 0
        for i in range(left_end, j):
            e += centre[i] != 1
        for i in range(j, right_start):
            e += centre[i] != 0
        if e > self.m:
            return 0
        left = self._pref.get((ag, dg, cg, t, fl))
        right = self._pref.get((ar, dr, cr, self.m - t, fr))
        if left is None or right is None:
            return 0
        total = 0
        limit = self.m - e
        for dl in range(limit + 1):
            lv = left[dl + 1] - left[dl]
            if lv:
                total += lv * right[limit - dl + 1]
        return total

    def row(self, centre: tuple[int, ...]) -> list[float]:
        return [self.count(centre, i) / self.sizes[i] for i in range(len(self.shells))]
