"""Resolved near-cut shells for the single-deletion union at budget m.

A shell is anchored at one pattern P_j = 1^j 0^{N-j}. The guard is the cut
window [max(0, j-1), min(N, j+1)). In fixed-guard mode that window matches
P_j. In free-guard mode each guard site is a match, a binary flip, or a
star in {2, 3}. The next `depth` sites on each side are likewise specified
one by one. The far arms, beyond that resolved band, carry a uniform edit
count, of which a capped number are binary flips and the rest are stars.

Every concrete word differs from P_j in exactly m coordinates, so the shell
lies in the radius-m ball of P_j and therefore in the target union.
"""

from __future__ import annotations

from math import comb


MATCH, FLIP, STAR = 0, 1, 2


def pattern(j: int, N: int) -> list[int]:
    return [1 if i < j else 0 for i in range(N)]


def guard_bounds(j: int, N: int) -> tuple[int, int]:
    return max(0, j - 1), min(N, j + 1)


def layout(j: int, N: int, depth: int) -> dict:
    g0, g1 = guard_bounds(j, N)
    guard = list(range(g0, g1))
    left = list(range(g0 - 1, -1, -1))
    right = list(range(g1, N))
    bL = min(depth, len(left))
    bR = min(depth, len(right))
    return {
        "guard": guard,
        "resL": left[:bL],
        "resR": right[:bR],
        "farL": left[bL:],
        "farR": right[bR:],
    }


def _edit_count(symbols: tuple[int, ...]) -> int:
    return sum(s != MATCH for s in symbols)


def list_shells(
    N: int,
    m: int,
    depth: int,
    flip_cap: int,
    r: int | None = None,
    free_guard: bool = False,
) -> list[tuple]:
    if r is None:
        r = N + 1
    shells = []
    for j in range(r):
        lay = layout(j, N, depth)
        g_opts = []
        if free_guard:
            gpos = lay["guard"]
            for mask in range(3 ** len(gpos) if gpos else 1):
                syms = []
                x = mask
                for _ in gpos:
                    syms.append(x % 3)
                    x //= 3
                g_opts.append(tuple(syms))
        else:
            g_opts.append(tuple(MATCH for _ in lay["guard"]))
        nL, nR = len(lay["resL"]), len(lay["resR"])
        for gsym in g_opts:
            eg = _edit_count(gsym)
            if eg > m:
                continue
            for maskL in range(3 ** nL if nL else 1):
                confL = []
                x = maskL
                for _ in range(nL):
                    confL.append(x % 3)
                    x //= 3
                eL = _edit_count(confL)
                if eg + eL > m:
                    continue
                for maskR in range(3 ** nR if nR else 1):
                    confR = []
                    x = maskR
                    for _ in range(nR):
                        confR.append(x % 3)
                        x //= 3
                    eR = _edit_count(confR)
                    used = eg + eL + eR
                    if used > m:
                        continue
                    rem = m - used
                    fL, fR = len(lay["farL"]), len(lay["farR"])
                    for t in range(rem + 1):
                        u = rem - t
                        if t > fL or u > fR:
                            continue
                        for fl in range(min(flip_cap, t) + 1):
                            for fr in range(min(flip_cap, u) + 1):
                                shells.append(
                                    (j, gsym, tuple(confL), tuple(confR), t, u, fl, fr)
                                )
    return shells


def shell_size(N: int, m: int, depth: int, shell: tuple) -> int:
    j, gsym, confL, confR, t, u, fl, fr = shell
    lay = layout(j, N, depth)
    if len(gsym) != len(lay["guard"]) or len(confL) != len(lay["resL"]) or len(confR) != len(lay["resR"]):
        return 0
    if t > len(lay["farL"]) or u > len(lay["farR"]) or fl > t or fr > u:
        return 0
    if _edit_count(gsym) + _edit_count(confL) + _edit_count(confR) + t + u != m:
        return 0
    ways = 1
    for s in gsym + confL + confR:
        if s == STAR:
            ways *= 2
    ways *= comb(len(lay["farL"]), t) * comb(t, fl)
    ways *= comb(len(lay["farR"]), u) * comb(u, fr)
    ways *= 1 << (t - fl + u - fr)
    return ways


def pattern_letter(j: int, i: int) -> int:
    return 1 if i < j else 0


def skeleton_distance(N: int, m: int, depth: int, shell: tuple, anchor: int) -> int | None:
    """Distance from every realization to P_anchor.

    Stars disagree with both binary letters, so the distance does not depend
    on the star filling. Returns None if the shell is not a legal m-edit of
    its own anchor.
    """
    j, gsym, confL, confR, t, u, fl, fr = shell
    if shell_size(N, m, depth, shell) <= 0:
        return None
    lay = layout(j, N, depth)
    # Build the disagreement set relative to P_anchor, counting a star as a disagreement
    # with either binary letter. Fixed matches contribute only where the written
    # letter differs from P_anchor.
    disagree = 0
    written = {}
    for pos, sym in zip(lay["guard"], gsym):
        written[pos] = sym
    for pos, sym in zip(lay["resL"], confL):
        written[pos] = sym
    for pos, sym in zip(lay["resR"], confR):
        written[pos] = sym
    # Far edits are not located. Distance to a foreign anchor is not constant
    # on the shell. Distance to the shell's own anchor is exactly m.
    if anchor == j:
        return m
    return None


def min_anchor_count(N: int, m: int, depth: int, shell: tuple, r: int) -> int:
    """How many library patterns meet every realization at distance <= m.

    Far-arm placements move the distance to a foreign pattern, so a shell is
    counted as multi-anchor only when every placement stays within m of a
    second pattern. The test uses the worst placement of the far edits.
    """
    j, gsym, confL, confR, t, u, fl, fr = shell
    lay = layout(j, N, depth)
    # Letter written, or None if the site is a far site that might be edited.
    # For the worst case against a foreign pattern, an editable far site can
    # be set to the foreign pattern's letter (decreasing distance) or not.
    # "Every realization" means the maximum, over realizations, of the distance
    # is still <= m. Stars are always disagreements with binary patterns.
    base = [pattern_letter(j, i) for i in range(N)]
    forced_star = set()
    forced_flip = set()
    for pos, sym in zip(lay["guard"], gsym):
        if sym == STAR:
            forced_star.add(pos)
        elif sym == FLIP:
            forced_flip.add(pos)
            base[pos] = 1 - pattern_letter(j, pos)
    for pos, sym in zip(lay["resL"], confL):
        if sym == STAR:
            forced_star.add(pos)
        elif sym == FLIP:
            forced_flip.add(pos)
            base[pos] = 1 - pattern_letter(j, pos)
    for pos, sym in zip(lay["resR"], confR):
        if sym == STAR:
            forced_star.add(pos)
        elif sym == FLIP:
            forced_flip.add(pos)
            base[pos] = 1 - pattern_letter(j, pos)
    far = set(lay["farL"]) | set(lay["farR"])
    # A far site that is edited disagrees with P_j. Relative to P_k it may
    # agree or disagree. The worst realization maximizes d(w, P_k).
    # On a far site, the adversary picks the letter in the shell's alphabet
    # at that site which disagrees with P_k, if any such letter is allowed.
    # Unedited far sites are fixed to P_j. Edited far sites are either a flip
    # or a star. The worst mix is at least as far as making every edited far
    # site disagree with P_k, which is always possible: a star disagrees with
    # both binary letters, and a flip is one specific letter. If the shell
    # allows a star on that edit (fl < t or the edit can be a star), the
    # adversary can choose a star. If ALL far edits on that side are flips
    # (fl == t and t > 0), those sites are the opposite of P_j, a fixed letter.
    # Conservative exact check: enumerate is impossible. Compute max distance
    # by giving the adversary, for each edited slot, the worse allowed letter.
    # Because the slots are chosen among the far positions, the adversary also
    # chooses WHICH positions are edited. Max over placements of the sum.
    hits = 0
    for k in range(r):
        Pk = pattern(k, N)
        # Fixed sites (not in far):
        dist = 0
        for i in range(N):
            if i in far:
                continue
            if i in forced_star:
                dist += 1
            else:
                dist += base[i] != Pk[i]
        # Far left: t edits. Each edited site contributes 1 if the adversary
        # can disagree with Pk, else 0. Unedited sites contribute 1 iff P_j disagrees with Pk.
        dist += _worst_far(lay["farL"], t, fl, j, Pk)
        dist += _worst_far(lay["farR"], u, fr, j, Pk)
        if dist <= m:
            hits += 1
    return hits


def _worst_far(positions: list[int], edits: int, flips: int, j: int, Pk: list[int]) -> int:
    """Maximum, over realizations, of the distance contribution of one far arm to Pk.

    Unedited sites stay at P_j. Exactly `edits` sites change. At least
    `edits - flips` of those are stars, and `flips` are binary flips. Stars
    always disagree with Pk. A flip writes 1 - P_j, which may agree with Pk.
    To maximize distance, put stars on sites where P_j already agrees with Pk
    (so the edit adds a disagreement), and put unavoidable flips on sites
    where the flip also disagrees, if possible.
    """
    if edits == 0:
        return sum(pattern_letter(j, i) != Pk[i] for i in positions)
    stars = edits - flips
    # Contribution if unedited, and the best added cost of a star or a flip.
    unedited = []
    star_gain = []
    flip_gain = []
    for i in positions:
        pj = pattern_letter(j, i)
        base = pj != Pk[i]
        unedited.append(base)
        star_gain.append(1 - base)  # star always disagrees: new value 1, gain 1-base
        flip_letter = 1 - pj
        flip_gain.append((flip_letter != Pk[i]) - base)
    # Start from all-unedited distance, then choose `stars` star-edits and
    # `flips` flip-edits on distinct sites to maximize the gain.
    # Greedy: assign each site its best allowed role. Exact for this additive
    # objective by taking the top gains.
    options = []
    for i in range(len(positions)):
        options.append((star_gain[i], 0, i))
        options.append((flip_gain[i], 1, i))
    # Knapsack on <= 26 sites is tiny: DP.
    n = len(positions)
    # dp[a][b][used mask] too big if n=20. n can be ~20, 2^20 no.
    # Gains are small integers. Use DP on (site, stars left, flips left).
    neg = -10**6
    dp = [[neg] * (flips + 1) for _ in range(stars + 1)]
    dp[0][0] = 0
    for i in range(n):
        nxt = [[neg] * (flips + 1) for _ in range(stars + 1)]
        for s in range(stars + 1):
            for f in range(flips + 1):
                cur = dp[s][f]
                if cur <= neg // 2:
                    continue
                # leave unedited
                if cur > nxt[s][f]:
                    nxt[s][f] = cur
                if s < stars and cur + star_gain[i] > nxt[s + 1][f]:
                    nxt[s + 1][f] = cur + star_gain[i]
                if f < flips and cur + flip_gain[i] > nxt[s][f + 1]:
                    nxt[s][f + 1] = cur + flip_gain[i]
        dp = nxt
    gain = dp[stars][flips]
    return sum(unedited) + gain


def build_arm_poly(max_len: int, m: int, flip_cap: int) -> dict:
    """(agree, disagree, nonbinary, edits, flips) -> ways at each distance.

    One canonical order of position types, matching the guarded-shell polynomial:
    the count depends only on how many arm sites agree with the pattern, how
    many are the opposite binary letter, and how many are nonbinary.
    """
    dist: dict[tuple[int, int, int, int, int], list[int]] = {}

    def put(a, b, c, s, f, d, v):
        if v <= 0 or s > m or f > flip_cap or d > m or a < 0:
            return
        key = (a, b, c, s, f)
        row = dist.get(key)
        if row is None:
            row = [0] * (m + 1)
            dist[key] = row
        row[d] += v

    put(0, 0, 0, 0, 0, 0, 1)
    have = {(0, 0, 0)}
    for length in range(1, max_len + 1):
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
                    raise RuntimeError("missing arm parent")
                for s in range(m + 1):
                    for f in range(flip_cap + 1):
                        row = dist.get((pa, pb, pc, s, f))
                        if not row:
                            continue
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
    pref = {}
    for key, row in dist.items():
        acc = [0]
        for v in row:
            acc.append(acc[-1] + v)
        pref[key] = acc
    return pref


def _site_letter_cost(symbol: int, centre_letter: int, pat: int, m: int) -> tuple[list[int], int]:
    """Return (ways_at_extra_distance of length m+1, unused)."""
    ways = [0] * (m + 1)
    if symbol == MATCH:
        ways[0 if centre_letter == pat else 1] = 1
        return ways
    if symbol == FLIP:
        flip = 1 - pat
        ways[0 if centre_letter == flip else 1] = 1
        return ways
    # star: two letters, 2 and 3
    if centre_letter >= 2:
        ways[0] += 1
        ways[1] += 1
    else:
        ways[1] += 2
    return ways


def count_in_ball(centre: list[int] | tuple[int, ...], N: int, m: int, depth: int, shell: tuple, pref: dict) -> int:
    j, gsym, confL, confR, t, u, fl, fr = shell
    lay = layout(j, N, depth)
    # Generating function over the fixed window, then the two far arms.
    cur = [0] * (m + 1)
    cur[0] = 1

    def convolve(gen, extra):
        out = [0] * (m + 1)
        for d, v in enumerate(gen):
            if not v:
                continue
            for e, w in enumerate(extra):
                if w and d + e <= m:
                    out[d + e] += v * w
        return out

    for pos, sym in zip(lay["guard"], gsym):
        cur = convolve(cur, _site_letter_cost(sym, centre[pos], pattern_letter(j, pos), m))
        if sum(cur) == 0:
            return 0
    for pos, sym in zip(lay["resL"], confL):
        cur = convolve(cur, _site_letter_cost(sym, centre[pos], 1, m))
    for pos, sym in zip(lay["resR"], confR):
        cur = convolve(cur, _site_letter_cost(sym, centre[pos], 0, m))

    def arm_types(positions, agree):
        ag = dg = cg = 0
        for i in positions:
            ch = centre[i]
            if ch == agree:
                ag += 1
            elif ch < 2:
                dg += 1
            else:
                cg += 1
        return ag, dg, cg

    ag, dg, cg = arm_types(lay["farL"], 1)
    ar, dr, cr = arm_types(lay["farR"], 0)
    left = pref.get((ag, dg, cg, t, fl))
    right = pref.get((ar, dr, cr, u, fr))
    if left is None or right is None:
        return 0
    total = 0
    for d, v in enumerate(cur):
        if not v:
            continue
        limit = m - d
        for dl in range(limit + 1):
            lv = left[dl + 1] - left[dl]
            if lv:
                total += v * lv * right[limit - dl + 1]
    return total


def enumerate_shell(N: int, m: int, depth: int, shell: tuple):
    j, gsym, confL, confR, t, u, fl, fr = shell
    lay = layout(j, N, depth)
    import itertools

    def fill_fixed(pos, sym, pat):
        if sym == MATCH:
            return [pat]
        if sym == FLIP:
            return [1 - pat]
        return [2, 3]

    fixed_choices = [[]]
    fixed_pos = []
    for pos, sym in zip(lay["guard"], gsym):
        fixed_pos.append(pos)
        letters = fill_fixed(pos, sym, pattern_letter(j, pos))
        fixed_choices = [prev + [let] for prev in fixed_choices for let in letters]
    for pos, sym in zip(lay["resL"], confL):
        fixed_pos.append(pos)
        letters = fill_fixed(pos, sym, 1)
        fixed_choices = [prev + [let] for prev in fixed_choices for let in letters]
    for pos, sym in zip(lay["resR"], confR):
        fixed_pos.append(pos)
        letters = fill_fixed(pos, sym, 0)
        fixed_choices = [prev + [let] for prev in fixed_choices for let in letters]

    def arm_words(positions, edits, flips, pat_of):
        out = []
        for chosen in itertools.combinations(positions, edits):
            for flip_set in itertools.combinations(chosen, flips):
                fset = set(flip_set)
                stars = [i for i in chosen if i not in fset]
                for bits in itertools.product((2, 3), repeat=len(stars)):
                    assign = {i: pattern_letter(j, i) for i in positions}
                    for i in fset:
                        assign[i] = 1 - pattern_letter(j, i)
                    for i, bit in zip(stars, bits):
                        assign[i] = bit
                    out.append(assign)
        if edits == 0:
            return [{i: pattern_letter(j, i) for i in positions}]
        return out

    left_assigns = arm_words(lay["farL"], t, fl, 1)
    right_assigns = arm_words(lay["farR"], u, fr, 0)
    for fixed in fixed_choices:
        base = {pos: let for pos, let in zip(fixed_pos, fixed)}
        for L in left_assigns:
            for R in right_assigns:
                word = pattern(j, N)
                for i, let in base.items():
                    word[i] = let
                for i, let in L.items():
                    word[i] = let
                for i, let in R.items():
                    word[i] = let
                yield tuple(word)


def hamming(a, b) -> int:
    return sum(x != y for x, y in zip(a, b))
