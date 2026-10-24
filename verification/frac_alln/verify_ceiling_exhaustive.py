"""Exhaustive, independent check of a released ceiling certificate (data/floor_certificates schema).

verify_ceiling_certificate.py checks a certificate on the r library centres and 90 sampled centres.
That shows consistency, not the ceiling. The exhaustive sweep that sets the ceiling lives in the
solver, certify_cell.py, and runs through the solver's own load tables. This file repeats that
sweep over EVERY centre class with a different algorithm, so a solver bug has to occur twice, the
same way, to survive.

It never calls F(p,w,s), load_table or anything from the solver. Distances are counted by a
dynamic program over explicit strings, as in verify_ceiling_certificate.py. All arithmetic is exact.

What a centre class is. Coordinates 0..AP-1 (AP = r-1) are the active prefix, where the library
patterns 1^k 0^(N-k), k < r, differ; the other L = N-AP coordinates are free, and every library
pattern holds 0 there. A support block is every word with a given active prefix over {0,1,X} (X =
either letter 2 or 3) and a given free weight w (w letters from {1,2,3} placed in every way), its
mass spread equally over its words. Such a measure is unchanged by swapping 2 and 3 at one active
coordinate, by permuting the free coordinates and by relabelling 1,2,3 on them. So a centre's load
depends only on its prefix read over {0,1,X} and its free weight p: 3^AP x (L+1) classes. Check 4
tests that reduction on explicit strings rather than assuming it.

Checks, each exits 1 on failure:
  1. every block mass is a nonnegative rational, and the masses sum to exact_total_mass;
  2. every block lies in the target: its distance to the nearest library pattern is at most m;
  3. EVERY class carries load at most 1, the maximum is reported with where it occurs, and it
     must equal exact_max_ball_mass;
  4. on sampled classes, a random explicit centre (2 and 3 mixed, free letters and positions
     random) carries exactly the class representative's load;
  5. total mass > r-1, so integrality gives K >= r.
Usage: python3 -I verify_ceiling_exhaustive.py cert_m9_r9.json [more certificates]
"""
import itertools, json, random, sys
from fractions import Fraction as Fr
from math import comb


def active_dist(prefix, u):
    """distance distribution over the active prefix: block prefix over 0/1/X against explicit u"""
    dist = {0: 1}
    for t, c in zip(prefix, u):
        letters = (2, 3) if t == 'X' else (int(t),)
        nxt = {}
        for d, n in dist.items():
            for a in letters:
                e = d + (a != c)
                nxt[e] = nxt.get(e, 0) + n
        dist = nxt
    return dist


def free_dist(w, v):
    """distance distribution: words of weight w over {0,1,2,3} on the free part against explicit v"""
    st = {(0, 0): 1}
    for c in v:
        nxt = {}
        for (k, d), n in st.items():
            e = d + (c != 0)
            nxt[(k, e)] = nxt.get((k, e), 0) + n
            if k < w:
                for a in (1, 2, 3):
                    e = d + (a != c)
                    nxt[(k + 1, e)] = nxt.get((k + 1, e), 0) + n
        st = nxt
    out = {}
    for (k, d), n in st.items():
        if k == w:
            out[d] = out.get(d, 0) + n
    return out


def cumulative(dist, top):
    """cum[s] = number of words at distance <= s, for s = 0..top"""
    cum, run = [], 0
    for s in range(top + 1):
        run += dist.get(s, 0)
        cum.append(run)
    return cum


def check(path, samples=40, seed=20261023):
    C = json.load(open(path))
    N, m, r, L = C['N'], C['m'], C['r'], C['free_length']
    AP = r - 1
    assert AP + L == N, 'free_length does not match N and r'
    lib = ['1' * k + '0' * (AP - k) for k in range(r)]
    blocks = [(d['prefix'], int(d['free_weight']), Fr(d['mass'])) for d in C['dual']]
    fails = []

    # 1. masses
    neg = [b for b in blocks if b[2] < 0]
    if neg:
        fails.append('%d negative masses' % len(neg))
    total = sum(b[2] for b in blocks)
    if total != Fr(C['exact_total_mass']):
        fails.append('masses sum to %s, certificate says %s' % (total, C['exact_total_mass']))

    # 2. target membership: library patterns hold 0 on the free part, so distance = prefix + w
    def pdist(y, t):
        return sum(1 for a, b in zip(y, t) if a != b)   # X differs from both 0 and 1
    far = [(y, w) for y, w, _ in blocks if min(pdist(y, t) for t in lib) + w > m]
    if far:
        fails.append('%d blocks outside the target, first %s' % (len(far), far[0]))

    # block sizes, counted by the same DPs at a centre of all zeros
    size = {}
    for y, w, _ in blocks:
        size[(y, w)] = sum(active_dist(y, [0] * AP).values()) * sum(free_dist(w, [0] * L).values())
        assert size[(y, w)] == 2 ** y.count('X') * comb(L, w) * 3 ** w

    # 3. every class: prefix over 0/1/X (X read as letter 2), free part 1^p 0^(L-p)
    ws = sorted({w for _, w, _ in blocks})
    fcum = {(w, p): cumulative(free_dist(w, [1] * p + [0] * (L - p)), m)
            for w in ws for p in range(L + 1)}
    prefixes = sorted({y for y, _, _ in blocks})
    worst, worst_at, n_classes, over = Fr(0), None, 0, 0
    for cp in itertools.product('01X', repeat=AP):
        u = [2 if c == 'X' else int(c) for c in cp]
        act = {y: active_dist(y, u) for y in prefixes}
        for p in range(L + 1):
            tot = Fr(0)
            for y, w, x in blocks:
                cf = fcum[(w, p)]
                inside = sum(n * cf[m - d1] for d1, n in act[y].items() if d1 <= m)
                if inside:
                    tot += x * Fr(inside, size[(y, w)])
            n_classes += 1
            if tot > 1:
                over += 1
            if tot > worst:
                worst, worst_at = tot, (''.join(cp), p)
    if n_classes != 3 ** AP * (L + 1):
        fails.append('swept %d classes, expected %d' % (n_classes, 3 ** AP * (L + 1)))
    if over:
        fails.append('%d classes carry load above 1' % over)
    if worst != Fr(C['exact_max_ball_mass']):
        fails.append('maximum load %s, certificate says %s' % (worst, C['exact_max_ball_mass']))

    # 4. the class reduction, on explicit strings
    rng = random.Random(seed)

    def explicit_load(u):
        tot = Fr(0)
        for y, w, x in blocks:
            a = active_dist(y, u[:AP])
            fc = cumulative(free_dist(w, u[AP:]), m)
            inside = sum(n * fc[m - d1] for d1, n in a.items() if d1 <= m)
            tot += x * Fr(inside, size[(y, w)])
        return tot

    bad = 0
    for _ in range(samples):
        cp = [rng.choice('01X') for _ in range(AP)]
        p = rng.randrange(L + 1)
        rep = [2 if c == 'X' else int(c) for c in cp] + [1] * p + [0] * (L - p)
        alt = [rng.choice((2, 3)) if c == 'X' else int(c) for c in cp]
        free = [0] * L
        for j in rng.sample(range(L), p):
            free[j] = rng.choice((1, 2, 3))
        if explicit_load(rep) != explicit_load(alt + free):
            bad += 1
    if bad:
        fails.append('%d of %d sampled classes: explicit centre load differs from the class load'
                     % (bad, samples))

    # 5. rounding
    if not total > r - 1:
        fails.append('total mass %s is not above r-1 = %d' % (total, r - 1))

    print('%s: N %d m %d r %d | %d blocks, mass %s = %.9f > %d %s | %d of %d classes swept, '
          'max load %s at prefix %s free weight %d (certificate: %s at %s) | reduction %d/%d equal'
          % (path.split('/')[-1], N, m, r, len(blocks), total, float(total), r - 1, total > r - 1,
             n_classes, 3 ** AP * (L + 1), worst, worst_at[0], worst_at[1],
             C['exact_max_ball_mass'], C.get('max_ball_at'), samples - bad, samples))
    for f in fails:
        print('  FAIL: ' + f)
    return not fails


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    ok = all([check(p) for p in sys.argv[1:]])
    print('ALL PASS' if ok else 'FAILED')
    sys.exit(0 if ok else 1)
