"""Exhaustive, independent check of a released ceiling certificate (data/floor_certificates schema).

verify_ceiling_certificate.py checks a certificate on the r library centres and 90 sampled centres.
That shows consistency, not the ceiling. The exhaustive sweep that sets the ceiling lives in the
solver, certify_cell.py, and runs through the solver's own load tables. This file repeats that
sweep over EVERY centre class without the solver's load tables. Agreement is an implementation
cross-check: both programs still rest on the same definitions, set out below.

It never calls F(p,w,s), load_table or anything from the solver. For each block it counts the
block's words by their distance to a given centre with a dynamic program over the coordinates of
that centre; it does not list the words one by one. All arithmetic is exact.

What a centre class is. Coordinates 0..AP-1 (AP = r-1) are the active prefix, where the library
patterns 1^k 0^(N-k), k < r, differ; the other L = N-AP coordinates are free, and every library
pattern holds 0 there. A support block is every word with a given active prefix over {0,1,X} (X =
either letter 2 or 3) and a given free weight w (w letters from {1,2,3} placed in every way), its
mass spread equally over its words. Such a measure is unchanged by swapping 2 and 3 at one active
coordinate, by permuting the free coordinates and by relabelling 1,2,3 independently at each
free coordinate. So a centre's load depends only on its prefix read over {0,1,X} and its free
weight p: 3^AP x (L+1) classes. Check 4 is a supplementary test of the class representatives; the reduction itself rests on the symmetry.

Checks, each exits 1 on failure:
  0. the file's format: integer N, m, r, free_length that fit together, every prefix of length r-1
     over 0, 1, X, every free weight in 0..free_length, no block listed twice, every mass and the
     two stated totals an exact rational (a JSON string or integer; decimal numbers are decoded
     exactly, never through a float; booleans are refused), max_ball_at, if present, a valid class
     (null is refused; the field may be left out), and no JSON object repeating a key;
  1. every block mass is a nonnegative rational, and the masses sum to exact_total_mass;
  2. every block lies in the target: its distance to the nearest library pattern is at most m;
  3. EVERY class carries load at most 1, the maximum is reported with where it occurs, and it
     must equal exact_max_ball_mass; the class named in max_ball_at, if any, must be a valid class
     carrying that maximum (ties are allowed). centre_classes_verified is printed, not tested:
     at m=12 the solver swept 136 binary-prefix classes by design;
  4. on sampled classes, a random centre of the class (2 and 3 mixed, free letters and positions
     random) carries exactly the load the sweep recorded for that class;
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


def no_repeated_keys(pairs):
    keys = [k for k, _ in pairs]
    if len(keys) != len(set(keys)):
        raise ValueError('a JSON object repeats a key: %s' % sorted({k for k in keys if keys.count(k) > 1}))
    return dict(pairs)


def exact(v):
    """a rational from a JSON string ("a/b", "1.5") or integer, never through a float; else None"""
    if type(v) is int or type(v) is Fr:
        return Fr(v)
    if type(v) is str:
        try:
            return Fr(v)
        except (ValueError, ZeroDivisionError):
            return None
    return None                           # bool, None, list, dict


def check(path, samples=40, seed=20261023):
    # Decimal tokens are decoded exactly (never through a binary float), and a repeated key is an error.
    C = json.load(open(path), parse_float=Fr, object_pairs_hook=no_repeated_keys)

    # 0. format: these are checks on the file, not mathematical assumptions
    bad_format = []
    if type(C) is not dict:
        print('%s: FORMAT' % path.split('/')[-1])
        print('  FAIL: the file is not a JSON object')
        return False
    N, m, r, L = (C.get(k) for k in ('N', 'm', 'r', 'free_length'))
    if type(C.get('dual')) is not list or not all(type(d) is dict for d in C['dual']):
        bad_format.append('dual must be a list of objects')
    elif not all(type(v) is int for v in (N, m, r, L)):
        bad_format.append('N, m, r and free_length must be integers')
    elif not (1 <= r <= N + 1 and 0 <= m <= N and (r - 1) + L == N):
        bad_format.append('need 1 <= r <= N+1, 0 <= m <= N and (r-1) + free_length = N')
    else:
        seen = set()
        for d in C.get('dual', []):
            y, w = d.get('prefix'), d.get('free_weight')
            if not (isinstance(y, str) and len(y) == r - 1 and set(y) <= set('01X')):
                bad_format.append('prefix %r is not a string of length r-1 over 0, 1, X' % (y,))
            if not (type(w) is int and 0 <= w <= L):
                bad_format.append('free weight %r is not an integer in 0..free_length' % (w,))
            if (y, w) in seen:
                bad_format.append('block (%s, %s) is listed twice' % (y, w))
            seen.add((y, w))
            if exact(d.get('mass')) is None:
                bad_format.append('mass %r is not an exact rational' % (d.get('mass'),))
        if not seen:
            bad_format.append('no dual blocks')
        for k in ('exact_total_mass', 'exact_max_ball_mass'):
            if exact(C.get(k)) is None:
                bad_format.append('%s %r is not an exact rational' % (k, C.get(k)))
        named = C.get('max_ball_at')
        if 'max_ball_at' in C and not (
                isinstance(named, list) and len(named) == 2 and isinstance(named[0], str)
                and len(named[0]) == r - 1 and set(named[0]) <= set('01X')
                and type(named[1]) is int and 0 <= named[1] <= L):
            bad_format.append('max_ball_at %r is not [prefix of length r-1 over 0, 1, X, '
                              'integer free weight in 0..free_length]' % (named,))
    if bad_format:
        print('%s: FORMAT' % path.split('/')[-1])
        for f in bad_format:
            print('  FAIL: ' + f)
        return False

    AP = r - 1
    lib = ['1' * k + '0' * (AP - k) for k in range(r)]
    blocks = [(d['prefix'], d['free_weight'], exact(d['mass'])) for d in C['dual']]
    fails = []

    # 1. masses
    neg = [b for b in blocks if b[2] < 0]
    if neg:
        fails.append('%d negative masses' % len(neg))
    total = sum(b[2] for b in blocks)
    if total != exact(C['exact_total_mass']):
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
    worst, worst_at, n_classes, over, loads = Fr(0), None, 0, 0, {}
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
            loads[(''.join(cp), p)] = tot
            if tot > 1:
                over += 1
            if tot > worst:
                worst, worst_at = tot, (''.join(cp), p)
    if n_classes != 3 ** AP * (L + 1):
        fails.append('swept %d classes, expected %d' % (n_classes, 3 ** AP * (L + 1)))
    if over:
        fails.append('%d classes carry load above 1' % over)
    if worst != exact(C['exact_max_ball_mass']):
        fails.append('maximum load %s, certificate says %s' % (worst, C['exact_max_ball_mass']))
    named = C.get('max_ball_at')
    if 'max_ball_at' in C:
        key = (named[0], named[1]) if isinstance(named, list) and len(named) == 2 else None
        if key not in loads:
            fails.append('max_ball_at %r is not a valid class' % (named,))
        elif loads[key] != worst:
            fails.append('max_ball_at %r carries %s, not the maximum %s' % (named, loads[key], worst))

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
        swept = loads[(''.join(cp), p)]
        if explicit_load(rep) != swept or explicit_load(alt + free) != swept:
            bad += 1
    if bad:
        fails.append('%d of %d sampled classes: explicit centre load differs from the class load'
                     % (bad, samples))

    # 5. rounding
    if not total > r - 1:
        fails.append('total mass %s is not above r-1 = %d' % (total, r - 1))

    print('%s: N %d m %d r %d | %d blocks, mass %s = %.9f > %d %s | %d of %d classes swept, '
          'max load %s at prefix %s free weight %d (certificate: %s at %s; declares %s classes) | '
          'reduction %d/%d equal'
          % (path.split('/')[-1], N, m, r, len(blocks), total, float(total), r - 1, total > r - 1,
             n_classes, 3 ** AP * (L + 1), worst, worst_at[0], worst_at[1],
             C['exact_max_ball_mass'], C.get('max_ball_at'), C.get('centre_classes_verified'),
             samples - bad, samples))
    for f in fails:
        print('  FAIL: ' + f)
    return not fails


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    results = []
    for p in sys.argv[1:]:
        try:
            results.append(check(p))
        except Exception as exc:          # a file this program cannot read is a failure, not a crash
            print('%s: FAIL: could not be checked: %s: %s' % (p.split('/')[-1], type(exc).__name__, exc))
            results.append(False)
    ok = all(results)
    print('ALL PASS' if ok else 'FAILED')
    sys.exit(0 if ok else 1)
