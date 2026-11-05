"""Check F(p,w,s), the closed form the whole sweep funnels through, over its ENTIRE domain.

The certificate's max-load sweep covers all (prefix type, free weight) classes, and every entry it
computes passes through F(p,w,s): the number of free-part words of weight w lying within distance s
of a centre whose free part has weight p. If F is wrong the sweep is wrong everywhere at once, and
the independent explicit-centre check would not catch it, because that check samples 101 centres out
of 590,490 classes.

So F is checked here against a dynamic program built from the explicit definition: walk the L free
positions, at each one decide whether the orbit point carries a zero or one of the three non-zero
letters, and track (non-zeros placed so far, Hamming distance to the centre so far). The centre's
free part is taken as p non-zeros in the FIRST p positions; by the symmetry that permutes free
coordinates only p can matter, and this file does not assume that -- it also re-runs the DP with the
centre's non-zeros in different positions and requires the same answer.

Nothing here calls F's recurrence, its comb/3**(w-o) sum, or its cache.
"""
import itertools, sys
from math import comb

import os
N = int(os.environ.get('CERT_N', 19))
M = int(os.environ.get('CERT_M', 8))

def dp_counts(L, w, centre_free):
    """exact distribution of Hamming distance between the orbit's weight-w words and centre_free"""
    st = {(0, 0): 1}
    for j in range(L):
        cj = centre_free[j]
        nxt = {}
        for (k, d), c in st.items():
            cost = 0 if cj == 0 else 1                      # the point carries 0 here
            nxt[(k, d+cost)] = nxt.get((k, d+cost), 0) + c
            if k < w:
                for letter in (1, 2, 3):                    # the point carries a non-zero here
                    cost = 0 if cj == letter else 1
                    nxt[(k+1, d+cost)] = nxt.get((k+1, d+cost), 0) + c
        st = nxt
    out = {}
    for (k, d), c in st.items():
        if k == w: out[d] = out.get(d, 0) + c
    return out

# the closed form under test, copied verbatim from certify_cell.py so a divergence is F's fault
_f = {}
def F(p, w, s, L):
    if s < 0: return 0
    k = (p, w, s)
    if k in _f: return _f[k]
    tt = 0
    for o in range(0, min(p, w)+1):
        b = comb(p, o)*comb(L-p, w-o)
        if b: tt += b*sum(comb(o, e)*2**(o-e) for e in range(o+1) if w+p-o-e <= s)*3**(w-o)
    _f[k] = tt
    return tt

fail = 0; checked = 0
import os
RS = [int(x) for x in os.environ.get('VF_R', '11,12,13,14,15,16').split(',')]
for R in RS:
    L = N - (R-1)
    _f.clear()
    for p in range(L+1):
        centre = [1]*p + [0]*(L-p)                      # p non-zeros, first p positions
        for w in range(L+1):
            cnt = dp_counts(L, w, centre)
            tot = sum(cnt.values())
            if tot != comb(L, w)*3**w:
                print('  r=%d p=%d w=%d: DP orbit size %d != comb(L,w)*3^w %d'
                      % (R, p, w, tot, comb(L, w)*3**w)); fail += 1
            for s in range(-1, L+2):
                want = sum(c for d, c in cnt.items() if d <= s)
                got = F(p, w, s, L)
                checked += 1
                if want != got:
                    print('  MISMATCH r=%d L=%d p=%d w=%d s=%d: DP %d, F %d'
                          % (R, L, p, w, s, want, got)); fail += 1
    # the symmetry claim itself: moving the centre's non-zeros must not change the distribution
    import random
    rng = random.Random(7)
    for _ in range(25):
        p = rng.randrange(L+1); w = rng.randrange(L+1)
        pos = rng.sample(range(L), p)
        c2 = [0]*L
        for j in pos: c2[j] = rng.choice([1, 2, 3])      # different positions AND letters
        a = dp_counts(L, w, [1]*p + [0]*(L-p))
        b = dp_counts(L, w, c2)
        if a != b:
            print('  SYMMETRY FAILS r=%d p=%d w=%d: free-part distance distribution depends on'
                  ' WHICH positions the centre uses, so free weight is not a sufficient label'
                  % (R, p, w)); fail += 1
print('%d values of F checked against the DP across r=%s; %d disagreements'
      % (checked, ','.join(map(str, RS)), fail))
if fail: sys.exit('F IS WRONG; every certificate built on it is void')
print('F agrees with the explicit dynamic program on its whole domain, and the free-part')
print('distribution depends only on the centre free WEIGHT, as the orbit reduction requires.')
