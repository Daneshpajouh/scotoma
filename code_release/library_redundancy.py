"""Is any library ball redundant? If one were, K <= r-1 would follow immediately.

For the staircase library P_k = 1^k 0^(N-k), k = 0..r-1, and a point z whose active-prefix type is
t over {0,1,X} with free weight w, the distance to P_k is a function of (t, w) alone:

    d(z, P_k) = #{i < k : t_i = 0} + #{i >= k : t_i = 1} + #{i : t_i = X} + w
              = A(k) + X + w,

because an X position disagrees with both library letters and every free position of P_k is 0. A(k)
is a +-1 walk in k: step +1 where t_i = 0, -1 where t_i = 1, 0 where t_i = X.

P_j is NOT redundant exactly when some (t, w) has A(j) + X + w <= m and A(k) + X + w > m for every
k != j. Setting c = X + w, that needs A(j) <= m - c and A(k) >= m - c + 1, and since A moves by at
most one per step the two together force c = m and a UNIQUE minimum of A at j.

CLAIMED WITNESS: t = X^{x1} 1^a 0^b X^{x2} with x1 + a + b + x2 = r-1 and x1 + x2 + w = m. Then A
falls to 0 at k = x1 + a and rises on both sides, so the minimum is unique there, and j = x1 + a can
be placed anywhere. This file checks that claim by computing A(k) directly from the definition for
every r, rather than trusting the derivation.
"""
import itertools, sys

N, M = 19, 8

def dists(t, w):
    """d(z,P_k) for k = 0..len(t), straight from the definition. No walk, no recurrence."""
    AP = len(t)
    X = sum(1 for c in t if c == 'X')
    out = []
    for k in range(AP+1):
        a = sum(1 for i in range(AP) if i < k and t[i] == '0')
        b = sum(1 for i in range(AP) if i >= k and t[i] == '1')
        out.append(a + b + X + w)
    return out

print('r   L   every ball has an exclusive point?   witness for j=0 / j=(r-1)/2 / j=r-1')
allgood = True
for r in range(3, N+1):
    AP, L = r-1, N-(r-1)
    ok = {}
    for j in range(r):
        found = None
        # the claimed family, searched over its free parameters
        for x1 in range(0, min(AP, M)+1):
            for x2 in range(0, min(AP-x1, M-x1)+1):
                w = M - x1 - x2
                if w < 0 or w > L: continue
                for a in range(0, AP-x1-x2+1):
                    b = AP - x1 - x2 - a
                    if x1 + a != j: continue
                    t = 'X'*x1 + '1'*a + '0'*b + 'X'*x2
                    d = dists(t, w)
                    inside = [k for k, v in enumerate(d) if v <= M]
                    if inside == [j]:
                        found = (t, w); break
                if found: break
            if found: break
        ok[j] = found
    good = all(ok[j] for j in range(r))
    allgood &= good
    reps = [ok[0], ok[(r-1)//2], ok[r-1]]
    print('%2d  %2d   %-34s %s' % (r, L, 'YES, none redundant' if good else
          'NO: balls %s have none' % [j for j in ok if not ok[j]],
          ' / '.join('%s w=%d' % v if v else 'none' for v in reps)))

print()
if allgood:
    print('No library ball is redundant at any run count from 3 to %d.' % N)
    print('CONSEQUENCE, and it is a negative: the cheapest route to an upper bound on K is closed.')
    print('Dropping a ball never covers U, so the boundary above cannot be found this way; an')
    print('(r-1)-cover, if one exists, must use centres that are not library patterns.')
else:
    sys.exit('some ball IS redundant; that r is an immediate K <= r-1 and must be checked by hand')
