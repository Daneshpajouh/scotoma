"""Exact staircase private-sphere loads. Standard library only; no repository imports."""
from functools import lru_cache
from fractions import Fraction
from math import comb


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(None)
def local_count(n, m, u, a, e, f):
    """u total free letters, a allowed free, e allowed binary misses, f forced binary misses."""
    same = n-a-e
    if min(n, m, u, a, e, f, same) < 0 or m > n:
        return 0
    ans = 0
    for t in range(min(a, m)+1):
        for b in range(max(0, m-t-same), min(e, m-t)+1):
            k = m-t-b
            cap = 2*t-u-e+b-f
            if cap < 0:
                continue
            good_letters = sum(choose(t, s) for s in range(min(t, cap)+1))
            ans += choose(a,t)*choose(e,b)*choose(same,k)*good_letters*2**(m-t)
    return ans


def profile(c, m):
    N = len(c)
    u = sum(x >= 2 for x in c)
    out = []
    for j in range(N+1):
        forced = {p for p in (j-1, j) if 0 <= p < N}
        n = N-len(forced)
        a = sum(c[p] >= 2 for p in range(N) if p not in forced)
        e = sum(c[p] < 2 and c[p] != int(p<j) for p in range(N) if p not in forced)
        f = sum(c[p] < 2 and c[p] != int(p<j) for p in forced)
        count = local_count(n,m,u,a,e,f)
        out.append(Fraction(count,choose(n,m)*2**m))
    return out


def load(c, m):
    return sum(profile(c,m), Fraction(0))


def L(m, N):
    return Fraction(m,N-2)+Fraction(m,N-m-1)


def nstar(m):
    N = max(m+2,3)
    while L(m,N) > 1:
        N += 1
    return N
