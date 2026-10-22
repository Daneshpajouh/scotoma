"""Does a separated private-witness system exist at the claimed sharp length?
Claim: at m=4, q=4, r>=4 the threshold is n>=16, i.e. N=n-1>=15.
Search cyclic templates mod N for N = 13..16 and report the smallest N that works."""
import itertools
def works(N, m, sup, Aset, R):
    Bset = set(sup) - set(Aset)
    def tau(k, s): return 1 + ((k + s - 1) % N)
    def P(k): return [1 if i <= k else 0 for i in range(1, N + 1)]
    def mk(k):
        w = P(k)
        for s in Aset: w[tau(k, s) - 1] = 2
        for s in Bset: w[tau(k, s) - 1] = 3
        return w
    d = lambda x, y: sum(a != b for a, b in zip(x, y))
    W = [mk(k) for k in range(R)]; Pc = [P(k) for k in range(R)]
    if any(d(W[k], Pc[k]) != m for k in range(R)): return None
    priv = min(d(W[k], Pc[j]) for k in range(R) for j in range(R) if j != k)
    if priv <= m: return None
    pair = min(d(W[i], W[j]) for i in range(R) for j in range(i + 1, R))
    if pair <= 2 * m: return None
    return priv, pair

m = 4
print("m = %d, searching cyclic templates; R = N+1 centres (the longest path at that length)" % m)
print("  N   n     result")
for N in range(12, 20):
    R = N + 1
    hit = None
    for sup in itertools.combinations(range(1, N + 1), m):
        for mask in range(1 << m):
            Aset = {sup[i] for i in range(m) if mask >> i & 1}
            r = works(N, m, sup, Aset, R)
            if r: hit = (sup, sorted(Aset), r); break
        if hit: break
    if hit:
        sup, Aset, (priv, pair) = hit
        print("  %2d  %2d    WORKS  support %s  letter-2 at %s  priv %d  pair %d" % (N, N + 1, sup, Aset, priv, pair))
    else:
        print("  %2d  %2d    none found among cyclic templates" % (N, N + 1))
