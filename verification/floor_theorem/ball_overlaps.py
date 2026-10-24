from math import comb, factorial
N, q, m = 19, 4, 4

def ball(N, q, m):
    return sum(comb(N, k) * (q - 1) ** k for k in range(m + 1))

def inter(delta, N=N, q=q, m=m):
    """|B_m(x) cap B_m(y)| for d(x,y)=delta, exact."""
    tot = 0
    # tail over the N-delta agreeing positions
    tail = [comb(N - delta, b) * (q - 1) ** b for b in range(0, N - delta + 1)]
    for u in range(delta + 1):
        for v in range(delta - u + 1):
            w = delta - u - v
            mult = factorial(delta) // (factorial(u) * factorial(v) * factorial(w)) * (q - 2) ** w
            if mult == 0:
                continue
            a, ap = v + w, u + w
            lim = m - max(a, ap)
            if lim < 0:
                continue
            tot += mult * sum(tail[b] for b in range(0, min(lim, N - delta) + 1))
    return tot

B = ball(N, q, m)
print("N=%d q=%d m=%d   |B_m| = %d" % (N, q, m, B))
print("\n delta   |B cap B|      as %% of |B|")
I = {}
for d_ in range(0, 10):
    I[d_] = inter(d_)
    print("  %2d   %12d   %8.3f%%" % (d_, I[d_], 100.0 * I[d_] / B))

# sanity: delta=0 must equal |B|
print("\ncontrol delta=0 equals |B|:", I[0] == B)
# two numbers quoted from an independent route, recomputed here as a control rather than taken on report
within3 = sum(comb(N, k) * (q - 1) ** k for k in range(4))
shell4_touching = comb(N - 1, 3) * (q - 1) ** 4
print("control delta=1: within-3 ball = %d, shell-4 touching the differing position = %d, sum = %d, computed = %d"
      % (within3, shell4_touching, within3 + shell4_touching, I[1]))

# union lower bound along the isometric path, r patterns at d(i,j)=|i-j|
for r in (14, 16, 20):
    lb = 0
    for i in range(r):
        sub = sum(I.get(i - j, 0) for j in range(i))
        lb += max(0, B - sub)
    print("\nr = %2d  union lower bound = %d   volume floor ceil(U/|B|) = %d"
          % (r, lb, -(-lb // B)))
