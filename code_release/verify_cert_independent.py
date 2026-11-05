"""Independent check of a floor certificate. Shares nothing with certify_cell.py but the JSON.

Derivation, from scratch. Library s_k = 1^k 0^(N-k), k < r. All s_k agree on positions AP..N-1,
so x = (xp, xt) with d(x, s_k) = dH(xp, P_k) + wt(xt). A support orbit is (prefix type t over
{0,1,2}, tail weight w): '2' means either of the two letters that are neither 0 nor 1.

For a centre c = (cp, ct) the ball count splits because prefix and tail are independent:
  |orbit & B(c,m)| = sum_{a+b<=m} Npre(a) * Ntail(w,b)
  Npre: per position, t_i=0 -> z^[cp_i!=0]; t_i=1 -> z^[cp_i!=1];
        t_i=2 -> (1+z) if cp_i is a non-{0,1} letter, else 2z.
  Ntail: per tail position, ct_j=0 -> (1 + 3yz); ct_j!=0 -> (z + y + 2yz);
        so Ntail(w,.) = [y^w] (1+3yz)^(L-g) (z+y+2yz)^g, depending on ct only through g=wt(ct).
Hence the load depends on cp only through its type in {0,1,2}^AP and on ct only through g.
Mass sits uniformly on each orbit, of size 2^c2 * C(L,w) * 3^w.
The mass is recomputed as the sum of the support numerators, never read from the file's total.
The scan over all centre classes is float64. Every term in it is non-negative, so each load carries
a relative rounding error far below 1e-12. Every class whose float load is within a relative 1e-9 of
the float maximum is then recomputed in exact rationals, and the exact maximum over that window is
the true maximum. The verdict uses the exact mass and the exact maximum only.
"""
import json, sys, itertools
import numpy as np
from math import comb
from fractions import Fraction as Fr

def tail_cum(L, w, g, m):
    """cum[a] = sum_{b <= m-a} Ntail(w,b), as exact ints, for a = 0..m."""
    # poly in y (deg L) and z (deg L)
    P = np.zeros((L+1, L+1), dtype=object); P[0,0] = 1
    for j in range(L):
        Q = np.zeros((L+1, L+1), dtype=object)
        for yy in range(L+1):
            for zz in range(L+1):
                v = P[yy,zz]
                if not v: continue
                if j < L-g:                      # ct_j = 0  -> 1 + 3yz
                    Q[yy,zz] += v
                    if yy+1<=L and zz+1<=L: Q[yy+1,zz+1] += 3*v
                else:                            # ct_j != 0 -> z + y + 2yz
                    if zz+1<=L: Q[yy,zz+1] += v
                    if yy+1<=L: Q[yy+1,zz] += v
                    if yy+1<=L and zz+1<=L: Q[yy+1,zz+1] += 2*v
        P = Q
    row = [int(P[w,b]) for b in range(L+1)]
    cum = []
    for a in range(m+1):
        k = m-a
        cum.append(sum(row[b] for b in range(0, min(k, L)+1)) if k >= 0 else 0)
    return cum

def main(path):
    d = json.load(open(path)); N,m,r,AP,L,Dx = d['N'],d['m'],d['r'],d['AP'],d['L'],d['Dx']
    sup = d['support']
    types = list(itertools.product([0,1,2], repeat=AP))
    nT = len(types)
    # per-orbit prefix polynomial over all centre types, built by DP (numpy float for the scan)
    print(f"cell N={N} m={m} r={r}  AP={AP} L={L}  orbits={len(sup)}  centre classes={nT*(L+1):,}")
    best = (-1.0, None)
    loads = np.zeros((nT, L+1))
    for (ti, w, num) in sup:
        t = types[ti]; c2 = t.count(2)
        size = (2**c2) * comb(L, w) * (3**w)
        # Npre over all centre types: array (nT, AP+1)
        Np = np.zeros((1, AP+1)); Np[0,0] = 1.0
        cur = Np
        for i in range(AP):
            new = np.zeros((cur.shape[0]*3, AP+1))
            for ci in range(3):                      # centre letter class 0,1,other
                blk = np.zeros_like(cur)
                if t[i] in (0,1):
                    dd = 0 if ci == t[i] else 1
                    blk[:, dd:] += cur[:, :AP+1-dd] if dd else cur
                else:
                    if ci == 2: # 1 + z
                        blk += cur; blk[:,1:] += cur[:, :AP]
                    else:       # 2z
                        blk[:,1:] += 2*cur[:, :AP]
                new[ci::3] = blk
            cur = new
        for g in range(L+1):
            cum = np.array(tail_cum(L, w, g, m), dtype=float)
            K = min(m+1, AP+1)
            loads[:, g] += (num/Dx) * (cur[:, :K] @ cum[:K]) / size
    gi = np.unravel_index(np.argmax(loads), loads.shape)
    mx = loads[gi]
    mass_x = Fr(sum(int(num) for (_, _, num) in sup), Dx)          # recomputed, not read
    mass_f = Fr(int(d['mass_numerator_over_Dx']), Dx)
    if mass_x != mass_f:
        print(f"  MASS DISAGREES: support sums to {float(mass_x):.9f}, file states {float(mass_f):.9f}")
    cand = np.argwhere(loads >= mx * (1 - 1e-9))
    mx_x = max(exact_load(types[int(i)], int(g), sup, types, AP, L, m, Dx) for i, g in cand)
    bound_x = mass_x / mx_x
    claim = d['bound_num']/d['bound_den']
    print(f"  mass {float(mass_x):.9f} (support sum)   float max load {mx:.9f} at type#{gi[0]} g={gi[1]}")
    print(f"  exact max load {float(mx_x):.12f} over {len(cand)} candidate class(es) within 1e-9 of it")
    print(f"  independent bound = exact mass / exact max load = {float(bound_x):.9f}")
    print(f"  certificate claims                                {claim:.9f}")
    ok = float(bound_x) >= claim - 1e-7 and mass_x == mass_f
    print(f"  -> {'AGREES' if ok else 'DISAGREES'}; needs > {r-1} for K={r}: "
          f"{'PROVED' if bound_x > r-1 else 'NOT PROVED'}")
    return loads, d, dict(mass=mass_x, maxload=mx_x, bound=bound_x, candidates=len(cand))

def exact_load(tc, g, sup, types, AP, L, m, Dx):
    """Load of the centre class (prefix type tc, tail weight g), in exact rationals."""
    tot = Fr(0)
    for (ti, w, num) in sup:
        t = types[ti]; c2 = t.count(2)
        size = (2**c2) * comb(L, w) * (3**w)
        poly = [1] + [0]*AP                       # Npre over z, exact ints
        for i in range(AP):
            new = [0]*(AP+1)
            for a, v in enumerate(poly):
                if not v: continue
                if t[i] in (0, 1):
                    dd = 0 if tc[i] == t[i] else 1
                    if a+dd <= AP: new[a+dd] += v
                elif tc[i] == 2:                  # 1 + z
                    new[a] += v
                    if a+1 <= AP: new[a+1] += v
                else:                             # 2z
                    if a+1 <= AP: new[a+1] += 2*v
            poly = new
        cum = tail_cum(L, w, g, m)
        K = min(m+1, AP+1)
        tot += Fr(int(num), Dx) * Fr(sum(poly[a]*cum[a] for a in range(K)), size)
    return tot

if __name__ == '__main__':
    main(sys.argv[1])
