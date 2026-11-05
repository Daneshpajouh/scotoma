"""My own exhaustive check of the two-deletion minimality counterexample, 2026-10-29.
Enumerates the whole space; no heuristics, no sampling. Calibrates on b=1, where
minimality is an established result, then exhibits the b=2 failure."""
from itertools import combinations, product
A = "ACGT"
def dele(g, k):
    return sorted({"".join(c for i, c in enumerate(g) if i not in idx)
                   for idx in combinations(range(len(g)), k)})
def ball(w, m):
    out = {w}
    for r in range(1, m + 1):
        for pos in combinations(range(len(w)), r):
            for subs in product(A, repeat=r):
                t = list(w); ok = True
                for p, c in zip(pos, subs):
                    if t[p] == c: ok = False; break
                    t[p] = c
                if ok: out.add("".join(t))
    return out
def min_cover(U, L, m):
    U = set(U)
    cands = [(c, ball(c, m) & U) for c in ("".join(p) for p in product(A, repeat=L))]
    cands = [(c, s) for c, s in cands if s]
    for k in range(1, len(U) + 1):
        for combo in combinations(range(len(cands)), k):
            cov = set()
            for i in combo: cov |= cands[i][1]
            if cov == U: return k, [cands[i][0] for i in combo]
    return None, None
if __name__ == "__main__":
    g, m, fails = "AACAA", 1, 0
    for b in (1, 2):
        lib = dele(g, b); L = len(lib[0])
        U = set().union(*(ball(w, m) for w in lib))
        k, wit = min_cover(U, L, m)
        print(f"{g} b={b} m={m}: {len(lib)} words, K={k}, minimality={'holds' if len(lib)==k else 'FAILS'}")
        if b == 1: assert len(lib) == k, "b=1 calibration failed"
        if b == 2:
            assert len(lib) != k, "expected the b=2 counterexample to fail minimality"
            print(f"  witness: {wit}  (one library word's ball is contained in the others')")
            fails += 1
    assert fails == 1
    print("calibration passed and the counterexample reproduces")
