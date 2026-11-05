import itertools
N=19; R=N+1  # centres k = 0..N, 20 patterns at N = 19, the largest library a 20-mer can have;
            # was R=19, which skipped the top case (see verify_construction.py).
def tau(k,s): return 1+((k+s-1)%19)
def P(k): return [1 if i<=k else 0 for i in range(1,N+1)]
def mk(k,Aset,Bset):
    w=P(k)
    for s in Aset: w[tau(k,s)-1]=2
    for s in Bset: w[tau(k,s)-1]=3
    return w
def d(x,y): return sum(a!=b for a,b in zip(x,y))

def test(Aset,Bset,m):
    Ws=[mk(k,Aset,Bset) for k in range(R)]
    Ps=[P(k) for k in range(R)]
    own=max(d(Ws[k],Ps[k]) for k in range(R))
    if own>m: return None
    priv=min(d(Ws[k],Ps[j]) for k in range(R) for j in range(R) if j!=k)
    pair=min(d(Ws[i],Ws[j]) for i in range(R) for j in range(i+1,R))
    return own,priv,pair,(priv>m and pair>2*m)

print("Searching size-m supports from the full coordinate set, letters 2 or 3 per position.")
print(" m  support            A        own priv pair  floor holds")
missing = []
for m in range(1,6):
    found=None
    # try subsets of the published S first, then any subset of 1..19
    pools=[list(itertools.combinations([4,6,9,12,16],m))] if m<=5 else []
    pools.append(list(itertools.combinations(range(1,20),m)))
    for pool in pools:
        for sup in pool:
            for amask in range(1<<m):
                Aset={sup[i] for i in range(m) if amask>>i&1}
                Bset=set(sup)-Aset
                r=test(Aset,Bset,m)
                if r and r[3]:
                    found=(sup,sorted(Aset),r); break
            if found: break
        if found: break
    if found:
        sup,Aset,r=found
        print(" %d  %-18s %-8s %3d %4d %4d  YES" % (m,str(sup),str(Aset),r[0],r[1],r[2]))
    else:
        print(" %d  none found" % m)
        missing.append(m)

# A verdict and an exit code. Until 2026-10-23 every script in this folder printed its checks
# and exited 0 whatever they said, so a reader running one in CI saw a pass even when the check
# had failed. A verifier that cannot fail verifies nothing.
print("\nVERDICT:", "PASS" if not missing else "FAIL",
      ("(a family at every budget 1 to 5)" if not missing
       else "(no family found at budget(s) %s)" % missing))
raise SystemExit(0 if not missing else 1)
