# EXPECTED OUTCOMES, written down so the run can be judged instead of read.
# The m = 4 family holds at every (N, r) row below. The m sweep holds to m = 5 and then FAILS
# at m = 6 and 7: that failure is the point of the sweep, not a fault, so a verdict here has to
# compare the pattern rather than demand every row say YES.
# Note what the failure does and does not show. The sweep extends the support with coordinates
# 17 and 18; m = 6 and 7 fail FOR THAT EXTENSION. It is not a statement that no support works
# at those budgets, and the README no longer says otherwise.
EXPECTED_ROWS = True            # every (N, r) row of the m = 4 table
EXPECTED_SWEEP = {1: True, 2: True, 3: True, 4: True, 5: True, 6: False, 7: False}
observed_rows, observed_sweep = [], {}

import itertools
def run(N,R,sup,Aset,m,mod=19):
    Bset=set(sup)-set(Aset)
    def tau(k,s): return 1+((k+s-1)%mod)
    def P(k): return [1 if i<=k else 0 for i in range(1,N+1)]
    def mk(k):
        w=P(k)
        for s in Aset: w[tau(k,s)-1]=2
        for s in Bset: w[tau(k,s)-1]=3
        return w
    d=lambda x,y: sum(a!=b for a,b in zip(x,y))
    Ws=[mk(k) for k in range(R)]; Ps=[P(k) for k in range(R)]
    own=max(d(Ws[k],Ps[k]) for k in range(R))
    priv=min(d(Ws[k],Ps[j]) for k in range(R) for j in range(R) if j!=k)
    pair=min(d(Ws[i],Ws[j]) for i in range(R) for j in range(i+1,R))
    return own,priv,pair,(own<=m and priv>m and pair>2*m)

print("m = 4, support (4,6,9,12), A = {4,9}")
print(" N   r   own priv pair  holds")
for N,R in ((19,19),(19,20),(19,14),(20,20),(25,20),(25,26),(30,31)):
    if R>N+1: continue
    try:
        o,p,q,ok=run(N,R,(4,6,9,12),{4,9},4)
        print(" %2d  %2d   %3d %4d %4d   %s" % (N,R,o,p,q,"YES" if ok else "NO"))
        observed_rows.append(ok)
    except Exception as e:
        print(" %2d  %2d   error %s" % (N,R,e))
print()
print("m sweep at N=19,r=19 with support = first m of (4,6,9,12,16):")
sup5=(4,6,9,12,16)
for m in range(1,8):
    s=tuple(sorted((sup5+tuple(range(17,30)))[:m]))
    try:
        o,p,q,ok=run(19,19,s,{x for i,x in enumerate(s) if i%2==0},m)
        print("  m=%d support=%s own=%d priv=%d pair=%d  %s" % (m,s,o,p,q,"YES" if ok else "NO"))
        observed_sweep[m]=ok
    except Exception as e:
        print("  m=%d error %s" % (m,e))

import sys
rows_ok  = bool(observed_rows) and all(r is EXPECTED_ROWS for r in observed_rows)
sweep_ok = observed_sweep == EXPECTED_SWEEP
print("\nVERDICT:", "PASS" if (rows_ok and sweep_ok) else "FAIL")
if not rows_ok:
    print(f"   m = 4 table: {sum(1 for r in observed_rows if not r)} of {len(observed_rows)} rows did not hold")
if not sweep_ok:
    print(f"   m sweep: observed {observed_sweep}, expected {EXPECTED_SWEEP}")
sys.exit(0 if (rows_ok and sweep_ok) else 1)
