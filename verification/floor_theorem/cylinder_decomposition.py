"""Is a Hamming ball a union of few cylinders? Settle it exactly, small case first."""
import itertools
from math import comb

def ball(c, m, A):
    N=len(c); S=set()
    for k in range(m+1):
        for pos in itertools.combinations(range(N),k):
            for sub in itertools.product(A,repeat=k):
                if any(sub[t]==c[pos[t]] for t in range(k)): continue
                x=list(c)
                for t in range(k): x[pos[t]]=sub[t]
                S.add("".join(x))
    return S

def cyl_union(c, m, A):
    """Union over every m-subset S of {words agreeing with c off S}."""
    N=len(c); U=set()
    for S in itertools.combinations(range(N),m):
        for sub in itertools.product(A,repeat=m):
            x=list(c)
            for t in range(m): x[S[t]]=sub[t]
            U.add("".join(x))
    return U

A="ACGT"
print(" N  m   |ball|   |union of C(N,m) cylinders|   equal   cylinders   points per cylinder")
for N in (5,6,7):
    for m in (1,2,3):
        if m>=N: continue
        c="A"*N
        B=ball(c,m,A); U=cyl_union(c,m,A)
        print(" %2d %2d  %7d   %24d   %5s   %9d   %d"
              % (N,m,len(B),len(U),B==U,comb(N,m),len(A)**m))

print()
N,q,m=19,4,4
Bsz=sum(comb(N,k)*(q-1)**k for k in range(m+1))
print("Deployed regime N=%d q=%d m=%d:" % (N,q,m))
print("  |ball| = %d" % Bsz)
print("  cylinders needed = C(%d,%d) = %d, each with %d points (q^m)" % (N,m,comb(N,m),q**m))
print("  cylinders x points = %d, which overcounts the ball by %.2fx" % (comb(N,m)*q**m, comb(N,m)*q**m/Bsz))
print("  largest cylinder that fits inside a radius-%d ball has at most m free positions, so at most q^m = %d points," % (m,q**m))
print("  hence at least ceil(|ball|/q^m) = %d cylinders are needed. The C(N,m) decomposition is within %.2fx of that floor."
      % (-(-Bsz//q**m), comb(N,m)/(-(-Bsz//q**m))))
