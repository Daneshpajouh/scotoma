"""Exact checks for the analytic one-free-letter lemma and threshold arithmetic."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
from exact_load import local_count, nstar, L, load
ROOT=Path(__file__).resolve().parent

def falling_ratio(k,n,d):
    if d>k:return F(0)
    return F(comb(n-d,k-d),comb(n,k))

def values(N,m):
    n=N-3;k=m-1;q=F(k,n)
    h=[falling_ratio(k,n,d) for d in range(k+4)]
    v=[]
    for d in range(k+3):
        A=h[d]+(F(d,2)*(h[d-1]-h[d]) if d else 0)
        v.append(A+sum(h[d:],F(0)))
    return h,v,q


def check():
    cells=0;identities=0;u1_words=0
    for m in range(2,61):
      Ns=nstar(m)
      for N in range(Ns-1,Ns+5):
        n=N-3;k=m-1;h,v,q=values(N,m)
        T=F(n+1,n-k+1)
        assert sum(h,F(0))==T
        assert v[:3]==[1+T,T-F(1,2)+q/2,T-1]
        Fpoly=N*N-(3*m+2)*N+m*m+6*m-1
        Qpoly=N*N-3*(m+1)*N+m*m+5*m+2
        assert Fpoly==Qpoly+N+m-3
        assert 2-T-q==F(Fpoly,(N-m-1)*(N-3))
        assert Fpoly>0 and q<F(1,2)
        assert all(v[d+1]<=v[d] for d in range(len(v)-1))
        assert all(v[d+1]<=q*F(d+1,d)*v[d] for d in range(2,len(v)-1))
        assert all(D*v[D-1]<=v[0] for D in range(1,len(v)+1))
        rho=F(m,N-2)
        assert rho*v[0]==L(m,N)
        if m<=20:
            for d in range(m+1):
                direct=F(local_count(N-2,m,1,1,d,0),comb(N-2,m)*2**m)
                formula=rho*(h[d]+(F(d,2)*(h[d-1]-h[d]) if d else 0))
                assert direct==formula
                if d:
                    direct=F(local_count(N-2,m,1,1,d-1,1),comb(N-2,m)*2**m)
                    assert direct==rho*h[d-1]/2
                identities+=1
        cells+=1
    for m,N in [(2,6),(3,9),(6,16),(9,24),(19,51)]:
        for a in range(N+1):
            b=[1]*a+[0]*(N-a)
            for p in range(N):
                c=b[:];c[p]=2
                assert load(c,m)<=L(m,N)
                u1_words+=1
    # Controls attack the new proof's necessary arithmetic and its nonzero term.
    h,v,q=values(22,10)
    assert not (2*v[1]<=v[0])
    N,m,d=16,6,1
    h,_,_=values(N,m)
    actual=F(local_count(N-2,m,1,1,d,0),comb(N-2,m)*2**m)
    wrong=F(m,N-2)*h[d]
    assert actual!=wrong
    result={'status':'VERIFIED-COMPUTATIONALLY','checked_parameter_cells':cells,'local_formula_checks':identities,'one_free_letter_centres_checked':u1_words,'controls':[{'control':'drop F(N,m)>=0 from valley comparison','outcome':'REFUTED_AS_REQUIRED','N':22,'m':10},{'control':'omit the one-uncovered-disagreement term','outcome':'REFUTED_AS_REQUIRED','N':16,'m':6,'actual':str(actual),'wrong':str(wrong)}]}
    (ROOT/'companion_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':check()
