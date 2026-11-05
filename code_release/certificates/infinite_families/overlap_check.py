"""Independent coordinate-polynomial recount of pair minima."""
from collections import Counter
from functools import lru_cache
from itertools import product, combinations_with_replacement
from math import comb
import numpy as np
from exact_tools import pair_minimum

@lru_cache(maxsize=15000)
def dense_rectangle(C,a,b,s,h,w):
    # Direct four-letter polynomial multiplication with exact signed int64.
    assert 4**C<2**63
    arr=np.zeros((w+1,C+1,C+1),dtype=np.int64)
    arr[0,0,0]=1
    av=(1,)*a+(0,)*(C-a)
    bv=(1,)*h+(2,)*(s-h)+(0,)*(a-s)+(1,)*(b-s)+(0,)*(C-a-b+s)
    assert len(av)==len(bv)==C
    for x,y in zip(av,bv):
        nxt=np.zeros_like(arr)
        for letter in range(4):
            dw,d1,d2=int(letter!=0),int(letter!=x),int(letter!=y)
            nxt[dw:,d1:,d2:]+=arr[:w+1-dw,:C+1-d1,:C+1-d2]
        arr=nxt
    hist=arr[w]
    assert int(hist.sum())==comb(C,w)*3**w
    return hist.cumsum(axis=0).cumsum(axis=1)


@lru_cache(None)
def concrete_prefix_histograms(a,b,p):
    # Enumerate all resolutions of B, rather than identifying equal-X counts.
    histograms=set()
    positions=[i for i,x in enumerate(b) if x==2]
    for values in product((2,3),repeat=len(positions)):
        bv=list(b)
        for i,value in zip(positions,values):
            bv[i]=value
        hist=Counter()
        for pv in product(*[(2,3) if x==2 else (x,) for x in p]):
            d1=sum(x!=y for x,y in zip(a,pv))
            d2=sum(x!=y for x,y in zip(bv,pv))
            hist[d1,d2]+=1
        histograms.add(tuple(sorted(hist.items())))
    return tuple(histograms)


@lru_cache(None)
def dense_minimum(C,m,ca,cb,point):
    a,aw=ca;b,bw=cb;p,w=point
    histograms=concrete_prefix_histograms(a,b,p)
    best=None
    for s in range(max(0,aw+bw-C),min(aw,bw)+1):
        for h in range(s+1):
            rect=dense_rectangle(C,aw,bw,s,h,w)
            for hist in histograms:
                value=sum(count*int(rect[min(C,m-d1),min(C,m-d2)])
                          for (d1,d2),count in hist if d1<=m and d2<=m)
                best=value if best is None else min(best,value)
    return best


def literal_calibration():
    tests=0
    for N,r in [(3,2),(3,3)]:
        C=N-r+1;L=r-1
        words=list(product(range(4),repeat=N))
        classes={}
        for i,w in enumerate(words):
            cls=(tuple(min(2,x) for x in w[:L]),sum(x!=0 for x in w[L:]))
            classes.setdefault(cls,[]).append(i)
        dist=[[sum(x!=y for x,y in zip(a,b)) for b in words] for a in words]
        masks={p:sum(1<<i for i in ids) for p,ids in classes.items()}
        for m in range(N+1):
            balls=[sum(1<<j for j,d in enumerate(row) if d<=m) for row in dist]
            for ca,cb in combinations_with_replacement(classes,2):
                for p in classes:
                    literal=min((balls[a]&balls[b]&masks[p]).bit_count()
                                for a in classes[ca] for b in classes[cb])
                    sparse=pair_minimum(C,m,ca,cb,p)
                    dense=dense_minimum(C,m,ca,cb,p)
                    assert literal==sparse==dense,(N,r,m,ca,cb,p,literal,sparse,dense)
                    tests+=1
    return tests


