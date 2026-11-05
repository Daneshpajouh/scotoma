"""Literal checks of the finite-grid reduction, independent of local load bounds."""
from fractions import Fraction as F
from itertools import product,combinations
from pathlib import Path
import json
import verify_w
ROOT=Path(__file__).resolve().parent


def vertices(c):
    z=0;o=c.count(1)
    out=[(z,o)]
    for a in c:
        z+=a==0;o-=a==1
        out.append((z,o))
    return out


def reward(c,i,table,K):
    z,o=vertices(c)[i]
    if z+o>=K:return F(0)
    prev=c[i-1] if i else 1
    nxt=c[i] if i<len(c) else 0
    s=int(prev==2)+int(nxt==2)
    f=int(prev==0)+int(nxt==1)
    return table[s,z+o-f,f]


def table_for(K):
    return {(s,e,f):F((7*(s+2)*(e+3)+11*f)%23+1,37) for s in range(3) for f in range(3-s) for e in range(K+1)}


def check():
    dp_cases=0;embeddings=0;level_checks=0
    for K,u in [(2,2),(3,2),(2,3)]:
        tab=table_for(K);N=2*K+u;best=F(-1)
        for stars in combinations(range(N),u):
            rem=[i for i in range(N) if i not in stars]
            for ones in combinations(rem,K):
                c=tuple(2 if i in stars else 1 if i in ones else 0 for i in range(N))
                v=sum((reward(c,j,tab,K) for j in range(N+1)),F(0))
                best=max(best,v)
                dp_cases+=1
        originalK=verify_w.K;verify_w.K=K
        got=verify_w.path_value(u,tab)
        verify_w.K=originalK
        assert got==best,(K,u,got,best)
    K=3;tab=table_for(K)
    for N in [5,6,7]:
      for c in product(range(3),repeat=N):
        walk=vertices(c);u=c.count(2)
        for h in range(N+1):
            assert sum(z+o<=h for z,o in walk)<=2*h+u+1
            level_checks+=1
        active=[i for i,(z,o) in enumerate(walk) if z+o<K]
        if not active:continue
        first,last=active[0],active[-1]
        zi,oi=walk[first];zj,oj=walk[last]
        block=c[first:last]
        extra=u-block.count(2)
        prefix=(2,)*extra+(0,)*zi+(1,)*(K-oi)
        suffix=(0,)*(K-zj)+(1,)*oj
        extended=prefix+block+suffix
        assert extended.count(0)==K and extended.count(1)==K and extended.count(2)==u
        for j in active:
            assert reward(c,j,tab,K)==reward(extended,len(prefix)+j-first,tab,K)
        assert sum((reward(c,j,tab,K) for j in active),F(0))<=sum((reward(extended,j,tab,K) for j in range(len(extended)+1)),F(0))
        embeddings+=1
    # Removing u from the level-count bound must fail on a free plateau.
    c=(2,)*7
    assert sum(z+o==0 for z,o in vertices(c))>1
    result={'status':'VERIFIED-COMPUTATIONALLY','literal_grid_paths':dp_cases,'embedded_centre_words':embeddings,'cumulative_level_checks':level_checks,'control':{'claim':'C_h<=2h+1, omitting free letters','centre':'2222222','h':0,'actual':8,'claimed_upper':1,'outcome':'REFUTED_AS_REQUIRED'}}
    (ROOT/'geometry_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':check()
