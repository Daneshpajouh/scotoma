"""Exact finite exceptions for the threshold selector. Standard library only."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
from functools import lru_cache
from itertools import product, combinations
import json
from exact_load import local_count, load, nstar, L
ROOT=Path(__file__).resolve().parent


def stratum(N,m,u,return_word=False,kill_free_step=False):
    best=F(-1);bestword='';bestO=-1
    for O in range(N-u+1):
        Z=N-u-O
        @lru_cache(None)
        def reward(z,o,x,previous,nxt):
            j=z+(O-o)+(u-x)
            n=N-1 if j in (0,N) else N-2
            s=int(previous==2)+int(nxt==2)
            f=int(previous==0)+int(nxt==1)
            e=z+o-f
            if e<0:return None
            return F(local_count(n,m,u,u-s,e,f),comb(n,m)*2**m)
        vals={}
        for z in range(Z,-1,-1):
          for o in range(O+1):
           for x in range(u+1):
            for previous in (0,1,2):
             if z==Z and o==0 and x==0:
                r=reward(z,o,x,previous,0)
                vals[z,o,x,previous]=r if r is not None else F(-1000)
                continue
             choices=[]
             for nxt in (0,1,2):
                if nxt==0:
                    if z==Z:continue
                    target=(z+1,o,x,0)
                elif nxt==1:
                    if o==0:continue
                    target=(z,o-1,x,1)
                else:
                    if x==0 or kill_free_step:continue
                    target=(z,o,x-1,2)
                r=reward(z,o,x,previous,nxt)
                if r is not None:choices.append(r+vals[target])
             vals[z,o,x,previous]=max(choices,default=F(-1000))
        value=vals[0,O,u,1]
        if value>best:
            best=value;bestO=O
            if return_word:
                z=0;o=O;x=u;prev=1;word=''
                while z+O-o+u-x<N:
                    found=False
                    for nxt in (0,1,2):
                        if nxt==0:
                            if z==Z:continue
                            target=(z+1,o,x,0)
                        elif nxt==1:
                            if o==0:continue
                            target=(z,o-1,x,1)
                        else:
                            if x==0:continue
                            target=(z,o,x-1,2)
                        r=reward(z,o,x,prev,nxt)
                        if r is not None and r+vals[target]==vals[z,o,x,prev]:
                            word+=str(nxt);z,o,x,prev=target;found=True;break
                    assert found
                bestword=word
    if return_word:
        assert load(tuple(map(int,bestword)),m)==best
    return best,bestword,bestO


def controls():
    out=[]
    for N,m,u in [(6,2,2),(7,2,3)]:
        brute=max(load(c,m) for c in product(range(3),repeat=N) if c.count(2)==u)
        dp,w,_=stratum(N,m,u,True)
        assert dp==brute
        bad,_,_=stratum(N,m,u,kill_free_step=True)
        assert bad!=brute
        out.append({'N':N,'m':m,'u':u,'brute':str(brute),'dp':str(dp),'missing_free_step':'REJECTED_AS_REQUIRED'})
        print('VERIFIED-COMPUTATIONALLY finite DP control',out[-1],flush=True)
    return out

if __name__=='__main__':
    ctrl=controls();rows=[]
    for m in [3,4,5,6,7,9]:
        N=nstar(m)-1;strata=[]
        for u in range(2,min(2*m,N)+1):
            v,w,O=stratum(N,m,u,True)
            assert v<1,(N,m,u,v)
            strata.append({'u':u,'value':str(v),'witness':w})
        value=max(F(r['value']) for r in strata)
        rows.append({'m':m,'N':N,'max_u_ge_2':str(value),'L':str(L(m,N)),'strata':strata})
        print('VERIFIED-COMPUTATIONALLY finite exception',m,N,'max_u_ge_2',value,'L',L(m,N),flush=True)
        (ROOT/'selector_finite_cells.json').write_text(json.dumps({'controls':ctrl,'rows':rows},indent=2)+'\n')
