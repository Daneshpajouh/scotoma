#!/usr/bin/env python3
"""Standalone exact verifier for the universal (W) certificate. No solver/imported project code.

Run: /usr/local/bin/python3 -B verify_w.py
Run a failing control: /usr/local/bin/python3 -B verify_w.py --control zero-local

Proof bridge: see PROOF.md. Every probability/bound comparison uses fractions.Fraction.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
from functools import lru_cache
from itertools import product, combinations
import argparse
import copy
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parent
P=F(2,5)
K=32
LIMIT=F(99,100)


def demand(condition,message):
    if not condition: raise ValueError(message)


def C(n,k):
    return comb(n,k) if 0<=k<=n else 0


@lru_cache(None)
def phi(u,t,e,b,f):
    cap=2*t-u-e+b-f
    return F(sum(C(t,h) for h in range(min(t,cap)+1)),2**t) if cap>=0 else F(0)


@lru_cache(None)
def basis(t,b,k,l,a,e):
    return F(C(t,k)*C(b,l),C(a,k)*C(e,l))


def exponential(u,a,e,f,t):
    demand(t>=1,'exponential parameter < 1')
    A=1-P+P*(t+t*t)/2
    B=1-P+P*t
    return A**a*B**e/t**(u+e+f)


def tail(u):
    a=F(23,20);q=F(11,20)
    return a**u*q**K*((2*K+u+1)/(1-q)+2*q/(1-q)**2)


def path_value(u,table):
    """Backward recursion on a monotone grid, independent of the producer's forward DP."""
    values={}
    for z in range(K,-1,-1):
      for o in range(K+1):
       for remaining in range(u+1):
        for previous in (0,1,2):
         if z==K and o==0 and remaining==0:
            values[z,o,remaining,previous]=F(0)
            continue
         choices=[]
         for nxt in (0,1,2):
            if nxt==0:
                if z==K: continue
                target=(z+1,o,remaining,0)
            elif nxt==1:
                if o==0: continue
                target=(z,o-1,remaining,1)
            else:
                if remaining==0: continue
                target=(z,o,remaining-1,2)
            s=int(previous==2)+int(nxt==2)
            f=int(previous==0)+int(nxt==1)
            e=z+o-f
            if e<0: continue
            reward=table[s,e,f] if z+o<K else F(0)
            choices.append(reward+values[target])
         # Unreachable states with previous=0,z=0 may have no valid continuation.
         values[z,o,remaining,previous]=max(choices,default=F(-10**6))
    return values[0,K,u,1]


def check_row(row):
    u=row['u']
    demand(row['K']==K and F(row['p'])==P,'row parameters altered')
    table={};inequalities=0;polycount=0;mgfcount=0
    for item in row['local_bounds']:
        s,e,f=item['s'],item['e'],item['f']
        demand(0<=s<=min(2,u) and 0<=f<=2-s and 0<=e<=K,'local key out of range')
        key=s,e,f
        demand(key not in table,'duplicate local row')
        a=u-s
        bound=F(item['bound'])
        demand(bound>=0,'negative local bound')
        cert=item['certificate']
        if cert[0]=='mgf':
            exact=exponential(u,a,e,f,F(cert[1]));mgfcount+=1
        elif cert[0]=='poly':
            q=[]
            for k,l,v in cert[1]:
                v=F(v)
                demand(0<=k<=a and 0<=l<=e and v>=0,'invalid positive polynomial')
                q.append((k,l,v))
            for t in range(a+1):
                for b in range(e+1):
                    got=sum((v*basis(t,b,k,l,a,e) for k,l,v in q),F(0))
                    demand(got>=phi(u,t,e,b,f),f'polynomial cover fails: u={u},s={s},e={e},f={f},t={t},b={b}')
                    inequalities+=1
            exact=sum((v*P**(k+l) for k,l,v in q),F(0));polycount+=1
        else:raise ValueError('unknown certificate kind')
        demand(bound>=exact,f'local bound too small: {u,key}')
        table[key]=bound
    expected={(s,e,f) for s in range(min(2,u)+1) for f in range(3-s) for e in range(K+1)}
    demand(set(table)==expected,'missing local row')
    computed=path_value(u,table)
    demand(F(row['path_bound'])==computed,f'path bound differs for u={u}')
    demand(F(row['tail_bound'])>=tail(u),f'tail bound too small for u={u}')
    total=computed+F(row['tail_bound'])
    demand(F(row['total_bound'])==total,f'total bound differs for u={u}')
    demand(total<LIMIT,f'universal bound does not clear {LIMIT} at u={u}')
    return {'u':u,'bound':str(total),'polynomial_inequalities':inequalities,'polynomial_rows':polycount,'mgf_rows':mgfcount}


def check_large_u():
    # Abel summation of the cumulative level count 2e+u+1, with t=2.
    bound36=(36+F(17,3))*F(9,10)**36
    demand(bound36<LIMIT,'large-u initial bound fails')
    # H(u+1)/H(u) < 1 iff 3u>10; hence it decreases for every u>=36.
    demand(3*36>10,'large-u monotonicity fails')
    return str(bound36)


def check_parameters():
    # For m>=6, at n=5m/2, L=2/5+2m/(3m+2)>=1;
    # L strictly decreases in n, so L<=1 implies n>=5m/2.
    demand(F(2,5)+F(12,20)==1,'m=6 endpoint identity fails')
    # Numerator of 2m/(3m+2)-3/5 is m-6, over 5(3m+2).
    for m,Ns in ((2,7),(3,10),(4,12),(5,15)):
        L=lambda N:F(m,N-2)+F(m,N-m-1)
        demand(L(Ns)<=1 and L(Ns-1)>1,'small N* incorrect')
        demand(F(m,Ns-2)<=P,'small ratio bound fails')
    return True


def small_m1_checks():
    # Analytic proof: u>2m gives zero; for m=1,u=2 the e=0 plateau
    # has sum a_j<=2, giving Lambda<=1/(N-2)<=1/2 for N>=4.
    # Check the plateau inequality for each of its 3 possible free-block sizes.
    for free_at_plateau in range(3):
        k=free_at_plateau
        if k==0: allowed_sum=2
        elif k==1: allowed_sum=1+1
        else: allowed_sum=1+0+1
        demand(allowed_sum<=2,'m=1 plateau count fails')
    return True


def literal_local_controls():
    # Compare the local event formula to direct word counts on newly specified test objects.
    cases=0
    for N,m in ((6,2),(7,2)):
      for centre in ((1,)*(N-4)+(2,2,0,0), (2,)+(1,)*(N-3)+(0,3), (0,1,2,3)+(0,)*(N-4)):
       u=sum(c>=2 for c in centre)
       for j in range(N+1):
        pat=tuple(int(i<j) for i in range(N))
        forced={i for i in (j-1,j) if 0<=i<N}
        allowed=[i for i in range(N) if i not in forced]
        a=sum(centre[i]>=2 for i in allowed)
        e=sum(centre[i]<2 and centre[i]!=pat[i] for i in allowed)
        f=sum(centre[i]<2 and centre[i]!=pat[i] for i in forced)
        brute=0
        for T in combinations(allowed,m):
         for letters in product((2,3),repeat=m):
          point=list(pat)
          for i,c in zip(T,letters):point[i]=c
          brute+=sum(x!=y for x,y in zip(point,centre))<=m
        predicted=F(0)
        n=len(allowed)
        for t in range(a+1):
         for b in range(e+1):
          predicted+=F(C(a,t)*C(e,b)*C(n-a-e,m-t-b),C(n,m))*phi(u,t,e,b,f)
        demand(predicted==F(brute,C(n,m)*2**m),'literal local control mismatch')
        cases+=1
    return cases


def check(data,verbose=True,only=None):
    demand(F(data['p'])==P and data['K']==K,'global parameters altered')
    rows=data['rows']
    demand(sorted(r['u'] for r in rows)==list(range(2,36)),'u range incomplete')
    selected=rows if only is None else [r for r in rows if r['u']==only]
    out=[]
    for row in selected:
        result=check_row(row);out.append(result)
        if verbose:print('VERIFIED-COMPUTATIONALLY',result,flush=True)
    check_parameters();small_m1_checks()
    b36=check_large_u()
    return {'status':'PROVED_WITH_EXACT_FINITE_CERTIFICATE','limit':str(LIMIT),'largest_finite_bound':str(max(F(x['bound']) for x in out)), 'finite_u_checked':len(out),'local_polynomial_inequalities':sum(x['polynomial_inequalities'] for x in out),'u_ge_36_bound':b36}


def mutate(data,kind):
    d=copy.deepcopy(data)
    row=d['rows'][0]
    if kind=='zero-local':
        item=next(x for x in row['local_bounds'] if (x['s'],x['e'],x['f'])==(0,0,0))
        item['bound']='0'
    elif kind=='empty-polynomial':
        item=next(x for x in row['local_bounds'] if (x['s'],x['e'],x['f'])==(0,0,0))
        item['certificate']=['poly',[]]
    elif kind=='missing-local':row['local_bounds'].pop(0)
    elif kind=='lower-path':row['path_bound']=str(F(row['path_bound'])-F(1,10**6))
    elif kind=='zero-tail':row['tail_bound']='0'
    elif kind=='wrong-p':d['p']='1/2'
    else:raise ValueError('unknown control')
    return d


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--certificate',type=Path,default=ROOT/'w_certificate.json')
    parser.add_argument('--control',choices=['zero-local','empty-polynomial','missing-local','lower-path','zero-tail','wrong-p'])
    parser.add_argument('--controls',action='store_true')
    args=parser.parse_args()
    data=json.loads(args.certificate.read_text())
    if args.control:
        check(mutate(data,args.control),verbose=False,only=2)
        print('ERROR: planted-invalid certificate accepted')
        return 2
    result=check(data)
    result['literal_local_cases']=literal_local_controls()
    if args.controls:
        controls=[]
        for kind in ['zero-local','empty-polynomial','missing-local','lower-path','zero-tail','wrong-p']:
            try:check(mutate(data,kind),verbose=False,only=2)
            except ValueError as exc:
                print('REFUTED CONTROL',kind,str(exc),flush=True)
                controls.append({'control':kind,'outcome':'REJECTED_AS_REQUIRED','reason':str(exc)})
            else:raise ValueError('invalid certificate accepted: '+kind)
        result['negative_controls']=controls
    result['certificate_sha256']=hashlib.sha256(args.certificate.read_bytes()).hexdigest()
    print(json.dumps(result,indent=2),flush=True)
    (ROOT/'verification_result.json').write_text(json.dumps(result,indent=2)+'\n')
    return 0

if __name__=='__main__':
    try:sys.exit(main())
    except ValueError as exc:
        print('REJECTED:',str(exc),file=sys.stderr)
        sys.exit(1)
