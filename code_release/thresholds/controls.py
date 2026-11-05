#!/usr/bin/env python3
"""Literal incidence cross-check and deliberately invalid certificates."""
import sys
sys.dont_write_bytecode=True
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import subprocess
from verify import Geometry,verify

ROOT=Path(__file__).resolve().parent


def literal():
    N,r=4,3
    C=N-r+1
    labels=[(a,p) for a in product(range(3),repeat=r-1) for p in range(C+1)]
    orbits={q:[] for q in labels}
    for x in product(range(4),repeat=N):
        a=tuple(2 if v>=2 else v for v in x[:r-1])
        p=sum(v!=0 for v in x[r-1:])
        orbits[a,p].append(x)
    count=0
    for m in range(N+1):
        G=Geometry(N,r,m)
        for c in labels:
            a,p=c
            centre=tuple(a)+(1,)*p+(0,)*(C-p)
            for q in labels:
                actual=F(sum(sum(v!=w for v,w in zip(centre,x))<=m
                             for x in orbits[q]),len(orbits[q]))
                assert actual==G.kernel(c,q),(m,c,q,actual,G.kernel(c,q))
                count+=1
    return {'status':'VERIFIED-COMPUTATIONALLY','literal_comparisons':count,'passes':True}


def main():
    results={'literal':literal(),'controls':[]}
    agreement=0
    for C in (3,7,11,20):
        G=Geometry(C+2,3,C)
        for p in range(3):
            cumulative=G.tails(p)
            for w in range(C+1):
                G.probability(w,p,0,0)
                direct=G._short_tails[w,p]
                for budget in range(C+1):
                    assert sum(v for distance,v in direct.items() if distance<=budget)==cumulative[w][budget]
                    agreement+=1
    results['short_tail_count_agreement']={'status':'VERIFIED-COMPUTATIONALLY',
                                           'exact_comparisons':agreement,'passes':True}
    # Preserve all private-region masses and move a packing to the first
    # disallowed radius. Rejection must occur at its actual centre load.
    full=json.loads((ROOT/'certificates/N16_r6_m7.json').read_text())
    wrong=dict(full,m=8)
    wrong['support']=[dict(t,tail_weight=t['tail_weight']+1) for t in full['support']]
    cases=[('radius_increase_preserving_region_masses',wrong,'packing overload')]
    dual=json.loads((ROOT/'certificates/N16_r6_m8.json').read_text())
    wrong=dict(dual)
    wrong['support']=[dict(t,mass=str(F(t['mass'])*F(999,1000))) for t in dual['support']]
    cases.append(('scale_private_dual_down',wrong,'private dual undercoverage'))
    wrong={'N':16,'r':6,'m':7,'kind':'EMPTY_REGION','support':[]}
    cases.append(('invent_empty_private_region',wrong,'claimed empty region'))
    for name,cert,reason in cases:
        path=ROOT/f'control_{name}.invalid.json'
        path.write_text(json.dumps(cert,indent=2)+'\n')
        proc=subprocess.run([sys.executable,'-B',str(ROOT/'verify.py'),str(path)],
                            text=True,capture_output=True)
        answer=json.loads(proc.stdout)
        assert proc.returncode==1 and answer['passes'] is False
        assert reason in answer['reason'],answer
        results['controls'].append({'name':name,'exit_code':proc.returncode,
                                    'outcome':'REJECTED-AS-REQUIRED','status':'REFUTED',
                                    'reason':answer['reason']})
    (ROOT/'control_results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':
    main()
