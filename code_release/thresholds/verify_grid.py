#!/usr/bin/env python3
"""Recheck every measured table boundary; standard library only, no repository."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse
import json
import time
from verify import verify,require

ROOT=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--r',type=int,nargs='+')
    ap.add_argument('--n',type=int,nargs='+')
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    rows=json.loads((ROOT/'threshold_table.json').read_text())
    require(len(rows)==432 and sum(x['status']=='VERIFIED-COMPUTATIONALLY' for x in rows)==412,
            'threshold table must retain 432 rows and 412 direct boundary pairs')
    checked=[]
    for x in rows:
        if x['status']!='VERIFIED-COMPUTATIONALLY' or args.r and x['r'] not in args.r:
            continue
        if args.n and x['n'] not in args.n:
            continue
        start=time.monotonic()
        lo=json.loads((ROOT/x['lower_certificate']).read_text())
        hi=json.loads((ROOT/x['upper_certificate']).read_text())
        require((lo['N'],lo['r'],lo['m'],lo['kind'])==(x['N'],x['r'],x['M'],'FULL'),
                'table lower-bound parameters mismatch')
        require((hi['N'],hi['r'],hi['m'])==(x['N'],x['r'],x['M']+1) and hi['kind']!='FULL',
                'table upper-bound parameters mismatch')
        lcheck=verify(lo)
        hcheck=verify(hi)
        answer={'n':x['n'],'r':x['r'],'M':x['M'],'status':'VERIFIED-COMPUTATIONALLY',
                'lower':lcheck,'upper':hcheck,'seconds':round(time.monotonic()-start,3)}
        checked.append(answer)
        print(json.dumps({'n':x['n'],'r':x['r'],'M':x['M'],'passes':True,
                          'seconds':answer['seconds']}),flush=True)
        if args.output:
            args.output.write_text(json.dumps(checked,indent=2)+'\n')
    require(bool(checked),'no matching direct boundary pairs were checked')
    if not args.r and not args.n:
        require(len(checked)==412,'incomplete full-grid replay')
    print(json.dumps({'all_passed':True,'cells_checked':len(checked)}),flush=True)


if __name__=='__main__':
    main()
