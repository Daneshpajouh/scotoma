#!/usr/bin/env python3
"""Portable exact proof verifier; optional full measurement-grid replay.

No repository imports and no numerical/scientific packages are required.
The supplied E4/E13/E16/E21 results remain explicit premises.
"""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse
import json
import time
from prove_family import prove
from verify import verify,require
from verify_additional import upper,constant_obstruction
import controls
import extra_controls
import verify_general

ROOT=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--full-grid',action='store_true',help='also replay every directly measured boundary pair')
    args=ap.parse_args()
    started=time.monotonic()
    results={'status':'VERIFIED-COMPUTATIONALLY','families':[],'upper_families':[],
             'supplied_premises':['E4 restriction','E13 six-run upper family',
                                  'E16 budget monotonicity','E21 five-run threshold']}
    verify_general.main()
    results['general_theorem_checks']=json.loads((ROOT/'general_verification.json').read_text())
    for name in ('family_r6.json','family_r7_lower14.json','family_r8.json',
                 'family_limit_r9.json','family_limit_r10.json'):
        cert=prove(ROOT/name)
        results['families'].append({'file':name,'r':cert['r'],'alpha':cert['alpha'],
                                    'min_n':cert['min_n'],'polynomials':cert['polynomial_rows'],
                                    'newton_coefficients':cert['checked_newton_coefficients']})
    for name in ('upper_r7.json','upper_r7_n46.json'):
        results['upper_families'].append(dict(file=name,**upper(json.loads((ROOT/name).read_text()))))
    extra=ROOT/'upper_r8_n83.json'
    eight_closed=False
    if extra.exists():
        candidate=json.loads(extra.read_text())
        if candidate.get('positive_from_z') is not None:
            results['upper_families'].append(dict(file=extra.name,**upper(candidate)))
            eight_closed=True
    results['constant_mass_obstruction']=constant_obstruction(json.loads((ROOT/'constant_mass_obstruction_r7.json').read_text()))
    rows=json.loads((ROOT/'threshold_table.json').read_text())
    critical6={10,11,12}
    critical7=set(range(10,16))|{16,21,26,31,36,41}
    critical8=set(range(10,51))|{53,58,63,68,73,78}
    available={(x['n'],x['r']) for x in rows if x['status']=='VERIFIED-COMPUTATIONALLY'}
    require({(n,6) for n in critical6}|{(n,7) for n in critical7}|{(n,8) for n in critical8} <= available,
            'missing required finite boundary certificate pair')
    count=0
    for row in rows:
        if row['status']!='VERIFIED-COMPUTATIONALLY':
            continue
        critical=(row['r']==6 and row['n'] in critical6 or
                  row['r']==7 and row['n'] in critical7 or
                  row['r']==8 and row['n'] in critical8)
        if args.full_grid or critical:
            for key,offset in (('lower_certificate',0),('upper_certificate',1)):
                cert=json.loads((ROOT/row[key]).read_text())
                require((cert['N'],cert['r'],cert['m'])==(row['N'],row['r'],row['M']+offset),
                        'table/certificate parameter mismatch')
                require((cert['kind']=='FULL')==(offset==0),'wrong boundary certificate kind')
                verify(cert)
            if row['r']==6:
                require(row['M']==(3*row['n']-12)//5+(row['n']==12),'six-run formula disagrees with a finite check')
            if row['r']==7:
                require(row['M']==(3*row['n']-14)//5+(row['n'] in {11,14,16,21,26,31,36,41}),
                        'seven-run formula disagrees with a finite check')
            if row['r']==8 and eight_closed:
                exceptional={11,16,21,26}|{13+5*j for j in range(14)}
                require(row['M']==(3*row['n']-15)//5+(row['n'] in exceptional),
                        'eight-run formula disagrees with a finite check')
            count+=1
    results['boundary_pairs_replayed']=count
    results['full_grid_requested']=args.full_grid
    controls.main()
    extra_controls.main()
    results['controls']=json.loads((ROOT/'control_results.json').read_text())['controls']+json.loads((ROOT/'extra_control_results.json').read_text())+results['general_theorem_checks']['controls']
    require(all(x['outcome']=='REJECTED-AS-REQUIRED' for x in results['controls']),
            'at least one invalid control was accepted')
    results['seconds']=round(time.monotonic()-started,3)
    results['passes']=True
    (ROOT/'verification_summary.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps({'passes':True,'families':len(results['families']),
                      'upper_families':len(results['upper_families']),
                      'boundary_pairs':count,'rejected_controls':len(results['controls']),
                      'seconds':results['seconds']}),flush=True)


if __name__=='__main__':
    main()
