#!/usr/bin/env python3
"""Executed failing controls for the infinite constructions and obstructions."""
import sys
sys.dont_write_bytecode=True
import json
from fractions import Fraction as F
from pathlib import Path
from verify import verify,require
from verify_additional import upper,constant_obstruction

ROOT=Path(__file__).resolve().parent


def finite_family(family,n,alpha=None):
    alpha=family['alpha'] if alpha is None else alpha
    r=family['r']
    m=(3*n-alpha)//5
    terms=[]
    for term in family['support']:
        if family.get('mass_type')=='linear_inverse_C':
            t=F(family['C0'],n-r)
            mass=(1-t)*F(term['mass_limit'])+t*F(term['mass_at_C0'])
        else:
            mass=F(term['mass'])
        if mass:
            terms.append({'prefix':term['prefix'],'tail_weight':m-term['distance'],'mass':str(mass)})
    return {'N':n-1,'r':r,'m':m,'kind':'FULL','support':terms}


def main():
    results=[]
    candidates=[('r6_use_five_run_boundary',finite_family(json.loads((ROOT/'family_r6.json').read_text()),17,11)),
                ('r7_finite_only_constant_candidate',finite_family(json.loads((ROOT/'family_r7_rejected_finite_candidate.json').read_text()),107)),
                ('r7_first_inverse_candidate',finite_family(json.loads((ROOT/'inverse_r7_rejected_first_candidate.json').read_text()),68))]
    if (ROOT/'family_r8.json').exists():
        candidates.append(('r8_ignore_late_transition',finite_family(json.loads((ROOT/'family_r8.json').read_text()),83,14)))
    for name,cert in candidates:
        path=ROOT/f'control_{name}.invalid.json'
        path.write_text(json.dumps(cert,indent=2)+'\n')
        try:
            verify(cert)
        except AssertionError as exc:
            require('packing overload' in str(exc),'wrong candidate rejection')
            results.append({'name':name,'status':'REFUTED','outcome':'REJECTED-AS-REQUIRED','reason':str(exc)})
        else:
            raise AssertionError('invalid candidate accepted: '+name)
    cert=json.loads((ROOT/'upper_r7.json').read_text())
    cert['weight_numerator_coefficients']=[[str(F(x)*F(999,1000)) for x in row]
                                          for row in cert['weight_numerator_coefficients']]
    try:
        upper(cert,check_bridges=False)
    except AssertionError as exc:
        require('upper coverage' in str(exc),'wrong upper-control rejection')
        results.append({'name':'r7_symbolic_upper_weights_scaled_down','status':'REFUTED',
                        'outcome':'REJECTED-AS-REQUIRED','reason':str(exc)})
    else:
        raise AssertionError('invalid symbolic upper control accepted')
    cert=json.loads((ROOT/'constant_mass_obstruction_r7.json').read_text())
    cert['support']=[dict(t,mass=str(F(t['mass'])*F(999,1000))) for t in cert['support']]
    try:
        constant_obstruction(cert)
    except AssertionError as exc:
        require('coefficient below one' in str(exc),'wrong constant-control rejection')
        results.append({'name':'r7_constant_obstruction_scaled_down','status':'REFUTED',
                        'outcome':'REJECTED-AS-REQUIRED','reason':str(exc)})
    else:
        raise AssertionError('invalid constant obstruction accepted')
    (ROOT/'extra_control_results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':
    main()
