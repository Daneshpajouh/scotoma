#!/usr/bin/env python3
"""Re-encode a cert_N_m_r certificate into the swept schema, so verify_ceiling_exhaustive.py can
sweep it over EVERY centre class instead of the hundred or so verify_ceiling_certificate.py samples.

Why this file exists. Two certificate schemas live in this directory and they carry the same object:
a list of (active-prefix type, free weight, mass). cert_N_m_r.json writes the prefix as an index into
itertools.product([0,1,2], repeat=AP) in lexicographic order; data/floor_certificates/cert_m*_r*.json
writes it as a string over {0, 1, X}, where X is the non-binary active prefix. Only the second schema
has an exhaustive checker, so two cells of the frontier table were readable only by the sampled one.
Nothing about them was weaker; the verifier could not read their encoding.

What it does and does not do. It changes the encoding and scales every mass by the certificate's own
claimed maximum load, so the maximum becomes exactly 1, which the swept schema requires. It computes
no loads and proves nothing. verify_ceiling_exhaustive.py then recomputes the maximum over all
(L+1)*3^AP classes: if the claimed maximum was right it finds exactly 1 and the total is the claimed
bound; if the true maximum is larger it finds a value above 1 and FAILS.

  usage:  convert_ceiling_certificate.py IN.json OUT.json
          convert_ceiling_certificate.py --self-check

--self-check converts cert_19_9_9.json and compares the result against the independently produced
data/floor_certificates/cert_m9_r9.json, which the exhaustive sweep already passes. That control
caught a real error while this file was being written: the prefix alphabet is {0, 1, X} and an
earlier version emitted {0, 1, 2}, which the verifier rejected. A converter aimed at an unchecked
cell must first reproduce a checked one.

Standard library only. Exit 0 on success, 1 on a failed self-check, 2 if it could not run.
"""
import itertools
import json
import sys
from fractions import Fraction as Fr

SYM = {0: '0', 1: '1', 2: 'X'}


def convert(src):
    d = json.load(open(src))
    for k in ('N', 'm', 'r', 'AP', 'L', 'Dx', 'Dglob', 'support',
              'max_load_numerator_over_Dx_times_Dglob', 'bound_num', 'bound_den'):
        if k not in d:
            raise SystemExit('NOT RUN: %s lacks the field %r; is it a cert_N_m_r certificate?'
                             % (src, k))
    N, m, r, AP, L = d['N'], d['m'], d['r'], d['AP'], d['L']
    max_load = Fr(d['max_load_numerator_over_Dx_times_Dglob'], d['Dx'] * d['Dglob'])
    prefixes = [''.join(SYM[t] for t in p) for p in itertools.product([0, 1, 2], repeat=AP)]
    dual, total = [], Fr(0)
    for idx, w, xnum in d['support']:
        if not 0 <= idx < len(prefixes):
            raise SystemExit('NOT RUN: prefix index %d out of range for AP=%d' % (idx, AP))
        mass = Fr(xnum, d['Dx']) / max_load
        total += mass
        dual.append({'prefix': prefixes[idx], 'free_weight': w,
                     'mass': '%d/%d' % (mass.numerator, mass.denominator)})
    cited = Fr(d['bound_num'], d['bound_den'])
    if total != cited:
        raise SystemExit('NOT RUN: rescaled total %s does not equal the certificate\'s own bound %s'
                         % (total, cited))
    return {'N': N, 'm': m, 'r': r, 'free_length': L,
            'status': 'RE_ENCODED_FOR_EXHAUSTIVE_SWEEP',
            'criterion': 'K >= ceil(LP*); a feasible dual of mass > r-1 proves K = r',
            'restriction': 'none; re-encoded from %s, masses scaled by its claimed max load'
                           % src.split('/')[-1],
            'exact_total_mass': '%d/%d' % (total.numerator, total.denominator),
            'exact_total_mass_float': float(total),
            'ceiling': r,
            'exact_max_ball_mass': '1',
            'max_ball_at': d.get('max_load_at'),
            'centre_classes_verified': (L + 1) * 3 ** AP,
            'support_orbits': len(dual),
            'provenance': 'encoding only; the bound is established by verify_ceiling_exhaustive.py',
            'dual': dual}


def self_check():
    got = convert('cert_19_9_9.json')
    want = json.load(open('data/floor_certificates/cert_m9_r9.json'))
    ok = True
    for k in ('N', 'm', 'r', 'free_length', 'exact_total_mass', 'centre_classes_verified'):
        if str(got[k]) != str(want[k]):
            print('SELF-CHECK FAIL %s: %r vs the published %r' % (k, got[k], want[k])); ok = False
    gp = sorted((e['prefix'], e['free_weight'], e['mass']) for e in got['dual'])
    wp = sorted((e['prefix'], e['free_weight'], e['mass']) for e in want['dual'])
    if gp != wp:
        print('SELF-CHECK FAIL: %d support orbits against the published %d; %d differ'
              % (len(gp), len(wp), sum(1 for a, b in zip(gp, wp) if a != b))); ok = False
    print('SELF-CHECK %s: cert_19_9_9 re-encodes to the published cert_m9_r9 '
          '(mass %s, %s classes)' % ('PASS' if ok else 'FAIL',
                                     got['exact_total_mass'], got['centre_classes_verified']))
    return 0 if ok else 1


def main(argv):
    if argv[:1] == ['--self-check']:
        return self_check()
    if len(argv) != 2:
        print('NOT RUN: usage: convert_ceiling_certificate.py IN.json OUT.json'
              '  |  --self-check')
        return 2
    out = convert(argv[0])
    json.dump(out, open(argv[1], 'w'), indent=1)
    print('%s -> %s: mass %s = %.9f, %s classes to sweep'
          % (argv[0].split('/')[-1], argv[1].split('/')[-1], out['exact_total_mass'],
             out['exact_total_mass_float'], out['centre_classes_verified']))
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
