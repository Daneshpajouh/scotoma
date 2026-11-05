"""Convert a certify_cell.py certificate into the released per-cell schema, and index it.

REFUSES TO WRITE UNLESS THE INDEPENDENT CHECKER PASSES FIRST. Five row-end certificates were
withdrawn on 2026-10-22 after being verified by the routine that produced them, so a certificate
reaches data/ here only by surviving verify_ceiling_certificate.py, which computes ball loads by a
dynamic program over explicit strings and shares no code with the solver, and
verify_F_independent.py, which checks the closed form the sweep funnels through over its whole
domain. Both are run as subprocesses and a non-zero exit stops the write.

The solver's measure has maximum load Mex over all centres, which may exceed 1; the released schema
records a dual whose maximum load is at most 1. The two differ only by the exact rational scaling
x -> x/Mex, which is applied here, so the released total mass IS the certified bound.

WHAT IS VERIFIED OVER WHAT. The released cells from the binary-supported dual verify their
constraint over the 2^(r-1) binary-prefix centre classes, because that dual's own support made the
reduction to binary prefixes available. These cells do not: their support carries non-binary
prefixes, so the constraint is verified over ALL 3^(r-1) x (free length + 1) centre classes. The
field centre_classes_verified says which, and restriction says none.
"""
import json, os, subprocess, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, 'data', 'floor_certificates')

def main(src):
    C0 = json.load(open(src))
    # THE CHECKER MUST RUN AT THE CERTIFICATE'S OWN BUDGET. verify_F_independent.py defaults to
    # m=8 and r=11..16; invoked bare it would validate an m=9 or m=12 certificate against the
    # closed form at m=8, which is a check running on the wrong object and worse than none. The
    # budget, the length and the run count all come from the certificate being emitted.
    env = dict(os.environ, CERT_N=str(C0['N']), CERT_M=str(C0['m']), VF_R=str(C0['r']))
    for checker, args in (('verify_F_independent.py', []),
                          ('verify_ceiling_certificate.py', [src])):
        p = os.path.join(HERE, checker)
        r = subprocess.run([sys.executable, '-I', p] + args, capture_output=True, text=True,
                           env=env)
        tail = (r.stdout or r.stderr).strip().splitlines()[-1:] or ['(no output)']
        print('  %-32s exit %d   %s' % (checker, r.returncode, tail[0][:90]))
        if r.returncode != 0:
            sys.exit('REFUSING TO WRITE: %s did not pass, so this certificate is not established.'
                     % checker)

    C = C0
    N, m, r, AP, L, Dx, DG = C['N'], C['m'], C['r'], C['AP'], C['L'], C['Dx'], C['Dglob']
    worst = C['max_load_numerator_over_Dx_times_Dglob']
    mass  = Fr(sum(n for _, _, n in C['support']), Dx)
    Mex   = Fr(worst, Dx*DG)
    bound = mass / Mex
    assert bound == Fr(C['bound_num'], C['bound_den'])

    import itertools
    TYPES = list(itertools.product([0, 1, 2], repeat=AP))
    dual = []
    for i, w, num in C['support']:
        scaled = Fr(num, Dx) / Mex                       # exact; max load becomes <= 1
        dual.append({'prefix': ''.join('01X'[t] for t in TYPES[i]),
                     'free_weight': w, 'mass': '%d/%d' % (scaled.numerator, scaled.denominator)})
    assert sum(Fr(d['mass']) for d in dual) == bound

    classes = (L + 1) * 3**AP
    rec = {
        'N': N, 'm': m, 'r': r, 'free_length': L, 'status': 'PROVED',
        'criterion': 'K >= ceil(LP*); a feasible dual of mass > r-1 proves K = r',
        'restriction': 'none; the dual support carries non-binary active prefixes',
        'exact_total_mass': '%d/%d' % (bound.numerator, bound.denominator),
        'exact_total_mass_float': float(bound),
        'ceiling': r,
        'exact_max_ball_mass': '1',
        'max_ball_at': C.get('max_load_at'),
        'centre_classes_verified': classes,
        'centre_classes_note': 'ALL (active-prefix type, free weight) classes, 3^%d x %d. Not the '
                               '2^%d binary-prefix classes: this dual has non-binary support, so '
                               'the binary reduction is unavailable and was not used.'
                               % (AP, L+1, AP),
        'support_orbits': len(dual),
        'rounding': 'solver masses rounded DOWN to denominator %d (preserves feasibility), then '
                    'scaled by the exact maximum load %s so the released dual is feasible '
                    'everywhere with maximum ball mass exactly 1'
                    % (Dx, '%d/%d' % (Mex.numerator, Mex.denominator)),
        'dual': dual,
        'DO_NOT_READ_THE_MASS_AS_A_MEASUREMENT':
            'exact_total_mass is a safe LOWER bound on tau*, not tau* itself: masses were rounded '
            'down and the column set offered to the LP was capped, both of which can only weaken '
            'it. ceil of it is what proves K = r.',
        'independently_checked_by': ['verify_ceiling_certificate.py', 'verify_F_independent.py'],
    }
    os.makedirs(OUT, exist_ok=True)
    # N IS PART OF THE CELL'S IDENTITY. The name used to be cert_m{m}_r{r}.json, so a certificate at
    # a different spacer length would have silently overwritten the panel's and corrupted the index,
    # which keys on (m, r) alone. N=19 keeps the established name so the existing files and every
    # figure quoting them stay valid; any other length carries its N.
    dest = os.path.join(OUT, ('cert_m%d_r%d.json' % (m, r)) if N == 19
                             else ('cert_N%d_m%d_r%d.json' % (N, m, r)))
    json.dump(rec, open(dest, 'w'), indent=1)
    print('wrote %s  mass %s = %.9f  ceil %d  r %d'
          % (os.path.basename(dest), rec['exact_total_mass'], float(bound),
             -(-bound.numerator // bound.denominator), r))

    ipath = os.path.join(OUT, 'INDEX.json')
    idx = json.load(open(ipath)) if os.path.exists(ipath) else {'criterion': rec['criterion'],
                                                                'cells': []}
    idx['cells'] = [c for c in idx['cells']
                    if not (c['m'] == m and c['r'] == r and c.get('N', 19) == N)]
    idx['cells'].append({'N': N, 'm': m, 'r': r, 'status': 'PROVED',
                         'exact_total_mass': rec['exact_total_mass'], 'ceiling': r,
                         'exact_max_ball_mass': '1',
                         'centre_classes_verified': classes,
                         'support_orbits': len(dual),
                         'restriction': 'none'})
    idx['cells'].sort(key=lambda c: (c.get('N', 19), c['m'], c['r']))
    json.dump(idx, open(ipath, 'w'), indent=1)
    print('indexed: %d cells, m=%d row now reaches r=%d'
          % (len(idx['cells']), m, max(c['r'] for c in idx['cells'] if c['m'] == m
                                       and c.get('N', 19) == N and c['status'] == 'PROVED')))

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'cert_19_8_11.json')
