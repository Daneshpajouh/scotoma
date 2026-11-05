#!/usr/bin/env python3
r"""Convert a cert_N_m_r certificate into the swept schema, so the release's EXHAUSTIVE verifier
can sweep it over every centre class.

WHY. The paper's m=9 r<=12 and m=10 r<=9 rows cite cert_19_9_12 and cert_19_10_9. Those belong to
the family whose only checker is SAMPLED, about a hundred centres. The other family has
verify_ceiling_exhaustive.py, which sweeps all (L+1)*3^AP classes. The two schemas carry the same
object: a list of (active-prefix type, free weight, mass). This converts the encoding and changes
nothing else, so the exhaustive verifier does the actual work and I reimplement no load computation.

The masses are scaled by the cited max load so that the maximum becomes exactly 1, which is what the
swept schema requires. If the cited max load is correct the verifier will recompute 1 and the total
will be the cited bound. If the true max load is LARGER, the verifier will see a maximum above 1 and
FAIL, which is precisely the check that does not currently exist for these two cells.
"""
import itertools, json, sys
from fractions import Fraction as Fr

src, out = sys.argv[1], sys.argv[2]
d = json.load(open(src))
N, m, r, AP, L, Dx, Dg = d['N'], d['m'], d['r'], d['AP'], d['L'], d['Dx'], d['Dglob']
max_load = Fr(d['max_load_numerator_over_Dx_times_Dglob'], Dx * Dg)
# The swept verifier's prefix alphabet is {0, 1, X}, not {0, 1, 2}: the third symbol is the
# non-binary active prefix, written X. The `max_load_at` fields of the cited certificates use
# X too, which I saw and did not follow up until the control rejected '22110201'.
SYM = {0: '0', 1: '1', 2: 'X'}
prefixes = [''.join(SYM[t] for t in p) for p in itertools.product([0, 1, 2], repeat=AP)]
assert len(prefixes) == 3 ** AP, "prefix enumeration"

dual, total = [], Fr(0)
for idx, w, xnum in d['support']:
    assert 0 <= idx < len(prefixes), f"index {idx} out of range for AP={AP}"
    mass = Fr(xnum, Dx) / max_load          # scale so the claimed maximum becomes exactly 1
    total += mass
    dual.append({'prefix': prefixes[idx], 'free_weight': w, 'mass': f"{mass.numerator}/{mass.denominator}"})

cited = Fr(d['bound_num'], d['bound_den'])
print(f"{src}: N={N} m={m} r={r} AP={AP} L={L}  support {len(dual)} orbits")
print(f"  cited max load      {max_load} = {float(max_load):.6f}")
print(f"  rescaled total mass {total} = {float(total):.6f}")
print(f"  cited bound         {cited} = {float(cited):.6f}   equal: {total == cited}")
print(f"  needs > {r-1} for K={r}: {total > r - 1}")

json.dump({'N': N, 'm': m, 'r': r, 'free_length': L,
           'status': 'CONVERTED_FOR_EXHAUSTIVE_CHECK',
           'criterion': 'K >= ceil(LP*); a feasible dual of mass > r-1 proves K = r',
           'restriction': 'none; converted from %s, masses scaled by its cited max load' % src.split('/')[-1],
           'exact_total_mass': f"{total.numerator}/{total.denominator}",
           'exact_total_mass_float': float(total),
           'ceiling': r,
           'exact_max_ball_mass': '1',
           'max_ball_at': d.get('max_load_at'),
           'centre_classes_verified': (L + 1) * 3 ** AP,
           'support_orbits': len(dual),
           'dual': dual}, open(out, 'w'), indent=1)
print(f"  wrote {out}")
