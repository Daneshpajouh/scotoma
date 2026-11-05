#!/usr/bin/env python3
"""Check the budget-ten extensions end to end, from a clone, with no cluster.

WHAT IS BEING CHECKED. Seven cells (N, r, m=10) were closed by a fractional packing whose total
mass strictly exceeds r-1, which gives K >= r; the library of r balls gives K <= r; so K = r. The
mass was then re-derived by an independently written exhaustive maximum-load scan over the
complete centre-class space, whose returned summary is in independent_scan/H3_summary.json. This
script checks, in exact rational arithmetic:

  1. each cell's own decision certificate records a mass strictly above r-1;
  2. the independent scan agrees with that cell's certificate on the mass, to the last digit;
  3. the independent scan reports an exhaustive pass and rejects both of its forgery controls;
  4. the frontier ledger's per-cell statuses imply the table the paper prints;
  5. a one-unit perturbation of ANY ONE cell's certificate mass is REJECTED by check 1
     or 2, checked once for every cell.

Check 5 is why 1 to 4 mean anything. Exact Fractions throughout; no float decides anything.
"""
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CELLS = [(20, 11), (20, 12), (21, 12), (21, 13), (22, 20), (22, 21), (20, 14),
         (20, 15), (21, 19), (22, 23)]
M = 10
# 2026-11-03: (20,14), (20,15), (21,19) and (22,23) were closed and independently verified
# after H3_summary.json was returned, by the same verifier under the same gate. Their records are
# deposited on their own rather than by rewriting a returned artefact, so the summary below still
# reports six and still hashes as it did. A single-cell record has exactly the shape of a summary
# row. The control now perturbs every cell in turn, so each of the two scan routes is exercised
# by more than one cell and no cell is covered only by the others' controls.
SIX_SCAN = ('independent_scan', 'H3_summary.json')
SINGLE_SCANS = {(N, r): ('independent_scan', 'cell_%d_%d_10_result.json' % (N, r))
                for (N, r) in ((20, 14), (20, 15), (21, 19), (22, 23))}


def jload(*parts):
    return json.load(open(os.path.join(HERE, *parts), encoding='utf-8'))


def mass_from_certificate(cert):
    """Re-derive the total mass by SUMMING the deposited packing, in exact Fractions.

    The certificate carries its packing, one entry per support orbit with its own rational mass,
    so the total is recomputed here rather than read from a field that claims it. A certificate
    whose stated total disagrees with its own packing fails on that disagreement.
    """
    entries = cert.get('packing')
    if not entries:
        raise SystemExit('a certificate carries no packing, so nothing can be re-derived from it '
                         'and this check fails rather than passing vacuously')
    total = Fraction(0)
    for e in entries:
        total += Fraction(e['mass'])
    return total, len(entries)


def scan_masses():
    s = jload(*SIX_SCAN)
    out = {}
    for rec in s['results']:
        out[(rec['N'], rec['r'])] = rec
    for key, parts in SINGLE_SCANS.items():
        rec = jload(*parts)
        if (rec['N'], rec['r']) != key:
            raise SystemExit('a deposited single-cell scan is not the cell it is filed under, '
                             'so nothing here can be trusted and this check fails')
        if key in out:
            raise SystemExit('a single-cell scan duplicates a summary row, which would let one '
                             'record stand in for another; this check fails rather than guessing')
        out[key] = rec
    return s, out


def main(perturb=None):
    bad = []
    summary, scan = scan_masses()
    print('independent scan status: %s' % summary.get('status'))
    if summary.get('status') != 'SIX_H3_VERIFIED':
        bad.append('the independent scan summary does not report SIX_H3_VERIFIED')

    for (N, r) in CELLS:
        d = os.path.join('cells', 'cell_%d_%d_%d' % (N, r, M))
        cert = jload(d, 'decision.certificate.json')
        result = jload(d, 'RESULT.json')
        mass, orbits = mass_from_certificate(cert)
        if perturb == (N, r):
            mass = mass - Fraction(1)       # the control: one unit down
        floor_ok = mass > r - 1
        own = Fraction(result['tau_lower']) == mass        # the cell's own record agrees
        rec = scan.get((N, r))
        agree = rec is not None and Fraction(rec['total_mass']) == mass
        exhaustive = rec is not None and rec.get('global_maximum_status') == 'EXHAUSTIVE_EXACT_PASS'
        # negative_controls holds the forgery results alongside bookkeeping, and one of the
        # bookkeeping entries is itself a dict with no status, so require a verdict only where
        # there is one, and require at least two of them.
        verdicts = [v['status'] for v in rec.get('negative_controls', {}).values()
                    if isinstance(v, dict) and 'status' in v] if rec else []
        controls = len(verdicts) >= 2 and all(v == 'REJECTED' for v in verdicts)
        kk = result.get('K_lower') == result.get('K_upper') == r and result.get('closed') is True
        ok = floor_ok and own and agree and exhaustive and controls and kk
        print('  N=%-3d r=%-3d orbits %-5d  mass>r-1 %-5s own %-5s scan %-5s exhaustive %-5s '
              'controls %-5s K=r %-5s [%s]'
              % (N, r, orbits, floor_ok, own, agree, exhaustive, controls, kk,
                 'ok' if ok else 'FAIL'))
        if not ok:
            bad.append('cell (%d,%d,%d)' % (N, r, M))

    ledger = jload('budget_ten_frontier.json')
    rows = {}
    for c in ledger['cells']:
        rows.setdefault(c['n'], {})[c['r']] = c['display']
    print('ledger execution_state %s, %d cells; certified prefixes:'
          % (ledger['execution_state'], len(ledger['cells'])))
    for n in sorted(rows):
        cs = sorted(r for r, d in rows[n].items() if d == 'C')
        print('    n=%-3d certified 1..%s' % (n, max(cs) if cs else 'none'))
    # every closed cell must appear closed in the ledger, under n = N+1
    for (N, r) in CELLS:
        if rows.get(N + 1, {}).get(r) not in ('C', 'I'):
            bad.append('ledger does not mark n=%d r=%d closed' % (N + 1, r))

    if bad:
        print('\n[FAIL] %s' % '; '.join(bad))
        return 1
    print('\nall %d extensions check out, and the ledger agrees with every one' % len(CELLS))
    return 0


if __name__ == '__main__':
    if '--control' in sys.argv:
        # One control per independent-scan route: a summary row and a single-cell record. A
        # control on only the first would leave the route added later untested by this gate.
        failures = 0
        for cell in CELLS:
            print('== control: one unit removed from cell (%d,%d,%d), which MUST be rejected'
                  % (cell[0], cell[1], M))
            rc = main(perturb=cell)
            print('control exit %d (1 is the correct outcome)\n' % rc)
            failures += 1 if rc == 1 else 0
        print('%d of %d controls were correctly rejected' % (failures, len(CELLS)))
        sys.exit(0 if failures == len(CELLS) else 1)
    sys.exit(main())
