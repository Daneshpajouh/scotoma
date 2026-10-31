#!/usr/bin/env python3
"""Check the five floor certificates for budgets m = 8..12 at N = 20.

WHY A CERTIFICATE AND NOT A SEARCH. The floor at run count r holds iff M*_r < r/(r-1), where

    M*_r = min over lam >= 0 of [ r*lam + max_c sum_j max(0, w_j(c) - lam) ].

EVERY lam gives a rigorous upper bound, so a single lam below the threshold PROVES admission and a
thin margin does not weaken the proof at all. What a thin margin threatens is REPRODUCTION: a
replicator searching a coarser lam set finds no witness, concludes the opposite, and their negative
arrives looking like a refutation of ours. At m = 10 the margin is 371/466752 = 0.000795, seven
times thinner than the next row, so that failure is a live possibility rather than a hypothetical.

The remedy is to publish the witness rather than to hedge the sentence. data/
floor_certificates_m8_m12.json records, per budget, the exact lam, the inner maximum the dual must
return there, and a centre attaining it. A checker then evaluates ONE dual at a KNOWN lam instead of
searching, and the attaining centre gives the matching lower bound on the inner max, so a checker
whose own DP returns less knows its DP is wrong rather than ours.

The lam values are implementation-independent: any correct dual evaluated at them reproduces the
bound. This script is our implementation of that check.
"""
import json, os, sys
from fractions import Fraction as F

# A referee reaching for the isolated interpreter, python3 -I, gets -P with it, which leaves this
# script's own folder off sys.path, so "import topr" fails and the package looks broken even when
# every file is sitting beside it. A packaging trap that can be removed in code should never be
# moved into prose: a README warning is read by the referee who already succeeded, not by the one
# staring at the traceback.
_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)  # -I leaves the script's folder off sys.path

try:
    import topr
except ImportError:
    # A referee who copies this one file out of curiosity, rather than the whole release, gets a
    # stack trace and concludes the package is broken. Three lines are cheaper than that.
    print("GATE NOT RUN: topr.py is not beside this script. It ships in the same directory;\n"
          "copy the whole release rather than this file alone.\n"
          "This is an incomplete release, not a finding about the claim.", file=sys.stderr)
    sys.exit(2)

# Resolved relative to THIS FILE, not to the working directory, and tried in both layouts, because
# this script exists in two places: at the project root, where the data sits under code_release/,
# and inside code_release itself, which is the tree that ships. A literal 'code_release/data/...'
# was correct in the first and resolved to code_release/code_release/data/... in the second, so the
# shipped copy could not find a file sitting beside it. Patching the two copies differently would
# have fixed it and broken verify_release_copies_match.py, which requires them byte-identical;
# resolving from __file__ satisfies both. Moving a file into the shipping tree does not make its
# paths shipping-tree paths.
_cert_default = next(
    (p for p in (os.path.join(_here, 'data', 'floor_certificates_m8_m12.json'),
                 os.path.join(_here, 'code_release', 'data', 'floor_certificates_m8_m12.json'))
     if os.path.exists(p)),
    os.path.join(_here, 'data', 'floor_certificates_m8_m12.json'))
# The path is overridable so that this checker is executable from a CLEAN COPY and not only from the
# project root. Six of eight blind round-5 lanes named an independently executable certificate and
# verifier package as the single change that would raise the paper most, and the editor lane's reject
# turned on verification being conditional on our build. Measured 2026-10-23: the r=20 certificate in
# certificates/r20_m8/ already passes from an empty directory under python3 -I, but THIS checker did
# not, because the constant below was relative to the project layout. A referee who copies the
# checker, its three stdlib-only modules and the certificate into one directory got a GATE FAILED
# naming a path that was right beside it. Overriding by argv or environment costs nothing here and
# removes that. Same pattern as SITES129 in code_release/verify_missed_actives.py.
CERT = (sys.argv[1] if len(sys.argv) > 1
        else os.environ.get('FLOOR_CERTIFICATES', _cert_default))


def main():
    try:
        rows = json.load(open(CERT))
    except FileNotFoundError:
        sys.exit('GATE FAILED: %s is missing. The appendix cites it as the witness for the m = 8 to '
                 '12 admitted rows, so its absence makes that paragraph unverifiable.' % CERT)
    if len(rows) != 5:
        sys.exit('GATE FAILED: %d certificates, expected 5 (m = 8..12)' % len(rows))

    fails = []
    for rec in sorted(rows, key=lambda r: r['m']):
        N, m, r = rec['N'], rec['m'], rec['r']
        lam = F(rec['lam'])
        thr = F(r, r - 1)

        inner, word = topr.dp_shifted(N, m, lam)
        bound = r * lam + inner

        # 1. the dual at the published lam must return the published inner maximum
        if inner != F(rec['inner_max']):
            fails.append('m=%d: dual at lam=%s returns %s, certificate says %s'
                         % (m, lam, inner, rec['inner_max']))
            continue
        # 2. the published centre must attain that inner maximum, which is the matching lower bound
        ws = topr.exhibit_weights([int(c) for c in rec['attaining_centre']], N, m)
        attained = sum(max(F(0), w - lam) for w in ws)
        if attained != inner:
            fails.append('m=%d: the published centre attains %s, not the inner maximum %s, so the '
                         'certificate does not bracket the dual' % (m, attained, inner))
            continue
        # 3. the bound must clear the threshold, which is what proves admission
        if not bound < thr:
            fails.append('m=%d: bound %s is not below the threshold %s, so run count %d is NOT '
                         'admitted by this certificate' % (m, bound, thr, r))
            continue
        print('  [ok  ] m=%2d r=%2d  lam=%-22s bound %s < %s   margin %s = %.8f'
              % (m, r, lam, bound, thr, thr - bound, float(thr - bound)))

    if fails:
        print('\nGATE FAILED:')
        for f in fails:
            print('  ' + f)
        return 2
    print('\nCLEAN: all five admitted rows verify from their published witnesses, one dual '
          'evaluation each, no search.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
