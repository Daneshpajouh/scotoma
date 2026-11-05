#!/usr/bin/env python3
"""One command that replays every exact maximum and sweeps all four dual certificates.

The paper's availability statement promises exactly that: "one documented command replays every
exact maximum and sweeps all four dual certificates, each with its reference output". Until now it
took three separate commands, each documented and each working, so the sentence was true in
substance and not in letter. This is the one command.

Standard library only. Run from the repository root:

    python3 -I verification/replay_all.py

Exit 0 means every stage passed. Any other exit names the stage that did not, and no stage's failure
is silent: each one's own exit code is checked and reported.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FRAC = os.path.join(HERE, 'frac_alln')
DUALS = ['cert_m9_r9.json', 'cert_m10_r7.json', 'cert_m11_r6.json', 'cert_m12_r4.json']

STAGES = [
    ('the four dual certificates, every centre class',
     ['verify_ceiling_exhaustive.py'] + [os.path.join('data', 'floor_certificates', d) for d in DUALS]),
    ('the exact maxima at the certified lengths',
     ['frac_alln_tail_check.py', '--finite', '12', '29', '30', '31', '32']),
    ('the neighbour lemmas behind the weighting',
     ['frac_alln_tail_check.py', '--lemmas']),
    ('the m = 6 covering number by dynamic program',
     ['step3_dp.py']),
]


def main(stages=None):
    stages = STAGES if stages is None else stages
    width = max(len(n) for n, _ in stages)
    bad = []
    for name, argv in stages:
        p = subprocess.run([sys.executable, '-I'] + argv, cwd=FRAC,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        tail = [l for l in p.stdout.strip().splitlines() if l.strip()][-1:] or ['(no output)']
        print('%-*s  exit %d  %s' % (width, name, p.returncode, tail[0][:96]))
        if p.returncode != 0:
            bad.append((name, p.returncode, p.stdout))
    print()
    if bad:
        for name, rc, out in bad:
            print('FAILED: %s (exit %d)' % (name, rc))
            print(out[-1500:])
        return 1
    print('ALL STAGES PASSED. Every exact maximum replayed and all four duals swept exhaustively.')
    return 0


def selftest():
    """A planted failure must be caught. A runner that cannot fail certifies nothing."""
    planted = STAGES + [('PLANTED FAILURE, must be caught',
                         ['verify_ceiling_exhaustive.py', 'no_such_certificate.json'])]
    rc = main(planted)
    print()
    if rc == 0:
        print('SELFTEST FAILED: the planted stage was not caught, so this runner cannot fail.')
        return 1
    print('SELFTEST PASSED: the planted stage was caught and the runner returned %d.' % rc)
    return 0


if __name__ == '__main__':
    if '--selftest' in sys.argv[1:]:
        sys.exit(selftest())
    sys.exit(main())
