"""Gate: the appendix's bulge-path completeness figures must equal the recorded producer outputs.

Round 44 added to the appendix the comparator-free check of the bulge path: all 23,692 measured
active bulge rows reached (17,267 RNA, 6,425 DNA, none unreachable, none out of class), plus the
non-vacuity control (35 of 36 flagged when the construction is shifted either way, 0 of 36 correct).
Those were the only quoted figures in the paper that nothing recomputed.

Two producers, both shipped because neither carries third-party data:
  data/l2_complete_66954118.out   the cluster run over the deposit
  data/l2_control_offbyone.out    the planted-shift control, reproducible locally

A missing or altered output is a FAILURE, not a skip. The appendix's own scope sentences are guarded
too: without them the counts read as a completeness claim the check does not support.
"""
import os, re, sys, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))

def _resolve(name):
    rel = os.path.join(HERE, 'data', name)
    sib = os.path.join(HERE, '..', 'producers_2026-10-20', name)
    return rel if os.path.exists(rel) else sib

RUN  = _resolve('l2_complete_66954118.out')
CTL  = _resolve('l2_control_offbyone.out')
TEX  = os.path.join(HERE, '..', 'paper_v7', 'appendix_moved.tex')

RUN_SHA = '49ab59579ffcfd136e91dc8b1da493e8ee6f264e7a2d45b1417e1cc6c24442ba'
CTL_SHA = 'fa93ca0b5bbaf40cd41ac3ca0fd662b780e2a9a938421cc4b0caf01e375c1b27'

def parse_run(t):
    g = {}
    for key, pat in (('rows',   r'active_bulge_rows\s+(\d+)'),
                     ('rna',    r'reached_RNA\s+(\d+)'),
                     ('dna',    r'reached_DNA\s+(\d+)')):
        m = re.search(pat, t)
        if not m: raise AssertionError('run output lacks %s' % key)
        g[key] = int(m.group(1))
    m = re.search(r'reached:\s*(\d+)\s+UNREACHED:\s*(\d+)', t)
    if not m: raise AssertionError('run output lacks the reached/UNREACHED summary')
    g['reached'], g['unreached'] = int(m.group(1)), int(m.group(2))
    # out-of-class buckets must be absent entirely or zero
    g['outofclass'] = sum(int(x) for x in re.findall(r'out_of_class\s+(\d+)', t))
    return g

def parse_ctl(t):
    m = re.search(r'correct construction\s*:\s*(\d+)', t)
    p = re.search(r'off by \+1:\s*(\d+)', t)
    n = re.search(r'off by -1:\s*(\d+)', t)
    if not (m and p and n): raise AssertionError('control output lacks its three counts')
    return int(m.group(1)), int(p.group(1)), int(n.group(1))

def check(run_t, ctl_t, tex):
    msgs, bad = [], []
    g = parse_run(run_t)
    c_ok, c_p1, c_m1 = parse_ctl(ctl_t)

    if g['rna'] + g['dna'] != g['rows']:
        bad.append('run is internally inconsistent: %d + %d != %d'
                   % (g['rna'], g['dna'], g['rows']))
    if g['reached'] != g['rows'] or g['unreached'] != 0:
        bad.append('run reports %d reached and %d unreached against %d rows'
                   % (g['reached'], g['unreached'], g['rows']))
    if g['outofclass'] != 0:
        bad.append('run reports %d out-of-class rows; the appendix says none' % g['outofclass'])

    for val, what in ((g['rows'], 'row count'), (g['rna'], 'RNA count'), (g['dna'], 'DNA count')):
        shown = '{:,}'.format(val).replace(',', '{,}')
        if shown not in tex and str(val) not in tex:
            bad.append('the appendix does not quote the %s %s' % (what, shown))
        else:
            msgs.append('  [ok  ] %-30s %s' % (what, shown))

    if c_ok != 0:
        bad.append('the control flags %d rows on a CORRECT construction; it is miscalibrated' % c_ok)
    if c_p1 < 1 or c_m1 < 1:
        bad.append('the control flags %d and %d on shifted constructions; it cannot fail, so the '
                   'completeness check is unsupported' % (c_p1, c_m1))
    else:
        msgs.append('  [ok  ] %-30s %d of 36 correct, %d and %d shifted'
                    % ('non-vacuity control', c_ok, c_p1, c_m1))
    for shown in ('35', '36'):
        if shown not in tex:
            bad.append('the appendix no longer quotes the control figure %s' % shown)

    # the scope sentences are load-bearing: without them the counts overclaim
    for probe, why in ((r'substitute for a gapped aligner',
                        'the appendix dropped the statement that this does not replace a gapped '
                        'aligner, so the counts read as full completeness'),
                       (r'not vacuous',
                        'the appendix dropped the non-vacuity statement'),
                       (r'property of the deposit.{0,3}s\s*\n?\s*annotation',
                        'the appendix dropped the note that in-class-ness is the annotation\'s '
                        'property, not evidence such actives do not exist')):
        if not re.search(probe, tex):
            bad.append(why)
    if not bad:
        msgs.append('  [ok  ] %-30s present' % 'three scope statements')
    return msgs, bad

def main():
    # 2026-10-23: BUILD gate, not a release artefact. It compares the paper's quoted figures with
    # the producers, so it needs the paper source, which is deliberately not in the released tree.
    # From a clean copy it must say so in one legible line, exit non-zero, and NOT report a pass it
    # did not earn. The clean-copy ruling exists because a check that quietly stops checking is not
    # a check; this one stops loudly. Caught by running the release gates from an empty directory.
    if not os.path.exists(TEX):
        print('%s: BUILD-ONLY GATE, NOT RUN. The paper source is not part of the release, so the '
              'paper-to-producer comparison cannot be made from a clean copy. This is NOT a pass.'
              % os.path.basename(__file__))
        return 2

    for path, sha, what in ((RUN, RUN_SHA, 'run'), (CTL, CTL_SHA, 'control')):
        if not os.path.exists(path):
            print('verify_l2_bulge_completeness: FAIL, %s output missing at %s' % (what, path))
            return 1
        got = hashlib.sha256(open(path, 'rb').read()).hexdigest()
        if got != sha:
            print('verify_l2_bulge_completeness: FAIL, %s sha256 %s does not match pinned %s'
                  % (what, got[:16], sha[:16]))
            return 1
    run_t = open(RUN, errors='replace').read()
    ctl_t = open(CTL, errors='replace').read()
    tex   = open(TEX, errors='replace').read()
    msgs, bad = check(run_t, ctl_t, tex)

    # self-control: the comparison must fire on a planted change
    _, cbad = check(run_t.replace('23692', '23691'), ctl_t, tex)
    if not cbad:
        print('verify_l2_bulge_completeness: FAIL, the comparison does not fire on a planted row '
              'count change, so it tests nothing')
        return 1
    msgs.append('  [ok  ] %-30s fires on a planted change' % 'self-control')

    print('\n'.join(msgs))
    if bad:
        print('\nL2 BULGE-COMPLETENESS GATE FAILURES')
        for b in bad: print('  [FAIL] %s' % b)
        return 1
    return 0

if __name__ == '__main__':
    sys.exit(main())
