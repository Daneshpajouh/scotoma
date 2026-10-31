"""Gate: the appendix's shipped-convention DNA figures must equal the recorded producer output.

Round 43 added three numbers to the appendix from `authors_arm.py`, the arm that scores at
CRISPRme's OWN convention (slice Bulge_Size from the front, then key within the trimmed string,
per PostProcess/new_simple_analysis.py read at source). They were the only quoted figures in the
paper that nothing recomputed, because round 42 measured them and never captured the output.

This reads the recorded run, parses its table, and compares every figure the appendix quotes from
it. A missing output is a FAILURE, not a skip: the lesson from the clean-copy ruling is that a
check which can quietly stop checking is not a check. The output carries no CRISPRme table, only
our own AUROCs, so unlike verify_cfd_gap_table.py it ships and needs no path argument.
"""
import os, re, sys, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
_REL = os.path.join(HERE, 'data', 'authors_arm_66950898.out')
_SIB = os.path.join(HERE, '..', 'producers_2026-10-20', 'authors_arm_66950898.out')
OUT  = _REL if os.path.exists(_REL) else _SIB
TEX  = os.path.join(HERE, '..', 'paper_v7', 'appendix_moved.tex')

# pinned so a silently regenerated or edited output cannot pass
OUT_SHA = '64360cf30c8f2bad35905d53774a31adc53c5610d7a3b1e9cee0aef08090abb2'

def parse(text):
    """-> (rows_scored, fail_open, {edits: gain})"""
    m = re.search(r'DNA-bulge rows scored:\s*(\d+)\s+skipped:\s*(\d+)\s+'
                  r'authors-branch table misses passed over:\s*(\d+)', text)
    if not m:
        raise AssertionError('producer output does not carry its own header line')
    gains = {}
    for line in text.splitlines():
        f = line.split()
        if len(f) == 6 and f[0].isdigit():
            gains[int(f[0])] = float(f[5])
    if not gains:
        raise AssertionError('producer output carries no scored strata')
    return int(m.group(1)), int(m.group(3)), gains

def check(text, tex):
    msgs, bad = [], []
    rows, failopen, gains = parse(text)

    def want(pat, what):
        if not re.search(pat, tex):
            bad.append('the appendix no longer quotes %s; this gate guards nothing for it' % what)
            return False
        return True

    # row count
    if want(r'331\{,\}995', 'the 331,995 row count'):
        if rows != 331995:
            bad.append('producer scored %d rows, appendix says 331,995' % rows)
        else:
            msgs.append('  [ok  ] %-34s %d' % ('rows scored', rows))
    # fail-open count
    if want(r'8\{,\}241', 'the 8,241 fail-open rows'):
        if failopen != 8241:
            bad.append('producer passed over %d table misses, appendix says 8,241' % failopen)
        else:
            msgs.append('  [ok  ] %-34s %d' % ('fail-open rows', failopen))
    # the three gains, as the appendix states them
    for edits, shown, sign in ((3, '0.0693', -1), (4, '0.0973', -1), (5, '0.0285', +1)):
        g = gains.get(edits)
        if g is None:
            bad.append('producer has no stratum at %d edits' % edits); continue
        if abs(abs(g) - float(shown)) > 5e-5:
            bad.append('at %d edits producer gain %+0.4f does not match the quoted %s'
                       % (edits, g, shown))
        elif (g < 0) != (sign < 0):
            bad.append('at %d edits producer gain %+0.4f has the opposite sign to the appendix'
                       % (edits, g))
        elif not re.search(re.escape(shown).replace(r'\.', r'\.'), tex):
            bad.append('the appendix no longer quotes %s at %d edits' % (shown, edits))
        else:
            msgs.append('  [ok  ] %-34s %+0.4f' % ('gain at %d edits' % edits, g))

    # the paper must still say these are sign-only, or it overclaims an untested magnitude
    if not re.search(r'not interval-tested', tex):
        bad.append('the appendix dropped the "not interval-tested" qualifier on these three '
                   'figures; the arm prints no DeLong interval, so the magnitude is unsupported')
    else:
        msgs.append('  [ok  ] %-34s present' % 'sign-only qualifier')
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

    if not os.path.exists(OUT):
        print('verify_authors_convention: FAIL, recorded producer output missing at %s' % OUT)
        return 1
    raw = open(OUT, 'rb').read()
    got = hashlib.sha256(raw).hexdigest()
    if got != OUT_SHA:
        print('verify_authors_convention: FAIL, output sha256 %s does not match the pinned %s'
              % (got[:16], OUT_SHA[:16]))
        return 1
    text = raw.decode('utf-8', 'replace')
    tex = open(TEX, errors='replace').read()
    msgs, bad = check(text, tex)

    # control: the comparison must be able to fail
    hurt = text.replace('331995', '331996')
    _, cbad = check(hurt, tex)
    if not cbad:
        print('verify_authors_convention: FAIL, the comparison does not fire on a planted '
              'row-count change, so it is not testing anything')
        return 1
    msgs.append('  [ok  ] %-34s fires on a planted change' % 'control')

    print('\n'.join(msgs))
    if bad:
        print('\nAUTHORS-CONVENTION GATE FAILURES')
        for b in bad: print('  [FAIL] %s' % b)
        return 1
    return 0

if __name__ == '__main__':
    sys.exit(main())
