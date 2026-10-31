"""Gate: every r=20 and pooled m=4 number the paper may state, recomputed from the released artifacts.

No number enters the paper that the build cannot recompute. The two exhaustive maxima (over 3^19
centre classes for r=20, over 4^19 centres for the pooled certificate) are too slow for a build gate
and are cited in the text as recorded runs with their hashes. Everything else is recomputed here,
including one leg that was expected to need a recorded run and does not: that no radius-4 ball holds
three of the pooled witnesses is decided EXACTLY by a dynamic program over coordinates, because three
words admit a common radius-4 centre iff some assignment keeps all three partial distances at most 4.

Artifacts, by hash, under certificates/:
  r20_m8/payload.bin            2ab76785...  564 templates, 46,100 words
  pooled_m4_1629/P_1629.txt     b749f086...  the pooled library
  pooled_m4_1629/W_3258.txt     7213bb41...  two private witnesses per library word
"""
import hashlib, itertools, os, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
# Prefer the copy INSIDE the release. This gate used to look only at os.path.dirname(HERE), the
# PARENT of the tree that ships, so a reader who cloned code_release got a gate pointing at a
# directory that was never in the download. The parent is kept as a fallback for running inside the
# working project, and CERT_DIR overrides both.
_CERT_INTREE = os.path.join(HERE, 'certificates')
_CERT_PARENT = os.path.join(os.path.dirname(HERE), 'certificates')
CERT = os.environ.get('CERT_DIR',
                      _CERT_INTREE if os.path.isdir(_CERT_INTREE) else _CERT_PARENT)
bad = []

def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

# ---------------------------------------------------------------- r = 20, m = 8
N, M = 19, 8
pay = os.path.join(CERT, 'r20_m8', 'payload.bin')
if not os.path.exists(pay):
    bad.append('r20 payload missing at %s' % pay)
else:
    print('  payload.bin sha256 %s...' % sha(pay)[:8])
    raw = open(pay, 'rb').read()
    recs = [(int.from_bytes(raw[i:i+4], 'little'), raw[i+4]) for i in range(0, len(raw), 5)]
    STAIR = ['1'*i + '0'*(N-i) for i in range(20)]
    words = 0; massnum = 0; bad_dist = 0
    for idx, b in recs:
        t = []
        v = idx
        for _ in range(N):
            t.append(v % 3); v //= 3
        if v: bad.append('r20: template index overflows %d base-3 digits' % N); break
        stars = [j for j in range(N) if t[j] == 2]
        k = len(stars)
        words += 1 << k
        massnum += b                                  # template mass is b/1573
        w0 = ''.join('2' if x == 2 else str(x) for x in t)
        if min(sum(a != c for a, c in zip(w0, s)) for s in STAIR) != M: bad_dist += 1
    mass = Fr(massnum, 1573)
    print('  r=20: %d templates, %d support words, mass %s = %.9f, templates off distance %d: %d'
          % (len(recs), words, mass, float(mass), M, bad_dist))
    if len(recs) != 564: bad.append('r20: %d templates, expected 564' % len(recs))
    if words != 46100: bad.append('r20: %d support words, expected 46,100' % words)
    if mass != Fr(30348, 1573): bad.append('r20: mass %s, paper states 30348/1573' % mass)
    if bad_dist: bad.append('r20: %d templates not at distance exactly %d from the staircase' % (bad_dist, M))

# ---------------------------------------------------------------- pooled, m = 4
P4 = 4
pdir = os.path.join(CERT, 'pooled_m4_1629')
if not os.path.isdir(pdir):
    bad.append('pooled certificate missing at %s' % pdir)
else:
    P = [l.strip() for l in open(os.path.join(pdir, 'P_1629.txt')) if l.strip()]
    W = [l.strip() for l in open(os.path.join(pdir, 'W_3258.txt')) if l.strip()]
    print('  P_1629 sha256 %s...   W_3258 sha256 %s...'
          % (sha(os.path.join(pdir, 'P_1629.txt'))[:8], sha(os.path.join(pdir, 'W_3258.txt'))[:8]))
    print('  pooled: library %d words, witnesses %d, mass %s'
          % (len(P), len(W), Fr(len(W), 2)))
    if len(P) != 1629: bad.append('pooled: library has %d words, paper states 1,629' % len(P))
    if len(W) != 3258: bad.append('pooled: %d witnesses, paper states 3,258' % len(W))
    if Fr(len(W), 2) != 1629: bad.append('pooled: mass %s, paper states 1,629' % Fr(len(W), 2))
    if len(set(P)) != len(P): bad.append('pooled: library words are not distinct')

    import numpy as np
    Pa = np.frombuffer(''.join(P).encode(), dtype='S1').reshape(len(P), N)
    Wa = np.frombuffer(''.join(W).encode(), dtype='S1').reshape(len(W), N)
    owner = []
    notpriv = 0
    for i in range(0, len(W), 256):
        blk = Wa[i:i+256]
        D = (blk[:, None, :] != Pa[None, :, :]).sum(axis=2)      # (blk, 1629)
        mn = D.min(axis=1); cnt = (D <= P4).sum(axis=1)
        owner.extend(D.argmin(axis=1).tolist())
        notpriv += int(((mn != P4) | (cnt != 1)).sum())
    print('  pooled privacy: witnesses at distance exactly %d from one library word and >%d from all'
          ' others: %d of %d' % (P4, P4, len(W) - notpriv, len(W)))
    if notpriv: bad.append('pooled: %d witnesses are not private at radius %d' % (notpriv, P4))
    from collections import Counter
    per = Counter(owner)
    if len(per) != len(P) or set(per.values()) != {2}:
        bad.append('pooled: witnesses are not exactly two per library word (%d words covered, counts %s)'
                   % (len(per), sorted(set(per.values()))[:4]))
    else:
        print('  pooled: exactly two witnesses per library word, all %d words represented' % len(per))

    # EXACT, not a recorded run: no radius-4 ball holds three witnesses. Three words share a ball
    # only if pairwise within 2*P4; for each such triple a DP over coordinates decides whether any
    # centre keeps all three distances <= P4.
    Wb = [w.encode() for w in W]
    close = [[] for _ in W]
    for i in range(0, len(W), 256):
        blk = Wa[i:i+256]
        D = (blk[:, None, :] != Wa[None, :, :]).sum(axis=2)
        for a in range(blk.shape[0]):
            close[i+a] = [j for j in np.nonzero(D[a] <= 2*P4)[0].tolist() if j > i+a]
    triples = 0; holds = 0
    LET = b'ACGT'
    for a in range(len(W)):
        ca = set(close[a])
        for b_ in close[a]:
            for c_ in close[b_]:
                if c_ not in ca: continue
                triples += 1
                st = {(0, 0, 0)}
                for j in range(N):
                    x, y, z = Wb[a][j], Wb[b_][j], Wb[c_][j]
                    nxt = set()
                    for L in LET:
                        dx, dy, dz = L != x, L != y, L != z
                        for (p, q, r) in st:
                            u, v, w_ = p+dx, q+dy, r+dz
                            if u <= P4 and v <= P4 and w_ <= P4: nxt.add((u, v, w_))
                    st = nxt
                    if not st: break
                if st: holds += 1
    print('  pooled: %d witness triples pairwise within %d; radius-%d balls holding all three: %d'
          % (triples, 2*P4, P4, holds))
    if holds: bad.append('pooled: %d triples DO share a radius-%d ball, so no ball collects at most'
                         ' one unit and the certificate fails' % (holds, P4))

# ---------------------------------------------------------------- control
# A gate that cannot fail proves nothing: corrupt one witness and the privacy check must catch it.
if not bad and os.path.isdir(pdir):
    import numpy as np
    Wc = list(W); Wc[0] = P[0]                       # a library word is never private to itself at distance 4
    Wca = np.frombuffer(''.join(Wc).encode(), dtype='S1').reshape(len(Wc), N)
    D = (Wca[:1][:, None, :] != Pa[None, :, :]).sum(axis=2)
    caught = not (D.min() == P4 and (D <= P4).sum() == 1)
    print('  control: a corrupted witness is rejected: %s' % caught)
    if not caught: bad.append('CONTROL FAILED: the privacy check accepts a corrupted witness')

print()
if bad:
    # See verify_maize_exclusive.py: a missing certificate and a refuted certificate are different
    # outcomes and shared one banner. 2 is NOT RUN, 1 is FAILED.
    absent = [b for b in bad if 'missing at' in b]
    if len(absent) == len(bad):
        print('GATE NOT RUN: %d certificate input(s) are not here, so nothing could be checked.'
              % len(absent))
        for b in absent: print('  [absent] %s' % b)
        print('This is an incomplete release, not a finding about the claim.')
        sys.exit(2)
    print('FLOOR NEW RESULTS GATE FAILURES:')
    for b in bad: print('  [FAIL] %s' % b)
    sys.exit(1)
print('every r=20 and pooled m=4 figure in the paper is recomputed from the released certificates')
