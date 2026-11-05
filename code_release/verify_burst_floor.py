# Does the floor theorem generalise to b-burst deletion balls, b >= 2?
# STRUCTURE (predicted, from the same 4-line derivation as b=1):
#   D_i deletes g_i..g_{i+b-1}, so D_i[p] = g_p for p<i and g_{p+b} for p>=i.
#   Hence d(D_i,D_j) = #{p in [i,j) : g_p != g_{p+b}}, patterns collide iff g_i = g_{i+b},
#   so DISTINCT patterns indexed by b-boundaries form an ISOMETRIC PATH, and every coordinate
#   carries at most the two letters {g_p, g_{p+b}} and switches once.
import json, itertools, os, sys
# The panel file is resolved relative to THIS script, or from argv[1], so the released copy
# runs for a reader. An absolute session path made it unrunnable for everyone but its author.
_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'data', 'panel_match.json')
_PANEL = sys.argv[1] if len(sys.argv) > 1 else _DEFAULT
if not os.path.exists(_PANEL):
    print('NOT RUN: match_changeseq.json not found. Pass its path as the first argument.')
    print('         looked for: %s' % os.path.normpath(_PANEL))
    raise SystemExit(2)          # 0 passed, 1 FAILED, 2 NOT RUN
d=json.load(open(_PANEL))
sp=[v['spacer'] for v in d.values()]
def bdels(g,b):
    out=[]
    for i in range(len(g)-b+1):
        w=g[:i]+g[i+b:]
        if w not in out: out.append(w)
    return out
def ham(a,b): return sum(1 for x,y in zip(a,b) if x!=y)

print('=== STRUCTURE, over the whole panel ===')
print('%-3s %-8s %-10s %-14s %-14s %s' % ('b','N=n-b','ball sizes','path violations','3+-letter cols','cols switching twice'))
for b in (1,2,3,4,5):
    pv=0; cols3=0; cols2=0; pairs=0; sizes=set()
    for g in sp:
        D=bdels(g,b); r=len(D); N=len(D[0]); sizes.add(r)
        for i in range(r):
            for j in range(i+1,r):
                pairs+=1
                if ham(D[i],D[j]) != j-i: pv+=1
        for p in range(N):
            col=[D[k][p] for k in range(r)]
            if len(set(col))>2: cols3+=1
            if sum(1 for q in range(len(col)-1) if col[q]!=col[q+1])>1: cols2+=1
    print('%-3d %-8d %-10s %-14d %-14d %d  (pairs %d)'
          % (b, 20-b, '%d..%d'%(min(sizes),max(sizes)), pv, cols3, cols2, pairs))

print()
print('=== the deficit-law ball size identity, independently ===')
ok=0; tot=0
for b in (1,2,3,4,5):
    for g in sp:
        tot+=1
        pred = 1 + sum(1 for i in range(len(g)-b) if g[i]!=g[i+b])
        if len(bdels(g,b))==pred: ok+=1
print('|B_b(g)| = 1 + #{i: g_i != g_{i+b}} holds for %d of %d guide-b pairs' % (ok,tot))

print()
print('=== CONSTRUCTION at the shorter pattern length: exhaustive support search ===')
def admissible(N, O, L, m, rmax):
    Os=list(O)
    for t in range(1, rmax):
        Si=set(Os); Sj={(o+t)%N for o in Os}; B={a%N for a in range(t)}
        if len(B - Si)==0: return False
        d=0
        for p in range(N):
            a,bb=p in Si,p in Sj
            if a and bb: d += 1 if L[p] != L[(p-t)%N] else 0
            elif a or bb: d+=1
            else: d += 1 if p in B else 0
        if d < 2*m+1: return False
    return True
def search(N,m,rmax):
    for O in itertools.combinations(range(N), m):
        for bits in itertools.product([2,3], repeat=m):
            L={O[k]:bits[k] for k in range(m)}
            if admissible(N,O,L,m,rmax): return O,bits
    return None
for b in (1,2,3):
    N=20-b
    for m in (4,5):
        if m>=N: continue
        hit=search(N,m,min(N+1,21))
        print('b=%d  N=%2d  m=%d : %s' % (b,N,m,
              ('support %s letters %s' % (list(hit[0]),list(hit[1]))) if hit else 'NO support of this form exists (exhaustive)'))
