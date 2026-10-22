import itertools
A="ACGT"
def runs(g):
    out=[];i=0
    while i<len(g):
        j=i
        while j+1<len(g) and g[j+1]==g[i]: j+=1
        out.append((i,j)); i=j+1
    return out
def patterns(g):
    # one distinct pattern per maximal run: delete one base from that run
    return [g[:s]+g[s+1:] for s,_e in runs(g)]
def ball(w,m,n):
    S=set()
    for k in range(m+1):
        for pos in itertools.combinations(range(n),k):
            for sub in itertools.product(A,repeat=k):
                if any(sub[t]==w[pos[t]] for t in range(k)): continue
                x=list(w)
                for t in range(k): x[pos[t]]=sub[t]
                S.add("".join(x))
    return S

g="ACGTA"; m=2
pats=patterns(g); n=len(pats[0])
words=["".join(p) for p in itertools.product(A,repeat=n)]
idx={w:i for i,w in enumerate(words)}
print("guide",g,"runs",len(runs(g)),"patterns",pats)
U=set()
for p in pats: U|=ball(p,m,n)
print("space",len(words)," |L| =",len(pats)," union size =",len(U))

Umask=0
for w in U: Umask|=1<<idx[w]
# every candidate centre anywhere in the space
cand=[]
for w in words:
    bm=0
    for x in ball(w,m,n): bm|=1<<idx[x]
    cand.append(bm)

def cover(k):
    """exact search for k centres whose balls contain U"""
    order=sorted(range(len(words)),key=lambda i:-bin(cand[i]&Umask).count("1"))
    best=None
    def rec(need,chosen,start):
        nonlocal best
        if best is not None: return
        if need==0:
            if chosen: best=list(chosen)
            return
        if not need: return
        # prune: remaining capacity
        rem=bin(need).count("1")
        if rem> (k-len(chosen))*67: return
        if len(chosen)==k: return
        # branch on the lowest uncovered word -- some chosen ball must cover it
        low=(need & -need).bit_length()-1
        w=words[low]
        for c in range(len(words)):
            if cand[c]>>low & 1:
                chosen.append(c); rec(need & ~cand[c],chosen,0); chosen.pop()
                if best is not None: return
    rec(Umask,[],0)
    return best

for k in (3,4,5):
    r=cover(k)
    print("k =",k,"->", "COVERED by "+str([words[i] for i in r]) if r else "no cover")
    if r:
        chk=set()
        for i in r: chk|=ball(words[i],m,n)
        print("   verify: union of chosen balls has",len(chk),"words; contains U:",U<=chk)
        break
