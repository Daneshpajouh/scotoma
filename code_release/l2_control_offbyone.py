# Is the L2 test vacuous? Plant a one-position error in our pattern construction, the exact
# class of bug that was live in CRISPRme's indexing, and see whether the test notices.
import csv
csv.field_size_limit(10**7)
rows=[]
# a handful of realistic single-bulge alignments, built by hand so this runs locally
cases=[('ACGTACGTACGTACGTACGT','ACGTACGTACGTACGTACGT'),   # 0 mm
       ('ACGTACGTACGTACGTACGT','ACGTTCGTACGTACGTACGT')]   # 1 mm
def mism(a,b):
    c=0
    for x,y in zip(a,b):
        if x=='*': continue
        if x!=y: c+=1
    return c
print("%-8s %-10s %-10s %-10s"%("case","correct","off_by_+1","off_by_-1"))
bad_correct=bad_p1=bad_m1=0
for gi,(g,t) in enumerate(cases):
    for j in range(1,19):
        # plant an RNA bulge at guide position j: target gets a gap there
        ts=t[:j]+'-'+t[j+1:]
        gs=g
        nm=mism(g[:j]+g[j+1:], t[:j]+t[j+1:])
        # correct construction
        d0=mism(gs[:j]+gs[j+1:], ts[:j]+ts[j+1:])
        # off-by-one constructions, the bug class
        d1=mism(gs[:j+1]+gs[j+2:], ts[:j]+ts[j+1:]) if j+2<=len(gs) else None
        d2=mism(gs[:j-1]+gs[j:],   ts[:j]+ts[j+1:]) if j>=1 else None
        if d0>nm: bad_correct+=1
        if d1 is not None and d1>nm: bad_p1+=1
        if d2 is not None and d2>nm: bad_m1+=1
print("%-8s %-10d %-10d %-10d"%("totals",bad_correct,bad_p1,bad_m1))
print()
print("rows the test would FLAG (ours > row's own count):")
print("  correct construction  : %d"%bad_correct)
print("  construction off by +1: %d"%bad_p1)
print("  construction off by -1: %d"%bad_m1)
print()
print("VERDICT:", "NOT VACUOUS - an off-by-one in our construction is flagged"
      if (bad_p1>0 or bad_m1>0) and bad_correct==0 else
      "VACUOUS - the test cannot distinguish a correct construction from an off-by-one")
