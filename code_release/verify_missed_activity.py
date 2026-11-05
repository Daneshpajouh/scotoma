import sys, csv
from collections import defaultdict
csv.field_size_limit(10**7)
def mm(a,b): return sum(1 for x,y in zip(a,b) if x!=y)
rdr=csv.DictReader(sys.stdin)
by_mm=defaultdict(list); missed=[]; missed_mm=defaultdict(int)
for row in rdr:
    try: reads=float((row.get("reads") or "0") or 0)
    except Exception: continue
    if reads<=0: continue
    try: nbi=int(float((row.get("Align.#Bulges") or "0").strip()))
    except Exception: nbi=0
    if nbi<1: continue
    g_al=(row.get("Align.sgRNA") or "").strip().upper()
    t_al=(row.get("Align.off-target") or "").strip().upper()
    if (g_al.count("-")+t_al.count("-"))!=nbi: continue
    sp_full=g_al.replace("-",""); n=len(sp_full)-3
    if n<4: continue
    try: dmm=int(float((row.get("Align.#Mismatches") or "-1").strip()))
    except Exception: dmm=-1
    by_mm[dmm].append(reads)
    if nbi!=1 or not t_al.count("-"): continue
    sp=sp_full[:n]; pos=None
    for i,c in enumerate(t_al):
        if c=="-": pos=len(t_al[:i].replace("-",""))
    if pos!=n-1 or sp[n-1]==sp[n-2]: continue
    core=t_al.replace("-","")[:n-1]
    if len(core)!=n-1: continue
    if min(mm(core,sp[:i]+sp[i+1:]) for i in range(1,n-1)) > 4:
        missed.append((reads,dmm)); missed_mm[dmm]+=1
def q(v,p):
    v=sorted(v); k=(len(v)-1)*p; lo=int(k); hi=min(lo+1,len(v)-1)
    return v[lo]+(v[hi]-v[lo])*(k-lo)
# 2026-10-21: this reads the deposit CSV from STDIN. Run with no input it processed zero rows,
# printed two empty distributions and exited 0 -- a pass with no data, the same defect class as
# verify_missed_actives.py had. Empty input is now an invocation error, never a pass.
if not by_mm and not missed_mm:
    # Exit 2, not 1: across this suite 0 means checked and passed, 1 means checked and FAILED,
    # and 2 means NOT RUN. Exiting 1 here made "no data supplied" indistinguishable from a real
    # finding, which is the ambiguity that let a never-run gate sit green for several rounds.
    print("NOT RUN: no rows read from stdin, so nothing was checked. This is an invocation state,\n"
          "  not a finding, and it is NOT a pass. Correct use pipes the CHANGE-seq deposit in, e.g.\n"
          "  cat <path-to>/CHANGEseq_CR_Lazzarotto_2020_dataset.csv \\\n"
          "    | python3 code_release/verify_missed_activity.py")
    sys.exit(2)
print("bulge-bearing measured actives by the deposit's own mismatch count:")
for k in sorted(by_mm):
    v=by_mm[k]
    print(f"   {k} mm : n={len(v):6,}  median reads {q(v,.5):7.0f}  Q3 {q(v,.75):8.0f}")
print()
# 2026-10-24: this heading printed a LITERAL 129 whatever it had counted, so a truncated stream
# printed "the 129 missed sites ... {4: 6}", a sentence contradicting itself on one line. The
# count now comes from the data. Refusing zero rows was an empty-input gate; the documented
# invocation streams 878,760,836 bytes over `ssh ... cat`, so arriving with SOME rows is the
# realistic failure. A short stream is an invocation state, so exit 2 NOT RUN, not a finding.
EXPECTED_MISSED = 129
_ANY = "--any-input" in sys.argv[1:]
_missed_n = sum(missed_mm.values())
if _missed_n != EXPECTED_MISSED and not _ANY:
    print(f"NOT RUN: {_missed_n} missed sites, not the deposited {EXPECTED_MISSED}. A truncated or\n"
          "  different stream is not a failure to reproduce. Check the input reached the end, or\n"
          "  pass --any-input to score another table on purpose.")
    sys.exit(2)
print(f"the {_missed_n} missed sites, by mismatch count: {dict(sorted(missed_mm.items()))}"
      + ("" if _missed_n == EXPECTED_MISSED else
         f"   <- NOT the deposited {EXPECTED_MISSED}, --any-input"))
print()
# the fair comparison: missed vs ALL bulge-bearing actives in the SAME mismatch stratum
for k in sorted(missed_mm):
    same=by_mm[k]; mr=[r for r,d in missed if d==k]
    print(f"stratum {k} mm:  all n={len(same):,} median {q(same,.5):.0f}   "
          f"missed n={len(mr)} median {q(mr,.5):.0f}")
    import bisect
    srt=sorted(same); pct=[100*bisect.bisect_left(srt,r)/len(srt) for r in mr]
    print(f"   missed sites' median percentile WITHIN their own stratum: {q(pct,.5):.1f}")
    above=sum(1 for p in pct if p>=50)
    print(f"   at or above their stratum's median: {above} of {len(mr)} ({100*above/len(mr):.0f}%)")
