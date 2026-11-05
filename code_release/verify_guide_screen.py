import json
from collections import defaultdict
import os
import sys
# SITES129 points at the omitted-position site records. Default is repo-relative so this runs
# anywhere; it previously held an absolute session path, which both leaked machine detail into a
# released file and made the script unrunnable for any reader.
SITES = os.environ.get("SITES129", os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "sites129.json"))
if not os.path.exists(SITES):
    sys.stderr.write(f"NOT RUN: set SITES129 to the site records; not found at {SITES}\n")
    raise SystemExit(2)
sites=json.load(open(SITES))

# 2026-10-24: this is a FALSIFICATION test, so the evidence it runs over has to be pinned. It was
# not: neither the record count, the file's hash, nor the panel size. Measured on the first 130 of
# the 260 records, the realised count fell from 40 to 18 and the script still printed CONFIRMED and
# exited 0. Removing evidence can only remove counterexamples, so a lost half reads as the same
# confirmation, and that is the one failure a falsification test must not have.
# A different corpus is an invocation error, so exit 2 NOT RUN, not 1; --any-corpus is the
# deliberate opt-out and labels its own output.
EXPECTED_SITES = 260          # records in the deposited data/sites129.json
EXPECTED_PANEL = 110          # guides in the panel index
EXPECTED_HITS  = 40           # guides realised as carrying an exclusively-missed site
ANY_CORPUS = "--any-corpus" in sys.argv[1:]
_n_sites = len(sites) if isinstance(sites, list) else len(sites.get("sites", sites))
if _n_sites != EXPECTED_SITES and not ANY_CORPUS:
    sys.stderr.write(
        "NOT RUN: %d site records; the deposited file holds %d.\n"
        "A falsification test over fewer records cannot find the counterexamples it is looking\n"
        "for, so a short corpus is a different corpus, not a confirmation. Pass --any-corpus to\n"
        "run another one on purpose.\n" % (_n_sites, EXPECTED_SITES))
    raise SystemExit(2)

comp={"A":"T","T":"A","G":"C","C":"G","N":"N"}
def rc(s): return "".join(comp.get(c,"N") for c in reversed(s))
def mm(a,b): return sum(1 for x,y in zip(a,b) if x!=y)
def bf_best(sp,win):
    n=len(sp); best=99
    for seq in (win, rc(win)):
        for i in range(0, max(0,len(seq)-(n+3))+1):
            w=seq[i:i+n+3]
            if len(w)<n+3 or "N" in w: continue
            if not (w[n+1]=="G" and w[n+2]=="G"): continue
            best=min(best, mm(sp,w[:n]))
    return best

# the panel's 110 spacers
IDX = os.environ.get("BULGE_IDX_DNA",
                     os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "bulge_census_full_idx_dna.json"))
if not os.path.exists(IDX):
    sys.stderr.write(f"NOT RUN: set BULGE_IDX_DNA to the DNA census index; not found at {IDX}\n")
    raise SystemExit(2)
idx=json.load(open(IDX))
panel=sorted({r[1] for rows in idx.values() for r in rows})
atrisk={g for g in panel if g[-1]!=g[-2]}        # law: exposed iff g_n != g_(n-1)
safe  ={g for g in panel if g[-1]==g[-2]}        # law: this omission cannot cost them anything
print(f"panel {len(panel)} guides -> law says EXPOSED {len(atrisk)}, PROTECTED {len(safe)}")
print(f"  (one character comparison, no genome, no scan)")
print()
# which guides actually carry an exclusively-missed site?
excl=defaultdict(int)
for h in sites:
    sp=h["spacer"]; n=len(sp); core=h["target"][:n-1]
    if len(core)!=n-1: continue
    emitted=[sp[:i]+sp[i+1:] for i in range(1,n-1)]
    if min(mm(core,e) for e in emitted) > 4 and bf_best(sp,h["window"]) > 4:
        excl[sp]+=1
hit=set(excl)
print(f"guides actually carrying an exclusively-missed site: {len(hit)}")
print(f"  of those, law said EXPOSED   : {len(hit & atrisk)}")
print(f"  of those, law said PROTECTED : {len(hit & safe)}   <- any of these falsifies the law")
print()
print(f"law's exposed set that carried no such site: {len(atrisk - hit)} of {len(atrisk)}")
print(f"  so the comparison is a NECESSARY condition, not a sufficient one:")
print(f"  it rules {len(safe)} of {len(panel)} guides out entirely, and {len(hit)} of the remaining "
      f"{len(atrisk)} are realised on this panel.")
print()
if hit & safe:
    # 2026-10-22: this branch printed and then fell off the end of the file, so the script
    # announced the proposition was falsified and still exited 0. The only SystemExit calls were
    # for missing inputs, which meant it failed loudly when it could not run and silently when the
    # science broke, and an exit-code sweep could not tell a falsification from a pass.
    print("  FALSIFIED: protected guides with missed sites:", sorted(hit & safe)[:4])
    # 2026-10-24: this exited 2, which the release convention reads as NOT RUN, an input the
    # checker could not find. A falsified law is the opposite: the checker ran and the claim is
    # refuted, which is exit 1. As written, the one outcome this script exists to detect was
    # reported with the code meaning "could not run", so an exit-code sweep still could not tell
    # a falsification from a missing file. That is the same confusion the 2026-10-22 note above
    # was written to end.
    raise SystemExit(1)
else:
    print("  CONFIRMED: no guide the law protects carries an exclusively-missed site.")
    print(f"  A designer can therefore exclude this exposure for {len(safe)} of {len(panel)} guides "
          f"({100*len(safe)/len(panel):.0f}%) before any scan.")
    # A confirmation is only worth the evidence behind it. The corpus size is pinned at load; the
    # panel and the realised count are pinned here, because a confirmation reached over a smaller
    # panel or with fewer realised cases is a weaker statement wearing the same words.
    _why = []
    if len(panel) != EXPECTED_PANEL:
        _why.append("panel is %d guides, not %d" % (len(panel), EXPECTED_PANEL))
    if len(hit) != EXPECTED_HITS:
        _why.append("%d guides realised, not %d" % (len(hit), EXPECTED_HITS))
    if _why and not ANY_CORPUS:
        sys.stderr.write("NOT RUN: " + "; ".join(_why) + ".\n"
                         "The law was not falsified, but this is not the deposited evidence, so\n"
                         "the confirmation above is not the paper's. Pass --any-corpus to run\n"
                         "another panel on purpose.\n")
        raise SystemExit(2)
    if _why:
        print("  --any-corpus: " + "; ".join(_why) + ". This confirmation is NOT the paper's.")
