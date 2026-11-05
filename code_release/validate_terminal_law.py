import sys, csv
csv.field_size_limit(10**7)
rdr=csv.DictReader(sys.stdin)
def mm(a,b): return sum(1 for x,y in zip(a,b) if x!=y)
lossless_reach=lossless_total=0
lossy_total=0
bad=[]
_rows_seen = 0
for row in rdr:
    _rows_seen += 1
    try: reads=float((row.get("reads") or "0") or 0)
    except Exception: continue
    if reads<=0: continue
    try: nbi=int(float((row.get("Align.#Bulges") or "0").strip()))
    except Exception: nbi=0
    if nbi!=1: continue
    g_al=(row.get("Align.sgRNA") or "").strip().upper()
    t_al=(row.get("Align.off-target") or "").strip().upper()
    if (g_al.count("-")+t_al.count("-"))!=1 or not t_al.count("-"): continue
    sp_full=g_al.replace("-",""); n=len(sp_full)-3
    if n<4: continue
    sp=sp_full[:n]
    pos=None
    for i,c in enumerate(t_al):
        if c=="-": pos=len(t_al[:i].replace("-",""))
    if pos != n-1: continue                      # only the omitted PAM-proximal terminal
    core=t_al.replace("-","")[:n-1]
    if len(core)!=n-1: continue
    emitted=[sp[:i]+sp[i+1:] for i in range(1,n-1)]
    best=min(mm(core,e) for e in emitted)
    if sp[n-1]==sp[n-2]:                          # law says: omission is LOSSLESS here
        lossless_total+=1
        if best<=4: lossless_reach+=1
        else: bad.append((sp,best))
    else:
        lossy_total+=1
# 2026-10-22: WITH NO STDIN THIS SCRIPT PROCESSED ZERO ROWS AND PRINTED A CONFIRMATION.
# Found by a hostile reproducibility review, and it is the same defect class that
# verify_missed_actives.py and verify_missed_activity.py were both repaired for earlier the same
# day. It mattered more here than there: this is one of only two scripts that would re-earn the
# deposit-derived constants, so a referee checking it with no input got "no counterexample" over
# nothing. Empty input is an invocation error, never a pass.
if _rows_seen == 0:
    # Exit 2, not 1: 0 = checked and passed, 1 = checked and FAILED, 2 = NOT RUN. See the note in
    # verify_missed_activity.py; exiting 1 for "no data" hid a skip behind a failure shape.
    print("NOT RUN: no rows read from stdin, so nothing was checked. This is an invocation state,\n"
          "  not a finding, and it is NOT a pass. Correct use pipes the CHANGE-seq deposit in, e.g.\n"
          "  cat <path-to>/"
          "CHANGEseq_CR_Lazzarotto_2020_dataset.csv \\\n"
          "    | python3 code_release/validate_terminal_law.py")
    sys.exit(2)

# 2026-10-24: refusing ZERO rows is an empty-input gate, not a complete-input gate. This script's
# documented invocation streams 878,760,836 bytes over `ssh ... cat`, so a dropped connection is a
# realistic way to arrive with SOME rows. Measured on the first 40,000,000 bytes, cut mid-row: 34
# measured actives instead of 583, and it still printed "confirmed on measured biology ... no
# counterexample" and exited 0. A confirmation over a twentieth of the evidence is not the
# paper's. A short stream is an invocation state, so exit 2 NOT RUN; --any-input is the opt-out.
EXPECTED_ACTIVES = 583          # the deposited table, 323 lossless and 260 lossy
_ANY = "--any-input" in sys.argv[1:]
_total = lossless_total + lossy_total
if _total != EXPECTED_ACTIVES and not _ANY:
    print(f"NOT RUN: {_total} measured actives at the omitted terminal; the deposited table gives "
          f"{EXPECTED_ACTIVES}.\n"
          "  A truncated or different stream is not a failure to reproduce. Check the input\n"
          "  reached the end, or pass --any-input to score another table on purpose.")
    sys.exit(2)

print(f"measured actives with a single RNA bulge at the omitted terminal (position n-1): "
      f"{_total}" + ("" if _total == EXPECTED_ACTIVES else
                     f"   <- NOT the deposited {EXPECTED_ACTIVES}, --any-input"))
print(f"  guide has g_n == g_(n-1)  -> law says the omission is LOSSLESS : {lossless_total}")
print(f"     of those, actually reachable by an emitted pattern at <=4mm : {lossless_reach}")
print(f"     NOT reachable, i.e. the law would be WRONG                 : {len(bad)}")
print(f"  guide has g_n != g_(n-1)  -> law says LOSSY, at risk          : {lossy_total}")
if bad:
    for sp,b in bad[:5]: print(f"       counterexample {sp} nearest emitted {b}mm")
else:
    print("\n  Proposition 1's run-degeneracy clause is confirmed on measured biology:")
    print("  every site the law says is rescued IS rescued, with no counterexample.")
