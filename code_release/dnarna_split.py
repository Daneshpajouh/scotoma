import sys, csv
from collections import Counter
csv.field_size_limit(10**7)
rdr=csv.DictReader(sys.stdin)
cls=Counter(); dna_pos=Counter(); rna_pos=Counter()
_rows_seen = 0
for row in rdr:
    _rows_seen += 1
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
    if g_al.count("-"):
        cls["DNA bulge (gap in guide)"]+=1
        for i,c in enumerate(g_al):
            if c=="-": dna_pos[len(g_al[:i].replace("-",""))]+=1
    if t_al.count("-"):
        cls["RNA bulge (gap in target)"]+=1
        for i,c in enumerate(t_al):
            if c=="-": rna_pos[len(t_al[:i].replace("-",""))]+=1
# 2026-10-22: WITH NO STDIN THIS SCRIPT PROCESSED ZERO ROWS AND PRINTED A CONFIRMATION.
# Found by a hostile reproducibility review, and it is the same defect class that
# verify_missed_actives.py and verify_missed_activity.py were both repaired for earlier the same
# day. It mattered more here than there: this is one of only two scripts that would re-earn the
# deposit-derived constants, so a referee checking it with no input got "no counterexample" over
# nothing. Empty input is an invocation error, never a pass.
if _rows_seen == 0:
    sys.exit("GATE FAILED: no rows read from stdin, so nothing was checked. This is an invocation\n"
             "  error, not a finding. Correct use pipes the CHANGE-seq deposit in, e.g.\n"
             "  cat <path-to>/"
             "CHANGEseq_CR_Lazzarotto_2020_dataset.csv \\\n"
             "    | python3 code_release/dnarna_split.py")

print("bulge-bearing measured-active sites by class:")
for k,v in cls.most_common(): print(f"   {k:28s} {v:,}")
print()
print("DNA-bulge wildcard slot distribution (slot >= 20 would be the omitted PAM-proximal one):")
print("   ", dict(sorted(dna_pos.items())))
print()
print("RNA-bulge deleted-position distribution (0 and 19 are the omitted terminals):")
print("   ", dict(sorted(rna_pos.items())))
