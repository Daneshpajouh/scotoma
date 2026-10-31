"""My own independent recheck of W1's five cellular witnesses, run 2026-10-29.
Uses the repo's filed pattern producer for the libraries and recomputes every distance here.
Checks BOTH halves of the claim: the omitted pattern is the only route within four
substitutions, and no ungapped NGG register in the supplied window gets close."""
import json, sys, importlib.util, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("mps", ROOT / "code_release/make_pattern_sets.py")
mps = importlib.util.module_from_spec(spec); spec.loader.exec_module(mps)
rc = lambda s: s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
hd = lambda a, b: sum(1 for x, y in zip(a, b) if x != y) if len(a) == len(b) else 10**6
recs = json.load(open(pathlib.Path(__file__).with_name("Tracking2_verified_witnesses.json")))
bad = 0
for r in recs:
    g, core, win, claim = r["guide"], r["canonical_core"], r["window"], r["locus_min"]
    complete = mps.rna_words(g, 1)
    emitted = mps.sets_for(g, 1)["rna_emitted"]
    full = min(hd(core, p) for p in complete)
    wrap = min(hd(core, p) for p in emitted)
    nob = min((hd(g, s[i:i+23][:20]) for s in (win, rc(win))
               for i in range(len(s) - 22) if s[i:i+23][21:23] == "GG"), default=10**6)
    ok = (full == claim["full"] and wrap == claim["wrapper"] and nob == claim["no_bulge"])
    bad += not ok
    print(f"{core} full={full} wrapper={wrap} no_bulge={nob} "
          f"claim={claim['full']}/{claim['wrapper']}/{claim['no_bulge']} {'OK' if ok else 'MISMATCH'}")
print(f"\n{len(recs)-bad} of {len(recs)} witnesses reproduced; mismatches: {bad}")
sys.exit(1 if bad else 0)
