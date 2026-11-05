#!/usr/bin/env python3
"""Generate the frontier manifest at the panel's own length from the certificate files themselves.

Every row is recomputed here, not copied from the file. For each cited cell this script loads the
certificate, recomputes its total mass and its maximum ball load over EVERY centre class with
verify_cert_independent.py (which shares nothing with the solver but the JSON), and reports the
centre-class count it scanned against the full count 3^AP x (L+1). A row is PROVED only if the
recomputed bound exceeds r - 1 AND the scan covered the full class set. A certificate whose scan
covered only the binary-prefix classes 2^AP x (L+1) is reported as such and not counted.

Usage: python3 make_frontier_manifest.py            prints the table
       python3 make_frontier_manifest.py --latex    prints a LaTeX tabular
"""
import json, sys, os, importlib.util
from fractions import Fraction
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("vc", os.path.join(HERE, "verify_cert_independent.py"))
vc = importlib.util.module_from_spec(spec); spec.loader.exec_module(vc)
D = os.path.join(HERE, "data", "floor_certificates")

# the cells the body states at n = 20 (library word length N = 19), top of each proved prefix
CITED = [(9, 12, "cert_19_9_12.json"), (10, 9, "cert_19_10_9.json"),
         (11, 6, "cert_19_11_6.json"), (12, 4, "cert_19_12_4.json")]

def row(m, r, fn):
    d = json.load(open(os.path.join(D, fn)))
    N = d["N"]; AP = r - 1; L = N - AP
    full = 3**AP * (L + 1); binary = 2**AP * (L + 1)
    if "support" in d:                              # solver schema: recompute independently
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            loads, _, ex = vc.main(os.path.join(D, fn))
        stated = Fraction(int(d["mass_numerator_over_Dx"]), int(d["Dx"]))
        mass = ex["mass"]; mx = ex["maxload"]; bound = ex["bound"]   # exact, recomputed
        scanned = loads.size if mass == stated else None              # a mass mismatch is not counted
    else:                                           # no support to recompute from: never counted
        mass = float(eval(str(d["exact_total_mass"]))); mx = float(eval(str(d.get("exact_max_ball_mass", 1))))
        bound = mass / mx; scanned = None
    cover = ("not recomputed" if scanned is None else "full" if scanned == full
             else "binary only" if scanned == binary else "partial")
    ok = bound > r - 1 and cover == "full"
    mass, mx, bound = float(mass), float(mx), float(bound)            # exact above, printed below
    return dict(N=N, m=m, r=r, mass=mass, maxload=mx, bound=bound, scanned=scanned, full=full,
                cover=cover, verdict="PROVED" if ok else "NOT ESTABLISHED", file=fn)

rows = [row(*c) for c in CITED]
if "--latex" in sys.argv:
    print(r"\begin{tabular}{rrrrrrrl}")
    print(r"\toprule")
    print(r"$N$ & $m$ & $r$ & mass & max load & bound & classes scanned & verdict \\")
    print(r"\midrule")
    for x in rows:
        print(f"{x['N']} & {x['m']} & {x['r']} & {x['mass']:.6f} & {x['maxload']:.6f} & "
              f"{x['bound']:.6f} & {(format(x['scanned'], ',') if x['scanned'] is not None else 'none')} of {x['full']:,} & {x['verdict']} \\\\".replace(",", "{,}"))
    print(r"\bottomrule")
    print(r"\end{tabular}")
else:
    print(f"{'N':>3} {'m':>3} {'r':>3} {'mass':>10} {'max load':>10} {'bound':>12} {'needs >':>8} {'classes scanned':>22}  verdict")
    for x in rows:
        print(f"{x['N']:3d} {x['m']:3d} {x['r']:3d} {x['mass']:10.6f} {x['maxload']:10.6f} {x['bound']:12.9f} "
              f"{x['r']-1:8d} {str(x['scanned'])+' of '+str(x['full']):>22}  {x['verdict']} ({x['cover']})")
    sys.exit(0 if all(x["verdict"] == "PROVED" for x in rows) else 1)
