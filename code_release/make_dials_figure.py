"""Producer for the five-dials figure: per-guide distribution of every choice that sets N.

Replaces the orphaned `figures/levers.pdf`, which had no producer anywhere in the repository and
plotted three choices, with a reproducible figure over the five dials the paper actually argues
after the two-base bulge term was priced.

Every input is a held per-guide producer in `producers_2026-10-20/`. The script verifies each dial's
median against the figure the paper states and refuses to write the figure if any disagrees, so a
silent data swap cannot produce a plausible-looking plot.

Run from the RECOMB_2027 directory:  python3 -I code_release/make_dials_figure.py
"""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 2026-10-22: this was a bare relative directory, resolved against the WORKING DIRECTORY,
# so the script passed from the paper directory and raised FileNotFoundError when run from
# its own folder, which is where the README runs everything else. Four scripts shared the
# defect and the release did not carry the inputs at all, so a referee could not run them
# from anywhere. The release now ships them under data/ and the path resolves against this
# file, preferring the in-release copy and falling back to the paper tree for development.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REL = os.path.join(_HERE, "data", "producers_2026-10-20")
D = _REL if os.path.isdir(_REL) else os.path.join(os.path.dirname(_HERE), "producers_2026-10-20")
OUT = "paper_v7/figures/dials.pdf"

def med(xs):
    xs = sorted(xs); n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2

G = {r["guide_index"]: r["sites"] for r in json.load(open(f"{D}/PANELREF_GRCh38_NGG.json"))["rows"]}
C = {r["guide_index"]: r["sites"] for r in json.load(open(f"{D}/PANELREF_CHM13_NGG.json"))["rows"]}
S1 = {r["guide_index"]: r for r in json.load(open(f"{D}/BULGE_CENSUS_size1_FULL_mm4.json"))["per_guide"]}
S2 = {r["guide_index"]: r for r in json.load(open(f"{D}/BULGE_CENSUS_size2_mm4.json"))["per_guide"]}
PM = {r["guide_index"]: r for r in json.load(open(f"{D}/PAM_PANEL_FACTOR_PERGUIDE.json"))["4"]["per_guide"]}
guides = sorted(set(G) & set(C) & set(S1) & set(S2) & set(PM))
print(f"guides with all five dials: {len(guides)}")

dials = [
    ("reference assembly\nGRCh38 to CHM13", [C[i]["4"] / G[i]["4"] for i in guides if G[i]["4"]], 1.0057),
    ("PAM set\nNGG to NGG+NAG",             [PM[i]["factor"] for i in guides],                    1.9898),
    ("bulge policy\none base",              [S1[i]["factor_merged"] for i in guides],             40.1),
    ("mismatch threshold\nfour to six",     [G[i]["6"] / G[i]["4"] for i in guides if G[i]["4"]], 69.9931),
    ("bulge policy\nup to two bases",       [S2[i]["factor_merged"] for i in guides],             98.44),
]

ok = True
for name, vals, expect in dials:
    m = med(vals)
    good = abs(m - expect) / expect < 0.01
    print(f"  {name.replace(chr(10),' '):42} median {m:9.4f}x  expected {expect:9.4f}x  "
          f"{'OK' if good else 'MISMATCH'}")
    ok = ok and good
if not ok:
    print("\nGATE FAILED: a dial median does not match the published figure. Figure not written.")
    sys.exit(2)
print("GATE PASSED: all five dial medians match the published figures to within 1%")

fig, ax = plt.subplots(figsize=(6.5, 3.8))
pos = list(range(1, len(dials) + 1))
ax.boxplot([d[1] for d in dials], positions=pos, widths=0.55,
           orientation="horizontal", showfliers=False,
           medianprops=dict(color="black", linewidth=1.4),
           boxprops=dict(color="black"), whiskerprops=dict(color="black"),
           capprops=dict(color="black"))
for p, (_, vals, _) in zip(pos, dials):
    ax.plot(vals, [p] * len(vals), marker=".", linestyle="none", markersize=2.5,
            color="0.45", alpha=0.55, zorder=1)
ax.axvline(1.0, linestyle=":", color="0.3", linewidth=1.0)
ax.set_xscale("log")
ax.set_xlabel(r"per-guide factor on $N$", fontsize=11)
ax.set_yticks(pos)
ax.set_yticklabels([d[0] for d in dials], fontsize=10.5)
ax.tick_params(axis="x", labelsize=10.5)
ax.invert_yaxis()
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT)
print(f"\nwrote {OUT}")
