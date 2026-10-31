#!/usr/bin/env python3
"""Validation-budget figure and tables for the paper. Every number is computed here.

The requirement is the exact hypergeometric zero-defect sample size: the smallest m such
that drawing m of N sites without replacement misses all k bad sites with probability at
most delta, that is C(N-k, m) / C(N, m) <= delta. It is computed by exact integer
comparison through m0_search.m0_hypergeom, not by a floating-point formula. The i.i.d.
column is the sample size of Theorem 2(i) at gamma = k/N, from m0_search.m0_iid.

The familiar closed form 1 - delta**(1/k) is the leading-order expression for the same
quantity with the continuity correction dropped; it is reported as the amount by which it
overshoots the exact requirement.

N is the number of candidate OFF-TARGET sites: the genomic windows returned by the scan
(data/exacel_enumeration_summary.json) minus the zero-mismatch window, which is the
on-target site and is cleaved by design. The unit is the window (a site), not the distinct
23-mer sequence: cleavage happens at a site, and the site count is the conservative
denominator.

Outputs, next to this script: budget_curve.pdf, budget_table.json, table_main.tex,
table_appendix.tex (the last two are \\input by the paper and its appendix).
"""
import json
import math
import pathlib
import sys
from fractions import Fraction

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from m0_search import m0_hypergeom, m0_iid  # noqa: E402

SUMMARY = json.loads((HERE / "data" / "exacel_enumeration_summary.json").read_text())
N = SUMMARY["off_target_windows"]
DELTAS = ("0.05", "0.025", "0.01")
KS_TABLE = (1, 2, 5, 10, 25)
KMAX = 40


def exact_sites(k, delta):
    return m0_hypergeom(N, k, delta)


def closed_form_sites(k, delta):
    return math.ceil((1.0 - float(delta) ** (1.0 / k)) * N)


def sci(x):
    """3.55e-05 -> LaTeX $3.55\\times10^{-5}$."""
    mant, exp = f"{x:.2e}".split("e")
    return f"${mant}\\times10^{{{int(exp)}}}$"


def main():
    plt.rcParams.update({
        "font.size": 10, "axes.labelsize": 10, "legend.fontsize": 10,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
        "font.family": "serif", "mathtext.fontset": "cm", "pdf.fonttype": 42,
    })
    ks = list(range(1, KMAX + 1))
    fig, ax = plt.subplots(figsize=(4.2, 2.7))
    for d, style in zip(DELTAS, ["-", "--", ":"]):
        ax.plot(ks, [100 * exact_sites(k, d) / N for k in ks], style, color="black",
                linewidth=1.3, label=rf"$\delta_1={d}$")
    for k in KS_TABLE[:4]:
        ax.plot([k], [100 * exact_sites(k, "0.05") / N], "o", color="black", markersize=4)
    ax.set_xlabel(r"$k$, sites an unclean guide cleaves")
    ax.set_ylabel("sites to validate (%)")
    ax.set_xlim(0, KMAX)
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.25, linewidth=0.5)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(HERE / "budget_curve.pdf", metadata={"CreationDate": None, "ModDate": None})

    rows = []
    for d in DELTAS:
        for k in KS_TABLE:
            ex, cf = exact_sites(k, d), closed_form_sites(k, d)
            rows.append({"k": k, "delta1": d, "gamma_lower_bound": k / N,
                         "exact_sites": ex, "exact_fraction": ex / N,
                         "iid_draws": m0_iid(Fraction(k, N), d),
                         "closed_form_sites": cf, "overshoot": cf - ex})
    (HERE / "budget_table.json").write_text(json.dumps(
        {"N_sites": N, "windows_total": SUMMARY["windows_total"],
         "method": "exact hypergeometric zero-defect via m0_search.m0_hypergeom",
         "source": "data/exacel_enumeration_summary.json",
         "rows": rows}, indent=2) + "\n")

    main_rows = [r for r in rows if r["delta1"] == "0.05"]
    lines = [r"\begin{tabular}{@{}rrrrr@{}}", r"\toprule",
             r"$k$ & $\gamma\ge k/N$ & sites to validate & fraction of $N$ & i.i.d.\ draws \\",
             r"\midrule"]
    for r in main_rows:
        lines.append(f"{r['k']} & {sci(r['gamma_lower_bound'])} & {r['exact_sites']:,} & "
                     f"{100 * r['exact_fraction']:.2f}\\% & {r['iid_draws']:,} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (HERE / "table_main.tex").write_text("\n".join(lines) + "\n")

    lines = [r"\begin{tabular}{@{}rrrrrrr@{}}", r"\toprule",
             r"& \multicolumn{2}{c}{$\delta_1=0.05$} & \multicolumn{2}{c}{$\delta_1=0.025$} "
             r"& \multicolumn{2}{c}{$\delta_1=0.01$} \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
             r"$k$ & exact sites & overshoot & exact sites & overshoot & exact sites & overshoot \\",
             r"\midrule"]
    for k in KS_TABLE:
        cells = []
        for d in DELTAS:
            r = next(x for x in rows if x["k"] == k and x["delta1"] == d)
            cells += [f"{r['exact_sites']:,}", f"$+{r['overshoot']}$"]
        lines.append(f"{k} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (HERE / "table_appendix.tex").write_text("\n".join(lines) + "\n")

    print(f"N = {N:,} candidate off-target sites ({SUMMARY['windows_total']:,} windows "
          f"minus {SUMMARY['on_target_windows']} on-target)")
    for d in DELTAS:
        print(f"delta1={d}: " + "  ".join(
            f"k={k}: {100 * exact_sites(k, d) / N:.2f}% ({exact_sites(k, d):,}"
            f"{'' if closed_form_sites(k, d) == exact_sites(k, d) else f', closed form +{closed_form_sites(k, d) - exact_sites(k, d)}'})"
            for k in KS_TABLE))


if __name__ == "__main__":
    main()
