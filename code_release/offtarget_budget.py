#!/usr/bin/env python3
"""Validation budget for an off-target screen: how many candidate sites must be tested.

Answers one question for your own guide. Given a candidate off-target set of N sites, and the
mechanistic assumption (A_k) that an unclean guide cleaves at least k of them, how many sites
must a validation campaign test so that it clears an unclean guide with probability at most
delta1?

    m*(N, k, delta1) = min { m : C(N-k, m) / C(N, m) <= delta1 }

Every number is an exact integer search in rational arithmetic, delegated to m0_search, which
is the authority for every sample size in the paper. Nothing here is estimated or fitted.

Two ways in:

  1. You already know N (from your own scan, or any nomination tool):

       python3 offtarget_budget.py --n 28175
       python3 offtarget_budget.py --n 28175 --k 1 10 25 --delta1 0.05

  2. You have a spacer and a reference genome, and want N computed first:

       python3 offtarget_budget.py --spacer CTAACAGTTGCTTTTATCAC --ref-dir /path/to/contigs

     This runs the same exhaustive scan as enumerate_sites.py (20-nt protospacer, NGG or NAG
     PAM, at most --max-mismatches substitutions, both strands, no bulges), excludes the
     zero-mismatch on-target window, and uses the remaining count as N.

Output is a table, and optionally JSON (--out-json) or a LaTeX tabular (--out-tex).

WHAT THE NUMBER MEANS, AND WHAT IT DOES NOT

  * delta1 is a repeated-sampling level: the probability of clearing an unclean guide is at
    most delta1. It is NOT the probability that a cleared guide is clean.
  * Without (A_k) the same row certifies only "fewer than k cleaved sites", never "none".
    (A_k) is what upgrades it, and (A_k) is an assumption about the enzyme and the assay that
    you must state, not something this tool or your data can supply.
  * (A_k) needs a FLOOR on the number of cleaved sites among unclean guides, not a typical
    value. A median is the wrong statistic here and will understate the budget.
  * The budget is conditional on the candidate set being the right search class. Sites outside
    it, from alignment bulges, from PAMs beyond the searched set, or created by variants absent
    from the reference, are not covered at any m.
  * It assumes the assay reports cleavage at every site it tests. With per-site sensitivity s,
    a guide with one cleaved site is cleared with probability 1 - s*m/N, so level delta1 needs
    s >= 1 - delta1; no m, including a census, repairs an assay below that.

Part of the code accompanying "When can a finite screen certify zero false alarms?".
"""
import argparse
import json
import os
import sys
from fractions import Fraction as F

# Work when run from any directory, and under `python3 -I`, which drops the script's own
# directory from sys.path.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from m0_search import as_fraction, m0_hypergeom, m0_iid

DEFAULT_K = (1, 2, 5, 10, 25)


def budget_rows(n_sites, k_values, delta1):
    """One row per k: the exact hypergeometric requirement and the i.i.d. comparison.

    n_sites is N, the number of candidate off-target sites with the on-target excluded.
    delta1 is passed through as_fraction, so give it as a string or a Fraction: '0.05' is
    exact, the float 0.05 is not 1/20 in binary and is refused upstream.
    """
    d = as_fraction(delta1)
    if not (0 < d < 1):
        raise ValueError(f"delta1 must lie in (0, 1), got {d}")
    if n_sites < 1:
        raise ValueError(f"N must be at least 1, got {n_sites}")
    rows = []
    for k in k_values:
        if not (1 <= k <= n_sites):
            raise ValueError(f"k must satisfy 1 <= k <= N = {n_sites}, got {k}")
        m_exact = m0_hypergeom(n_sites, k, d)
        rows.append({
            "k": k,
            "gamma": f"{k}/{n_sites}",
            "gamma_float": k / n_sites,
            "sites_to_validate": m_exact,
            "fraction_of_N": m_exact / n_sites,
            "iid_draws": m0_iid(F(k, n_sites), d),
        })
    return rows


def enumerate_n(spacer, ref_dir, max_mismatches, quiet=False):
    """Scan a reference for candidate sites and return (N, summary dict).

    Composed exactly as enumerate_sites.main does, so the N used here and the N behind the
    paper's worked example come from the same code path. enumerate_sites is imported lazily
    because --n needs neither a reference genome nor numpy.
    """
    import enumerate_sites as es

    spacer_enc = es.encode(spacer.encode())
    if spacer_enc.size != es.PROTOSPACER_LEN or (spacer_enc > 3).any():
        raise SystemExit(f"the spacer must be {es.PROTOSPACER_LEN} nt over ACGT")
    records, contigs = [], 0
    for name, raw in es.iter_contigs(ref_dir):
        contigs += 1
        records.extend(es.scan_encoded(es.encode(raw), spacer_enc, max_mismatches))
        if not quiet:
            print(f"[scan] {name:<12} windows so far: {len(records):>8,}",
                  file=sys.stderr, flush=True)
    if contigs == 0:
        raise SystemExit(
            f"no sequence files in {ref_dir!r}: expected one .fa, .fasta or .seq per contig"
        )
    summary = es.summarise(records, spacer, max_mismatches)
    summary["contigs_scanned"] = contigs
    if summary["on_target_windows"] == 0:
        print("warning: the scan found no zero-mismatch window, so this spacer does not occur "
              "exactly in the reference; N below is the full candidate count with nothing "
              "excluded.", file=sys.stderr)
    return summary["off_target_windows"], summary


def render_table(rows, n_sites, delta1):
    out = [
        f"Validation budget, N = {n_sites:,} candidate off-target sites, delta1 = {delta1}",
        "",
        f"  {'k':>6}  {'gamma >= k/N':>14}  {'sites to test':>14}  {'of N':>8}  {'i.i.d. draws':>13}",
    ]
    for r in rows:
        out.append(
            f"  {r['k']:>6}  {r['gamma_float']:>14.3e}  {r['sites_to_validate']:>14,}"
            f"  {100 * r['fraction_of_N']:>7.2f}%  {r['iid_draws']:>13,}"
        )
    out += [
        "",
        "Each row certifies \"fewer than k cleaved sites\" on its own; it certifies \"none\" only",
        "under (A_k), which asserts an unclean guide cleaves at least k enumerated sites.",
        "(A_k) needs a floor, not a typical value: a median understates the budget.",
    ]
    return "\n".join(out)


def render_tex(rows, n_sites, delta1):
    lines = [
        r"\begin{tabular}{@{}rrrrr@{}}",
        r"\toprule",
        r"$k$ & $\gamma\ge k/N$ & sites to validate & fraction of $N$ & i.i.d.\ draws \\",
        r"\midrule",
    ]
    for r in rows:
        g = r["gamma_float"]
        mant, exp = f"{g:.2e}".split("e")
        lines.append(
            f"{r['k']} & ${mant}\\times10^{{{int(exp)}}}$ & {r['sites_to_validate']:,} & "
            f"{100 * r['fraction_of_N']:.2f}\\% & {r['iid_draws']:,} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Exact validation budget for an off-target screen.",
        epilog="Give either --n, or both --spacer and --ref-dir.",
    )
    src = ap.add_argument_group("the candidate set")
    src.add_argument("--n", type=int, help="N, candidate off-target sites, on-target excluded")
    src.add_argument("--spacer", help="20-nt spacer, to compute N by scanning a reference")
    src.add_argument("--ref-dir", help="directory of one sequence file per contig")
    src.add_argument("--max-mismatches", type=int, default=6,
                     help="substitutions allowed in the scan (default 6)")
    ap.add_argument("--k", type=int, nargs="+", default=list(DEFAULT_K),
                    help=f"values of k in (A_k) (default {' '.join(map(str, DEFAULT_K))})")
    ap.add_argument("--delta1", default="0.05",
                    help="level, as an exact string or fraction (default 0.05)")
    ap.add_argument("--out-json", help="write the rows as JSON")
    ap.add_argument("--out-tex", help="write a LaTeX tabular")
    a = ap.parse_args(argv)

    if (a.n is None) == (a.spacer is None):
        ap.error("give exactly one of --n, or --spacer with --ref-dir")
    summary = None
    if a.n is not None:
        n_sites = a.n
    else:
        if not a.ref_dir:
            ap.error("--spacer requires --ref-dir")
        n_sites, summary = enumerate_n(a.spacer, a.ref_dir, a.max_mismatches)
        print(f"scan: {n_sites:,} candidate off-target sites "
              f"(on-target window excluded)\n", file=sys.stderr)

    try:
        rows = budget_rows(n_sites, sorted(set(a.k)), a.delta1)
    except ValueError as e:
        raise SystemExit(f"offtarget_budget: {e}")
    print(render_table(rows, n_sites, a.delta1))

    if a.out_json:
        payload = {"N": n_sites, "delta1": str(a.delta1), "rows": rows}
        if summary:
            payload["enumeration_summary"] = summary
        with open(a.out_json, "w") as fh:
            json.dump(payload, fh, indent=1, sort_keys=True)
            fh.write("\n")
    if a.out_tex:
        with open(a.out_tex, "w") as fh:
            fh.write(render_tex(rows, n_sites, a.delta1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
