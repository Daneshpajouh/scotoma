#!/usr/bin/env python3
"""Reproduce the in vitro bulge shares, 55.7% and 45.4%, from the CHANGE-seq measurement table.

WHAT IS CHECKED. The appendix section "Bulge policy, priced on both axes" reports two shares of
measured active off-target sites that carry a bulge, taken at a matched budget of four:

  55.7% matched on mismatches: bulged active sites at <= 4 mismatches over all active sites at
        <= 4 mismatches, 23,692 of 42,572.
  45.4% matched on total edits: a bulge spends one edit, so bulged sites at <= 3 mismatches plus
        bulge-free sites at <= 4, 15,698 of 34,578.

It also says the published analysis of the same deposit counts the same 23,692 and reports it as
11.7% of 202,040. This script recomputes every one of those numbers, the bulged rows by mismatch
count behind them (9, 246, 2,715, 12,728 and 7,994 at 0 to 4 mismatches), and the 1,749,342
single-bulge rows the deposit holds.

THE ON-TARGET ROWS. An active row is a row with reads > 0. The zero-edit active rows (no bulge, no
mismatch) are the on-targets, one per guide that has a measured on-target: 109 rows over 109
guides. Both shares exclude them. Kept in, the shares are 55.5% and 45.3%, so the exclusion is a
condition of reproducing the printed figures, and the script checks both versions. The cited
202,040 is every active row, on-targets included.

THE SOURCE. CHANGEseq_CR_Lazzarotto_2020_dataset.csv, the CHANGE-seq table of Lazzarotto et al.
(2020) as redistributed with CRISPR-Bulge: datasets/CHANGEseq/include_on_targets/ inside
files/datasets.zip of the OrensteinLab/CRISPR-Bulge repository on GitHub. The zip is stored in Git
LFS (524,400,344 bytes, SHA-256 f892f70ba4ac3b05b03b2171b4ad38746630de08ad630e650f355dd61203eab0);
the Zenodo record 10.5281/zenodo.10902445 holds only the LFS pointer, so fetch it with Git LFS. The
same zip holds GUIDE-seq tables with near-identical names; the include_on_targets CHANGE-seq table
is the one. It is a third-party deposit and is not redistributed here.

THE CONTROLS. Two run every time. Keeping the on-targets must move both shares off their printed
values. And the 45.4% must pin its counting rule: on this table, each neighbouring rule (bulged
sites to 4 mismatches, or bulge-free sites to 3 or to 5) must give a count other than 15,698 of
34,578, or matching that count would not show which rule was used.

Usage:  share_strata.py DEPOSIT.csv
Reads the 879 MB table once, streaming: about 20 seconds and 20 MB of memory when last checked.
Exit 0 only when every check above holds and both controls fire; 1 otherwise; 2 if the file could
not be read or parsed as this table (no digest is printed then).
"""
import argparse
import csv
import hashlib
import io
import sys
from collections import Counter

EXPECTED_SHA256 = "68f3e1d7eb9e9e18e3ddc859752be8d156378a3e681b0fff329fe679366adc11"
N_ROWS, N_SINGLE_BULGE, N_ON_TARGETS = 4936278, 1749342, 109
BULGED_BY_MISMATCH = {0: 9, 1: 246, 2: 2715, 3: 12728, 4: 7994}
SHARE_MM = (23692, 42572, "55.7")       # matched on mismatches
SHARE_TOTAL = (15698, 34578, "45.4")    # matched on total edits
WITH_ON_TARGETS = ("55.5", "45.3")
CITED = (202040, "11.7")


class HashingReader(io.RawIOBase):
    """Pass the file through unchanged while hashing every byte read."""

    def __init__(self, fh):
        self.fh, self.h = fh, hashlib.sha256()

    def readable(self):
        return True

    def readinto(self, buf):
        n = self.fh.readinto(buf)
        self.h.update(memoryview(buf)[:n])
        return n


def pct(num, den):
    return f"{100 * num / den:.1f}"


def bulge_share(tab, bulged_max, free_max):
    """Bulged rows to bulged_max mismatches over those plus bulge-free rows to free_max."""
    num = sum(v for (b, k), v in tab.items() if b >= 1 and k <= bulged_max)
    den = num + sum(v for (b, k), v in tab.items() if b == 0 and k <= free_max)
    return num, den, pct(num, den)


def shares(tab):
    """Both shares from active rows counted by (bulges, mismatches)."""
    bulged = Counter()
    for (b, k), v in tab.items():
        if b >= 1:
            bulged[k] += v
    return {"bulged": {k: bulged[k] for k in range(5)},
            "mm": bulge_share(tab, 4, 4), "total": bulge_share(tab, 3, 4)}


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("deposit")
    args = parser.parse_args(argv[1:])      # a usage error exits 2, as documented

    tab, on_targets, bulge_values, rows = Counter(), Counter(), Counter(), 0
    try:
        with open(args.deposit, "rb") as raw:
            reader = HashingReader(raw)
            text = io.TextIOWrapper(io.BufferedReader(reader), encoding="utf-8", newline="")
            for row in csv.DictReader(text):
                bulges, mismatches = int(row["Align.#Bulges"]), int(row["Align.#Mismatches"])
                bulge_values[bulges] += 1
                if float(row["reads"]) > 0:
                    if (bulges, mismatches) == (0, 0):
                        on_targets[row["sgRNA"]] += 1
                    else:
                        tab[bulges, mismatches] += 1
                rows += 1
            digest = reader.h.hexdigest()
    except OSError as exc:
        print(f"NOT RUN: {exc}")
        return 2
    except KeyError as exc:
        print(f"NOT RUN: {args.deposit} has no column {exc}")
        return 2
    except (ValueError, TypeError, csv.Error) as exc:     # not text, or a row that is not this table's
        print(f"NOT RUN: {args.deposit} at data row {rows + 1:,}: cannot parse: {exc}")
        return 2
    if not rows:
        print(f"NOT RUN: {args.deposit} has no rows")
        return 2

    print(f"deposit  {args.deposit}")
    print(f"         SHA-256 {digest}")
    print(f"active off-target rows by (bulges, mismatches): {dict(sorted(tab.items()))}")
    print()

    got = shares(tab)
    n_on = sum(on_targets.values())
    all_active = sum(tab.values()) + n_on

    # Control 1: put the on-targets back. Both shares must move, and to the values the paper gives.
    kept = Counter(tab)
    kept[0, 0] = n_on
    c1 = shares(kept)
    # Control 2: the neighbouring counting rules for the total-edit share, as (bulged, bulge-free)
    # mismatch limits. Each must miss the printed count, or the match would not identify the rule.
    c2 = {rule: bulge_share(tab, *rule) for rule in ((4, 4), (3, 3), (3, 5))}

    checks = [
        ("input is the CHANGE-seq include_on_targets table", digest == EXPECTED_SHA256,
         digest[:8]),
        (f"{N_ROWS:,} rows, {N_SINGLE_BULGE:,} of them single-bulge",
         rows == N_ROWS and bulge_values[1] == N_SINGLE_BULGE and set(bulge_values) == {0, 1},
         f"{rows:,} rows, bulge counts {dict(sorted(bulge_values.items()))}"),
        (f"{N_ON_TARGETS} zero-edit active rows, one per guide",
         n_on == len(on_targets) == N_ON_TARGETS, f"{n_on} rows over {len(on_targets)} guides"),
        ("bulged active rows at 0 to 4 mismatches", got["bulged"] == BULGED_BY_MISMATCH,
         ", ".join(f"{got['bulged'][k]:,}" for k in range(5))),
        (f"matched on mismatches: {SHARE_MM[0]:,} of {SHARE_MM[1]:,}, {SHARE_MM[2]}%",
         got["mm"] == SHARE_MM, f"{got['mm'][0]:,} of {got['mm'][1]:,}, {got['mm'][2]}%"),
        (f"matched on total edits: {SHARE_TOTAL[0]:,} of {SHARE_TOTAL[1]:,}, {SHARE_TOTAL[2]}%",
         got["total"] == SHARE_TOTAL,
         f"{got['total'][0]:,} of {got['total'][1]:,}, {got['total'][2]}%"),
        (f"cited: {SHARE_MM[0]:,} is {CITED[1]}% of {CITED[0]:,} active rows",
         (all_active, pct(got["mm"][0], all_active)) == CITED,
         f"{pct(got['mm'][0], all_active)}% of {all_active:,}, on-targets included"),
        (f"control: on-targets kept give {WITH_ON_TARGETS[0]}% and {WITH_ON_TARGETS[1]}%",
         (c1["mm"][2], c1["total"][2]) == WITH_ON_TARGETS
         and c1["mm"][2] != SHARE_MM[2] and c1["total"][2] != SHARE_TOTAL[2],
         f"{c1['mm'][2]}% and {c1['total'][2]}%"),
        (f"control: neighbouring rules miss {SHARE_TOTAL[0]:,} of {SHARE_TOTAL[1]:,}",
         all(c[:2] != SHARE_TOTAL[:2] for c in c2.values()),
         ", ".join(f"({b},{f}) {c[0]:,} of {c[1]:,}" for (b, f), c in c2.items())),
    ]

    width = max(len(name) for name, _, _ in checks)
    for name, ok, detail in checks:
        print(f"  {name:<{width}}  {'ok ' if ok else 'NO '}  {detail}")

    good = all(ok for _, ok, _ in checks)
    print()
    print("VERDICT:", "PASS" if good else "FAIL")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
