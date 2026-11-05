#!/usr/bin/env python3
"""Reproduce the three bulge shares, 55.7% in vitro and 25.8% and 27.6% in cells, from the deposits.

WHAT IS CHECKED. The appendix section "Bulge policy, priced on both axes" reports shares of
measured active off-target sites that carry a bulge, taken at a matched budget of four. In vitro,
from the CHANGE-seq table:

  55.7% matched on mismatches: bulged active sites at <= 4 mismatches over all active sites at
        <= 4 mismatches, 23,692 of 42,572.
  45.4% matched on total edits: a bulge spends one edit, so bulged sites at <= 3 mismatches plus
        bulge-free sites at <= 4, 15,698 of 34,578.

It also says the published analysis of the same deposit counts the same 23,692 and reports it as
11.7% of 202,040. For this table the script recomputes every one of those numbers, the bulged rows
by mismatch count behind them (9, 246, 2,715, 12,728 and 7,994 at 0 to 4 mismatches), and the
1,749,342 single-bulge rows the deposit holds.

In cells, from the GUIDE-seq and NewGUIDE-seq tables, the same paragraph gives the share matched on
mismatches, on the same rule as the 55.7%: 312 of 1,207 (25.8%) and 147 of 533 (27.6%). For each
the script checks that share, the bulged rows at 0 to 4 mismatches, the row count and the
single-bulge rows.

THE ON-TARGET ROWS. An active row is a row with reads > 0. The zero-edit active rows (no bulge, no
mismatch) are the on-targets, one per guide that has a measured on-target: 109 rows over 109
guides in vitro, 58 in GUIDE-seq, 20 in NewGUIDE-seq. Every share excludes them. Kept in, the
shares move (55.5% and 45.3% in vitro, 24.7% and 26.6% in cells), so the exclusion is a condition
of reproducing the printed figures, and the script checks both versions. The cited 202,040 is every
active in vitro row, on-targets included.

THE BULGED ROWS AT ZERO MISMATCHES. Every share counts the active sites that carry a bulge and no
mismatch: 9 in vitro, 2 and 3 in cells. Dropped, the shares are 55.6%, 25.7% and 27.2%. The 23,692
the published analysis counts includes them, so they stay in, and the script checks that dropping
them moves each share.

THE SOURCE. The CHANGE-seq, GUIDE-seq and NewGUIDE-seq tables of Lazzarotto et al. (2020) as
redistributed with CRISPR-Bulge: datasets/<assay>/include_on_targets/ inside files/datasets.zip of
the OrensteinLab/CRISPR-Bulge repository on GitHub. The zip is stored in Git LFS (524,400,344
bytes, SHA-256 f892f70ba4ac3b05b03b2171b4ad38746630de08ad630e650f355dd61203eab0); the Zenodo record
10.5281/zenodo.10902445 holds only the LFS pointer, so fetch it with Git LFS. The script knows each
table by its SHA-256. The tables are third-party deposits and are not redistributed here.

THE CONTROLS. Three run every time. Keeping the on-targets must move every share. Dropping the
bulged rows at zero mismatches must move every share. And in vitro the 45.4% must pin its counting
rule: on that table, each neighbouring rule (bulged sites to 4 mismatches, or bulge-free sites to 3
or to 5) must give a count other than 15,698 of 34,578, or matching that count would not show which
rule was used.

Usage:  share_strata.py [--table {changeseq,guideseq,newguideseq}] DEPOSIT.csv
The table is named by the digest of what was read; --table names it instead, for a file that is
not one of the three. Reads the table once, streaming: about 20 seconds for the 879 MB CHANGE-seq
table and 20 MB of memory when last checked. A table can be piped straight out of the archive:
unzip -p datasets.zip <member> | share_strata.py /dev/stdin
Exit 0 only when every check above holds and every control fires; 1 otherwise; 2 if the file could
not be read or parsed as such a table (no digest is printed then).
"""
import argparse
import csv
import hashlib
import io
import sys
from collections import Counter

TABLES = {
    "changeseq": dict(
        name="CHANGE-seq include_on_targets table (in vitro)",
        sha256="68f3e1d7eb9e9e18e3ddc859752be8d156378a3e681b0fff329fe679366adc11",
        rows=4936278, single_bulge=1749342, on_targets=109,
        bulged={0: 9, 1: 246, 2: 2715, 3: 12728, 4: 7994},
        share_mm=(23692, 42572, "55.7"), kept_mm="55.5", dropped_mm=(23683, 42563, "55.6"),
        share_total=(15698, 34578, "45.4"), kept_total="45.3", cited=(202040, "11.7")),
    "guideseq": dict(
        name="GUIDE-seq include_on_targets table (cells)",
        sha256="a3398cba2ef452ec37ebebda38a44aa62615c66165a43699691ca210f467bcd2",
        rows=2120340, single_bulge=739570, on_targets=58,
        bulged={0: 2, 1: 47, 2: 139, 3: 95, 4: 29},
        share_mm=(312, 1207, "25.8"), kept_mm="24.7", dropped_mm=(310, 1205, "25.7")),
    "newguideseq": dict(
        name="NewGUIDE-seq include_on_targets table (cells)",
        sha256="ac9ca55a384a9af5ba37d667766e6723a7c211b5fa1d54875e450db23a13a23b",
        rows=1150709, single_bulge=438375, on_targets=20,
        bulged={0: 3, 1: 14, 2: 57, 3: 57, 4: 16},
        share_mm=(147, 533, "27.6"), kept_mm="26.6", dropped_mm=(144, 530, "27.2")),
}


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
    parser.add_argument("--table", choices=sorted(TABLES))
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

    key = args.table or next((k for k, t in TABLES.items() if t["sha256"] == digest), None)
    if key is None:
        print(f"  input is one of the three tables  NO   {digest[:8]} is none of "
              + ", ".join(t["sha256"][:8] for t in TABLES.values()) + "; name one with --table")
        print()
        print("VERDICT: FAIL")
        return 1
    t = TABLES[key]

    got = shares(tab)
    n_on = sum(on_targets.values())
    all_active = sum(tab.values()) + n_on

    # Control 1: put the on-targets back. Every share must move, and to the values given above.
    kept = Counter(tab)
    kept[0, 0] = n_on
    c1 = shares(kept)
    # Control 2: drop the bulged rows at zero mismatches. The share must move to the other rule's value.
    dropped = Counter({bk: v for bk, v in tab.items() if not (bk[0] >= 1 and bk[1] == 0)})
    c2 = shares(dropped)

    checks = [
        (f"input is the {t['name']}", digest == t["sha256"], digest[:8]),
        (f"{t['rows']:,} rows, {t['single_bulge']:,} of them single-bulge",
         rows == t["rows"] and bulge_values[1] == t["single_bulge"] and set(bulge_values) == {0, 1},
         f"{rows:,} rows, bulge counts {dict(sorted(bulge_values.items()))}"),
        (f"{t['on_targets']} zero-edit active rows, one per guide",
         n_on == len(on_targets) == t["on_targets"], f"{n_on} rows over {len(on_targets)} guides"),
        ("bulged active rows at 0 to 4 mismatches", got["bulged"] == t["bulged"],
         ", ".join(f"{got['bulged'][k]:,}" for k in range(5))),
        (f"matched on mismatches: {t['share_mm'][0]:,} of {t['share_mm'][1]:,}, {t['share_mm'][2]}%",
         got["mm"] == t["share_mm"], f"{got['mm'][0]:,} of {got['mm'][1]:,}, {got['mm'][2]}%"),
    ]
    if "share_total" in t:
        total, cited = t["share_total"], t["cited"]
        checks += [
            (f"matched on total edits: {total[0]:,} of {total[1]:,}, {total[2]}%",
             got["total"] == total, f"{got['total'][0]:,} of {got['total'][1]:,}, {got['total'][2]}%"),
            (f"cited: {t['share_mm'][0]:,} is {cited[1]}% of {cited[0]:,} active rows",
             (all_active, pct(got["mm"][0], all_active)) == cited,
             f"{pct(got['mm'][0], all_active)}% of {all_active:,}, on-targets included"),
            (f"control: on-targets kept give {t['kept_mm']}% and {t['kept_total']}%",
             (c1["mm"][2], c1["total"][2]) == (t["kept_mm"], t["kept_total"])
             and c1["mm"][2] != t["share_mm"][2] and c1["total"][2] != total[2],
             f"{c1['mm'][2]}% and {c1['total'][2]}%"),
        ]
        # The neighbouring rules for the total-edit share, as (bulged, bulge-free) mismatch limits.
        # Each must miss the printed count, or the match would not identify the rule.
        c3 = {rule: bulge_share(tab, *rule) for rule in ((4, 4), (3, 3), (3, 5))}
        checks.append((f"control: neighbouring rules miss {total[0]:,} of {total[1]:,}",
                       all(c[:2] != total[:2] for c in c3.values()),
                       ", ".join(f"({b},{f}) {c[0]:,} of {c[1]:,}" for (b, f), c in c3.items())))
    else:
        checks.append((f"control: on-targets kept give {t['kept_mm']}%",
                       c1["mm"][2] == t["kept_mm"] and c1["mm"][2] != t["share_mm"][2],
                       f"{c1['mm'][2]}%"))
    dm = t["dropped_mm"]
    checks.append((f"control: zero-mismatch bulged rows dropped give {dm[0]:,} of {dm[1]:,}, {dm[2]}%",
                   c2["mm"] == dm and c2["mm"][2] != t["share_mm"][2],
                   f"{c2['mm'][0]:,} of {c2['mm'][1]:,}, {c2['mm'][2]}%"))

    width = max(len(name) for name, _, _ in checks)
    for name, ok, detail in checks:
        print(f"  {name:<{width}}  {'ok ' if ok else 'NO '}  {detail}")

    good = all(ok for _, ok, _ in checks)
    print()
    print("VERDICT:", "PASS" if good else "FAIL")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
