#!/usr/bin/env python3
"""Look up bare figures as whole tokens in data files and notes.

usage: figure_lookup.py FIGURES.txt KIT.tsv ROOT...

FIGURES.txt holds one figure per line, as printed (commas allowed), and no
other text. KIT.tsv is superseded_figures.tsv from the stale-figures kit.
Each ROOT is a file or folder; .md files count as NOTES, other text files
(.tsv .csv .txt .json, plain or .gz) count as DATA. A figure matches with or
without its thousands commas or space groups, and never inside a longer
number or word, so 978,103 does not match 51978103, chr978103 or a hash. Occurrence is not provenance: DATA means the number
is present in a deposited data file, not that the file derives it.

One TSV row per figure: figure, status, data hits, note hits, first places.
Status is STALE (old value in the kit), DATA, NOTE_ONLY or NONE.
Exit 0 on success, 2 on a bad input.
"""
import csv, gzip, os, re, sys

DATA_EXT = (".tsv", ".csv", ".txt", ".json")


SEPS = (" ", "\u2009", "\u202f")


def forms(fig):
    s = fig.replace(",", "")
    for sep in SEPS:
        s = s.replace(sep, "")
    ip, dot, fp = s.partition(".")
    g = "{:,}".format(int(ip))
    plain = {ip + dot + fp, g + dot + fp}
    spaced = {g.replace(",", sep) + dot + fp for sep in SEPS} - plain
    return plain, spaced


def pattern(fig):
    # Whole token: no letter, digit or underscore on either side, so a figure
    # inside a hash (fae46854a39) or a name (chr978103) is not a hit; printed
    # with or without thousands commas or space groups (202 040).
    edge_l, edge_r = r"(?<![\w.,])", r"(?![\w]|[.,]\d)"
    plain, spaced = forms(fig)
    alts = ["(?:%s)" % "|".join(re.escape(a) for a in sorted(plain))]
    if spaced:
        sp = "[ \u2009\u202f]"
        alts.append("(?<!\\d%s)(?:%s)(?!%s\\d)" % (sp, "|".join(re.escape(a) for a in sorted(spaced)), sp))
    return re.compile(edge_l + "(?:" + "|".join(alts) + ")" + edge_r)


def kind(path):
    name = path[:-3] if path.endswith(".gz") else path
    if name.endswith(".md"):
        return "NOTE"
    if name.endswith(DATA_EXT):
        return "DATA"
    return None


def files(roots, skip):
    for root in roots:
        if os.path.isfile(root):
            walk = [(os.path.dirname(root), [], [os.path.basename(root)])]
        else:
            walk = os.walk(root)
        for d, _, names in walk:
            for n in sorted(names):
                p = os.path.join(d, n)
                if kind(p) and os.path.realpath(p) not in skip:
                    yield p


def lines(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        return f.read().splitlines()


def main(argv):
    if len(argv) < 4:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    figs_path, kit_path, roots = argv[1], argv[2], argv[3:]
    with open(figs_path) as f:
        figs = [x.strip() for x in f if x.strip()]
    bad = [x for x in figs if not re.fullmatch(r"\d[\d,]*(\.\d+)?", x)]
    if bad:
        print("NOT A BARE FIGURE: %s" % bad[0], file=sys.stderr)
        return 2
    with open(kit_path, newline="") as f:
        stale = {r["old"].replace(",", ""): r["current"]
                 for r in csv.DictReader(f, delimiter="\t") if r["kind"] == "number"}
    skip = {os.path.realpath(p) for p in (figs_path, kit_path)}
    texts = []
    for p in files(roots, skip):
        try:
            texts.append((p, kind(p), lines(p)))
        except (OSError, UnicodeDecodeError, EOFError) as e:
            print("UNREADABLE %s: %s" % (p, e), file=sys.stderr)
            return 2
    print("figure\tstatus\tdata_hits\tnote_hits\tfirst_places")
    for fig in figs:
        pat = pattern(fig)
        hits = {"DATA": [], "NOTE": []}
        for p, k, ls in texts:
            for n, line in enumerate(ls, 1):
                if pat.search(line):
                    hits[k].append("%s:%d" % (p, n))
        if fig.replace(",", "") in stale:
            status = "STALE->" + stale[fig.replace(",", "")]
        elif hits["DATA"]:
            status = "DATA"
        elif hits["NOTE"]:
            status = "NOTE_ONLY"
        else:
            status = "NONE"
        places = (hits["DATA"] + hits["NOTE"])[:3]
        print("%s\t%s\t%d\t%d\t%s" % (fig, status, len(hits["DATA"]), len(hits["NOTE"]), "; ".join(places)))
    print("FILES %d (%d data, %d notes)" % (len(texts), sum(k == "DATA" for _, k, _ in texts),
                                           sum(k == "NOTE" for _, k, _ in texts)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
