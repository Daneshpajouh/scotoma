#!/usr/bin/env python3
"""Derive two benchmark figures, and the counts they rest on, from saved outputs.

usage: derive_figures.py PAPER_DIR CRZ_DIR

PAPER_DIR holds the paper-guide SWOffinder outputs listed in paper_inputs.sha256:
the forward run's CSV and paper_keys.tsv (PAM keys from the plus-strand rows of a
forward and a reverse-complement run, milestone 3). CRZ_DIR holds crz_dna/ and
crz_rna/, the two CRISPRitz passes of scope A, as rerun_benchmark.sh (kit v5)
writes them; expected_crispritz.sha256 is that kit's list. Every input is hashed
first; a missing or different input stops the run (exit 1). Writes one TSV row
per figure to stdout: name, value, value with thousands commas, how it is computed.
"""
import collections, csv, hashlib, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FWD = "paper_fwd_ATTGAGATAGTGTGGGGAAGNGG.csv"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check(listing, root):
    bad = 0
    with open(os.path.join(HERE, listing)) as f:
        for line in f:
            want, name = line.split(None, 1)
            path = os.path.join(root, name.strip())
            if not os.path.isfile(path):
                print("MISSING %s" % path, file=sys.stderr); bad += 1
            elif sha256(path) != want:
                print("DIFFERENT %s" % path, file=sys.stderr); bad += 1
    return bad


def paper(root):
    fwd = list(csv.DictReader(open(os.path.join(root, FWD))))
    keys = list(csv.DictReader(open(os.path.join(root, "paper_keys.tsv")), delimiter="\t"))
    plus_rows = [r for r in fwd if r["Strand"] == "+"]
    minus_rows = [r for r in fwd if r["Strand"] == "-"]
    # PAM of a forward minus row = EndPosition - site length (milestone 3, strand_audit.py)
    fwd_minus = {(r["Chromosome"], int(r["EndPosition"]) - len(r["AlignedText"].replace("-", "")))
                 for r in minus_rows}
    key_plus = {(r["chrom"], int(r["pam_coord"])) for r in keys if r["strand"] == "+"}
    key_minus = {(r["chrom"], int(r["pam_coord"])) for r in keys if r["strand"] == "-"}
    return [
        ("swo_paper_rows_forward_run", len(fwd), "rows in the forward run CSV (the paper's count)"),
        ("swo_paper_rows_plus", len(plus_rows), "forward rows on the plus strand"),
        ("swo_paper_rows_minus", len(minus_rows), "forward rows on the minus strand"),
        ("swo_paper_pams_in_forward_rows", len(key_plus) + len(fwd_minus),
         "distinct PAMs among the forward rows: plus keys + distinct minus PAMs"),
        ("swo_paper_minus_pams_in_forward_rows", len(fwd_minus), "distinct minus PAMs among the forward rows"),
        ("swo_paper_minus_pams_rc_run", len(key_minus), "minus PAMs from the reverse-complement run"),
        ("swo_paper_pam_loci_both_strands", len(key_plus) + len(key_minus),
         "distinct PAM loci, both strands (plus keys + reverse-complement minus keys)"),
        ("swo_paper_extra_minus_rows", len(minus_rows) - len(fwd_minus),
         "forward minus rows beyond one per PAM"),
        ("swo_paper_minus_pams_lost", len(key_minus - fwd_minus),
         "reverse-complement minus PAMs with no forward row"),
        ("swo_paper_plus_rows_equal_keys", int({(r["Chromosome"], int(r["EndPosition"])) for r in plus_rows} == key_plus),
         "1 if the forward plus rows and the plus keys are the same set"),
    ]


def crispritz(root):
    out = []
    total = 0
    for half in ("crz_dna", "crz_rna"):
        kinds = collections.Counter()
        d = os.path.join(root, half)
        for name in sorted(os.listdir(d)):
            if name.endswith(".targets.txt"):
                with open(os.path.join(d, name)) as f:
                    for line in f:
                        if not line.startswith("#"):
                            kinds[line.split("\t", 1)[0]] += 1
        total += sum(kinds.values())
        out.append(("crz_%s_rows" % half[4:], sum(kinds.values()), "all rows in %s/*.targets.txt" % half))
        for k in sorted(kinds):
            out.append(("crz_%s_rows_bulge_type_%s" % (half[4:], k), kinds[k],
                        "rows in %s with bulge type %s (X = no bulge)" % (half, k)))
    out.append(("crz_rows_both_passes", total, "all rows, both passes"))
    return out


def main(argv):
    if len(argv) != 3:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    if check("paper_inputs.sha256", argv[1]) + check("expected_crispritz.sha256", argv[2]):
        return 1
    rows = paper(argv[1]) + crispritz(argv[2])
    v = dict((n, x) for n, x, _ in rows)
    ident = (v["swo_paper_pam_loci_both_strands"] + v["swo_paper_extra_minus_rows"]
             - v["swo_paper_minus_pams_lost"] == v["swo_paper_rows_forward_run"])
    rows.append(("identity_loci_plus_extra_minus_lost_equals_rows", int(ident),
                 "1 if pam_loci_both_strands + extra_minus_rows - minus_pams_lost = rows_forward_run"))
    print("name\tvalue\tprinted\thow")
    for n, x, how in rows:
        print("%s\t%d\t%s\t%s" % (n, x, format(x, ","), how))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
