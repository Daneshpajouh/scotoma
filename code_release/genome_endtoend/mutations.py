#!/usr/bin/env python3
"""Build deliberately corrupted copies of this directory's inputs and run the checker on each.

    python3 -I mutations.py            # build, run, print one line per corruption

Every corruption below is invalid for a reason visible in the deposited files alone, with one
exception that is marked and expected to pass: changing one mismatching reference base to
another mismatching base cannot be detected without the reference, which this deposit does not
carry. A corruption that the checker accepts is a defect in the checker, so this script exits
non-zero if any unmarked corruption passes, or if the marked one fails.
"""
import csv, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
AL = "corrected_only_alignments.part000.tsv"
PC, PG, ACT = "per_chromosome.tsv", "per_guide.tsv", "active_sites.part000.tsv"


def load(d, n):
    with open(os.path.join(d, n), newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows, list(rows[0])


def save(d, n, rows, fields):
    with open(os.path.join(d, n), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def one_locus_repeated(d):
    rows, f = load(d, AL)
    save(d, AL, [dict(rows[0]) for _ in rows], f)


def coordinates_outside_assembly(d):
    rows, f = load(d, AL)
    for r in rows:
        r["start"] = str(int(r["start"]) + 2 ** 40)
        r["pam_coord"] = str(int(r["pam_coord"]) + 2 ** 40)
    save(d, AL, rows, f)


def minus_footprint_outside_record(d):
    rows, f = load(d, AL)
    r = next(x for x in rows if x["strand"] == "-")
    # place the start so that the 22 nt footprint ends one base past chr21
    r["chrom"], r["strand"], r["start"] = "chr21", "-", str(46709983 - 21)
    r["pam_coord"] = r["start"]
    save(d, AL, rows, f)


def nonexistent_chromosome(d):
    rows, f = load(d, AL)
    rows[0]["chrom"] = "chrDOES_NOT_EXIST"
    save(d, AL, rows, f)


def invalid_strand(d):
    rows, f = load(d, AL)
    rows[0]["strand"] = "?"
    save(d, AL, rows, f)


def overlong_crrna(d):
    rows, f = load(d, AL)
    rows[0]["crRNA"] = rows[0]["crRNA"] + "A"
    save(d, AL, rows, f)


def relocated_gap(d):
    rows, f = load(d, AL)
    r = rows[0]
    i = r["DNA"].find("-")
    j = 3 if i != 3 else 5
    s = list(r["DNA"].replace("-", ""))
    s.insert(j, "-")
    r["DNA"] = "".join(s)
    save(d, AL, rows, f)


def negative_chromosome_counts(d):
    rows, f = load(d, PC)
    rows[0]["corrected_only"] = str(-1)
    save(d, PC, rows, f)


def negative_guide_counts(d):
    rows, f = load(d, PG)
    rows[0]["corrected_only"] = str(-1)
    save(d, PG, rows, f)


def redistribute_chromosome_losses(d):
    rows, f = load(d, PC)
    a, b = rows[0], rows[1]
    a["corrected_only"] = str(int(a["corrected_only"]) + 5)
    b["corrected_only"] = str(int(b["corrected_only"]) - 5)
    save(d, PC, rows, f)


def redistribute_guide_losses(d):
    rows, f = load(d, PG)
    a = next(r for r in rows if int(r["corrected_only"]) > 0)
    b = next(r for r in rows if r is not a)
    a["corrected_only"] = str(int(a["corrected_only"]) - 1)
    b["corrected_only"] = str(int(b["corrected_only"]) + 1)
    save(d, PG, rows, f)


def per_record_count_algebra(d):
    rows, f = load(d, PC)
    a, b = rows[0], rows[1]
    a["corrected_rna1"] = str(int(a["corrected_rna1"]) + 1)
    b["corrected_rna1"] = str(int(b["corrected_rna1"]) - 1)
    save(d, PC, rows, f)


def per_guide_count_algebra(d):
    rows, f = load(d, PG)
    a, b = rows[0], rows[1]
    a["corrected_rna1"] = str(int(a["corrected_rna1"]) + 1)
    b["corrected_rna1"] = str(int(b["corrected_rna1"]) - 1)
    save(d, PG, rows, f)


def active_identity_destroyed(d):
    rows, f = load(d, ACT)
    for r in rows:
        r["guide"] = "NO_GUIDE"
        r["start"] = r["pam_coord"] = str(-(2 ** 40))
        r["deleted_position_0based"] = "0"
    save(d, ACT, rows, f)


def active_status_permutation(d):
    rows, f = load(d, ACT)
    a = next(r for r in rows if r["status"] == "VERIFIED_NGG")
    b = next((r for r in rows if r["status"] != "VERIFIED_NGG"), None)
    if b is None:
        a["status"] = "NON_NGG"
    else:
        a["status"], b["status"] = b["status"], a["status"]
    save(d, ACT, rows, f)


def active_return_bits_swapped(d):
    rows, f = load(d, ACT)
    a = next(r for r in rows if r["released_returned"] == "1")
    b = next(r for r in rows if r["released_returned"] == "0")
    a["released_returned"], b["released_returned"] = "0", "1"
    save(d, ACT, rows, f)


def renamed_record(d):
    rows, f = load(d, PC)
    rows[0]["chrom"] = "chr99"
    save(d, PC, rows, f)


def changed_difference_total(d):
    rows, f = load(d, "difference_classes.tsv")
    rows[0]["distinct_loci"] = str(int(rows[0]["distinct_loci"]) + 1)
    save(d, "difference_classes.tsv", rows, f)


def wrong_analysis_flag(d):
    rows, f = load(d, PC)
    r = next((x for x in rows if x["analysis_status"] != "PASS"), rows[0])
    r["analysis_status"] = "PASS"
    save(d, PC, rows, f)


def target_the_released_build_reaches(d):
    """A corrected-only target rewritten so an interior deletion word reaches it exactly, which
    contradicts the export's own definition: these are the loci only the terminal word finds."""
    rows, f = load(d, AL)
    pg, _ = load(d, PG)
    r = rows[0]
    spacer = next(x["spacer"] for x in pg if x["guide"] == r["guide"])
    interior = spacer[:5] + spacer[6:]          # the word that deletes spacer position 6
    dna = list(r["DNA"])
    gap = dna.index("-")
    r["DNA"] = interior[:gap] + "-" + interior[gap:] + r["DNA"].replace("-", "")[19:22]
    save(d, AL, rows, f)


def shift_assay_pam_boundary(d):
    """Every assay PAM coordinate moved by one, with all other identities untouched. The
    relation to start and strand is fixed, so this is detectable without the reference."""
    rows, f = load(d, ACT)
    for r in rows:
        r["pam_coord"] = str(int(r["pam_coord"]) + 1)
    save(d, ACT, rows, f)


def break_lost_assay_join(d):
    """One exported row for a claimed-lost assay locus moved by a single base, consistently in
    start and PAM coordinate. Every distance, domain and subgroup check still holds; what breaks
    is the join between the assay table and the export, which is the claim itself."""
    act, _ = load(d, ACT)
    lost = {(r["guide"], r["chrom"], r["strand"], r["pam_coord"]) for r in act
            if r["released_returned"] == "0" and r["corrected_returned"] == "1"}
    rows, f = load(d, AL)
    r = next(x for x in rows
             if (x["guide"], x["chrom"], x["strand"], x["pam_coord"]) in lost)
    r["start"] = str(int(r["start"]) + 1)
    r["pam_coord"] = str(int(r["pam_coord"]) + 1)
    save(d, AL, rows, f)


def invent_reference_base(d):
    """EXPECTED TO PASS. One mismatching base becomes another mismatching base. The counts,
    coordinates, PAM, gap and distances are all unchanged, and only the reference can say the
    new base is not what GRCh38 holds there."""
    rows, f = load(d, AL)
    r = rows[0]
    cr, dna = r["crRNA"], list(r["DNA"])
    for i, (a, b) in enumerate(zip(cr, r["DNA"])):
        if b != "-" and a != "N" and a != b:
            dna[i] = next(c for c in "ACGT" if c != b and c != a)
            break
    r["DNA"] = "".join(dna)
    save(d, AL, rows, f)


CASES = [
    (one_locus_repeated, False), (coordinates_outside_assembly, False),
    (minus_footprint_outside_record, False), (nonexistent_chromosome, False),
    (invalid_strand, False), (overlong_crrna, False), (relocated_gap, False),
    (negative_chromosome_counts, False), (negative_guide_counts, False),
    (redistribute_chromosome_losses, False), (redistribute_guide_losses, False),
    (per_record_count_algebra, False), (per_guide_count_algebra, False),
    (active_identity_destroyed, False), (active_status_permutation, False),
    (active_return_bits_swapped, False), (renamed_record, False),
    (changed_difference_total, False), (wrong_analysis_flag, False),
    (target_the_released_build_reaches, False), (shift_assay_pam_boundary, False),
    (break_lost_assay_join, False),
    (invent_reference_base, True),
]

rejected = accepted = wrong = 0
print("%-38s %-8s %s" % ("corruption", "verdict", "expected"))
for fn, expect_pass in CASES:
    # the checker reads ../data/sites129.json, so the sibling tree is staged too
    root = tempfile.mkdtemp()
    d = os.path.join(root, os.path.basename(HERE))
    os.mkdir(d)
    for n in os.listdir(HERE):
        if n.endswith((".tsv", ".py", ".json", ".md")):
            shutil.copy(os.path.join(HERE, n), d)
    src = os.path.join(HERE, os.pardir, "data")
    if os.path.isdir(src):
        shutil.copytree(src, os.path.join(root, "data"))
    try:
        fn(d)
        proc = subprocess.run([sys.executable, "-I", "-B",
                               os.path.join(d, "verify_genome_endtoend.py")],
                              capture_output=True, text=True, cwd=d)
        passed = proc.returncode == 0
    finally:
        shutil.rmtree(root, ignore_errors=True)
    ok = passed == expect_pass
    wrong += not ok
    rejected += not passed
    accepted += passed
    print("%-38s %-8s %s%s" % (fn.__name__, "accepted" if passed else "rejected",
                               "accepted" if expect_pass else "rejected",
                               "" if ok else "   <-- UNEXPECTED"))

print()
print("%d corruptions built, %d rejected, %d accepted, %d unexpected"
      % (len(CASES), rejected, accepted, wrong))
print("The one accepted corruption is the disclosed limit: an invented reference base cannot be")
print("detected without the reference, which this deposit does not carry.")
sys.exit(1 if wrong else 0)
