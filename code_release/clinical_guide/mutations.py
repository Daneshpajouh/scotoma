#!/usr/bin/env python3
"""Plant corrupted clinical_guide inputs; exit non-zero if an unmarked one is accepted."""
import csv, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
VERIFIER = os.path.join(HERE, "verify_clinical_guide.py")


def copy_tree(dst):
    # copy only the TSV inputs the verifier reads, plus the script itself via cwd
    for name in os.listdir(HERE):
        if name.endswith(".tsv") or name.endswith(".json"):
            shutil.copy(os.path.join(HERE, name), os.path.join(dst, name))
    shutil.copy(VERIFIER, os.path.join(dst, "verify_clinical_guide.py"))


def run(dst):
    return subprocess.run([sys.executable, "-I", "-B", "verify_clinical_guide.py"],
                          cwd=dst, capture_output=True, text=True).returncode


def rewrite(path, fn):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
        fields = list(rows[0].keys())
    rows = fn(rows)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def corrupt_sweep_mm(dst):
    # r5 attack: budget-2 corrected RNA sweep row made perfect, label kept at 2
    p = os.path.join(dst, "sweep_mm2_corrected_rna1.tsv")
    def go(rows):
        for r in rows:
            if r["Bulge type"] == "RNA" and int(r["Mismatches"]) == 2:
                cr, dna = list(r["crRNA"]), list(r["DNA"])
                # force spacer columns to match where DNA has letters
                for i, (a, b) in enumerate(zip(cr[:-3], dna[:-3])):
                    if a != "-" and b != "-":
                        dna[i] = a
                r["DNA"] = "".join(dna)
                # leave Mismatches=2
                break
        return rows
    rewrite(p, go)


def corrupt_shared_pam(dst):
    # r5 attack: shared budget-4 PAM -> AA on both builds
    for name in ("released_rna1.tsv", "corrected_rna1.tsv"):
        p = os.path.join(dst, name)
        def go(rows):
            for r in rows:
                dna = r["DNA"]
                r["DNA"] = dna[:-2] + "AA"
            return rows
        rewrite(p, go)


CASES = [
    ("clinical_sweep_mm", corrupt_sweep_mm, False),
    ("clinical_non_NGG", corrupt_shared_pam, False),
]

print("%-28s %-8s %s" % ("corruption", "verdict", "expected"))
unexpected = 0
for name, fn, expect_accept in CASES:
    d = tempfile.mkdtemp()
    try:
        copy_tree(d)
        fn(d)
        code = run(d)
        accepted = code == 0
        verdict = "accepted" if accepted else "rejected"
        expected = "accepted" if expect_accept else "rejected"
        bad = accepted != expect_accept
        unexpected += bad
        print("%-28s %-8s %s%s" % (name, verdict, expected, "   <-- UNEXPECTED" if bad else ""))
    finally:
        shutil.rmtree(d)
print("%d corruptions, %d unexpected" % (len(CASES), unexpected))
sys.exit(1 if unexpected else 0)
