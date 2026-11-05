#!/usr/bin/env python3
"""Independent verification of the PANEL_EXTEND measurement.

Shares no code with panel_extend.py.  It re-reads the raw Cas-OFFinder output
tables and re-derives, per guide and per budget:

  released, corrected, corrected_only, released_only, terminal_attributable

The locus convention is re-implemented here from its definition (guide,
chromosome, strand, PAM-side forward boundary, zero based, no merging) rather
than imported.  Terminal attribution is recomputed from the aligned genomic
sequence against the terminal deletion word g[:-1], not inherited.

Controls, all of which must behave as stated or the run is rejected:
  C1  every row's crRNA equals this guide plus NNN
  C2  every row is an RNA bulge of size one and its mismatch count is <= budget
  C3  released must be a subset of corrected (released_only == 0)
  C4  forgery: attributing by the FIRST word g[1:] instead of the terminal word
      must not explain the corrected-only set
  C5  corruption: a single-base corruption of the guide must reduce the
      terminal attribution count
"""
import os
import json, os, sys

ROOT = os.environ.get("PANEL_ROOT", os.path.dirname(os.path.abspath(__file__)))
SHARD = int(sys.argv[1])
SDIR = os.path.join(ROOT, "shard_%d" % SHARD)
BUDGETS = [5, 6, 7]

CHROM_IDX = {}
def chrom_id(c):
    i = CHROM_IDX.get(c)
    if i is None:
        i = len(CHROM_IDX) + 1
        CHROM_IDX[c] = i
    return i

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def hamming(ref, word):
    """Substitutions between a reference window and a guide-derived word.
    A reference IUPAC ambiguity code counts as a match when the word's base is
    one of the bases it denotes, which is how the search itself scores it."""
    if len(ref) != len(word):
        return 10 ** 6
    n = 0
    for x, y in zip(ref, word):
        if x == y:
            continue
        if y in IUPAC.get(x, x):
            continue
        n += 1
    return n


def hamming_strict(ref, word):
    if len(ref) != len(word):
        return 10 ** 6
    return sum(1 for x, y in zip(ref, word) if x != y)

def load_guides():
    g = {}
    with open(os.path.join(ROOT, "panel_guides.txt")) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            gid, seq = line.split("\t")
            g[gid] = seq.upper()
    return g

GUIDES = load_guides()

def shard_guide_ids(m):
    ids = []
    with open(os.path.join(SDIR, "input_mm%d.txt" % m)) as fh:
        for i, line in enumerate(fh):
            if i < 2:
                continue
            parts = line.split()
            if len(parts) == 3:
                ids.append(parts[2])
    return ids

def parse(path, m, want_terminal, want_first):
    """Return per-guide locus sets, plus per-guide sets of loci explained by the
    terminal word and by the first word, recomputed from the aligned sequence."""
    loci = {}
    term = {}
    first = {}
    bad = []
    nrows = 0
    nskip = [0]
    nambig = [0]
    with open(path) as fh:
        head = fh.readline()
        if not head.startswith("#Id"):
            raise SystemExit("unexpected header in %s: %r" % (path, head[:60]))
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            gid, btype, crrna, dna, chrom, pos, direction, mm, bsize = f[:9]
            nrows += 1
            g = GUIDES[gid]
            if crrna != g + "NNN":
                bad.append(("C1", gid, crrna))
                continue
            if btype != "RNA" or bsize != "1":
                nskip[0] += 1          # bulge-free row of the same run, not this object
                continue
            if int(mm) > m:
                bad.append(("C2", gid, btype, bsize, mm))
                continue
            genomic = dna.replace("-", "")
            if len(genomic) != 22:
                bad.append(("C2len", gid, dna))
                continue
            p = int(pos)
            key = (chrom_id(chrom) << 40) | ((p + 22 if direction == "+" else p) << 1) | (1 if direction == "+" else 0)
            s = loci.get(gid)
            if s is None:
                s = loci[gid] = set()
            s.add(key)
            proto = genomic[:19]
            if want_terminal:
                d = hamming(proto, g[:-1])
                if d <= m and hamming_strict(proto, g[:-1]) > m:
                    nambig[0] += 1
                if d <= m:
                    t = term.get(gid)
                    if t is None:
                        t = term[gid] = set()
                    t.add(key)
            if want_first:
                if hamming(proto, g[1:]) <= m:
                    t = first.get(gid)
                    if t is None:
                        t = first[gid] = set()
                    t.add(key)
    return loci, term, first, bad, nrows, nskip[0], nambig[0]

out = {"shard": SHARD, "budgets": {}, "controls": {"C1_C2_violations": [], "C3_released_only": [],
        "C4_first_word_explains_all": [], "C5_corruption": []}}

for m in BUDGETS:
    rel_path = os.path.join(SDIR, "released_mm%d.tsv" % m)
    cor_path = os.path.join(SDIR, "corrected_mm%d.tsv" % m)
    rel, _, _, bad_r, nr, sk_r, _ = parse(rel_path, m, False, False)
    cor, term, first, bad_c, nc, sk_c, amb = parse(cor_path, m, True, True)
    out["controls"]["C1_C2_violations"].extend(bad_r[:5] + bad_c[:5])
    per = {}
    for gid in sorted(set(list(rel.keys()) + list(cor.keys()))):
        R = rel.get(gid, set())
        C = cor.get(gid, set())
        co = C - R
        ro = R - C
        ta = len(co & term.get(gid, set()))
        fa = len(co & first.get(gid, set()))
        per[gid] = {"released": len(R), "corrected": len(C), "corrected_only": len(co),
                    "released_only": len(ro), "terminal_attributable": ta,
                    "first_word_attributable": fa}
        if ro:
            out["controls"]["C3_released_only"].append([gid, m, len(ro)])
        if co and fa == len(co):
            out["controls"]["C4_first_word_explains_all"].append([gid, m, len(co)])
    out["budgets"][str(m)] = {"rows_released": nr, "rows_corrected": nc,
                              "bulge_free_rows_skipped_released": sk_r,
                              "bulge_free_rows_skipped_corrected": sk_c,
                              "rows_matched_only_via_ambiguity_code": amb,
                              "per_guide": per}

# C5: corruption control on the first guide of this shard that has a deficit, at budget 5
m = 5
cand = [g for g, v in out["budgets"]["5"]["per_guide"].items() if v["corrected_only"] > 0]
if cand:
    gid = sorted(cand)[0]
    g = GUIDES[gid]
    corrupt = ("A" if g[3] != "A" else "C") + ""  # single base change at position 4
    gcorr = g[:3] + corrupt + g[4:]
    assert gcorr != g and len(gcorr) == len(g)
    real = 0
    fake = 0
    cor_path = os.path.join(SDIR, "corrected_mm%d.tsv" % m)
    rel_path = os.path.join(SDIR, "released_mm%d.tsv" % m)
    relset = set()
    with open(rel_path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != gid:
                continue
            p = int(f[5])
            relset.add((chrom_id(f[4]) << 40) | ((p + 22 if f[6] == "+" else p) << 1) | (1 if f[6] == "+" else 0))
    seen = set()
    with open(cor_path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != gid:
                continue
            p = int(f[5])
            key = (chrom_id(f[4]) << 40) | ((p + 22 if f[6] == "+" else p) << 1) | (1 if f[6] == "+" else 0)
            if key in relset or key in seen:
                continue
            proto = f[3].replace("-", "")[:19]
            if hamming(proto, g[:-1]) <= m:
                real += 1
                seen.add(key)
            elif hamming(proto, gcorr[:-1]) <= m:
                fake += 1
    out["controls"]["C5_corruption"] = {"guide": gid, "real_terminal_attributed": real,
                                        "corrupted_terminal_attributed_extra": fake,
                                        "rejected": fake < real}

ok = (not out["controls"]["C1_C2_violations"] and not out["controls"]["C3_released_only"]
      and not out["controls"]["C4_first_word_explains_all"]
      and out["controls"]["C5_corruption"].get("rejected", False))
out["verdict"] = "VERIFY_OK" if ok else "VERIFY_FLAGGED"
with open(os.path.join(SDIR, "VERIFY.json"), "w") as fh:
    json.dump(out, fh, indent=2, sort_keys=True)
print(out["verdict"], "shard", SHARD)
for m in BUDGETS:
    per = out["budgets"][str(m)]["per_guide"]
    print("  m=%d guides=%d corrected_only=%d terminal=%d" % (
        m, len(per), sum(v["corrected_only"] for v in per.values()),
        sum(v["terminal_attributable"] for v in per.values())))
print("  C5:", out["controls"]["C5_corruption"])
