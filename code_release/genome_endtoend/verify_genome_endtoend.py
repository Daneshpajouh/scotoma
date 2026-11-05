#!/usr/bin/env python3
"""Replay every genome-scale claim in the paper from the exported tables.

Recomputes each load-bearing number from the tables in this directory and prints a
claim-by-claim verdict. Exit status 0 only when every claim is reproduced.

    python3 -I verify_genome_endtoend.py

No arguments, no network, no genome download: the tables in this directory are the
complete evidence for the numbers the paper states.
"""
import csv
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
VERDICTS = []


def claim(ok, label, got, want):
    VERDICTS.append((ok, label, got, want))
    print("  %-5s %-58s got %-22s want %s" % ("PASS" if ok else "FAIL", label, got, want))


def rows(name):
    with open(os.path.join(HERE, name), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- per chromosome
pc = rows("per_chromosome.tsv")
claim(len(pc) == 25, "records covered (chr1-22, X, Y, M)", len(pc), 25)
# 24 records pass every control. chr22 is flagged for one reason only: it carries
# alignments whose matched PAM is not literally GG, accepted through an IUPAC-ambiguous
# reference base. The flag is retained deliberately rather than suppressed.
passed = [r["chrom"] for r in pc if r["analysis_status"] == "PASS"]
flagged = [r["chrom"] for r in pc if r["analysis_status"] != "PASS"]
claim(len(passed) == 24 and flagged == ["chr22"],
      "24 records PASS every control; chr22 alone is flagged",
      (len(passed), flagged), (24, "['chr22']"))

released = sum(int(r["released_rna1"]) for r in pc)
corrected = sum(int(r["corrected_rna1"]) for r in pc)
conly = sum(int(r["corrected_only"]) for r in pc)
ronly = sum(int(r["released_only"]) for r in pc)
claim(released == 937934, "loci found by the released build", released, 937934)
claim(corrected == 972875, "loci found by the corrected build", corrected, 972875)
claim(conly == 34941, "loci found only by the corrected build", conly, 34941)
claim(ronly == 0, "loci found only by the released build", ronly, 0)

per_record_ronly = [r["chrom"] for r in pc if int(r["released_only"]) != 0]
claim(not per_record_ronly,
      "released-only is zero on each record separately, not only in total",
      per_record_ronly or "none", "none")

r0a = sum(int(r["released_rna0"]) for r in pc)
r0b = sum(int(r["corrected_rna0"]) for r in pc)
claim(r0a == r0b == 34573,
      "bulge-free control identical for the two builds",
      (r0a, r0b), (34573, 34573))

# candidate overhead, as an exact fraction then as a percentage
claim(abs(conly / released - 0.037253) < 5e-7,
      "candidate overhead of the repair",
      "%d/%d = %.6f" % (conly, released, conly / released), "0.037253")

# ---------------------------------------------------------------- the mechanism
diff = rows("difference_classes.tsv")
claim(len(diff) == 1, "difference classes present in the whole assembly", len(diff), 1)
d = diff[0]
claim((d["direction"], d["class"], d["rna_position_1based"], d["mismatches"]) ==
      ("corrected_only", "RNA", "20", "4"),
      "the single class is a terminal RNA deletion at four mismatches",
      (d["direction"], d["class"], d["rna_position_1based"], d["mismatches"]),
      "('corrected_only', 'RNA', '20', '4')")
claim(int(d["distinct_loci"]) == conly,
      "class table total agrees with the per-record total", d["distinct_loci"], conly)

al = rows("corrected_only_alignments.part000.tsv")
claim(len(al) == conly, "corrected-only alignments exported", len(al), conly)
for field, want in (("class", "RNA"), ("rna_position_1based", "20"), ("mismatches", "4")):
    seen = Counter(r[field] for r in al)
    claim(list(seen) == [want],
          "every corrected-only alignment has %s = %s" % (field, want),
          dict(seen), "{%r: %d}" % (want, conly))
claim(all(r["reference_verified"] == "1" for r in al),
      "every corrected-only alignment verified against the reference",
      len({r["reference_verified"] for r in al}), 1)

# The decisive check on the chr22 flag: the ambiguous-PAM alignments are in the set
# BOTH builds find, not in the difference. If even one of the 34,941 had a non-GG PAM,
# the headline difference would inherit the ambiguity.
nongg = [r for r in al if r["DNA"][-2:] != "GG"]
claim(not nongg,
      "corrected-only alignments with a PAM that is not literally GG",
      len(nongg), 0)

ro = rows("released_only_alignments.part000.tsv")
claim(len(ro) == 0, "released-only alignments exported", len(ro), 0)

# ---------------------------------------------------------------- the prediction
act = rows("active_sites.part000.tsv")
claim(len(act) == 260, "deposited measured-active rows examined", len(act), 260)
status = Counter(r["status"] for r in act)
claim(status.get("VERIFIED_NGG") == 237 and status.get("NON_NGG") == 23,
      "rows in and out of the paper's PAM scope",
      dict(status), "{'VERIFIED_NGG': 237, 'NON_NGG': 23}")
rel_ret = sum(int(r["released_returned"]) for r in act)
cor_ret = sum(int(r["corrected_returned"]) for r in act)
claim(rel_ret == 108, "measured-active rows the released build returns", rel_ret, 108)
claim(cor_ret == 237, "measured-active rows the corrected build returns", cor_ret, 237)
claim(cor_ret - rel_ret == 129,
      "difference, against the combinatorial census value 129",
      cor_ret - rel_ret, 129)
claim(status.get("VERIFIED_NGG", 0) + status.get("NON_NGG", 0) == len(act),
      "scope accounting closes", status.get("VERIFIED_NGG", 0) + status.get("NON_NGG", 0), 260)

# ---------------------------------------------------------------- per guide cross-check
pg = rows("per_guide.tsv")
claim(len(pg) == 110, "guides in the panel", len(pg), 110)
claim(sum(int(r["corrected_only"]) for r in pg) == conly,
      "per-guide total agrees with the per-record total",
      sum(int(r["corrected_only"]) for r in pg), conly)
claim(all(r["released_rna1_ungapped"] == r["corrected_rna1_ungapped"] for r in pg),
      "bulge-free counts equal for the two builds on every guide",
      sum(1 for r in pg if r["released_rna1_ungapped"] != r["corrected_rna1_ungapped"]), 0)
claim(all(int(r["released_only"]) == 0 for r in pg),
      "no guide loses a locus under the repair",
      sum(1 for r in pg if int(r["released_only"]) != 0), 0)

# ------------------------------------------------- recomputation, not label-reading
# Everything above this line compares reported totals. A reported total can be
# self-consistent and still wrong, so each check below RECONSTRUCTS the quantity from
# the sequences and coordinates and compares the result with the label beside it. A
# review pass showed that without these, several deliberately corrupted inputs, each
# preserving the totals, passed the whole file.

EXPECTED_RECORDS = (["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY", "chrM"])
names = [r["chrom"] for r in pc]
claim(sorted(names) == sorted(EXPECTED_RECORDS) and len(set(names)) == 25,
      "the 25 records are exactly the primary assembly's, each once",
      sorted(set(EXPECTED_RECORDS) ^ set(names)) or "exact match", "exact match")

# The flag on chr22 carries its reason. Any other status would be a different finding.
claim([r["analysis_status"] for r in pc if r["chrom"] == "chr22"] == ["REVIEW_REQUIRED"],
      "chr22's flag is the ambiguous-PAM review, not another failure",
      [r["analysis_status"] for r in pc if r["chrom"] == "chr22"], "['REVIEW_REQUIRED']")

# The bulge-free control, per guide as well as per record, and reconciled across the
# two groupings. The per-record and per-guide tables are independent aggregations of
# the same run, so their totals must agree or one of them is wrong.
unequal0 = [r["guide"] for r in pg if r["released_rna0"] != r["corrected_rna0"]]
claim(not unequal0,
      "bulge-free control equal for the two builds on every guide",
      unequal0[:4] or "none", "none")
pg_r0 = sum(int(r["released_rna0"]) for r in pg)
pg_c0 = sum(int(r["corrected_rna0"]) for r in pg)
claim(pg_r0 == pg_c0 == r0a,
      "per-guide bulge-free totals reconcile with the per-record total",
      (pg_r0, pg_c0, r0a), (34573, 34573, 34573))
pg_r1 = sum(int(r["released_rna1"]) for r in pg)
pg_c1 = sum(int(r["corrected_rna1"]) for r in pg)
claim((pg_r1, pg_c1) == (released, corrected),
      "per-guide bulge-aware totals reconcile with the per-record totals",
      (pg_r1, pg_c1), (released, corrected))

spacer_of = {r["guide"]: r["spacer"] for r in pg}
claim(len(spacer_of) == 110 and all(len(s) == 20 for s in spacer_of.values()),
      "every guide carries a 20 nt spacer", len(spacer_of), 110)


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


# Every exported alignment, reconstructed from its own strings. The loop recomputes the
# PAM-side boundary, the gap placement and the mismatch count; none of the three is read
# from the annotation it is compared against.
bad_coord = bad_gap = bad_crrna = bad_mm = bad_alpha = 0
first = {}
for r in al:
    dna = r["DNA"]
    ungapped = dna.replace("-", "")
    start, pam_coord = int(r["start"]), int(r["pam_coord"])
    want = start + len(ungapped) if r["strand"] == "+" else start
    if pam_coord != want:
        bad_coord += 1
        first.setdefault("coord", (r["guide"], r["chrom"], start, pam_coord, want))
    if len(dna) != 23 or dna.count("-") != 1 or dna.index("-") != 19:
        bad_gap += 1
        first.setdefault("gap", (r["guide"], dna))
    if set(ungapped) - set("ACGT"):
        bad_alpha += 1
        first.setdefault("alpha", (r["guide"], dna))
    cr = r["crRNA"]
    if not (cr.endswith("NNN") and cr[:20] == spacer_of.get(r["guide"])):
        bad_crrna += 1
        first.setdefault("crrna", (r["guide"], cr))
    # the aligned guide is the spacer with its terminal base withheld: 19 real bases
    if hamming(cr[:19], dna[:19]) != int(r["mismatches"]):
        bad_mm += 1
        first.setdefault("mm", (r["guide"], cr, dna, r["mismatches"],
                                hamming(cr[:19], dna[:19])))

claim(bad_coord == 0,
      "PAM boundary recomputed from start and alignment length",
      first.get("coord", bad_coord), 0)
claim(bad_gap == 0,
      "deletion found at position 20 of the string, not of the label",
      first.get("gap", bad_gap), 0)
claim(bad_alpha == 0,
      "no alignment carries an ambiguous DNA base",
      first.get("alpha", bad_alpha), 0)
claim(bad_crrna == 0,
      "every alignment's crRNA is its own guide's spacer",
      first.get("crrna", bad_crrna), 0)
claim(bad_mm == 0,
      "mismatch count recomputed by Hamming distance on all %d" % conly,
      first.get("mm", bad_mm), 0)

# ------------------------------------------- the prediction, row by row, not in total
# The difference 972,875 - 937,934 is a total, and totals survive permutation. The
# census route predicts WHICH rows each build returns, from the assay sequences alone.
# That prediction is computed here and compared row by row, so exchanging two rows'
# return bits while preserving both totals fails this check.
import json

SITES = os.path.join(HERE, os.pardir, "data", "sites129.json")
if not os.path.exists(SITES):
    claim(False, "data/sites129.json present for the row-by-row prediction",
          "missing", "present")
else:
    sites = json.load(open(SITES))
    claim(len(sites) == len(act),
          "selected assay rows match the deposited genome rows", len(sites), len(act))
    by_index = {int(r["row_index"]): r for r in act}
    claim(set(by_index) == set(range(len(sites))),
          "genome rows indexed 0..%d without gaps" % (len(sites) - 1),
          len(by_index), len(sites))
    agree = disagree = 0
    pred_rel = pred_cor = 0
    for i, s in enumerate(sites):
        at, spacer = s["aligned_target"], s["spacer"]
        dna19, pam = at[:19], at[20:23]
        pam_ok = pam[1:3] == "GG"
        costs = [hamming(spacer[:j] + spacer[j + 1:], dna19) for j in range(20)]
        rel_bit = int(pam_ok and min(costs[:19]) <= 4)
        cor_bit = int(pam_ok and min(costs) <= 4)
        pred_rel += rel_bit
        pred_cor += cor_bit
        g = by_index[i]
        for mine, theirs in ((rel_bit, "released_returned"), (cor_bit, "corrected_returned")):
            if mine == int(g[theirs]):
                agree += 1
            else:
                disagree += 1
    claim((agree, disagree) == (2 * len(sites), 0),
          "every return bit predicted from sequence alone",
          "%d agree, %d disagree" % (agree, disagree),
          "%d agree, 0 disagree" % (2 * len(sites)))
    claim((pred_rel, pred_cor) == (rel_ret, cor_ret),
          "predicted totals equal the run's own totals",
          (pred_rel, pred_cor), (rel_ret, cor_ret))
    # and the per-row control: a permutation that keeps the totals must be caught
    hi = next(k for k in sorted(by_index) if by_index[k]["released_returned"] == "1")
    lo = next(k for k in sorted(by_index)
              if by_index[k]["released_returned"] == "0"
              and by_index[k]["status"] == "VERIFIED_NGG")
    claim(hi != lo, "the permutation control has two distinct rows to exchange",
          (hi, lo), "distinct")

# ------------------------------------------------- domains, identities and reachability
# A reviewer built sixteen corruptions that the checks above all accepted: a single locus
# repeated 34,941 times, coordinates translated out of the assembly, invented chromosomes and
# strands, counts moved between subgroups until one went negative, and a target that the
# released build can in fact reach. Each check below closes one of those, by deriving the
# quantity rather than reading a label.

GRCh38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    "chrM": 16569,
}
claim(sum(GRCh38.values()) == 3088286401,
      "the record lengths sum to the assembled bases scanned",
      sum(GRCh38.values()), 3088286401)
claim(set(GRCh38) == set(EXPECTED_RECORDS), "record lengths cover exactly those records",
      sorted(set(GRCh38) ^ set(EXPECTED_RECORDS)) or "exact", "exact")

# every exported alignment is a distinct locus, in a real record, on a real strand, inside
# that record, with a 23-character crRNA
keys = Counter((r["guide"], r["chrom"], r["strand"], r["pam_coord"]) for r in al)
claim(len(keys) == conly, "exported alignments are distinct loci", len(keys), conly)
claim(max(keys.values()) == 1, "no exported locus appears twice", max(keys.values()), 1)
off = [r["chrom"] for r in al if r["chrom"] not in GRCh38]
claim(not off, "every alignment sits on a declared record", sorted(set(off))[:2] or "all", "all")
bad_strand = {r["strand"] for r in al} - {"+", "-"}
claim(not bad_strand, "every alignment has a real strand", sorted(bad_strand) or "ok", "ok")
# A reviewer placed a minus-strand hit whose start was inside chr21 but whose 22 nt footprint
# ended past the record, which the earlier one-ended test could not see.
oob = []
for r in al:
    L, s, pm = GRCh38.get(r["chrom"]), int(r["start"]), int(r["pam_coord"])
    span = len(r["DNA"].replace("-", ""))
    if L is None or not (0 <= s and s + span <= L and 0 <= pm <= L):
        oob.append((r["chrom"], s, s + span, L))
claim(not oob, "every alignment's whole footprint lies inside its record", oob[:1] or 0, 0)
badlen = [r["crRNA"] for r in al if len(r["crRNA"]) != 23]
claim(not badlen, "every crRNA is 23 characters", badlen[:1] or 0, 0)

# the released build must NOT be able to reach a corrected-only target
reachable = 0
for r in al:
    sp = spacer_of[r["guide"]]
    dna19 = r["DNA"].replace("-", "")[:19]
    if min(hamming(sp[:j] + sp[j + 1:], dna19) for j in range(19)) <= 4:
        reachable += 1
claim(reachable == 0,
      "corrected-only targets the released patterns could reach", reachable, 0)

# subgroup totals recomputed from the exports, not read from the summaries
by_chrom = Counter(r["chrom"] for r in al)
stated_chrom = {r["chrom"]: int(r["corrected_only"]) for r in pc}
claim(all(by_chrom.get(k, 0) == v for k, v in stated_chrom.items()),
      "per-record corrected-only recomputed from the exports",
      sum(1 for k, v in stated_chrom.items() if by_chrom.get(k, 0) != v), 0)
by_guide = Counter(r["guide"] for r in al)
stated_guide = {r["guide"]: int(r["corrected_only"]) for r in pg}
claim(all(by_guide.get(k, 0) == v for k, v in stated_guide.items()),
      "per-guide corrected-only recomputed from the exports",
      sum(1 for k, v in stated_guide.items() if by_guide.get(k, 0) != v), 0)

COUNT_FIELDS = ("released_rna1", "corrected_rna1", "released_rna0", "corrected_rna0",
                "corrected_only", "released_only")
neg = [(r.get("chrom") or r.get("guide"), f) for r in pc + pg for f in COUNT_FIELDS
       if f in r and int(r[f]) < 0]
claim(not neg, "every reported count is nonnegative", neg[:1] or 0, 0)

# Counts moved between subgroups by equal and opposite amounts leave every total intact, so
# totals cannot catch them. The set identity can: corrected minus released is exactly
# corrected-only minus released-only, on every record and every guide.
D = lambda r, a, b: int(r[a]) - int(r[b])
vio_c = [r["chrom"] for r in pc
         if D(r, "corrected_rna1", "released_rna1") != D(r, "corrected_only", "released_only")]
claim(not vio_c, "per record, the rna1 difference equals the subgroup difference",
      vio_c[:2] or 0, 0)
vio_g = [r["guide"] for r in pg
         if D(r, "corrected_rna1", "released_rna1") != D(r, "corrected_only", "released_only")]
claim(not vio_g, "per guide, the rna1 difference equals the subgroup difference",
      vio_g[:2] or 0, 0)
r0_split = [r["chrom"] for r in pc if r["released_rna0"] != r["corrected_rna0"]]
claim(not r0_split, "bulge-free control equal on every record separately",
      r0_split[:2] or "none", "none")
ung = [r["guide"] for r in pg
       if int(r["released_rna1_ungapped"]) > int(r["released_rna1"])
       or int(r["corrected_rna1_ungapped"]) > int(r["corrected_rna1"])]
claim(not ung, "bulge-free subset never exceeds the bulge-aware total",
      ung[:2] or "none", "none")
# the bulge-free subset of the bulge-aware run IS the bulge-free run, per guide and per build
split = [r["guide"] for r in pg
         if not (r["released_rna1_ungapped"] == r["released_rna0"]
                 == r["corrected_rna0"] == r["corrected_rna1_ungapped"])]
claim(not split, "bulge-free subset equals the bulge-free run on every guide",
      split[:2] or "none", "none")

# the guide table as a whole is bound, so a guide with no alignments cannot be rewritten
import hashlib
digest = hashlib.sha256("\n".join(sorted(r["spacer"].upper() for r in pg)).encode()).hexdigest()
claim(digest == "0b9bce3d297d36848c825a4b6ce222bca73fc1cbef4e8cbc4be8337322c109ad",
      "the panel's spacer set matches its pinned digest", digest[:16] + "...", "0b9bce3d297d3684...")
alpha = [r["guide"] for r in pg
         if len(r["spacer"]) != 20 or set(r["spacer"].upper()) - set("ACGT")]
claim(not alpha, "every spacer is 20 nt over the four letters", alpha[:2] or "none", "none")

# the measured-active rows keep their assay identity, and their PAM status is recomputed
if os.path.exists(SITES):
    mismatched_id = 0
    for i, srow in enumerate(sites):
        g = by_index[i]
        if g["chrom"] != srow["chrom"] or g["strand"] != srow["strand"]:
            mismatched_id += 1
    claim(mismatched_id == 0,
          "each genome row keeps its assay row's record and strand", mismatched_id, 0)

    # Record and strand alone left the rest of the identity free: a reviewer replaced every
    # guide with NO_GUIDE and every coordinate with a large negative number and nothing saw it.
    spacer_to_guide = {r["spacer"]: r["guide"] for r in pg}
    wrong_guide = sum(1 for i, srow in enumerate(sites)
                      if by_index[i]["guide"] != spacer_to_guide.get(srow["spacer"]))
    claim(wrong_guide == 0, "  and its guide, by the spacer that guide carries", wrong_guide, 0)
    wrong_coord = sum(1 for i, srow in enumerate(sites)
                      if by_index[i]["start"] != str(srow["start"]))
    claim(wrong_coord == 0, "  and its start coordinate", wrong_coord, 0)
    wrong_del = sum(1 for i, srow in enumerate(sites)
                    if int(by_index[i]["deleted_position_0based"]) != int(srow["deleted_position"]))
    claim(wrong_del == 0, "  and its deleted position", wrong_del, 0)
    outside, pamrel = [], 0
    for i, srow in enumerate(sites):
        g = by_index[i]
        L = GRCh38.get(g["chrom"])
        st, pm = int(g["start"]), int(g["pam_coord"])
        if L is None or not (0 <= st <= L and 0 <= pm <= L) \
                or not (0 <= int(g["deleted_position_0based"]) < 20):
            outside.append((g["chrom"], st, pm))
        # the PAM coordinate is not a free field: it is the footprint away from the start on
        # the plus strand and the start itself on the minus strand, the same relation the
        # exported alignments satisfy. A reviewer shifted all 260 by one and nothing saw it.
        want = st + len(srow["aligned_target"].replace("-", "")) if g["strand"] == "+" else st
        if pm != want:
            pamrel += 1
    claim(not outside, "  and lies inside its record at a real deleted position",
          outside[:1] or 0, 0)
    claim(pamrel == 0, "  and its PAM coordinate is the one its start and strand imply",
          pamrel, 0)

    # The 129 are claimed lost BY the corrected-only export, so the two tables must join. A
    # reviewer moved one export row by a single base, consistently, and every distance, domain,
    # uniqueness and subgroup check still passed while the export held only 128 of the 129
    # assay keys. Nothing had ever compared them.
    lost = [r for r in act if r["released_returned"] == "0"
            and r["corrected_returned"] == "1"]
    key = lambda r: (r["guide"], r["chrom"], r["strand"], r["pam_coord"])
    exported = {key(r) for r in al}
    joined = sum(1 for r in lost if key(r) in exported)
    claim(len(lost) == 129, "assay rows the released build misses", len(lost), 129)
    claim(joined == len(lost),
          "  every one of them is a row of the corrected-only export",
          "%d of %d" % (joined, len(lost)), len(lost))
    wrong_status = 0
    for i, srow in enumerate(sites):
        want = "VERIFIED_NGG" if srow["aligned_target"][20:23][1:3] == "GG" else "NON_NGG"
        if by_index[i]["status"] != want:
            wrong_status += 1
    claim(wrong_status == 0,
          "each row's PAM status recomputed from its own PAM bases", wrong_status, 0)
    assay_spacers = {r["spacer"] for r in sites}
    unbound = [g for g in spacer_of.values() if False]
    claim(assay_spacers <= set(spacer_of.values()),
          "every assay spacer appears in the guide table",
          len(assay_spacers - set(spacer_of.values())), 0)
    # the named permutation control, actually executed rather than merely selected
    hi = next(k for k in sorted(by_index) if by_index[k]["released_returned"] == "1")
    lo = next(k for k in sorted(by_index)
              if by_index[k]["released_returned"] == "0"
              and by_index[k]["status"] == "VERIFIED_NGG")
    swapped = dict(by_index)
    swapped[hi] = dict(by_index[hi], released_returned="0")
    swapped[lo] = dict(by_index[lo], released_returned="1")
    dis = 0
    for i, srow in enumerate(sites):
        at, sp2 = srow["aligned_target"], srow["spacer"]
        d19, pm = at[:19], at[20:23]
        ok = pm[1:3] == "GG"
        cost = [hamming(sp2[:j] + sp2[j + 1:], d19) for j in range(20)]
        mine = int(ok and min(cost[:19]) <= 4)
        if mine != int(swapped[i]["released_returned"]):
            dis += 1
    claim(dis == 2,
          "  control: exchanging two return bits is caught, totals preserved", dis, 2)

# ---------------------------------------------------------------- summary
bad = [v for v in VERDICTS if not v[0]]
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - len(bad), len(bad)))
if bad:
    sys.exit(1)
print("Every genome-scale claim in the paper is reproduced from these tables.")
