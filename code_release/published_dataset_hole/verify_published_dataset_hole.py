#!/usr/bin/env python3
"""A published off-target dataset is missing two whole alignment classes.

    python3 -I verify_published_dataset_hole.py [DIR]

With no argument it checks the deposited census, `gap_position_census.json`. Given a directory
holding the two public CSVs it recomputes that census from them first, checks their SHA-256
against the pinned values, and then runs the same assertions. Nothing else is read: no genome,
no executable, no network.

The claim. CRISPR-Net is a published model for off-target activity with indels and mismatches;
its paper states that Cas-OFFinder was used to enumerate candidate off-targets in hg38,
including with indels, and its two indel datasets are deposited. In every one of the 238,705
enumerated RNA-bulge rows across both datasets and all sixteen guides, the deletion sits at
spacer position 2 through 19. Position 1 and position 20 never occur.

That range is not a convention. The 2016 Cas-OFFinder bulge wrapper enumerates deletion starts
with `range(1, len(target) - bulge_size)`, which for a 20 nt spacer and one bulge is exactly
1-based positions 2 through 19. The deposited data carry the loop's signature.

Two mechanisms fit that support and the rows do not separate them: a loop that never visits
the endpoints, and an alignment policy that re-anchors an endpoint gap inward. An earlier
version of this file argued that no canonicalisation could produce the observed support. That
argument was wrong and was withdrawn; executed realignment policies reproduce it. Nothing here
depends on which mechanism applies.
"""
import csv
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
VERDICTS = []
SPACER, BULGE = 20, 1

PINNED = {
    "CIRCLE_seq_10gRNA_wholeDataset.csv":
        "86d88947ef9c0a1b545814c4e4241e4f6103d15a4d056d697a4ea5b6017648a8",
    "elevation_6gRNA_wholeDataset.csv":
        "56f1c95a0a3d61b873af984f7a2110e27538540b49904b913982d63764dcbb55",
}
COLUMNS = {
    "CIRCLE_seq_10gRNA_wholeDataset.csv": ("sgRNA_seq", "off_seq", "label"),
    "elevation_6gRNA_wholeDataset.csv": ("crRNA", "DNA", "label"),
}


def claim(ok, label, got, want):
    VERDICTS.append(ok)
    print("  %-5s %-52s got %-20s want %s" % ("PASS" if ok else "FAIL", label, got, want))


def recompute(d):
    import hashlib
    out = {}
    for fn, (gc, oc, lc) in COLUMNS.items():
        p = os.path.join(d, fn)
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        rna, dna, spacers, n = defaultdict(Counter), defaultdict(Counter), Counter(), 0
        widths = Counter()
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                n += 1
                g, o, lab = r[gc], r[oc], r[lc]
                widths[(len(g), len(o))] += 1
                spacers[g.replace("_", "").replace("-", "")[:SPACER]] += 1
                if "_" in o:
                    rna[lab][o.index("_")] += 1
                elif "_" in g:
                    dna[lab][g.index("_")] += 1
        out[fn] = dict(rows=n, sha256=h, distinct_spacers=len(spacers),
                       widths={str(k): v for k, v in widths.items()},
                       rna={k: {str(i): c for i, c in sorted(v.items())} for k, v in rna.items()},
                       dna={k: {str(i): c for i, c in sorted(v.items())} for k, v in dna.items()})
    return out


BOUND = len(sys.argv) > 1
if BOUND:
    census = recompute(sys.argv[1])
    for fn, want in PINNED.items():
        claim(census[fn]["sha256"] == want, "%.30s: public file digest" % fn,
              census[fn]["sha256"][:12] + "...", want[:12] + "...")
    stored = json.load(open(os.path.join(HERE, "gap_position_census.json")))
    claim(census == stored, "recomputation equals the deposited census",
          "identical" if census == stored else "differs", "identical")
else:
    census = json.load(open(os.path.join(HERE, "gap_position_census.json")))
    # The census is a derived summary and cannot certify itself. A reviewer altered it six ways
    # and every check still passed, because nothing bound it to the files it came from. The
    # digests it records are checked here, which catches a census claiming the wrong source, but
    # only the recomputing form below binds the counts themselves.
    for fn, want in PINNED.items():
        claim(census.get(fn, {}).get("sha256") == want,
              "%.30s: census names the pinned file" % fn,
              (census.get(fn, {}).get("sha256") or "absent")[:12] + "...", want[:12] + "...")

# the 2016 wrapper's own loop bound, written out rather than asserted
WRAPPER = [i + 1 for i in range(1, SPACER - BULGE)]     # range(1, len - bulge), 1-based
claim(WRAPPER == list(range(2, 20)),
      "the wrapper's loop covers 1-based positions", "%d..%d" % (WRAPPER[0], WRAPPER[-1]), "2..19")

TOTAL = 0
for fn, exp_rows, exp_guides in (
        ("CIRCLE_seq_10gRNA_wholeDataset.csv", 584949, 10),
        ("elevation_6gRNA_wholeDataset.csv", 213943, 6)):
    d = census[fn]
    print()
    print("  %s" % fn)
    claim(d["rows"] == exp_rows, "  rows", d["rows"], exp_rows)
    claim(d["distinct_spacers"] == exp_guides, "  distinct 20 nt spacers",
          d["distinct_spacers"], exp_guides)
    claim(list(d["widths"]) == ["(24, 24)"], "  every row is a 24-column alignment",
          list(d["widths"]), "['(24, 24)']")

    neg = {int(k): v for k, v in d["rna"].get("0.0", {}).items()}
    TOTAL += sum(neg.values())
    claim(sorted(neg) == WRAPPER, "  enumerated RNA-bulge positions",
          "%d..%d" % (min(neg), max(neg)), "2..19")
    claim(1 not in neg and SPACER not in neg,
          "  positions 1 and 20 among %d rows" % sum(neg.values()),
          sorted(set(neg) & {1, SPACER}) or "absent", "absent")

    dslots = {int(k): v for k, v in d["dna"].get("0.0", {}).items()}
    claim(sorted(dslots) == list(range(1, SPACER)), "  enumerated DNA-bulge slots",
          "%d..%d" % (min(dslots), max(dslots)), "1..19")
    claim(0 not in dslots and SPACER not in dslots,
          "  slots 0 and 20 among %d rows" % sum(dslots.values()),
          sorted(set(dslots) & {0, SPACER}) or "absent", "absent")

    # control: the MEASURED rows are a different, sparser set, so the assertion above is
    # about the enumerated rows and is not a property every row trivially has
    pos = {int(k): v for k, v in d["rna"].get("1.0", {}).items()}
    claim(set(pos) < set(neg) and len(pos) < len(neg),
          "  control: measured rows occupy a strictly smaller set",
          "%d positions, %d rows" % (len(pos), sum(pos.values())),
          "fewer than %d positions" % len(neg))

print()
claim(TOTAL == 238705, "enumerated RNA-bulge rows over both datasets", TOTAL, 238705)

# ------------------------------------------------- what the two omitted positions cost
# A whole-assembly scan of GRCh38, 3,088,286,401 bases, for every distinct deletion word of
# these sixteen guides, at the settings the datasets' own rows exhibit: three mismatches, PAM
# unconstrained, bulge sizes zero so the enumerator under audit is never exercised.
print()
print("  what the two omitted positions cost, from a whole-assembly scan")
# ------------------------------------------------- the tables, derived rather than read
# A reviewer collapsed unreachable_sequences.tsv to 10,888 copies of one row and swapped two
# absence flags, and every check still passed, because the counts above were read off the file
# instead of rebuilt. Both deposited tables are now regenerated here by the deposited producer
# from the scan statistics and, when bound, the two public CSVs, and compared byte for byte.
import subprocess, tempfile

TABLES = ("per_guide_hole.tsv", "unreachable_sequences.tsv")
if BOUND:
    tmp = tempfile.mkdtemp()
    proc = subprocess.run([sys.executable, "-I", "-B", os.path.join(HERE, "derive_hole.py"),
                           HERE, sys.argv[1], tmp], capture_output=True, text=True)
    claim(proc.returncode == 0, "the producer rebuilds the tables", proc.returncode, 0)
    for t in TABLES:
        a = open(os.path.join(HERE, t), "rb").read()
        b = open(os.path.join(tmp, t), "rb").read() if proc.returncode == 0 else b""
        claim(a == b, "  %s is what the producer makes" % t,
              "identical" if a == b else "differs", "identical")
        os.path.exists(os.path.join(tmp, t)) and os.remove(os.path.join(tmp, t))
    os.rmdir(tmp)

# Bind the deposited raw scan to the reachability statistics. A reviewer replaced
# scan_raw_outputs.tgz by an empty archive and left reached_*.tsv untouched; every check
# still passed because those two files were trusted as inputs. They are now rebuilt from
# the archive (plus the deposited queries.tsv derived from guides_meta) and compared.
import tempfile, tarfile, shutil, subprocess, glob as _glob
ARCHIVE = os.path.join(HERE, "scan_raw_outputs.tgz")
claim(os.path.isfile(ARCHIVE), "raw scan archive is deposited",
      "present" if os.path.isfile(ARCHIVE) else "absent", "present")
claim(os.path.isfile(os.path.join(HERE, "queries.tsv")),
      "queries.tsv that names each scanned word is deposited",
      "present" if os.path.isfile(os.path.join(HERE, "queries.tsv")) else "absent", "present")
scan_tmp = tempfile.mkdtemp()
try:
    with tarfile.open(ARCHIVE, "r:gz") as tf:
        members = tf.getmembers()
        claim(len(members) > 0, "raw scan archive is non-empty", len(members), ">0")
        tf.extractall(scan_tmp)
    n_tsv = len(_glob.glob(os.path.join(scan_tmp, "out", "chr*.tsv")))
    claim(n_tsv == 25, "raw scan holds one TSV per primary-assembly record", n_tsv, 25)
    # row count before reduction
    raw_rows = 0
    for p in sorted(_glob.glob(os.path.join(scan_tmp, "out", "chr*.tsv"))):
        with open(p, newline="") as fh:
            raw_rows += sum(1 for _ in csv.DictReader(fh, delimiter="\t"))
    claim(raw_rows == 850463, "raw scan alignment rows", raw_rows, 850463)
    shutil.copy(os.path.join(HERE, "queries.tsv"), os.path.join(scan_tmp, "queries.tsv"))
    shutil.copy(os.path.join(HERE, "guides_meta.json"),
                os.path.join(scan_tmp, "guides_meta.json"))
    out_tmp = tempfile.mkdtemp()
    proc = subprocess.run(
        [sys.executable, "-I", "-B", os.path.join(HERE, "make_reached.py"), scan_tmp, out_tmp],
        capture_output=True, text=True)
    claim(proc.returncode == 0, "make_reached rebuilds statistics from the archive",
          proc.returncode, 0)
    for t in ("reached_positions.tsv", "reached_loci.tsv"):
        a = open(os.path.join(HERE, t), "rb").read()
        b = open(os.path.join(out_tmp, t), "rb").read() if proc.returncode == 0 else b""
        claim(a == b, "  deposited %s equals the archive rebuild" % t,
              "identical" if a == b else "differs", "identical")
    shutil.rmtree(out_tmp)
finally:
    shutil.rmtree(scan_tmp)

# whatever the tables say, the quantities the paper states are recomputed from the scan
EM = set(range(2, SPACER))
reach = defaultdict(dict)
with open(os.path.join(HERE, "reached_positions.tsv"), newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        reach[r["guide_id"]][r["off_target_19mer"]] = \
            {int(x) for x in r["positions"].split(",")}
meta = json.load(open(os.path.join(HERE, "guides_meta.json")))
lrows = list(csv.DictReader(open(os.path.join(HERE, "reached_loci.tsv")), delimiter="\t"))

# The query words are part of the specification, not decoration. A reviewer replaced all 241
# with nineteen As and every assertion stayed green, because the reducer reads only query ids
# and position labels. Each word is now rebuilt from its guide and its stated positions.
qpath = os.path.join(HERE, "queries.tsv")
if os.path.exists(qpath):
    qrows = list(csv.DictReader(open(qpath), delimiter="\t"))
    claim(len(qrows) == 241, "deposited query words", len(qrows), 241)
    claim(len({r["word"] for r in qrows}) == 241,
          "  all distinct", len({r["word"] for r in qrows}), 241)
    wrong = 0
    for r in qrows:
        g = r["query_id"].split("_")[0]
        sp = meta[g]["spacer"]
        pos = [int(x) for x in r["positions"].split(",")]
        if any(sp[:i - 1] + sp[i:] != r["word"].upper() for i in pos):
            wrong += 1
    claim(wrong == 0,
          "  each is its guide's spacer with the stated position deleted", wrong, 0)
    claim(all(r["word"].upper() in meta[r["query_id"].split("_")[0]]["words"] for r in qrows),
          "  and is a word of that guide's library", "all", "all")

claim(len(reach) == 16, "guides in the scan statistic", len(reach), 16)
claim(sum(len(meta[g]["words"]) for g in reach) == 241,
      "distinct deletion words required", sum(len(meta[g]["words"]) for g in reach), 241)
claim(sum(1 for g in reach for _, ps in meta[g]["words"].items() if not (set(ps) & EM)) == 21,
      "words sitting only at position 1 or 20",
      sum(1 for g in reach for _, ps in meta[g]["words"].items() if not (set(ps) & EM)), 21)
claim(sum(int(r["loci"]) for r in lrows) == 769330, "loci the complete library reaches",
      sum(int(r["loci"]) for r in lrows), 769330)
claim(sum(int(r["unreachable_loci"]) for r in lrows) == 35272,
      "loci reachable only through the omitted positions",
      sum(int(r["unreachable_loci"]) for r in lrows), 35272)
claim(sum(len(v) for v in reach.values()) == 122170,
      "distinct off-target sequences, complete library",
      sum(len(v) for v in reach.values()), 122170)

UN = {(g, s) for g in reach for s, ps in reach[g].items() if not (ps & EM)}
claim(len(UN) == 10888, "sequences reachable only through those positions", len(UN), 10888)

# the deposited sequence table must be exactly that set, as a set of distinct pairs
seqs = list(csv.DictReader(open(os.path.join(HERE, "unreachable_sequences.tsv")),
                           delimiter="\t"))
pairs = {(r["guide_id"], r["off_target_19mer"]) for r in seqs}
claim(len(seqs) == 10888, "deposited rows", len(seqs), 10888)
claim(len(pairs) == 10888, "  deposited DISTINCT guide and sequence pairs", len(pairs), 10888)
claim(pairs == UN, "  the deposited set is the recomputed set",
      "identical" if pairs == UN else "differs", "identical")

# every deposited sequence is bound to its guide: reachable at an endpoint, out of reach inside
spacer_of = {g: meta[g]["spacer"] for g in meta}


def ham(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


worst_in, bad_in, bad_end = 99, 0, 0
for g, t in sorted(UN):
    sp = spacer_of[g]
    interior = min(ham(sp[:j] + sp[j + 1:], t) for j in range(1, SPACER - 1))
    endpoint = min(ham(sp[1:], t), ham(sp[:-1], t))
    worst_in = min(worst_in, interior)
    bad_in += interior <= 3
    bad_end += endpoint > 3
claim(bad_end == 0, "  each is within budget at position 1 or 20", bad_end, 0)
claim(bad_in == 0, "  and beyond budget at every position 2..19", bad_in, 0)
claim(worst_in == 4, "  the closest interior word to any of them", worst_in, 4)

# the rescue law is the control: a guide whose boundary bases agree loses nothing
pg = list(csv.DictReader(open(os.path.join(HERE, "per_guide_hole.tsv")), delimiter="\t"))
I = lambda r, k: int(r[k])
rescued = [r for r in pg if r["spacer"][-1] == r["spacer"][-2] and r["spacer"][0] == r["spacer"][1]]
atrisk = [r for r in pg if r not in rescued]
claim(all(I(r, "unreachable_sequences") == 0 for r in rescued),
      "  control: guides the rescue law exempts lose nothing",
      [I(r, "unreachable_sequences") for r in rescued], "all zero")
claim(all(I(r, "unreachable_sequences") > 0 for r in atrisk),
      "  control: every at-risk guide loses something",
      min(I(r, "unreachable_sequences") for r in atrisk), "> 0")
claim(len(rescued) == 3 and len(atrisk) == 13, "rescued and at-risk split",
      (len(rescued), len(atrisk)), (3, 13))

# ------------------------------------------------- the complete search policy, not one class
# The decisive objection to a deletion-only audit is that the datasets enumerated more than one
# class. They also enumerated ungapped candidates, and the arithmetic settles the rest without
# the reference: deleting a spacer endpoint and restoring any flanking base changes one column,
# so a locus an endpoint word reaches within three is within four ungapped whatever base sits
# beside it. The bound is therefore universal over the unknown flank, and it is checked here
# against the substitution budget these datasets' own ungapped rows exhibit.
by_spacer = {meta[g]["spacer"]: g for g in meta}
if BOUND:
    # The budget is a property of the ENUMERATED rows, which carry label 0. A measured positive
    # is in the file because it was observed, not because the enumerator reached it, so mixing
    # the two infers seven where the policy is six.
    ung_mm, ung_pos = Counter(), Counter()
    for fn, (gc, oc, lc) in COLUMNS.items():
        with open(os.path.join(sys.argv[1], fn), newline="") as fh:
            for r in csv.DictReader(fh):
                g, o = r[gc], r[oc]
                sp = g.replace("_", "").replace("-", "")[:SPACER]
                if sp not in by_spacer or "_" in g or "_" in o:
                    continue
                proto = o.replace("-", "").upper()[:-3]
                if len(proto) != SPACER:
                    continue
                d = sum(1 for x, y in zip(sp, proto) if x != y)
                (ung_mm if str(r[lc]).startswith("0") else ung_pos)[d] += 1
    claim(max(ung_mm) == 6,
          "the enumerated ungapped rows reach exactly six substitutions", max(ung_mm), 6)
    claim(max(ung_pos) > 6,
          "  control: a measured positive goes past it, so the two must not be mixed",
          max(ung_pos), "> 6")

    # The two endpoints are not alike, and the PAM decides it. At a fixed PAM boundary the
    # ungapped protospacer is the 20 nt immediately 5' of it, that is b + t for an unknown base
    # one position further out. Restoring a position-1 deletion prepends that base at the far
    # end and leaves the boundary alone. A terminal deletion sits against the PAM, so restoring
    # it would move the boundary; holding the boundary fixed instead compares the guide against
    # a shifted window, and nothing makes that small.
    P1 = {(g, t) for (g, t) in UN if reach[g][t] == {1}}
    P20 = {(g, t) for (g, t) in UN if reach[g][t] == {SPACER}}
    claim(len(P1) == 5235 and len(P20) == 5653 and len(P1) + len(P20) == len(UN),
          "endpoint-only sequences split by which endpoint",
          (len(P1), len(P20)), (5235, 5653))

    w1 = max(max(ham(spacer_of[g], b + t) for b in "ACGT") for g, t in P1)
    claim(w1 <= 4, "  position 1: within four ungapped at the fixed PAM, every flank", w1, "<= 4")
    b20 = min(min(ham(spacer_of[g], b + t) for b in "ACGT") for g, t in P20)
    claim(b20 > 6, "  position 20: beyond six at the fixed PAM, every flank", b20, "> 6")

    # and the datasets' DNA-bulge class cannot stand in either: at the same PAM a DNA bulge
    # reads the 21 nt immediately 5' of it, two unknown context bases then t, with one DNA
    # position unpaired against the 20 nt guide
    def dna_bulge_best(g, t):
        sp, best = spacer_of[g], 99
        for c in "ACGT":
            for b in "ACGT":
                s21 = c + b + t
                for j in range(21):
                    d = ham(sp, s21[:j] + s21[j + 1:])
                    if d < best:
                        best = d
        return best

    closest = min(dna_bulge_best(g, t) for g, t in P20)
    claim(closest == 4,
          "  position 20: no DNA-bulge placement reaches three, over every context",
          closest, 4)

    # The limit that must travel with the number. Moving the PAM boundary one base outward makes
    # the ungapped twenty-mer t + b instead of b + t, and there every terminal pair is close. So
    # what the policy cannot align is the register, not the neighbourhood, and saying otherwise
    # would repeat the conflation that made an earlier version of this file claim 10,887.
    shifted = [min(ham(spacer_of[g], t + b) for b in "ACGT") for g, t in P20]
    claim(max(shifted) == 3,
          "  but one base out every terminal pair is ungapped at three",
          "min %d max %d" % (min(shifted), max(shifted)), "all 3")
    claim(sum(1 for d in shifted if d <= 6) == len(P20),
          "    so all of them are eligible at the shifted boundary",
          "%d of %d" % (sum(1 for d in shifted if d <= 6), len(P20)), len(P20))

    # the same split at the locus level, from the scan statistic
    l1 = sum(int(r["loci_position_1_only"]) for r in lrows)
    l20 = sum(int(r["loci_position_20_only"]) for r in lrows)
    claim(l1 == 6407, "  loci reachable only at position 1, which are rescued", l1, 6407)
    claim(l20 == 28865,
          "  loci reachable only at position 20, which are not", l20, 28865)
    claim(l1 + l20 == sum(int(r["unreachable_loci"]) for r in lrows),
          "    and the two account for every endpoint-only locus",
          "%d + %d" % (l1, l20), sum(int(r["unreachable_loci"]) for r in lrows))


# ------------------------------------------------- what the files hold anyway
# The sites are mostly not missing. They are in the published files under a different
# alignment, so the omission costs the deposited data an alignment class and not the sites.
# This is the conservative reading and it is recomputed from the CSVs, not read from a flag.
if BOUND:
    held, kinds = defaultdict(set), defaultdict(set)
    for fn, (gc, oc, _) in COLUMNS.items():
        with open(os.path.join(sys.argv[1], fn), newline="") as fh:
            for r in csv.DictReader(fh):
                sp = r[gc].replace("_", "").replace("-", "")[:SPACER]
                if sp not in by_spacer:
                    continue
                gid, o = by_spacer[sp], r[oc]
                proto = o.replace("-", "").replace("_", "").upper()[:-3]
                kind = {19: "gapped19", 20: "ungapped20", 21: "dnabulge21"}.get(len(proto))
                if kind is None:
                    continue
                for i in range(len(proto) - 19 + 1):
                    held[gid].add(proto[i:i + 19])
                    kinds[(gid, proto[i:i + 19])].add(kind)
    inf = {(g, s) for (g, s) in UN if s in held[g]}
    claim(len(inf) == 10854, "unreachable sequences the files hold anyway", len(inf), 10854)
    claim(len(UN) - len(inf) == 34, "  and the ones absent outright",
          len(UN) - len(inf), 34)
    ung = sum(1 for k in inf if kinds[k] == {"ungapped20"})
    dna = sum(1 for k in inf if kinds[k] == {"dnabulge21", "ungapped20"})
    claim(ung == 6267, "    held as an ungapped twenty-mer only", ung, 6267)
    claim(dna == 4587, "    held inside a DNA-bulge twenty-one-mer too", dna, 4587)
    claim(ung + dna == len(inf), "    the two routes account for all of them",
          ung + dna, len(inf))
    # control: no unreachable sequence is held by a GAPPED row, which is the whole point
    gap = sum(1 for k in inf if "gapped19" in kinds[k])
    claim(gap == 0, "    control: none is held by a deletion row", gap, 0)
    # control: the deposited flag column agrees with this recomputation row by row
    flag = {(r["guide_id"], r["off_target_19mer"]): r["in_published_files"] for r in seqs}
    dis = sum(1 for k in UN if flag.get(k) != ("yes" if k in inf else "no"))
    claim(dis == 0, "    control: every deposited flag matches the recomputation", dis, 0)

bad = VERDICTS.count(False)
print()
print("%d claims checked, %d reproduced, %d failed" % (len(VERDICTS), len(VERDICTS) - bad, bad))
if bad:
    sys.exit(1)
print("Across %d enumerated RNA-bulge rows and sixteen guides, the deletion never sits at" % TOTAL)
print("spacer position 1 or 20, so those two alignment classes are absent from the published")
print("data. The two endpoints then part company at the PAM: the 6,407 loci reachable only at")
print("position 1 are rescued ungapped within four, while the 28,865 reachable only at position")
print("20 stay beyond six for every flank and admit no DNA-bulge placement within three.")
if not BOUND:
    print()
    print("UNBOUND: the counts above were read from the deposited census, which is a derived")
    print("summary and cannot certify itself. Re-run with the directory holding the two public")
    print("CSVs to recompute them from the files and bind this result. Exiting 2, not 0.")
    sys.exit(2)
