"""How many candidate sites does admitting one single-base bulge add?

A candidate is a genomic locus for a guide, keyed (guide, chrom, position, strand). The
ungapped set comes from the validated Cas-OFFinder run that reproduced our own enumerator
exactly through five mismatches, so both sides of the comparison are produced by the same tool and
the same reference bytes; nothing here depends on our enumerator. At six it is higher by four, the
three genomic-N windows and the one both-strand window that our universe excludes by convention.
"""
import json, os, sys, collections

# The path scrub that prepared this release wrote the shell form "${SCRATCH:?...}" into these
# three strings. Python does not expand it, so the script stopped on a literal path and could not
# be run as shipped. Read the variable the way the .sh files in this directory do, and fail with
# their message when it is unset.
SCRATCH = os.environ.get("SCRATCH")
if not SCRATCH:
    sys.exit("set SCRATCH to your own scratch directory")

MM = int(sys.argv[1]) if len(sys.argv) > 1 else 4
W = sys.argv[2] if len(sys.argv) > 2 else SCRATCH + "/bulge_census"
UNGAPPED = os.environ.get("UNGAPPED", SCRATCH + "/cof_h2h/out_ngg_mm6.txt")  # <=6; filtered to <=MM below
PANEL = os.environ.get("PANEL", SCRATCH + "/umbra_run/PANELREF_GRCh38_NGG.json")

rows = json.load(open(PANEL))["rows"]
gi_of_spacer = {r["spacer"]: r["guide_index"] for r in rows}
ours = {r["guide_index"]: r["sites"][str(MM)] for r in rows}

# ungapped loci, from the validated run, at this threshold
U = collections.defaultdict(set)
with open(UNGAPPED) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6 or int(p[5]) > MM:
            continue
        U[gi_of_spacer[p[0][:20]]].add((p[1], p[2], p[4]))
print(f"ungapped loci at <=​{MM}: {sum(len(v) for v in U.values()):,} "
      f"(our enumerator reports {sum(ours.values()):,} sites)")

B = {k: collections.defaultdict(set) for k in ("dna", "rna")}
for kind in ("dna", "rna"):
    owners = json.load(open(f"{W}/idx_{kind}.json"))
    n = 0
    with open(f"{W}/out_{kind}.txt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6 or int(p[5]) > MM:
                continue
            n += 1
            for gi, sp, pos in owners[p[0]]:
                B[kind][gi].add((p[1], p[2], p[4]))
    print(f"{kind}: {n:,} rows -> {sum(len(v) for v in B[kind].values()):,} distinct (guide, locus)")

# A bulged 24mer at position p and an ungapped 23mer at p+1 are different coordinates but
# heavily overlapping loci, and one amplicon covers both. So report two figures: STRICT, where
# a new coordinate is a new candidate, and MERGED, where a bulged hit within TOL bp of an
# ungapped hit for the same guide and strand counts as already covered. The truth for assay
# burden lies between them, and they bound it.
TOL = 3

def covered(pos_by_strand, chrom, pos, strand):
    """Is there an ungapped hit for this guide within TOL bp, same chromosome and strand?"""
    xs = pos_by_strand.get((chrom, strand))
    if not xs:
        return False
    import bisect
    i = bisect.bisect_left(xs, pos)
    for j in (i - 1, i):
        if 0 <= j < len(xs) and abs(xs[j] - pos) <= TOL:
            return True
    return False

out = {"threshold": MM, "tolerance_bp": TOL, "per_guide": [], "totals": {}}
tot_u = tot_add_d = tot_add_r = tot_add_any = tot_add_merged = 0
for gi in sorted(ours):
    u = U[gi]
    d = B["dna"][gi] - u
    r = B["rna"][gi] - u
    both = (B["dna"][gi] | B["rna"][gi]) - u
    pos_by_strand = {}
    for c, pp, st in u:
        pos_by_strand.setdefault((c, st), []).append(int(pp))
    for k in pos_by_strand:
        pos_by_strand[k].sort()
    merged = {x for x in both if not covered(pos_by_strand, x[0], int(x[1]), x[2])}
    tot_u += len(u); tot_add_d += len(d); tot_add_r += len(r)
    tot_add_any += len(both); tot_add_merged += len(merged)
    out["per_guide"].append({"guide_index": gi, "ungapped": len(u), "add_dna": len(d),
                             "add_rna": len(r), "add_any": len(both), "add_merged": len(merged),
                             "factor": (len(u) + len(both)) / len(u) if u else None,
                             "factor_merged": (len(u) + len(merged)) / len(u) if u else None})
def q(xs):
    xs = sorted(xs)
    n = len(xs)
    return (xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2, xs[0], xs[-1])

f_med, f_min, f_max = q([g["factor"] for g in out["per_guide"] if g["factor"]])
m_med, m_min, m_max = q([g["factor_merged"] for g in out["per_guide"] if g["factor_merged"]])
out["totals"] = {"ungapped": tot_u, "add_dna": tot_add_d, "add_rna": tot_add_r,
                 "add_any": tot_add_any, "add_merged": tot_add_merged,
                 "panel_factor_strict": (tot_u + tot_add_any) / tot_u,
                 "panel_factor_merged": (tot_u + tot_add_merged) / tot_u,
                 "per_guide_factor_strict": [f_med, f_min, f_max],
                 "per_guide_factor_merged": [m_med, m_min, m_max]}
print()
print(f"ungapped candidates            {tot_u:,}")
print(f"added by DNA bulge             {tot_add_d:,}")
print(f"added by RNA bulge             {tot_add_r:,}")
print(f"added by either, STRICT        {tot_add_any:,}")
print(f"added by either, MERGED +-{TOL}bp  {tot_add_merged:,}")
print()
print(f"panel factor  STRICT {out['totals']['panel_factor_strict']:.3f}x   "
      f"MERGED {out['totals']['panel_factor_merged']:.3f}x")
print(f"per-guide     STRICT median {f_med:.3f}x (min {f_min:.3f}, max {f_max:.3f})")
print(f"per-guide     MERGED median {m_med:.3f}x (min {m_min:.3f}, max {m_max:.3f})")
json.dump(out, open(f"{W}/BULGE_CENSUS_mm{MM}.json", "w"), indent=2)
print(f"\nwrote {W}/BULGE_CENSUS_mm{MM}.json")
