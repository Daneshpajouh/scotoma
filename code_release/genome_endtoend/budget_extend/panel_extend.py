#!/usr/bin/env python3
"""The 110-guide panel at budgets 5, 6 and 7, sharded by guide across a Slurm array.

Why sharded. One guide at budget seven produced 403 MB and 1.89 M loci, so the whole panel in one
input would be about 44 GB per build and an unmeasured memory footprint for the pattern set. Shards
bound both and let a single failure cost one shard instead of the run.

Gates before any search, same as the clinical run: both binaries must hash to the values the build
gate recorded, and the staged reference must be the 25 primary records at the recorded byte count.

Counts loci under the paper's convention and, for each guide, reports how many corrected-only loci
are attributable to the terminal deletion word, because that mechanism claim has to be retested at
every budget rather than inherited from budget four.
"""
import os
import csv, hashlib, json, os, subprocess, sys, time
from pathlib import Path

ROOT = Path(os.environ.get("PANEL_ROOT", Path(__file__).resolve().parent))
FASTA = Path(os.environ["PANEL_FASTA"])      # staged GRCh38 primary assembly, 25 records
BIN = Path(os.environ["PANEL_BIN"])          # directory holding both pinned binaries
PATTERN = "N" * 21 + "GG"
EXPECT_BIN = {"released": "b8991a89143763809f1684517613eadc6e179a362044001da5a152ce9754849c",
              "corrected": "cee9231169ec193bd7e7756981dd96411e786033b75922d6a9c4fb2df50350fa"}
EXPECT_FASTA_BYTES = 3132404952
BUDGETS = [5, 6, 7]


def sha256(p, chunk=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def gates():
    for name, want in EXPECT_BIN.items():
        got = sha256(BIN / name)
        if got != want:
            raise SystemExit(f"GATE FAILED: {name} is not the published build")
    fa = sorted(f for f in os.listdir(FASTA) if f.endswith((".fa", ".fasta")))
    tot = sum(os.path.getsize(FASTA / f) for f in fa)
    if len(fa) != 25 or tot != EXPECT_FASTA_BYTES:
        raise SystemExit("GATE FAILED: staged reference is not the published universe")
    print(f"gates ok: binaries match, {len(fa)} records, {tot:,} bytes", flush=True)


def locus(row):
    p = int(row["Position"])
    span = len(row["DNA"].replace("-", ""))
    return (row["#Id"], row["Chromosome"], row["Direction"], p + span if row["Direction"] == "+" else p)


def by_guide(path):
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out.setdefault(row["#Id"], {}).setdefault(locus(row), []).append(row)
    return out


def terminal_attributable(rows, spacer, budget):
    """True if NO word from positions 1..n-1 reaches this locus within budget."""
    n = len(spacer)
    native = {spacer[:i] + spacer[i + 1:] for i in range(n - 1)}
    for r in rows:
        core = r["DNA"].replace("-", "")[:-3]
        if r["Bulge type"] == "X":
            if len(core) == n and sum(a != b for a, b in zip(core, spacer)) <= budget:
                return False
            continue
        if len(core) != n - 1:
            continue
        if min(sum(a != b for a, b in zip(core, w)) for w in native) <= budget:
            return False
    return True


def main():
    idx = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    nsh = int(os.environ.get("NSHARDS", "8"))
    guides = [l.split("\t") for l in (ROOT / "panel_guides.txt").read_text().split("\n") if l.strip()]
    mine = [(g, s) for i, (g, s) in enumerate(guides) if i % nsh == idx]
    out = ROOT / f"shard_{idx}"
    out.mkdir(parents=True, exist_ok=True)
    gates()
    print(f"shard {idx} of {nsh}: {len(mine)} guides", flush=True)

    rec = {"shard": idx, "nshards": nsh, "host": os.uname().nodename,
           "job": os.environ.get("SLURM_JOB_ID"), "guides": {}}
    for mm in BUDGETS:
        inp = out / f"input_mm{mm}.txt"
        with inp.open("w") as fh:
            fh.write(f"{FASTA}\n{PATTERN} 0 1\n")
            for g, s in mine:
                fh.write(f"{s}NNN {mm} {g}\n")
        res = {}
        for build in ("released", "corrected"):
            tsv = out / f"{build}_mm{mm}.tsv"
            t0 = time.time()
            r = subprocess.run([str(BIN / build), str(inp), "C", str(tsv)],
                               capture_output=True, text=True)
            if r.returncode != 0:
                print(f"  mm={mm} {build} FAILED rc={r.returncode}", flush=True)
                res[build] = None
                continue
            res[build] = by_guide(tsv)
            print(f"  mm={mm} {build:9s} {time.time()-t0:6.1f}s "
                  f"{os.path.getsize(tsv)/1e6:8.1f} MB", flush=True)
        if not res.get("released") or not res.get("corrected"):
            continue
        for g, s in mine:
            ra = set(res["released"].get(g, {}))
            rb = res["corrected"].get(g, {})
            only = set(rb) - ra
            attr = sum(1 for k in only if terminal_attributable(rb[k], s, mm))
            rec["guides"].setdefault(g, {})[f"mm{mm}"] = {
                "released": len(ra), "corrected": len(rb),
                "corrected_only": len(only), "released_only": len(ra - set(rb)),
                "terminal_attributable": attr}
    (out / "record.json").write_text(json.dumps(rec, indent=2, sort_keys=True))
    print(f"SHARD_{idx}_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
