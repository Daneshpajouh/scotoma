#!/usr/bin/env python3
"""Per-haplotype reference-absent fraction. The partition step alone, no guide scoring.

Counts, for one haplotype assembly, how many SpCas9 target windows it carries and how many of
those the reference universe does not contain. This is the quantity the coverage result reports,
computed here on many more haplotypes so its ancestry ordering can be tested on the campaign's own
measure rather than on a different one.

Codes are counted per record and released, so peak memory is set by the largest record rather than
by the genome. Absence is defined from the haplotype's side, windows the haplotype carries that the
reference lacks, so a haplotype lacking a chromosome cannot inflate its own rate.
"""
import argparse
import csv
import json
import os
import time

import numpy as np

from scotoma.coding import iter_fasta, oriented_sites
from scotoma.universe import contains, load_universe


def read_manifest_row(path, task_id):
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if int(row["task_index"]) == task_id:
                return row
    raise SystemExit(f"no manifest row {task_id} in {path}")


def count_absence(fasta_path, universe_codes, pam):
    """Return (records, assembled_bp, windows_total, windows_absent)."""
    records = 0
    assembled_bp = 0
    windows_total = 0
    windows_absent = 0
    for _name, bases in iter_fasta(fasta_path):
        sites = oriented_sites(bases, pam=pam)
        assembled_bp += int(bases.size)
        if sites.size == 0:
            continue
        records += 1
        windows_total += int(sites.size)
        windows_absent += int((~contains(universe_codes, sites)).sum())
        del sites
    return records, assembled_bp, windows_total, windows_absent


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--task-id", type=int, default=None)
    ap.add_argument("--asm-dir", required=True)
    ap.add_argument("--universe", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pam", default="NGG")
    a = ap.parse_args(argv)

    task_id = a.task_id
    if task_id is None:
        task_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
    row = read_manifest_row(a.manifest, task_id)

    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, f"ABS_{row['donor']}_{row['hap']}.json")
    if os.path.exists(out) and os.path.getsize(out) > 50:
        print(f"[skip] already done: {out}", flush=True)
        return 0

    fasta = os.path.join(a.asm_dir, row["file"])
    if not os.path.exists(fasta):
        raise SystemExit(f"GATE FAIL missing assembly {fasta}")

    started = time.time()
    universe_codes, meta = load_universe(a.universe)
    stored_pam = meta.get("pam")
    if stored_pam is not None and stored_pam != a.pam:
        raise SystemExit(
            f"GATE FAIL universe PAM scope {stored_pam!r} does not match requested {a.pam!r}"
        )
    print(f"universe {meta.get('label')}: {universe_codes.size:,} codes, "
          f"recorded PAM {stored_pam or '<unrecorded, asserted ' + a.pam + '>'}, "
          f"loaded in {time.time() - started:.0f}s", flush=True)

    records, assembled_bp, total, absent = count_absence(fasta, universe_codes, a.pam)
    if total == 0:
        raise SystemExit(f"GATE FAIL no target windows found in {fasta}")

    fraction = absent / total
    present_fraction = 1.0 - fraction
    record = {
        "donor": row["donor"],
        "hap": row["hap"],
        "label": row["label"],
        "file": row["file"],
        "pam": a.pam,
        "records_with_sites": records,
        "assembled_bp": assembled_bp,
        "windows_total": total,
        "windows_absent": absent,
        "absent_fraction": fraction,
        "absent_percent": 100.0 * fraction,
        "present_fraction": present_fraction,
        "universe": meta.get("label"),
        "universe_codes": int(universe_codes.size),
        "wall_seconds": round(time.time() - started, 1),
    }
    # A human haplotype against a human reference must be overwhelmingly present. Report rather
    # than refuse, so an odd assembly is visible in the aggregate instead of silently missing.
    record["present_fraction_above_0.90"] = bool(present_fraction > 0.90)

    with open(out, "w") as fh:
        json.dump(record, fh, indent=2)
    print(json.dumps(record), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
