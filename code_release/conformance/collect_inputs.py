#!/usr/local/bin/python3 -I
"""Snapshot local read-only evidence into this deposit; never execute source-tree scripts."""
import csv
import hashlib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    records = []

    def copy(source, dest):
        source = ROOT / source
        target = HERE / dest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        records.append({"source": str(source.relative_to(ROOT)), "snapshot": dest,
                        "bytes": target.stat().st_size,
                        "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})

    audit = "results_2026-10-29/T1_cross_tool_audit/cas_baselines"
    copy("results_2026-10-29/HOSTILE_REVIEW_2026-11-01/r3_published_hole/upstream/cas-offinder-bulge",
         "upstream/cas-offinder-bulge")
    copy("paper_v7/sections/03_method.tex", "inputs/authority/03_method.tex")
    copy("paper_v7/main.tex", "inputs/authority/main.tex")
    for name in ("main.bbl", "sections/01_introduction.tex", "sections/07_discussion.tex"):
        copy("paper_v7/" + name, "inputs/authority/" + Path(name).name)
    copy("code_release/tool_audit/selftest/guides.tsv", "inputs/panel_original.tsv")
    copy("code_release/tool_audit/selftest/fixture_rna.fa", "inputs/fixture_rna.fa")
    copy("code_release/tool_audit/audit_selected_counts.tsv", "inputs/audit_selected_counts.tsv")
    copy("code_release/tool_audit/audit_all110_native_counts.tsv", "inputs/audit_all110_native_counts.tsv")
    copy("code_release/tool_conformance/per_guide_counts.tsv", "inputs/per_guide_counts.tsv")
    for version in ("2.7.1", "2.9.0"):
        for source, dest in (
            (f"code_release/tool_conformance/crispritz_{version}_rna1.targets.tsv", f"inputs/crispritz_{version}.tsv"),
            (f"code_release/tool_audit/crispritz_audit/{'run' if version == '2.9.0' else 'run_2_7_1'}/rna1.targets.txt",
             f"inputs/crispritz_{version}_audit.tsv")):
            copy(source, dest)
    for name in ("input.txt", "input_bulge.txt"):
        copy(f"code_release/crisprseek_executed/{name}", f"inputs/crisprseek/{name}")
    for name in ("rna.tsv", "rna_reverse.tsv"):
        copy(f"results_2026-10-29/A5_casoffinder2_and_usage/run/{name}", f"inputs/rust/{name}")
    copy("results_2026-10-29/A5_casoffinder2_and_usage/lead_cli_check/out_rna_local.tsv",
         "inputs/rust/cli.tsv")
    for name in ("rna_loop_3.0.0b.txt", "rna_loop_3.0.0b3.txt"):
        copy(f"results_2026-10-29/A2_casoffinder_v3_tags/{name}", f"inputs/static/{name}")
    guide_lines = (HERE / "inputs/panel_original.tsv").read_text().splitlines()
    guides = dict(line.split("\t") for line in guide_lines)
    (HERE / "inputs/guides.tsv").write_text("guide_id\tguide\n" + "\n".join(guide_lines) + "\n")
    for gid, guide in guides.items():
        for suffix in ("bulge.txt", "native_release.tsv", "native_develop.tsv"):
            name = f"{gid}_{guide}_{suffix}"
            copy(f"{audit}/raw/{name}", f"inputs/cas_raw/{name}")

    inventory = []
    dirs = ("tool_audit", "tool_conformance", "crisprseek_executed", "genome_endtoend", "clinical_guide")
    for directory in dirs:
        for source in sorted((ROOT / "code_release" / directory).rglob("*")):
            if source.suffix not in (".tsv", ".txt") or not source.is_file():
                continue
            if directory in ("genome_endtoend", "clinical_guide"):
                reason = "not a pattern library: selected genomic output or summary; target universe not exhausted"
            elif source.name in ("input.txt", "input_bulge.txt"):
                reason = "checked as CRISPRseek input / pooled emitted pattern file"
            elif "targets" in source.name:
                reason = "checked from raw CRISPRitz alignments, including duplicate deposited copies"
            elif source.name == "guides.tsv":
                reason = "panel input, checked against all per-guide evidence"
            else:
                reason = "count/source summary or input; no emitted word set to certify"
            with source.open(newline="") as handle:
                count = sum(1 for row in csv.reader(handle, delimiter="\t") if row)
            inventory.append({"source": str(source.relative_to(ROOT)), "nonempty_rows_including_header": count,
                              "sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "disposition": reason})
    (HERE / "inputs/manifest.json").write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
    (HERE / "inputs/output_inventory.json").write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"snapshot_files": len(records), "inventoried_files": len(inventory), "guides": len(guides)}, sort_keys=True))


if __name__ == "__main__":
    main()
