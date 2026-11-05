#!/usr/local/bin/python3 -I
"""Classify deposited word projections; Table 1 alone determines execution labels."""
import csv
import hashlib
import json
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent
SUITE = runpy.run_path(str(HERE / "suite.py"))
GUIDES = SUITE["read_guides"](HERE / "inputs/guides.tsv")


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def raw_cas(kind):
    emitted = {}
    for gid, guide in GUIDES.items():
        path = HERE / "inputs/cas_raw" / f"{gid}_{guide}_{kind}"
        lines = path.read_text().splitlines()
        patterns = [line.split()[0] for line in (lines[2:] if kind == "bulge.txt" else lines)]
        # Exact layout of single RNA deletion with DNA padding budget one, trailing PAM.
        emitted[gid] = [p[2:-3] for p in patterns if len(p) == len(guide) + 4
                        and p.startswith("NN") and p.endswith("NNN")
                        and set(p[2:-3]) <= set("ACGT")]
        if kind != "bulge.txt":
            # Parser dumps serialize one row per placement metadata entry, while
            # the enumerator's pattern map contains each string key only once.
            emitted[gid] = sorted(set(emitted[gid]))
    return emitted


def alignment_words(path):
    result = {gid: [] for gid in GUIDES}
    by_spacer = {g: gid for gid, g in GUIDES.items()}
    for row in read_tsv(path):
        if row.get("#Bulge type", row.get("Bulge type")) != "RNA":
            raise ValueError(f"unexpected non-RNA row in {path}")
        guide, dna = row["crRNA"], row["DNA"]
        gid = by_spacer.get(guide[:-3])
        if gid is None or len(guide) != 23 or not guide.endswith("NNN"):
            raise ValueError("unknown or malformed alignment guide")
        if row.get("#Id", gid) != gid or row["Direction"] not in ("+", "-"):
            raise ValueError("alignment guide ID or strand inconsistent")
        if len(dna) != len(guide) or dna.count("-") != 1 or dna.index("-") >= 20:
            raise ValueError("invalid single-deletion alignment shape")
        if row["Mismatches"] != "0" or row["Bulge Size"] != "1" or dna[-2:] != "GG":
            raise ValueError("alignment budget or PAM inconsistent")
        i = dna.index("-")
        word = dna.replace("-", "")[:-3]
        if guide[:i] + guide[i + 1:20] != word:
            raise ValueError("gap placement does not spell the reported exact word")
        result[gid].append(word)
    return result


def main():
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    authority = (HERE / "inputs/authority/03_method.tex").read_text()
    table = authority.split("\\midrule", 1)[1].split("\\bottomrule", 1)[0]
    evidence = [line.strip() for line in table.splitlines() if "\\tabularnewline" in line]
    if not all(text in table for text in ("Executed wrapper dry-run", "Static source and loop-bound emulation",
                                         "Executed parser;", "Executed CPU CLI;", "Executed in R, default route",
                                         "Native RNA search", "pipelines untested", "native runtime untested",
                                         "Static alignment source")):
        raise ValueError("Table 1 authority changed; review execution labels")
    (out / "table1_authority.json").write_text(json.dumps(evidence, indent=2) + "\n")
    scoreboard = []
    reports = {}

    def add(route, name, emitted, state, scope):
        report = SUITE["panel"](GUIDES, emitted, budget=0)
        reports[route] = report
        record = {k: v for k, v in report.items() if k != "per_guide"}
        record.update(route=route, tool=name, evidence_state=state, scope=scope,
                      word_set_status="invalid" if report["invalid_words"] else
                      "incomplete" if report["missing_words"] else "exact")
        scoreboard.append(record)
        (out / f"{route}.emitted.json").write_text(json.dumps(emitted, indent=2, sort_keys=True) + "\n")
        (out / f"{route}.report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    wrapper = raw_cas("bulge.txt")
    add("wrapper", "Cas-OFFinder wrapper 2016", wrapper, "executed", "deposited official dry-run patterns")
    for suffix, route, name in (("native_release.tsv", "native_b3", "Cas-OFFinder 3.0.0b3"),
                                ("native_develop.tsv", "native_develop", "Cas-OFFinder develop")):
        add(route, name, raw_cas(suffix), "executed",
            "deposited parser map keys (placement metadata collapsed); full native search untested")
    static_loop = (HERE / "inputs/static/rna_loop_3.0.0b.txt").read_text()
    if "j = 1; j < sline[0].find_last_not_of('N') + 1 - i" not in static_loop:
        raise ValueError("unexpected early-beta static loop")
    static_emitted = {gid: sorted({g[:i] + g[i + 1:] for i in range(1, len(g) - 1)})
                      for gid, g in GUIDES.items()}
    add("native_early_static", "Cas-OFFinder 3.0.0b / 3.0.0b2", static_emitted, "static",
        "loop-bound emulation; native parsers not executed")
    for version in ("2.7.1", "2.9.0"):
        emitted = alignment_words(HERE / f"inputs/crispritz_{version}.tsv")
        other = alignment_words(HERE / f"inputs/crispritz_{version}_audit.tsv")
        if emitted != other:
            raise ValueError("CRISPRitz duplicate deposits disagree")
        add("crispritz_" + version, "CRISPRitz " + version, emitted, "executed",
            "deposited native RNA fixture; projected row duplicates reported separately")
    for route, path, scope in (("rust_cli", "cli.tsv", "deposited CPU CLI RNA fixture"),
                               ("rust_forward", "rna.tsv", "deposited serial CPU harness, forward fixture"),
                               ("rust_reverse", "rna_reverse.tsv", "deposited serial CPU harness, reverse fixture")):
        add(route, "Cas-OFFinder 2 Rust 1.0.0", alignment_words(HERE / "inputs/rust" / path), "executed", scope)

    # Pooled pattern files discard guide ownership. Never credit a shared word to all guides.
    seek_input = (HERE / "inputs/crisprseek/input.txt").read_text().splitlines()
    if {line.split()[0][:-3] for line in seek_input[2:]} != set(GUIDES.values()):
        raise ValueError("CRISPRseek input is not exactly the audited panel")
    if seek_input[1].split()[1:] != ["0", "1"]:
        raise ValueError("CRISPRseek input does not request the audited bulge budgets")
    seek_lines = (HERE / "inputs/crisprseek/input_bulge.txt").read_text().splitlines()
    seek_patterns = [line.split()[0] for line in seek_lines[2:]]
    seek_words = {p[1:-3] for p in seek_patterns if p.startswith("N")}
    union = {word for words in wrapper.values() for word in words}
    if seek_words != union or any(line.split()[1] != "4" for line in seek_lines[2:]):
        raise ValueError("CRISPRseek pooled written file differs from wrapper per-guide union/budget")
    if len(seek_patterns) != len(set(seek_patterns)):
        raise ValueError("CRISPRseek pooled patterns duplicated")
    if {p[:-3] for p in seek_patterns if not p.startswith("N")} != set(GUIDES.values()):
        raise ValueError("CRISPRseek ungapped patterns disagree with the panel")
    add("crisprseek", "CRISPRseek 1.52.0 default bulge route", wrapper, "executed",
        "Table 1: executed in R; pooled file equals wrapper union; per-guide ownership from wrapper dumps")
    (out / "crisprseek_pool.json").write_text(json.dumps({"pooled_deletion_words": len(seek_words),
        "per_guide_attribution": "wrapper snapshots, not membership credit to every possible owner"}, indent=2) + "\n")

    # Validate every deposited panel count against derived words; summaries are cross-checks only.
    checks = 0
    tool_routes = {"Cas-OFFinder-bulge v1.2 beta (2016 script)": "wrapper",
                   "Cas-OFFinder 3.0.0 develop": "native_develop", "Cas-OFFinder tag 3.0.0b3": "native_b3",
                   "CRISPRitz 2.7.1 native TST": "crispritz_2.7.1", "CRISPRitz 2.9.0 native TST": "crispritz_2.9.0"}
    for filename in ("audit_all110_native_counts.tsv", "per_guide_counts.tsv", "audit_selected_counts.tsv"):
        for row in read_tsv(HERE / "inputs" / filename):
            gid = row["guide_id"]
            if row["guide"] != GUIDES[gid]:
                raise ValueError("count table guide disagrees")
            if row["tool"] in tool_routes:
                report = reports[tool_routes[row["tool"]]]["per_guide"][gid]
                field = "patterns_emitted" if row["patterns_emitted"] != "N/A" else "distinct_deletion_words_recovered"
                if int(row[field]) != report["valid_distinct"] or int(row["deficit"]) != len(report["missing_words"]):
                    raise ValueError("deposited count table disagrees with raw word projection")
            elif row["tool"] == "CHOPOFF symbolic RNA-gap component":
                r = SUITE["REF"]["run_count"](GUIDES[gid])
                if int(row["source_derived_deletion_words"]) != r:
                    raise ValueError("CHOPOFF source count differs from run count")
            checks += 1

    for name, scope in (("CRISPRme 2.1.14 / CRISPRme+ 2.8.0", "delegated CRISPRitz engines checked separately; pipelines untested"),
                        ("CHOPOFF 1.0.3", "symbolic source-component counts match; no native full-word output deposited"),
                        ("SWOffinder 1.0.1", "alignment architecture; native runtime untested; count criterion inapplicable"),
                        ("CALITAS 7a3456c", "alignment architecture; native runtime untested; count criterion inapplicable")):
        scoreboard.append({"tool": name, "evidence_state": "static", "status": "not_certifiable",
                           "word_set_status": "N/A", "scope": scope})
    fixture = "".join(line.strip() for line in (HERE / "inputs/fixture_rna.fa").read_text().splitlines()
                      if not line.startswith(">"))
    expected = {w for g in GUIDES.values() for w in SUITE["REF"]["required_words"](g)}
    absent = [word for word in expected if word + "AGG" not in fixture]
    if absent:
        raise ValueError("fixture is missing required planted words")
    summary = {"guides": len(GUIDES), "required_guide_specific_words": sum(SUITE["REF"]["run_count"](g) for g in GUIDES.values()),
               "pooled_required_words": len(expected), "fixture_bases": len(fixture),
               "fixture_planted_words_absent": len(absent), "summary_rows_crosschecked": checks,
               "reference_source_lines": len((HERE / "reference.py").read_text().splitlines())}
    (out / "panel_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (out / "scoreboard.json").write_text(json.dumps(scoreboard, indent=2, sort_keys=True) + "\n")
    lines = ["| Route | Evidence from Table 1 | Distinct / required | Missing (guides) | Repeated word rows | Verdict |",
             "|---|---|---:|---:|---:|---|"]
    for row in scoreboard:
        if "required_words" not in row:
            cells = ("N/A", "N/A", "N/A", "not certifiable")
        else:
            cells = (f"{row['valid_distinct_words']} / {row['required_words']}",
                     f"{row['missing_words']} ({row['guides_incomplete']})", str(row["redundant_rows"]),
                     row["status"] + ("; exact word set" if row["word_set_status"] == "exact" else ""))
        lines.append(f"| {row['tool']} ({row.get('route', 'source')}) | {row['evidence_state']} | " + " | ".join(cells) + " |")
    (out / "SCOREBOARD.md").write_text("\n".join(lines) + "\n\nScopes and per-guide missing words: `scoreboard.json` and `*.report.json`.\n")
    print(json.dumps(summary, sort_keys=True))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
