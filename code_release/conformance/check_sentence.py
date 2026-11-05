#!/usr/local/bin/python3 -I
"""Compile local copies and measure the proposed sentence in the paper's typography."""
import hashlib
import json
from pathlib import Path
import re
import runpy
import shutil
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SUPPORT = runpy.run_path(str(HERE / "run_support.py"))
run = SUPPORT["run"]
OLD = ("No audited route emits the minimum\nfamily and nothing else: two fall below its word set and two attain it but serialise it in more\nrows.")
NEW = ("The emitter, suite and patch are in\n"
       "\\path{HOSTILE_REVIEW_2026-11-01/R19_conformance_standard/conformance/}.")


def main():
    area = HERE / "build/manuscript"
    source_dir = HERE / "inputs/authority"
    files = ("main.tex", "main.bbl", "sections/01_introduction.tex", "sections/03_method.tex", "sections/07_discussion.tex")
    sources = {}
    reports = {}
    for version in ("before", "after"):
        directory = area / version
        directory.mkdir(parents=True, exist_ok=True)
        for name in files:
            source = source_dir / Path(name).name
            target = directory / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            sources[name] = hashlib.sha256(source.read_bytes()).hexdigest()
        method = directory / "sections/03_method.tex"
        if method.read_text().count(OLD) != 1:
            raise ValueError("replacement sentence not uniquely located")
        if version == "after":
            method.write_text(method.read_text().replace(OLD, NEW))
        for pass_number in (1, 2):
            run(f"paper-{version}-pass{pass_number}", ["/Library/TeX/texbin/pdflatex", "-no-shell-escape",
                 "-halt-on-error", "-interaction=nonstopmode", "main.tex"], cwd=directory)
        pdf = directory / "main.pdf"
        info = run(f"paper-{version}-pdfinfo", ["/opt/homebrew/bin/pdfinfo", pdf])
        run(f"paper-{version}-text", ["/opt/homebrew/bin/pdftotext", "-layout", pdf, directory / "main.txt"])
        run(f"paper-{version}-bbox", ["/opt/homebrew/bin/pdftotext", "-bbox", pdf, directory / "main.html"])
        pages = (directory / "main.txt").read_text().split("\f")
        if not pages[-1].strip():
            pages.pop()
        references_page = next(i + 1 for i, page in enumerate(pages) if re.search(r"^References\s*$", page, re.M))
        reports[version] = {"total_pdf_pages": len(pages), "body_pages": references_page - 2,
                            "references_physical_page": references_page,
                            "pdfinfo_pages": int(re.search(r"^Pages:\s+(\d+)", info.stdout, re.M).group(1))}
        if reports[version]["total_pdf_pages"] != reports[version]["pdfinfo_pages"]:
            raise ValueError("PDF page counts disagree")

    # Measure the complete replacement paragraph in the actual width, font and margin settings.
    measure = area / "measure"
    measure.mkdir(exist_ok=True)
    preamble = re.split(r"(?m)^\\begin\{document\}", (source_dir / "main.tex").read_text(), maxsplit=1)[0]
    probe = (preamble + "\n\\begin{document}\n"
             "\\setbox0=\\vbox{\\hsize=\\textwidth\\noindent " + OLD.replace("\n", " ") +
             "\\par\\global\\count0=\\prevgraf}\n\\typeout{OLDLINES=\\the\\count0}\n"
             "\\setbox1=\\vbox{\\hsize=\\textwidth\\noindent " + NEW.replace("\n", " ") +
             "\\par\\global\\count1=\\prevgraf}\n\\typeout{NEWLINES=\\the\\count1}\n"
             "\\noindent Existing sentence:\\par\\box0\\bigskip\\noindent Proposed sentence:\\par\\box1\n"
             "\\end{document}\n")
    (measure / "measure.tex").write_text(probe)
    result = run("sentence-line-measure", ["/Library/TeX/texbin/pdflatex", "-no-shell-escape", "-halt-on-error",
                                          "-interaction=nonstopmode", "measure.tex"], cwd=measure)
    old_lines = int(re.search(r"OLDLINES=(\d+)", result.stdout).group(1))
    new_lines = int(re.search(r"NEWLINES=(\d+)", result.stdout).group(1))
    if new_lines > 3 or new_lines > old_lines or reports["after"]["body_pages"] > reports["before"]["body_pages"]:
        raise ValueError("replacement does not meet line/page budget")
    run("sentence-render", ["/opt/homebrew/bin/pdftoppm", "-scale-to", "1200", "-singlefile", "-png",
                            measure / "measure.pdf", measure / "measure"])
    document = ET.parse(area / "after/main.html")
    pages = document.findall(".//{*}page")
    selected_page = next(i + 1 for i, page in enumerate(pages)
                         if "The emitter, suite and patch are in" in " ".join(word.text or "" for word in page.findall(".//{*}word")))
    run("sentence-context-render", ["/opt/homebrew/bin/pdftoppm", "-f", selected_page, "-l", selected_page,
                                    "-scale-to", "1600", "-singlefile", "-png", area / "after/main.pdf", area / "after/sentence-page"])
    report = {"old_sentence_tex": OLD, "new_sentence_tex": NEW, "old_printed_lines": old_lines,
              "new_printed_lines": new_lines, "paper_builds": reports, "source_sha256": sources,
              "rendered_context_physical_page": selected_page,
              "caveat": "local replacement proposal only; paper_v7 was not edited"}
    (HERE / "results/sentence_check.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
