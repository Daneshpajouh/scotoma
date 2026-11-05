#!/usr/bin/env python3
"""Typeset the proposed sentence and exact removal at the current body geometry."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from integrity import sha,require
from render_findings import proposed_sentence,removal_clause
import json
import os
import re
import subprocess


def main():
    out=HERE/'layout'
    out.mkdir(exist_ok=True)
    summary=json.loads((HERE/'results/summary.json').read_text())
    manuscript=HERE.parents[2]/'paper_v7/main.tex'
    frozen=out/'body_format_source.tex'
    original=frozen.read_bytes() if frozen.exists() else manuscript.read_bytes()
    (out/'body_format_source.tex').write_bytes(original)
    # Copy the exact font, width and margin setup preceding the unrelated hyperlink setup.
    preamble=original.decode().split('\\usepackage[colorlinks=true,allcolors=blue]{hyperref}',1)[0]
    preamble+='\\usepackage{url}\n\\newcount\\proposedlines\n\\newcount\\removedlines\n'
    def box(content,register,label):
        return ('\\setbox0=\\vbox{\\hsize=\\textwidth\\noindent '+content+
                '\\par\\global\\'+register+'=\\prevgraf}\n'+
                '\\typeout{'+label+'=\\the\\'+register+'}\n\\box0\n')
    tex=preamble+'\\begin{document}\n'+box(proposed_sentence(summary),'proposedlines','PROPOSED_LINES')
    tex+='\\bigskip\n'+box(removal_clause(),'removedlines','REMOVED_LINES')+'\\end{document}\n'
    (out/'sentence_layout.tex').write_text(tex)
    env=os.environ.copy()
    for name in ('TEXMFVAR','TEXMFCONFIG','TEXMFCACHE'):
        dest=out/name.lower(); dest.mkdir(exist_ok=True); env[name]=str(dest)
    with (out/'build.log').open('w') as log:
        completed=subprocess.run(['nice','-n','10','pdflatex','-interaction=nonstopmode','-halt-on-error',
                                  '-output-directory',str(out),str(out/'sentence_layout.tex')],
                                 stdout=log,stderr=subprocess.STDOUT,env=env,cwd=out,timeout=60)
    require(completed.returncode==0,'Layout build failed; see layout/build.log')
    log=(out/'sentence_layout.log').read_text()
    p=int(re.search(r'PROPOSED_LINES=(\d+)',log).group(1))
    r=int(re.search(r'REMOVED_LINES=(\d+)',log).group(1))
    require(p<=4 and r>=p,'Proposed sentence exceeds the line budget')
    require('Overfull \\hbox' not in log,'Layout contains an overfull line')
    bbox=subprocess.run(['pdftotext','-bbox-layout',str(out/'sentence_layout.pdf'),'-'],
                        capture_output=True,text=True,check=True,timeout=30)
    (out/'sentence_layout.html').write_text(bbox.stdout)
    subprocess.run(['pdftoppm','-scale-to','1500','-png','-singlefile',str(out/'sentence_layout.pdf'),
                    str(out/'sentence_layout')],check=True,timeout=30)
    receipt=dict(proposed_lines=p,removed_lines=r,maximum_proposed_lines=4,build_exit_code=completed.returncode,
                 format_source_sha256=sha(out/'body_format_source.tex'),
                 tex_sha256=sha(out/'sentence_layout.tex'),pdf_sha256=sha(out/'sentence_layout.pdf'),
                 scope='Isolated sentence/removal at the current 11pt body font and text width; no full manuscript rebuild.')
    (out/'layout_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print(json.dumps(receipt,sort_keys=True))


if __name__=='__main__':
    main()
