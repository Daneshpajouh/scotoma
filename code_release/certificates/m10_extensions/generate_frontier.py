#!/usr/bin/env python3
"""Regenerate the frontier table and vector figure from pinned exact records.

This reads small result records and checks their H3 provenance; it never repeats
large maximum-load sweeps. All input paths are explicit.
"""
from __future__ import annotations
import hashlib
import json
import os
import subprocess
from fractions import Fraction
from pathlib import Path

DRAFT = Path(__file__).resolve().parents[1]
# 2026-11-02: ROOT was DRAFT.parents[2], which lands two directories above the project, so every
# input read failed with FileNotFoundError. The load guard below runs first, so on a busy machine
# the script refused before it ever reached a read and the path bug stayed invisible; the build
# never ran this producer at all, which is why the rendered table and the figure disagreed for a
# day. Find the root by a marker instead of a hop count, so a staged copy at any depth also works.
def _find_root(start):
    override = os.environ.get('FRONTIER_ROOT')
    if override:
        return Path(override).resolve()
    for candidate in (start, *start.parents):
        if (candidate / 'build_all.sh').is_file() and (candidate / 'code_release').is_dir():
            return candidate
    raise SystemExit('cannot locate the project root above ' + str(start)
                     + '; set FRONTIER_ROOT to it')


ROOT = _find_root(DRAFT)
OUT = DRAFT / 'generated'
RESULT_ROOT = 'results_2026-10-29/ZR3_exact_tau_n21_23/cluster_returned/full_1576400'
H3_PATH = 'results_2026-10-29/ZX2_verify_ZR3_closures/cluster_returned/full_1665238/runs/full/H3_summary.json'
# 2026-11-03: the N=20 ladder (fir job 2105334) closed r=14 after the six cells above were
# returned, so that cell's producer record and its independent verification sit in their own
# trees rather than in the two returned directories the six share. Its verification is a
# single-cell result.json, not a summary, because the summary artefact above is the 31 October
# evidence and is left byte-identical.
RESULT_ROOTS = {
    (20, 14): 'results_2026-10-29/ZR3_v2_ladder_n20_2105334/cluster_returned/ladder_N20',
    (20, 15): 'results_2026-10-29/ZR3_v2_solve_n20r15_2150456/cluster_returned',
    (21, 19): 'results_2026-10-29/ZR3_v2_ladders_2078168/cluster_returned/ladder_N21',
    (22, 23): 'results_2026-10-29/ZR3_v2_ladders_2078168/cluster_returned/ladder_N22',
}
CELL_DIRS = {(20, 15): 'solver'}
H3_SINGLES = {
    (20, 14): 'results_2026-10-29/ZX3_n20_r14/cluster_returned/runs/full/cell_20_14_10/result.json',
    (20, 15): 'results_2026-10-29/ZX4_n20_r15/cluster_returned/runs/full/cell_20_15_10/result.json',
    (21, 19): 'results_2026-10-29/ZX5_ladders/cluster_returned/runs/full/cell_21_19_10/result.json',
    (22, 23): 'results_2026-10-29/ZX5_ladders/cluster_returned/runs/full/cell_22_23_10/result.json',
}
SEEDS = ((20, 11), (20, 12), (21, 12), (21, 13), (22, 20), (22, 21), (20, 14),
         (20, 15), (21, 19), (22, 23))
CLASS_COUNTS = (649539, 1771470, 1948617, 5314410, 4649045868, 10460353203, 12754584,
                33480783, 1549681956, 31381059609)
SOURCES = {}


def read(relative, parse=True):
    path = ROOT / relative
    if path.stat().st_size > 10_000_000:
        raise ValueError('Unexpectedly large input: ' + relative)
    data = path.read_bytes()
    SOURCES[relative] = hashlib.sha256(data).hexdigest()
    print(f'Read {relative} ({len(data)} bytes)', flush=True)
    return json.loads(data) if parse else data


def digit_block(value):
    digits = str(value)
    return '\n'.join(r'\noindent\texttt{' + digits[i:i+64] + r'}\par'
                     for i in range(0, len(digits), 64))


def plot(cells):
    from reportlab.pdfgen import canvas
    from reportlab.lib.colors import HexColor, white
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    import reportlab
    font_file = Path(reportlab.__file__).parent / 'fonts/Vera.ttf'
    pdfmetrics.registerFont(TTFont('FrontierSans', str(font_file)))
    pdf = OUT / 'budget_ten_frontier.pdf'
    chart = canvas.Canvas(str(pdf), pagesize=(464.4, 140.4), invariant=1,
                          initialFontName='FrontierSans', initialFontSize=10)
    chart.setTitle('Budget-ten single-deletion covering frontier')
    chart.setAuthor('Amirhossein Daneshpajouh, Megan Fowler, Kay C. Wiese')
    chart.setCreator('scripts/generate_frontier.py')
    chart.setFont('FrontierSans', 10)
    chart.drawCentredString(250, 128, 'Run count r')
    for r in range(1, 25):
        chart.drawCentredString(50 + (r-1)*17.6, 107, str(r))
    chart.saveState()
    chart.translate(10, 60)
    chart.rotate(90)
    chart.drawCentredString(0, 0, 'Guide length n')
    chart.restoreState()
    for n in range(20, 25):
        y = 90 - (n-20)*15
        chart.drawRightString(36, y-3.4, str(n))
        chart.setStrokeColor(HexColor('#dddddd'))
        chart.setLineWidth(0.4)
        chart.line(44, y, 460, y)
    def mark(code, x, y):
        chart.setLineWidth(0.8)
        if code == 'I':
            chart.setStrokeColor(HexColor('#ac591f'))
            chart.setFillColor(HexColor('#ac591f'))
            p = chart.beginPath()
            p.moveTo(x, y+3.1)
            p.lineTo(x-3.1, y-2.6)
            p.lineTo(x+3.1, y-2.6)
            p.close()
            chart.drawPath(p, stroke=1, fill=1)
        else:
            chart.setStrokeColor(HexColor('#176b80' if code == 'C' else '#444444'))
            chart.setFillColor(HexColor('#176b80') if code == 'C' else white)
            chart.circle(x, y, 2.7, stroke=1, fill=1)
    for c in cells:
        if c['display'] != '-':
            mark(c['display'], 50+(c['r']-1)*17.6, 90-(c['n']-20)*15)
    chart.setFillColor(HexColor('#222222'))
    for code, label, x in [('C', 'Certificate', 104), ('I', 'Integer argument', 218), ('O', 'Open', 356)]:
        mark(code, x, 10)
        chart.setFillColor(HexColor('#222222'))
        chart.drawString(x+9, 6.6, label)
    chart.showPage()
    chart.save()
    subprocess.run(['pdftoppm', '-singlefile', '-r', '150', '-png', str(pdf),
                    str(OUT / 'budget_ten_frontier')], check=True, timeout=20,
                   capture_output=True)
    print('Generated vector frontier figure and preview', flush=True)


def main():
    print(f'Frontier: nice={os.getpriority(os.PRIO_PROCESS, 0)}, load={os.getloadavg()[0]:.2f}', flush=True)
    if os.getpriority(os.PRIO_PROCESS, 0) < 15 or os.getloadavg()[0] > 35:
        raise SystemExit('Require nice >=15 and one-minute load <=35; wait and retry.')
    OUT.mkdir(exist_ok=True)
    gaps = read('code_release/gap_census/exact_gaps.json')
    upper = read('code_release/upper_rows/extended_exact_table.json')
    records = {x['r']: x for x in upper}
    for rec in gaps:
        if rec['m'] == 10:
            assert rec['r'] - 1 < Fraction(rec['tau_star']) < rec['r']
            assert rec['minimum_coverage'] == rec['maximum_load'] == '1'
    assert {x['r'] for x in gaps if x['m'] == 10} == set(range(5, 11))
    for r, rec in records.items():
        assert Fraction(rec['tau_upper']) <= r - Fraction(2277, 2000)
        if r <= 15:
            assert rec['tau_star'] == rec['tau_lower'] == rec['tau_upper']
    closures = []
    for filename in ('packing_N22_r16_m10_d2_c2_all.json', 'packing_N23_r24_m10_d2_c1_all.json'):
        rec = read('code_release/certificates/m10_closure/' + filename)
        assert rec['m'] == 10 and Fraction(rec['total_mass']) > rec['r'] - 1
        assert Fraction(rec['max_load']) == 1 and rec['status'] == 'VERIFIED-COMPUTATIONALLY'
        closures.append(rec)
    h3 = read(H3_PATH)
    assert h3['status'] == 'SIX_H3_VERIFIED'
    verified = {(v['N'], v['r']): v for v in h3['results']}
    for key, single_path in H3_SINGLES.items():
        single = read(single_path)
        assert (single['N'], single['r']) == key and key not in verified
        verified[key] = single
    assert set(verified) == set(SEEDS)
    seeds = []
    for (N, r), count in zip(SEEDS, CLASS_COUNTS):
        root = RESULT_ROOTS.get((N, r), RESULT_ROOT)
        # the direct solve of (20,15,10) wrote its records straight into solver/, with no
        # enclosing cell directory, because it was not run through the ladder driver
        leaf = CELL_DIRS.get((N, r), f'cell_{N}_{r}_10')
        result_path = f'{root}/{leaf}/RESULT.json'
        certificate_path = f'{root}/{leaf}/decision.certificate.json'
        rec = read(result_path)
        read(certificate_path, parse=False)
        check = verified[(N, r)]
        assert (rec['N'], rec['r'], rec['m']) == (N, r, 10)
        assert rec['status'] == 'CLOSED' and rec['K'] == r
        assert Fraction(rec['tau_lower']) == Fraction(check['total_mass']) > r - 1
        assert check['support_inside_U'] and check['closure_verified']
        assert check['execution_state'] == 'H3'
        assert check['global_maximum_status'] == 'EXHAUSTIVE_EXACT_PASS'
        assert Fraction(check['exact_maximum']) == 1
        assert check['classes_checked'] == count
        assert check['certificate_sha256'] == SOURCES[certificate_path]
        for key in ('scaled_one_class', 'moved_one_class'):
            assert check['negative_controls'][key]['status'] == 'REJECTED'
        seeds.append({'N': N, 'n': N+1, 'r': r, 'm': 10,
                      'tau_lower': rec['tau_lower'], 'tau_upper': rec['tau_upper'],
                      'centre_classes': count, 'source': result_path,
                      'certificate': certificate_path,
                      'certificate_sha256': check['certificate_sha256'],
                      'state': 'H3', 'maximum_load': '1',
                      'native_elapsed_seconds': check['elapsed_seconds'],
                      'principal_estimate_bytes': check['resource_estimate']['principal_estimate_bytes']})
    assert Fraction('10.7983263360') < Fraction(seeds[0]['tau_lower'])
    assert Fraction(seeds[0]['tau_upper']) < Fraction('10.8053096925')
    cells = []
    base_reach = max(rec['r'] for rec in gaps if rec['m'] == 10)
    reach = {n: base_reach for n in range(20, 25)}
    for rec in closures:
        n = rec['N'] + 1
        reach[n] = max(reach[n], rec['r'])
    for rec in seeds:
        reach[rec['n']] = max(reach[rec['n']], rec['r'])
    for n in range(21, 25):
        reach[n] = max(reach[n], reach[n-1])
    for n in range(20, 25):
        for r in range(1, 25):
            code, reason = 'O', 'No established closure claimed.'
            if r > n:
                code, reason = '-', 'Inadmissible: r > n.'
            elif r <= reach[n]:
                code = 'C'
                reason = ('Exact full-weight threshold or matched fractional pair.' if n == 20 else
                          'H3 packing and certificate-preserving run-count step-down.' if n < 24 else
                          'Released diagonal closure and run-count step-down.')
            elif n == 20 and r == 11:
                code, reason = 'I', 'Five-case integer obstruction, K(19,11,10)=11.'
            c = {'n': n, 'N': n-1, 'r': r, 'm': 10, 'display': code,
                 'reason': reason, 'fractional_rounding_refuted': n == 20 and 11 <= r <= 20}
            if n == 20 and r in records:
                c['fractional_record'] = {k: records[r].get(k) for k in
                                         ('tau_star', 'tau_lower', 'tau_upper', 'resolution')}
            cells.append(c)
    table = r'''% Regenerate with scripts/generate_frontier.py; do not edit by hand.
\begingroup\small
\setlength{\tabcolsep}{9pt}
\newcommand{\fc}{\textcolor{teal!70!black}{$\bullet$}}
\newcommand{\fiarg}{\textcolor{orange!70!black}{$\blacktriangle$}}
\newcommand{\fo}{\ensuremath{\circ}}
\begin{tabular}{rlll}
\toprule
$n$ & \fc\ Certificate & \fiarg\ Integer argument & \fo\ Open \\
\midrule
'''
    def range_text(selected):
        if not selected:
            return '--'
        runs = [c['r'] for c in selected]
        assert runs == list(range(min(runs), max(runs)+1))
        expression = (str(runs[0]) if len(runs) == 1 else
                      ','.join(map(str, runs)) if len(runs) == 2 else
                      str(runs[0]) + r'\!:\!' + str(runs[-1]))
        if all(c['fractional_rounding_refuted'] for c in selected):
            expression += r'^{\rm R}'
        return '$' + expression + '$'
    for n in range(20, 25):
        groups = [range_text([c for c in cells if c['n'] == n and c['display'] == status])
                  for status in ('C', 'I', 'O')]
        table += str(n) + ' & ' + ' & '.join(groups) + r' \\' + '\n'
    table += r'''\bottomrule
\end{tabular}
\endgroup
'''
    (OUT / 'budget_ten_frontier.tex').write_text(table)
    values = r'''% Regenerate with scripts/generate_frontier.py.
\subsection{Exact fractional values behind the budget-ten frontier}
\label{app:frontier-values}
\verified For $N=19,m=10$, there are exact optima $t_r$ at
$11\le r\le15$ and certified rational intervals at $16\le r\le20$.
Every exact upper bound is at most $r-2277/2000<r-1$, refuting the
fractional route to $K=r$. The integer theorem closes $r=11$; $K$ remains
open at $r=12,\ldots,20$. Define $t_r=A_r/B_r$ using the digit strings
below. Successive lines are concatenated with no padding or omitted digits.
The source is \path{code_release/upper_rows/extended_exact_table.json}.
Complete fractions, input hashes and H3 seed provenance are preserved in
\path{generated/budget_ten_frontier.json} beside the manuscript source.
'''
    for r in range(11, 16):
        frac = Fraction(records[r]['tau_star'])
        values += f'\n\\paragraph{{Run count {r}.}} $A_{{{r}}}=$\n{{\\small\n'
        values += digit_block(frac.numerator) + '\n}\n'
        values += f'\\noindent $B_{{{r}}}=$\n{{\\small\n'
        values += digit_block(frac.denominator) + '\n}\n'
    values += r'''
For $16\le r\le20$ only the certified bound
$\tau^*\le r-2277/2000$ is asserted, not an exact optimum.
\subsection{Exact endpoints of the length-twenty-one bracket}
\verified The bracket for $\tau^*(20,11,10)$ is an interval, not an
exact optimum. Its exact lower and upper endpoints $A/B$ are:
'''
    for key in ('tau_lower', 'tau_upper'):
        f = Fraction(seeds[0][key])
        values += '\n\\paragraph{' + key.replace('_', ' ') + r'.} $A=$' + '\n{\\small\n'
        values += digit_block(f.numerator) + '\n}\n\\noindent $B=$\n{\\small\n'
        values += digit_block(f.denominator) + '\n}\n'
    (OUT / 'frontier_tau_values.tex').write_text(values)
    plot(cells)
    manifest = {'budget': 10, 'execution_state': 'H3',
                'verification_scope': 'Pinned returned exhaustive checks; no local numerical replay.',
                'sources_sha256': SOURCES, 'verified_seed_records': seeds, 'cells': cells}
    (OUT / 'budget_ten_frontier.json').write_text(json.dumps(manifest, indent=2) + '\n')
    tsv = ['n\tr\tstatus\tfractional_rounding_refuted\treason']
    tsv += [f"{c['n']}\t{c['r']}\t{c['display']}\t{c['fractional_rounding_refuted']}\t{c['reason']}" for c in cells]
    (OUT / 'budget_ten_frontier.tsv').write_text('\n'.join(tsv) + '\n')
    print('Generated table, figure, exact-value supplement and 120-cell manifest', flush=True)


if __name__ == '__main__':
    main()
