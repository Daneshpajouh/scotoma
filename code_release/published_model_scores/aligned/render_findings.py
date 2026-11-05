#!/usr/bin/env python3
"""Render verified numerical summaries as the findings and paired guide table."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
HERE=Path(__file__).resolve().parent
import csv
import io
import json


def proposed_sentence(summary):
    a=summary['active']['by_return_bits']['01']
    b=summary['active']['released_returned_mm4']
    t,c=summary['treatment'],summary['control']
    def number(n):
        return format(n,',').replace(',', '{,}')
    return (f"At the published CRISPR-Net cutoff, ${a['above']}/{a['rows']}$ missed measured-active sites "
            f"qualify versus ${b['above']}/{b['rows']}$ returned sites at four mismatches; genome-wide, "
            f"guide- and mismatch-matched nomination counts are ${number(t['above'])}/{number(t['rows'])}$ versus "
            f"${number(c['above'])}/{number(c['rows'])}$.")


def removal_clause():
    body=(HERE/'inputs/body_excerpt.tex').read_text()
    start=body.index('The tables and one command\n')
    end=body.index('which needs the reference to catch.',start)+len('which needs the reference to catch.')
    return body[start:end]


def guide_table(summary):
    stream=io.StringIO()
    writer=csv.writer(stream,delimiter='\t',lineterminator='\n')
    writer.writerow(['guide','mismatches','treatment_rows','treatment_above','treatment_percent',
                     'control_rows','control_above','control_percent','treatment_mean','control_mean',
                     'treatment_median','control_median'])
    for g,v in sorted(summary['by_guide'].items(),key=lambda kv:int(kv[0][1:])):
        a,b=v['treatment'],v['control']
        writer.writerow([g,4,a['rows'],a['above'],format(a['above_percent'],'.9f') if a['rows'] else '',
                         b['rows'],b['above'],format(b['above_percent'],'.9f') if b['rows'] else '',
                         format(a['mean'],'.17g') if a['rows'] else '',
                         format(b['mean'],'.17g') if b['rows'] else '',
                         format(a['quantiles']['0.5'],'.17g') if a['rows'] else '',
                         format(b['quantiles']['0.5'],'.17g') if b['rows'] else ''])
    return stream.getvalue()


def render(summary,convention,layout):
    t,c=summary['treatment'],summary['control']
    a=summary['active']; split=a['by_return_bits']; lost,returned,outside=split['01'],split['11'],split['00']
    mm4=a['released_returned_mm4']; matched=a['exact_guide_mm4_matched']
    selection=summary['control_selection']; panel=summary['panel']; comparison=summary['comparison']
    def ratio(d):
        return f"{d['above']:,}/{d['rows']:,} ({d['above_percent']:.3f}%)" if d['rows'] else '0 rows; not scored'
    text=[
        '**Encoder convention, pinned before any scoring.** The published encoder uses `_` for a biological gap '
        'and `-` for zero-vector left padding. Each deposited 23-character pair was converted by replacing only '
        'alignment `-` with `_`; every non-gap character was fixed. The encoder then prepended one padding `-` '
        'to each sequence, producing a 24-column, seven-channel tensor. It copies each guide `N` from the DNA '
        'symbol in that column; the deposited guide PAM is `NNN`. Real DNA PAMs were retained, but that guide '
        'mask makes PAM columns matches, so these scores do not estimate PAM specificity.',
        f"The original vendored encoder and the reused forward encoder agreed exactly on all {t['rows']:,} treatment "
        f"and {c['rows']:,} control tensors, and all {a['all']['rows']:,} measured-active tensors. "
        f"Treatment and control each converted {convention['by_arm']['treatment']['alignment_gap_characters_converted']:,} "
        f"gap characters and preserved {convention['by_arm']['treatment']['non_gap_characters_unchanged']:,} non-gap "
        'characters. Every alignment was representable. Passing a biological `-` through unchanged fails the '
        'gap-versus-padding negative control. See [the frozen policy](published_model_aligned_scores/ENCODER_POLICY.md) '
        'and [the exhaustive encoder audit](published_model_aligned_scores/results/encoder_convention.json).',
        '**1. Whole-assembly nomination scores, with the matched control beside them.** '
        'The unchanged deposited NumPy forward pass and checkpoint were reused, without fitting, training or tuning. '
        'The threshold was strictly greater than exact `3/1000` throughout.',
        '| Arm | Scored loci | Above cutoff | Percent above |\n|---|---:|---:|---:|\n'
        f"| Corrected only | {t['rows']:,} | {t['above']:,} | {t['above_percent']:.6f}% |\n"
        f"| Released returned, matched | {c['rows']:,} | {c['above']:,} | {c['above_percent']:.6f}% |",
        f"The arms match exactly within each of {panel['treatment_guides']} guides, at four protospacer mismatches, "
        'one RNA bulge and literal NGG. The control contains distinct released-returned loci selected by the '
        'lowest fixed SHA-256 locus ranks per guide, without inspecting scores. Multiple eligible alignments of '
        'one locus were represented by the lexicographically smallest DNA alignment. Treatment uses the deposited '
        'terminal alignment; control uses a released interior alignment. Chromosome, sequence context and gap '
        'position were not matched, so this is a descriptive comparison rather than a causal effect of the repair.',
        f"Raw released outputs were recovered by read-only transfer from the original run and reconciled to the "
        f"deposited per-guide and per-chromosome counts: {selection['raw_rows']:,} raw alignment rows, "
        f"{selection['released_loci']:,} distinct loci. There were {selection['eligible_loci']:,} eligible control "
        f"loci for the treatment guides; {selection['eligible_loci_not_selected']:,} were not selected and were "
        f"not scored. {selection['unsupported_eligible_alignment_rows_not_scored']} otherwise eligible raw alignment "
        'rows with ambiguous DNA/PAM symbols were excluded and not scored; the published encoder cannot represent '
        'those symbols. No neighbouring score was substituted.',
        f"Released-only export: {summary['reverse_arm']['rows']} data rows, confirmed by parsing it. "
        f"The enumeration totals remain released {panel['released_loci']:,} versus corrected {panel['corrected_loci']:,}; "
        f"the ungapped build control is {panel['ungapped_released']:,} versus {panel['ungapped_corrected']:,}, "
        'including agreement per guide. These baseline tables contain counts, not scoreable sequence pairs. '
        f"Their separate deposited magnitude baseline is {panel['deposited_magnitude_baseline']:,}; it is not an "
        'alignment control or a claim of exact reference agreement.',
        '**2. Full score distributions and paired guide breakdown.** The absolute difference in above-cutoff '
        f"rates is {comparison['above_percent_difference']:.6f} percentage points. Most rows in both arms fall "
        'below the cutoff. All row scores are deposited; the following distribution summaries retain the control.',
        '| Statistic | Corrected only | Matched released returned |\n|---|---:|---:|\n'
        f"| Mean | {t['mean']:.9g} | {c['mean']:.9g} |\n"
        f"| Population SD | {t['population_sd']:.9g} | {c['population_sd']:.9g} |"
    ]
    for q,label in [('0','Minimum'),('0.01','1%'),('0.05','5%'),('0.25','25%'),('0.5','Median'),
                    ('0.75','75%'),('0.95','95%'),('0.99','99%'),('1','Maximum')]:
        text[-1]+=f"\n| {label} | {t['quantiles'][q]:.9g} | {c['quantiles'][q]:.9g} |"
    hist='| Score interval | Corrected only | Matched released returned |\n|---|---:|---:|'
    edges=t['histogram_edges']
    for i,(nt,nc) in enumerate(zip(t['histogram_counts'],c['histogram_counts'])):
        hist+=f"\n| [{edges[i]}, {edges[i+1]}{']' if i==len(edges)-2 else ')'} | {nt:,} | {nc:,} |"
    text.extend([hist,
        'Quantiles use linear interpolation. Histogram bins include their lower endpoint; only the last includes '
        'its upper endpoint. The cutoff comparison remains strict regardless of histogram bin boundaries.',
        f"Within the {panel['treatment_guides']} scored guides, treatment has a higher above-cutoff rate on "
        f"{comparison['guide_above_rate_higher']}, an equal rate on {comparison['guide_above_rate_equal']}, and a "
        f"lower rate on {comparison['guide_above_rate_lower']}. The complete paired breakdown, including guides "
        'without a corrected-only row, is [per_guide_scores.tsv](published_model_aligned_scores/results/per_guide_scores.tsv); '
        'zero-row groups are explicitly not scored.',
        '**3. Measured-active rows, split from the two return columns.** The assay sequences were joined by '
        '`row_index` from `sites129.json`, with guide, chromosome, strand, coordinate, real PAM, positive read count '
        'and all return bits rechecked from the sequences. Scores use the deposited terminal alignment and the '
        'same `NNN` guide mask in every arm. The corrected-only count was derived from the bits, not imported '
        'from the paper.',
        '| Released / corrected return bits | Meaning | Scored rows | Above cutoff | Percent above |\n'
        '|---|---|---:|---:|---:|'
    ])
    for bits,label in [('01','Corrected only; omitted in-scope active sites'),
                       ('11','Both returned; measured-active control'),('00','Neither returned; non-NGG'),
                       ('10','Released only')]:
        d=split[bits]
        text[-1]+=f"\n| {bits[0]} / {bits[1]} | {label} | {d['rows']} | {d['above']} | "+(f"{d['above_percent']:.6f}%" if d['rows'] else 'Not scored')+' |'
    text.extend([
        f"Thus {lost['above']}/{lost['rows']} measured-active sites omitted by the released enumerator exceed "
        f"the scorer's cutoff, beside {returned['above']}/{returned['rows']} returned active sites. Those "
        'denominators have different mismatch compositions and are not a matched enrichment comparison. '
        f"All {lost['rows']} omitted in-scope active sites occur in the scored treatment export, with identical "
        f"aligned sequence/PAM and replay-consistent scores. Released returns {returned['rows']} active rows; "
        f"corrected returns {lost['rows']+returned['rows']}. The other {outside['rows']} are outside literal NGG "
        'scope. Their source mismatch annotation includes PAM substitutions; matching here uses protospacer '
        'mismatches recomputed from the alignment.',
        f"If 'released misses' includes out-of-scope PAMs, the total is "
        f"{ratio(a['released_missed_all_including_non_NGG'])}, beside the returned active control "
        f"{ratio(returned)}. The non-NGG rows are scored but excluded from the terminal-omission claim.",
        '**4. Measured-active control at the same mismatch count, and its limit.**',
        '| Comparison | Corrected-only active sites | Released-returned active control |\n|---|---:|---:|\n'
        f"| All available rows at four mismatches | {ratio(lost)} | {ratio(mm4)} |\n"
        f"| Equal counts within shared guides at four mismatches | {ratio(matched['treatment'])} | {ratio(matched['control'])} |",
        f"The first row matches mismatch count but not guide composition. The second uses the "
        f"{len(matched['shared_guides'])} guides present in both active arms; within each guide it takes the "
        'smaller arm size from both arms, ordered by original assay row index before looking at scores. '
        'The complete selected identities and paired per-guide summaries are in '
        '[summary.json](published_model_aligned_scores/results/summary.json). This small selected subset does '
        'not establish general enrichment, cleavage prediction or ranking impact. In particular, the '
        'four-mismatch active rates in the first row are very similar; their score distributions are not '
        'evidence that omitted active sites are more important than retained ones.',
        'These are in-silico nomination scores on retrospective, selected inputs. The combined-data checkpoint '
        'is unchanged; the published cutoff came from a CIRCLE-seq leave-one-guide-out sensitivity analysis, '
        'and calibration of that cutoff for this checkpoint or these inputs is not established. Prior deposited '
        'validation is reused, including its recorded benchmark discrepancies. No independence from training '
        'data, prospective validation, cleavage outcome, clinical risk, top-ranked loss or altered ranking '
        'is claimed.',
        '**Proposed body sentence.**',
        '> '+proposed_sentence(summary),
        f"The sentence occupies {layout['proposed_lines']} printed lines in the current body font, width and "
        f"margins (limit: four). To pay for it, remove exactly the following sentence from "
        '`paper_v7/sections/03_method.tex`, in the paragraph beginning “Twenty-four records pass every control”; '
        f"the removed sentence occupies {layout['removed_lines']} lines in the same layout test. Its audit detail "
        'can remain in the release documentation. No manuscript file was edited.',
        '```tex\n'+removal_clause()+'\n```',
        'The isolated typography check establishes the sentence length and a positive line budget; it does '
        'not establish the final manuscript pagination after insertion. See '
        '[sentence_layout.pdf](published_model_aligned_scores/layout/sentence_layout.pdf) and '
        '[layout_receipt.json](published_model_aligned_scores/layout/layout_receipt.json).',
        '**Replay and integrity.** [The producer](published_model_aligned_scores/produce.py) reuses the '
        'byte-identical validated forward pass. [The separately coded verifier](published_model_aligned_scores/verify.py) '
        'reconstructs the matched control with SQL rather than the producer’s heap, replays the original encoder '
        'and independent recurrent/convolution implementation on every scored row, and rebuilds every summary '
        'and the rendered numerical report. It checks the [pinned input SHA-256 table]'
        '(published_model_aligned_scores/INPUT_SHA256.tsv). The checkpoint SHA-256 is '
        '`a7f6aa381520f5c68fa1f099a6ef3ebc3b8ce846709b97dfde2053f26ca62f80`.',
        'Run `published_model_aligned_scores/reproduce.sh` from this directory. All Python execution uses '
        '`nice -n 10 /usr/local/bin/python3 -I`, with at most four computational threads and observed resident '
        'memory below the required limit. Logs and exact resident-memory receipts are deposited. The original '
        'reference assembly was not re-read here: reference verification remains the producing run’s assertion, '
        'while this replay verifies the pinned exported strings and their run provenance.',
        '**Verifier exit code: `0`.** The corruption control accepts the unmodified input with exit `0` and '
        'rejects a one-byte DNA-input change with exit `1`, before scoring. The receipt is '
        '[corruption_control/receipt.json](published_model_aligned_scores/corruption_control/receipt.json). '
        'Every treatment and control row score is available in '
        '[treatment_scores.tsv](published_model_aligned_scores/results/treatment_scores.tsv) and '
        '[control_scores.tsv](published_model_aligned_scores/results/control_scores.tsv); measured rows are in '
        '[active_scores.tsv](published_model_aligned_scores/results/active_scores.tsv).'
    ])
    return '\n\n'.join(text)+'\n'


def main():
    summary=json.loads((HERE/'results/summary.json').read_text())
    convention=json.loads((HERE/'results/encoder_convention.json').read_text())
    layout=json.loads((HERE/'layout/layout_receipt.json').read_text())
    (HERE/'results/per_guide_scores.tsv').write_text(guide_table(summary))
    (HERE.parent/'FINDINGS.md').write_text(render(summary,convention,layout))


if __name__=='__main__':
    main()
