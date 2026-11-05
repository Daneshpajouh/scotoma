#!/usr/bin/env python3
"""Replay the deposited scorer on aligned treatment, matched control and assay rows."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from integrity import check_pins, sha, require, peak_rss
import csv
import hashlib
import heapq
import importlib.util
import io
import json
import tarfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction

sys.path.insert(0, str(HERE/'vendor'))
from model_numpy import NumpyModel, batches, encode, np


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def table(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def dump(path, data):
    path.write_text(json.dumps(data, sort_keys=True, indent=2)+'\n')


def write_table(path, fields, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows({k: r[k] for k in fields} for r in rows)


def cut(x):
    a, b = float(x).as_integer_ratio()
    return a*1000 > b*3


def distribution(values):
    p = np.asarray(values, dtype=np.float64)
    if not len(p):
        return {'rows': 0, 'above': 0}
    require(np.isfinite(p).all() and ((0 <= p) & (p <= 1)).all(), 'Invalid scores')
    passed = sum(map(cut, p))
    edges = ('0', '0.000001', '0.0001', '0.001', '0.003', '0.01', '0.1', '0.5', '1')
    hist = [0]*(len(edges)-1)
    rational_edges = list(map(Fraction, edges))
    for v in p:
        x = Fraction.from_float(float(v))
        for i in range(len(hist)):
            if rational_edges[i] <= x and (x < rational_edges[i+1] or
                                          i == len(hist)-1 and x <= rational_edges[i+1]):
                hist[i] += 1
                break
    require(sum(hist) == len(p), 'Histogram accounting differs')
    return dict(rows=len(p), above=passed, above_fraction=str(Fraction(passed, len(p))),
                above_percent=100*passed/len(p), mean=float(p.mean()),
                population_sd=float(p.std()),
                quantiles={str(q): float(np.quantile(p, q, method='linear'))
                           for q in (0, .01, .05, .25, .5, .75, .95, .99, 1)},
                histogram_edges=list(edges), histogram_counts=hist,
                nearest_cutoff_distance=float(min(abs(Fraction.from_float(float(v))-Fraction(3,1000))
                                                  for v in p)))


def prepare_treatment():
    rows = table(HERE/'inputs/corrected_only_alignments.part000.tsv')
    pg = table(HERE/'inputs/per_guide.tsv')
    counts = Counter(r['guide'] for r in rows)
    spacers = {r['guide']: r['spacer'] for r in pg}
    require(len(rows) == sum(int(r['corrected_only']) for r in pg), 'Treatment totals differ')
    require(all(counts[r['guide']] == int(r['corrected_only']) for r in pg), 'Treatment strata differ')
    seen = set()
    for i, r in enumerate(rows, 1):
        r['data_row_1based'] = i
        key = (r['guide'], r['chrom'], r['strand'], int(r['pam_coord']))
        require(key not in seen, 'Duplicate treatment locus')
        seen.add(key)
        require(r['crRNA'] == spacers[r['guide']]+'NNN', 'Treatment guide differs')
        dna = r['DNA']
        require(len(dna) == 23 and dna.count('-') == 1 and dna[19] == '-', 'Unexpected gap')
        mm = sum(a != b for a,b in zip(r['crRNA'][:20], dna[:20]) if b != '-')
        require(mm == int(r['mismatches']) == 4 and r['class'] == 'RNA' and
                r['rna_position_1based'] == '20' and dna[-2:] == 'GG', 'Unexpected treatment class')
        require(set(dna) <= set('ACGT-') and r['reference_verified'] == '1', 'Unexpected DNA alphabet')
        require(int(r['pam_coord']) == int(r['start'])+(22 if r['strand'] == '+' else 0),
                'Treatment PAM boundary differs')
    reverse = table(HERE/'inputs/released_only_alignments.part000.tsv')
    require(not reverse and all(int(r['released_only']) == 0 for r in pg), 'Reverse arm is nonempty')
    baseline = table(HERE/'inputs/ungapped_baseline_total.tsv')
    require({r['guide'] for r in baseline} == set(spacers), 'Baseline guide keys differ')
    for r in pg:
        b = next(x for x in baseline if x['guide'] == r['guide'])
        require(len({int(r[k]) for k in ('released_rna0','corrected_rna0',
                         'released_rna1_ungapped','corrected_rna1_ungapped')}) == 1,
                'Ungapped build control differs')
        require(all(r[k] == b[k] for k in ('released_rna0','corrected_rna0',
                         'released_rna1_ungapped','corrected_rna1_ungapped')), 'Baseline reconciliation differs')
    return rows, pg, counts, spacers


def prepare_active(spacers, treatment):
    selected = table(HERE/'inputs/active_sites.part000.tsv')
    assay = json.loads((HERE/'inputs/sites129.json').read_text())
    require(len(selected) == len(assay), 'Assay join length differs')
    require({int(r['row_index']) for r in selected} == set(range(len(assay))), 'Assay index differs')
    treatment_keys = {(r['guide'],r['chrom'],r['strand'],int(r['pam_coord'])):r for r in treatment}
    out = []
    for r in sorted(selected, key=lambda r:int(r['row_index'])):
        s = assay[int(r['row_index'])]
        require(all(str(s[k]) == r[k] for k in ('chrom','strand','start')), 'Assay identity differs')
        require(s['spacer'] == spacers[r['guide']] and s['reads'] > 0, 'Assay guide/activity differs')
        dna, spacer = s['aligned_target'], s['spacer']
        require(len(dna) == 23 and dna[19] == '-' and dna.count('-') == 1 and
                set(dna) <= set('ACGT-') and s['target'] == dna.replace('-',''), 'Assay alignment differs')
        costs = [sum(a != b for a,b in zip(spacer[:j]+spacer[j+1:],dna[:19])) for j in range(20)]
        pam_substitutions = sum(a != b for a,b in zip(s['sgRNA'][-2:],dna[-2:]))
        require(costs[19]+pam_substitutions == int(s['mm']), 'Assay mismatch label differs')
        pam_ok = dna[-2:] == 'GG'
        require((int(r['released_returned']),int(r['corrected_returned'])) ==
                (int(pam_ok and min(costs[:19]) <= 4),int(pam_ok and min(costs) <= 4)),
                'Assay return bits differ from sequences')
        require(int(r['pam_coord']) == int(r['start'])+(22 if r['strand'] == '+' else 0),
                'Assay PAM boundary differs')
        if r['released_returned'] == '0' and r['corrected_returned'] == '1':
            key = r['guide'],r['chrom'],r['strand'],int(r['pam_coord'])
            require(key in treatment_keys and treatment_keys[key]['DNA'] == dna, 'Missed assay locus absent from treatment')
        out.append(dict(r, data_row_1based=int(r['row_index'])+1, crRNA=spacer+'NNN', DNA=dna,
                        mismatches=costs[19], deposited_mm_including_PAM=int(s['mm']),
                        PAM_substitutions=pam_substitutions, reads=s['reads']))
    return out


def control_sample(targets, spacers):
    """Lowest SHA-256 locus keys per guide; score-independent sampling without replacement."""
    heaps = {g:[] for g in targets}
    kept = {g:{} for g in targets}
    eligible = Counter()
    all_counts = Counter()
    ungapped_counts = Counter()
    raw_rows = 0
    unsupported_rows = 0
    chrom_audit = []
    chroms = [f'chr{i}' for i in range(1,23)]+['chrX','chrY','chrM']
    with tarfile.open(HERE/'inputs/released_run.tar.gz') as archive:
        require(archive.extractfile('runs/collected/per_guide.tsv').read() ==
                (HERE/'inputs/per_guide.tsv').read_bytes(), 'Retrieved guide export differs')
        require(archive.extractfile('runs/collected/per_chromosome.tsv').read() ==
                (HERE/'inputs/per_chromosome.tsv').read_bytes(), 'Retrieved chromosome export differs')
        for chrom in chroms:
            prefix = 'runs/'+('pilot' if chrom == 'chr21' else 'scaled')+'/'+chrom
            summary = json.load(archive.extractfile(prefix+'/analysis/summary.json'))
            run = json.load(archive.extractfile(prefix+'/run.json'))
            require(run['status'] == 'COMPLETE' and run['chrom'] == chrom, 'Incomplete source run')
            require(sha_bytes(archive.extractfile(prefix+'/analysis/summary.json').read()) == run['analysis_sha256'],
                    'Source summary seal differs')
            raw = archive.extractfile(prefix+'/released_rna1.tsv').read()
            require(sha_bytes(raw) == summary['input_sha256']['released_rna1'], 'Raw released output seal differs')
            seen, eligible_seen, x_seen = set(), set(), set()
            lines = csv.DictReader(io.StringIO(raw.decode()), delimiter='\t')
            for data_row, r in enumerate(lines,1):
                raw_rows += 1
                gid, dna = r['#Id'], r['DNA']
                require(gid in spacers and r['Chromosome'] == chrom and r['Direction'] in ('+','-'),
                        'Raw identity differs')
                require(r['crRNA'] == spacers[gid]+'NNN' and len(dna) == 23, 'Raw alignment differs')
                start = int(r['Position'])
                pam = start+len(dna.replace('-','')) if r['Direction'] == '+' else start
                key = (gid, chrom, r['Direction'], pam)
                seen.add(key)
                if r['Bulge type'] == 'X':
                    x_seen.add(key)
                if gid not in targets or r['Bulge type'] != 'RNA' or int(r['Mismatches']) != 4:
                    continue
                if set(dna) - set('ACGT-') or dna[-2:] != 'GG':
                    unsupported_rows += 1
                    continue
                require(dna.count('-') == 1 and dna.index('-') < 19 and r['Bulge Size'] == '1',
                        'Released gap placement differs')
                require(sum(a != b for a,b in zip(spacers[gid],dna[:20]) if b != '-') == 4,
                        'Raw mismatch count differs')
                eligible_seen.add(key)
                # Decimal coordinate strings and tab separators define the frozen rank key.
                identity = '\t'.join(map(str,key))
                rank = int(hashlib.sha256(('aligned-control-v1\t'+identity).encode()).hexdigest(),16)
                record = dict(guide=gid, chrom=chrom, strand=r['Direction'], pam_coord=pam,
                              start=start, **{'class':'RNA'}, rna_position_1based=dna.index('-')+1,
                              mismatches=4, crRNA=r['crRNA'], DNA=dna,
                              source_member=prefix+'/released_rna1.tsv', source_data_row_1based=data_row,
                              selection_sha256=f'{rank:064x}')
                if key in kept[gid]:
                    if (dna, start) < (kept[gid][key]['DNA'], kept[gid][key]['start']):
                        kept[gid][key] = record
                    continue
                heap = heaps[gid]
                if len(heap) < targets[gid]:
                    heapq.heappush(heap,(-rank,key)); kept[gid][key] = record
                elif (-rank,key) > heap[0]:
                    _, removed = heapq.heapreplace(heap,(-rank,key))
                    del kept[gid][removed]
                    kept[gid][key] = record
            all_counts.update(k[0] for k in seen)
            ungapped_counts.update(k[0] for k in x_seen)
            eligible.update(k[0] for k in eligible_seen)
            require(len(seen) == summary['loci']['released_rna1'], 'Source locus count differs')
            chrom_audit.append(dict(chrom=chrom, released_loci=len(seen),
                                    raw_sha256=sha_bytes(raw), job_id=run['job_id']))
            print('control_source_verified:'+chrom, flush=True)
            peak_rss()
    pg = table(HERE/'inputs/per_guide.tsv')
    require(all(all_counts[r['guide']] == int(r['released_rna1']) and
                ungapped_counts[r['guide']] == int(r['released_rna1_ungapped']) for r in pg),
            'Raw released populations do not reconcile with per-guide deposit')
    selected = []
    for gid in sorted(targets, key=lambda g:int(g[1:])):
        require(len(kept[gid]) == targets[gid], 'Insufficient matched controls for '+gid)
        selected.extend(sorted(kept[gid].values(), key=lambda r:r['selection_sha256']))
    for i,r in enumerate(selected,1):
        r['data_row_1based'] = i
    return selected, dict(method='Lowest SHA-256 of aligned-control-v1<TAB>guide<TAB>chrom<TAB>strand<TAB>pam_coord; '
                                  'one locus per key, smallest DNA string on alignment ties; no score used.',
                          matching=['guide','mismatches=4','RNA single-bulge class','strict NGG PAM'],
                          raw_rows=raw_rows, released_loci=sum(all_counts.values()),
                          unsupported_eligible_alignment_rows_not_scored=unsupported_rows,
                          eligible_loci=sum(eligible.values()), eligible_by_guide=dict(eligible),
                          selected_loci=len(selected), eligible_loci_not_selected=sum(eligible.values())-len(selected),
                          scored_reverse_arm_rows=0, by_chromosome=chrom_audit)


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def pin_encoder(arms):
    original = load_module('deposited_encoder', HERE/'vendor/Encoder_sgRNA_off.py')
    by_arm = {}
    for arm, rows in arms.items():
        stats = Counter()
        for start in range(0,len(rows),512):
            block = rows[start:start+512]
            pairs = []
            for r in block:
                before = (r['crRNA'],r['DNA'])
                after = tuple(s.replace('-','_') for s in before)
                require(all(len(a) == len(b) == 23 and all(x == y if x != '-' else y == '_'
                            for x,y in zip(a,b)) for a,b in zip(before,after)), 'Non-gap conversion changed')
                stats['alignment_gap_characters_converted'] += sum(s.count('-') for s in before)
                stats['non_gap_characters_unchanged'] += sum(len(s)-s.count('-') for s in before)
                pairs.append(after)
            enc = encode(pairs)
            for j,pair in enumerate(pairs):
                want = original.Encoder(*pair).on_off_code
                require(want.shape == (24,7) and np.array_equal(enc[j],want), 'Original encoder disagreement')
                require(not np.any(want[0]), 'Left padding column is not zero')
                stats['original_encoder_tensor_equal_rows'] += 1
        by_arm[arm] = dict(stats, rows=len(rows))
    # A real gap and zero padding must remain distinguishable.
    toy = ('A'*20+'NNN','A'*19+'-'+'TGG')
    raw = original.Encoder(*toy).on_off_code
    converted = original.Encoder(toy[0],toy[1].replace('-','_')).on_off_code
    require(not np.array_equal(raw,converted) and converted[20,4] == converted[20,5] == 1 and
            raw[20,4] == raw[20,5] == 0, 'Gap-vs-padding control failed')
    return dict(published_gap='_', published_padding='-', input_length=23, tensor_columns=24,
                tensor_channels=7, left_padding_columns=1,
                conversion="Within deposited alignment strings only, '-' -> '_'; all other characters fixed. "
                           "The encoder subsequently prepends its own '-' padding to both sequences.",
                guide_N_policy='Each guide N is replaced by the DNA symbol at that same column; guide PAM is NNN.',
                all_alignments_representable=True, gap_vs_padding_negative_control_pass=True, by_arm=by_arm)


def summary(treatment, control, active, score_t, score_c, score_a, pg, audit):
    def dist(rows, scores, predicate):
        return distribution([p for r,p in zip(rows,scores) if predicate(r)])
    per_guide = {}
    for g in sorted((r['guide'] for r in pg),key=lambda x:int(x[1:])):
        per_guide[g] = {'treatment':dist(treatment,score_t,lambda r:r['guide']==g),
                        'control':dist(control,score_c,lambda r:r['guide']==g)}
    split = {}
    for bits in ('01','11','00','10'):
        split[bits] = dist(active,score_a,lambda r:r['released_returned']+r['corrected_returned']==bits)
    missed_all = dist(active,score_a,lambda r:r['released_returned']=='0')
    control_mm4 = dist(active,score_a,lambda r:r['released_returned']=='1' and int(r['mismatches'])==4)
    shared = sorted({r['guide'] for r in active if r['released_returned']=='0' and r['corrected_returned']=='1'} &
                    {r['guide'] for r in active if r['released_returned']=='1' and int(r['mismatches'])==4},
                    key=lambda x:int(x[1:]))
    active_matched = {}
    selected_t, selected_c = [], []
    for g in shared:
        ti = [i for i,r in enumerate(active) if r['guide']==g and r['released_returned']=='0' and r['corrected_returned']=='1']
        ci = [i for i,r in enumerate(active) if r['guide']==g and r['released_returned']=='1' and int(r['mismatches'])==4]
        n = min(len(ti),len(ci))
        selected_t.extend(ti[:n]); selected_c.extend(ci[:n])
        active_matched[g] = dict(treatment_indices=[int(active[i]['row_index']) for i in ti[:n]],
                                control_indices=[int(active[i]['row_index']) for i in ci[:n]],
                                treatment=distribution(score_a[ti[:n]]), control=distribution(score_a[ci[:n]]))
    t, c = distribution(score_t), distribution(score_c)
    wins = ties = loss = 0
    for g in per_guide:
        a, b = per_guide[g]['treatment'],per_guide[g]['control']
        if a['rows']:
            wins += a['above'] > b['above']; ties += a['above'] == b['above']; loss += a['above'] < b['above']
    return dict(cutoff=dict(exact='3/1000', comparison='strictly greater',
                            provenance='Lin et al. 2020 sensitivity analysis; carried over unchanged.',
                            calibration_limit='Published CIRCLE-seq leave-one-guide-out cutoff; calibration for this combined-data checkpoint is not established.'),
                treatment=t, control=c, by_guide=per_guide,
                comparison=dict(above_percent_difference=t['above_percent']-c['above_percent'],
                                above_count_difference=t['above']-c['above'],
                                guide_above_rate_higher=int(wins), guide_above_rate_equal=int(ties),
                                guide_above_rate_lower=int(loss)),
                reverse_arm=distribution([]), control_selection=audit,
                active=dict(by_return_bits=split, released_missed_all_including_non_NGG=missed_all,
                            released_returned_mm4=control_mm4,
                            exact_guide_mm4_matched=dict(shared_guides=shared,by_guide=active_matched,
                                treatment=distribution(score_a[selected_t]),control=distribution(score_a[selected_c]),
                                policy='For each shared guide at four mismatches, first min(n_missed,n_returned) rows by original assay row_index in each arm.'),
                            all=distribution(score_a)),
                panel=dict(guides=len(pg), treatment_guides=sum(bool(per_guide[g]['treatment']['rows']) for g in per_guide),
                           released_loci=sum(int(r['released_rna1']) for r in pg),
                           corrected_loci=sum(int(r['corrected_rna1']) for r in pg),
                           ungapped_released=sum(int(r['released_rna0']) for r in pg),
                           ungapped_corrected=sum(int(r['corrected_rna0']) for r in pg),
                           deposited_magnitude_baseline=sum(int(r['genomewide_baseline']) for r in table(HERE/'inputs/ungapped_baseline_total.tsv'))))


def run():
    begin = time.monotonic()
    check_pins()
    out = HERE/'results'
    events = []
    def event(name):
        events.append(dict(event=name,utc=datetime.now(timezone.utc).isoformat(),peak_rss_bytes=peak_rss()))
        dump(out/'execution.json',events)
        print(name,flush=True)
    event('input_seals_verified')
    treatment,pg,counts,spacers = prepare_treatment()
    active = prepare_active(spacers,treatment)
    control,audit = control_sample(counts,spacers)
    event('score_independent_controls_frozen')
    convention = pin_encoder(dict(treatment=treatment,control=control,active=active))
    dump(out/'encoder_convention.json',convention)
    event('encoder_convention_pinned_all_rows_before_any_score')
    prior = json.loads((HERE/'vendor/prior_validation.json').read_text())
    require(prior['gate']['pass_for_scoring'], 'Previously validated checkpoint gate failed')
    event('deposited_validation_gate_confirmed')
    model = NumpyModel(HERE/'inputs/CRISPR_Net_CIRCLE_elevation_SITE_weights.h5')
    predictions = []
    for name,rows in (('treatment',treatment),('control',control),('active',active)):
        pairs = [(r['crRNA'].replace('-','_'),r['DNA'].replace('-','_')) for r in rows]
        scores = batches(model,pairs)
        predictions.append(scores)
        for r,p in zip(rows,scores):
            r.update(score=format(float(p),'.17g'), above_3_over_1000=int(cut(p)),
                     encoder_crRNA=r['crRNA'].replace('-','_'),encoder_DNA=r['DNA'].replace('-','_'))
        fields=list(rows[0])
        write_table(out/(name+'_scores.tsv'),fields,rows)
        event(name+'_scored')
    score_t,score_c,score_a = predictions
    # The assay endpoint is already present in the treatment export; compare exact inputs and scores.
    t_by_key = {(r['guide'],r['chrom'],r['strand'],int(r['pam_coord'])):float(r['score']) for r in treatment}
    for r in active:
        if r['released_returned'] == '0' and r['corrected_returned'] == '1':
            require(abs(float(r['score'])-t_by_key[r['guide'],r['chrom'],r['strand'],int(r['pam_coord'])]) < 2e-12,
                    'Overlapping assay/treatment scores differ')
    result = summary(treatment,control,active,score_t,score_c,score_a,pg,audit)
    dump(out/'summary.json',result)
    event('all_score_summaries_written')
    dump(out/'resources.json',dict(wall_seconds=time.monotonic()-begin,peak_rss_bytes=peak_rss(),threads=4,
                                    executable=sys.executable,python=sys.version))
    dump(out/'outputs.sha256.json',{p.name:sha(p) for p in sorted(out.iterdir())
                                   if p.is_file() and p.name not in
                                   ('outputs.sha256.json','verification.json','per_guide_scores.tsv')})
    print(json.dumps({k:result[k] for k in ('treatment','control','comparison')},sort_keys=True),flush=True)


if __name__ == '__main__':
    run()
