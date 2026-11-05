#!/usr/bin/env python3
"""Independent reconstruction of controls, encodings, scores and all reported summaries."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from integrity import check_pins, require, sha, peak_rss
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import sqlite3
import statistics
import subprocess
import tarfile
import time
from collections import Counter, defaultdict
from fractions import Fraction
import numpy as np


def module(name, path):
    spec = importlib.util.spec_from_file_location(name,path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def read_rows(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def compare(got, want, path='summary'):
    if isinstance(want,dict):
        require(isinstance(got,dict) and set(got)==set(want),path+': keys differ')
        for key in want:
            compare(got[key],want[key],path+'.'+key)
    elif isinstance(want,list):
        require(len(got)==len(want),path+': lengths differ')
        for i,(a,b) in enumerate(zip(got,want)):
            compare(a,b,path+f'[{i}]')
    elif isinstance(want,float):
        require(math.isclose(float(got),want,rel_tol=2e-12,abs_tol=1e-300),path+': numerical value differs')
    else:
        require(got==want,path+': value differs')


def decision(p):
    return Fraction.from_float(float(p)) > Fraction(3,1000)


def ref_distribution(values):
    ordered = sorted(map(float,values))
    n = len(ordered)
    if not n:
        return {'rows':0,'above':0}
    require(all(math.isfinite(x) and 0<=x<=1 for x in ordered),'Invalid score domain')
    passed = sum(decision(x) for x in ordered)
    quantiles = {}
    for q in (0,.01,.05,.25,.5,.75,.95,.99,1):
        loc=(n-1)*q
        a,b=math.floor(loc),math.ceil(loc)
        quantiles[str(q)] = ordered[a]+(ordered[b]-ordered[a])*(loc-a)
    boundaries=['0','0.000001','0.0001','0.001','0.003','0.01','0.1','0.5','1']
    exact=[Fraction(x) for x in boundaries]
    histogram=[]
    rational=[Fraction.from_float(x) for x in ordered]
    for i in range(len(exact)-1):
        histogram.append(sum(exact[i]<=x and (x<=exact[i+1] if i==len(exact)-2 else x<exact[i+1])
                             for x in rational))
    return dict(rows=n,above=passed,above_fraction=str(Fraction(passed,n)),above_percent=100*passed/n,
                mean=statistics.mean(ordered),population_sd=statistics.pstdev(ordered),
                quantiles=quantiles,histogram_edges=boundaries,histogram_counts=histogram,
                nearest_cutoff_distance=float(min(abs(x-Fraction(3,1000)) for x in rational)))


def ref_source_and_selection(root, treatment, pg):
    work=HERE/'verify_work'
    work.mkdir(exist_ok=True)
    database=work/'controls.sqlite'
    database.unlink(missing_ok=True)
    db=sqlite3.connect(database)
    db.execute('PRAGMA journal_mode=OFF')
    db.execute('PRAGMA cache_size=-16000')
    db.execute('CREATE TABLE loci (guide TEXT, chrom TEXT, strand TEXT, pam INTEGER, '
               'PRIMARY KEY(guide,chrom,strand,pam)) WITHOUT ROWID')
    db.execute('CREATE TABLE ungapped (guide TEXT, chrom TEXT, strand TEXT, pam INTEGER, '
               'PRIMARY KEY(guide,chrom,strand,pam)) WITHOUT ROWID')
    db.execute('CREATE TABLE eligible (guide TEXT, chrom TEXT, strand TEXT, pam INTEGER, '
               'dna TEXT, start INTEGER, member TEXT, rownum INTEGER, rank TEXT, '
               'PRIMARY KEY(guide,chrom,strand,pam)) WITHOUT ROWID')
    desired=Counter(r['guide'] for r in treatment)
    spacers={r['guide']:r['spacer'] for r in pg}
    raw_count=unsupported=0
    audit=[]
    with tarfile.open(root/'inputs/released_run.tar.gz') as t:
        for name in ('per_guide.tsv','per_chromosome.tsv'):
            require(t.extractfile('runs/collected/'+name).read()==(root/'inputs'/name).read_bytes(),
                    'Source/deposit export mismatch: '+name)
        for chrom in [f'chr{i}' for i in range(1,23)]+['chrX','chrY','chrM']:
            base='runs/'+('pilot' if chrom=='chr21' else 'scaled')+'/'+chrom
            summary_bytes=t.extractfile(base+'/analysis/summary.json').read()
            source_summary=json.loads(summary_bytes)
            run=json.load(t.extractfile(base+'/run.json'))
            require(run['chrom']==chrom and run['status']=='COMPLETE','Source run identity/state differs')
            require(hashlib.sha256(summary_bytes).hexdigest()==run['analysis_sha256'],'Source summary seal differs')
            raw=t.extractfile(base+'/released_rna1.tsv').read()
            digest=hashlib.sha256(raw).hexdigest()
            require(digest==source_summary['input_sha256']['released_rna1'],'Raw output source digest differs')
            # A distinct implementation uses SQL deduplication and sorting, not a bounded heap.
            for rownum,line in enumerate(raw.decode().splitlines()[1:],1):
                raw_count+=1
                fields=line.split('\t')
                require(len(fields)==9,'Raw field count differs')
                g,kind,on,off,ch,pos,strand,mm,bulge=fields
                require(g in spacers and ch==chrom and strand in ['+','-'],'Raw identity differs')
                require(on==spacers[g]+'NNN' and len(off)==23,'Raw guide differs')
                start=int(pos)
                pam=start+sum(x!='-' for x in off) if strand=='+' else start
                key=(g,ch,strand,pam)
                db.execute('INSERT OR IGNORE INTO loci VALUES (?,?,?,?)',key)
                if kind=='X':
                    db.execute('INSERT OR IGNORE INTO ungapped VALUES (?,?,?,?)',key)
                if g not in desired or kind!='RNA' or int(mm)!=4:
                    continue
                if any(x not in 'ACGT-' for x in off) or off[-2:]!='GG':
                    unsupported+=1
                    continue
                gaps=[i for i,x in enumerate(off) if x=='-']
                require(len(gaps)==1 and gaps[0]<19 and bulge=='1','Unsupported released gap')
                require(sum(a!=b for a,b in zip(spacers[g],off[:20]) if b!='-')==4,'Raw mismatch annotation differs')
                identity='aligned-control-v1\t'+g+'\t'+ch+'\t'+strand+'\t'+str(pam)
                rank=hashlib.sha256(identity.encode()).hexdigest()
                db.execute('INSERT INTO eligible VALUES (?,?,?,?,?,?,?,?,?) '
                           'ON CONFLICT(guide,chrom,strand,pam) DO UPDATE SET '
                           'dna=excluded.dna,start=excluded.start,member=excluded.member,rownum=excluded.rownum '
                           'WHERE excluded.dna<eligible.dna OR (excluded.dna=eligible.dna AND excluded.start<eligible.start)',
                           (*key,off,start,base+'/released_rna1.tsv',rownum,rank))
            db.commit()
            count=db.execute('SELECT count(*) FROM loci WHERE chrom=?',(chrom,)).fetchone()[0]
            require(count==source_summary['loci']['released_rna1'],'Chromosome released count differs')
            audit.append(dict(chrom=chrom,released_loci=count,raw_sha256=digest,job_id=run['job_id']))
            print('independent_source_replay:'+chrom,flush=True)
            peak_rss()
    allcounts=dict(db.execute('SELECT guide,count(*) FROM loci GROUP BY guide'))
    ungapped=dict(db.execute('SELECT guide,count(*) FROM ungapped GROUP BY guide'))
    eligible=dict(db.execute('SELECT guide,count(*) FROM eligible GROUP BY guide'))
    for r in pg:
        require(allcounts.get(r['guide'],0)==int(r['released_rna1']) and
                ungapped.get(r['guide'],0)==int(r['released_rna1_ungapped']),'Guide population differs')
    selected=[]
    for g in sorted(desired,key=lambda g:int(g[1:])):
        rows=list(db.execute('SELECT * FROM eligible WHERE guide=? ORDER BY rank,chrom,strand,pam LIMIT ?',
                             (g,desired[g])))
        require(len(rows)==desired[g],'Incomplete guide matching')
        for gid,ch,strand,pam,dna,start,member,rownum,rank in rows:
            selected.append(dict(guide=gid,chrom=ch,strand=strand,pam_coord=str(pam),DNA=dna,start=str(start),
                                 source_member=member,source_data_row_1based=str(rownum),selection_sha256=rank,
                                 crRNA=spacers[gid]+'NNN',mismatches='4',**{'class':'RNA'},
                                 rna_position_1based=str(dna.index('-')+1),data_row_1based=str(len(selected)+1)))
    output_audit=dict(method='Lowest SHA-256 of aligned-control-v1<TAB>guide<TAB>chrom<TAB>strand<TAB>pam_coord; '
                            'one locus per key, smallest DNA string on alignment ties; no score used.',
                     matching=['guide','mismatches=4','RNA single-bulge class','strict NGG PAM'],
                     raw_rows=raw_count,released_loci=sum(allcounts.values()),
                     unsupported_eligible_alignment_rows_not_scored=unsupported,
                     eligible_loci=sum(eligible.values()),eligible_by_guide=eligible,
                     selected_loci=len(selected),eligible_loci_not_selected=sum(eligible.values())-len(selected),
                     scored_reverse_arm_rows=0,by_chromosome=audit)
    return selected,output_audit,db


def verify(root,integrity_only=False,scores_only=False):
    begin=time.monotonic()
    check_pins(root)
    expected_hashes=read_rows(root/'INPUT_SHA256.tsv') if (root/'INPUT_SHA256.tsv').exists() else []
    if expected_hashes:
        compare({r['path']:r['sha256'] for r in expected_hashes},json.loads((root/'pins.json').read_text()),'input_hash_table')
    print('PASS input SHA-256 seals',flush=True)
    if integrity_only:
        return
    out=HERE/'results'
    hashes=json.loads((out/'outputs.sha256.json').read_text())
    for name,expected in hashes.items():
        require(sha(out/name)==expected,'Produced output digest differs: '+name)
    original=module('original_encoder_independent',root/'vendor/Encoder_sgRNA_off.py')
    reference=module('validated_independent_forward',root/'vendor/verify.py')
    forward=reference.ReferenceModel(root/'inputs/CRISPR_Net_CIRCLE_elevation_SITE_weights.h5')
    pg=read_rows(root/'inputs/per_guide.tsv')
    spacers={r['guide']:r['spacer'] for r in pg}
    raw_t=read_rows(root/'inputs/corrected_only_alignments.part000.tsv')
    reverse=read_rows(root/'inputs/released_only_alignments.part000.tsv')
    require(not reverse,'Nonempty reverse arm')
    require(len({(r['guide'],r['chrom'],r['strand'],r['pam_coord']) for r in raw_t})==len(raw_t),'Duplicate treatment locus')
    counts=Counter(r['guide'] for r in raw_t)
    for r in pg:
        require(counts[r['guide']]==int(r['corrected_only']) and int(r['released_only'])==0,'Treatment stratum differs')
    selected,audit,db=ref_source_and_selection(root,raw_t,pg)
    convention=json.loads((out/'encoder_convention.json').read_text())
    compare({k:convention[k] for k in ('published_gap','published_padding','input_length','tensor_columns',
                                      'tensor_channels','left_padding_columns','all_alignments_representable',
                                      'gap_vs_padding_negative_control_pass')},
            dict(published_gap='_',published_padding='-',input_length=23,tensor_columns=24,tensor_channels=7,
                 left_padding_columns=1,all_alignments_representable=True,gap_vs_padding_negative_control_pass=True),
            'encoder_constants')
    # Directly exhibit why passing '-' through unchanged is not a valid encoding.
    rawtoy=original.Encoder('A'*20+'NNN','A'*19+'-'+'TGG').on_off_code
    gappedtoy=original.Encoder('A'*20+'NNN','A'*19+'_'+'TGG').on_off_code
    require(not np.array_equal(rawtoy,gappedtoy) and gappedtoy[20,4]==gappedtoy[20,5]==1 and
            rawtoy[20,4]==rawtoy[20,5]==0,'Gap/padding falsification failed')
    scored={name:read_rows(out/(name+'_scores.tsv')) for name in ('treatment','control','active')}
    for i,(r,s) in enumerate(zip(raw_t,scored['treatment']),1):
        require(all(r[k]==s[k] for k in r) and s['data_row_1based']==str(i),'Treatment output identity differs')
        require(r['crRNA']==spacers[r['guide']]+'NNN' and r['DNA'][19]=='-' and
                r['DNA'].count('-')==1 and r['DNA'][-2:]=='GG','Treatment alignment class differs')
        require(sum(a!=b for a,b in zip(r['crRNA'][:20],r['DNA'][:20]) if b!='-')==4 and
                r['mismatches']=='4' and r['rna_position_1based']=='20' and r['class']=='RNA',
                'Treatment mismatch/class annotation differs')
        require(int(r['pam_coord'])==int(r['start'])+(22 if r['strand']=='+' else 0),'Treatment coordinates differ')
    require(len(raw_t)==len(scored['treatment']) and len(selected)==len(scored['control']),'Arm length differs')
    for r,s in zip(selected,scored['control']):
        require(all(r[k]==s[k] for k in r),'Independently selected control differs')
        require(not db.execute('SELECT 1 FROM eligible WHERE guide=? AND chrom=? AND strand=? AND pam=?',
                              (r['guide'],r['chrom'],r['strand'],int(r['pam_coord']))).fetchone() is None,
                'Control is not in released output')
    tkeys={(r['guide'],r['chrom'],r['strand'],r['pam_coord']) for r in raw_t}
    ckeys={(r['guide'],r['chrom'],r['strand'],r['pam_coord']) for r in selected}
    require(tkeys.isdisjoint(ckeys),'Treatment and control loci overlap')
    raw_a=read_rows(root/'inputs/active_sites.part000.tsv')
    assay=json.loads((root/'inputs/sites129.json').read_text())
    require({int(r['row_index']) for r in raw_a}==set(range(len(assay))) and len(raw_a)==len(assay),
            'Assay identities incomplete')
    require(len(scored['active'])==len(raw_a),'Assay score rows incomplete')
    bits=Counter()
    for row,r in zip(sorted(raw_a,key=lambda x:int(x['row_index'])),scored['active']):
        i=int(row['row_index']); s=assay[i]
        require(all(r[k]==row[k] for k in row),'Assay return columns changed')
        require(all(str(s[k])==r[k] for k in ('chrom','strand','start')) and
                s['spacer']==spacers[r['guide']] and float(r['reads'])==float(s['reads'])>0,
                'Assay join differs')
        require(r['DNA']==s['aligned_target'] and r['crRNA']==s['spacer']+'NNN','Assay scoring input differs')
        require(s['target']==s['aligned_target'].replace('-','') and s['aligned_target'][19]=='-',
                'Assay source alignment differs')
        cost=[sum(a!=b for a,b in zip(s['spacer'][:j]+s['spacer'][j+1:],s['aligned_target'][:19]))
              for j in range(20)]
        pam=int(s['aligned_target'][-2:]=='GG')
        require((int(r['released_returned']),int(r['corrected_returned']))==
                (int(pam and min(cost[:-1])<=4),int(pam and min(cost)<=4)),'Assay return-bit replay differs')
        pcost=sum(a!=b for a,b in zip(s['sgRNA'][-2:],r['DNA'][-2:]))
        require(int(r['mismatches'])==cost[-1] and int(r['PAM_substitutions'])==pcost and
                int(r['deposited_mm_including_PAM'])==cost[-1]+pcost==int(s['mm']),
                'Assay protospacer/PAM mismatch split differs')
        require(int(r['pam_coord'])==int(r['start'])+(22 if r['strand']=='+' else 0),'Assay PAM boundary differs')
        bits[r['released_returned']+r['corrected_returned']]+=1
        if r['released_returned']+r['corrected_returned']=='01':
            require((r['guide'],r['chrom'],r['strand'],r['pam_coord']) in tkeys,'Missed active absent from treatment')
    arrays={}
    replay_max=0.
    for name,rows in scored.items():
        stats=Counter()
        predictions=[]
        for start in range(0,len(rows),512):
            block=rows[start:start+512]
            pairs=[]
            for r in block:
                a,b=r['crRNA'],r['DNA']
                require(len(a)==len(b)==23,'Wrong alignment width')
                converted=(a.translate(str.maketrans({'-':'_'})),b.translate(str.maketrans({'-':'_'})))
                require(r['encoder_crRNA']==converted[0] and r['encoder_DNA']==converted[1],
                        'Gap conversion output differs')
                for old,new in zip((a,b),converted):
                    require(len(old)==len(new) and all((x=='-' and y=='_') or (x!='-' and x==y)
                                                      for x,y in zip(old,new)),'Non-gap column changed')
                    stats['alignment_gap_characters_converted']+=old.count('-')
                    stats['non_gap_characters_unchanged']+=len(old)-old.count('-')
                pairs.append(converted)
            encoded=reference.encode_reference(pairs)
            for j,pair in enumerate(pairs):
                original_tensor=original.Encoder(*pair).on_off_code
                require(original_tensor.shape==(24,7) and np.array_equal(encoded[j],original_tensor),
                        'Independent/original encoder differs')
                require(not np.any(original_tensor[0]),'Nonzero left-padding column')
                stats['original_encoder_tensor_equal_rows']+=1
            got=forward.predict(encoded)
            want=np.array([float(r['score']) for r in block])
            replay_max=max(replay_max,float(np.max(np.abs(got-want))))
            require(np.allclose(got,want,rtol=2e-11,atol=2e-12),'Independent forward score differs')
            require([decision(x) for x in got]==[bool(int(r['above_3_over_1000'])) for r in block] and
                    [decision(x) for x in want]==[decision(x) for x in got],'Exact threshold decision differs')
            # Stored score decimals are independently replayed; summaries refer to these exact decimals.
            predictions.extend(want)
        arrays[name]=np.array(predictions)
        compare(dict(stats,rows=len(rows)),convention['by_arm'][name],'encoder_audit.'+name)
        print('PASS independent tensor and forward replay: '+name,flush=True)
        peak_rss()
    events=json.loads((out/'execution.json').read_text())
    labels=[x['event'] for x in events]
    require(labels.index('encoder_convention_pinned_all_rows_before_any_score') <
            min(labels.index(n+'_scored') for n in scored),'Encoding was not pinned before scoring')
    require(json.loads((root/'vendor/prior_validation.json').read_text())['gate']['pass_for_scoring'],
            'Prior checkpoint gate failed')
    result=json.loads((out/'summary.json').read_text())
    ref={}
    for name in ('treatment','control'):
        ref[name]=ref_distribution(arrays[name])
        compare(result[name],ref[name],name)
    perguide={}
    for row in pg:
        g=row['guide']; perguide[g]={}
        for name in ('treatment','control'):
            perguide[g][name]=ref_distribution([p for r,p in zip(scored[name],arrays[name]) if r['guide']==g])
        require(perguide[g]['treatment']['rows']==perguide[g]['control']['rows'],'Guide arms are not matched')
    compare(result['by_guide'],perguide,'by_guide')
    compare(result['control_selection'],audit,'control_selection')
    compare(result['reverse_arm'],ref_distribution([]),'reverse_arm')
    outcomes=Counter()
    for arms in perguide.values():
        a,b=arms['treatment'],arms['control']
        if a['rows']:
            outcomes[(a['above']>b['above'])-(a['above']<b['above'])]+=1
    contrast=dict(above_percent_difference=ref['treatment']['above_percent']-ref['control']['above_percent'],
                  above_count_difference=ref['treatment']['above']-ref['control']['above'],
                  guide_above_rate_higher=outcomes[1],guide_above_rate_equal=outcomes[0],guide_above_rate_lower=outcomes[-1])
    compare(result['comparison'],contrast,'comparison')
    arows=scored['active']; ap=arrays['active']
    def group(predicate):
        return ref_distribution([p for r,p in zip(arows,ap) if predicate(r)])
    split={b:group(lambda r:r['released_returned']+r['corrected_returned']==b) for b in ('01','11','00','10')}
    for b in split:
        require(split[b]['rows']==bits[b],'Active bit count differs')
    active_result=dict(by_return_bits=split,
                       released_missed_all_including_non_NGG=group(lambda r:r['released_returned']=='0'),
                       released_returned_mm4=group(lambda r:r['released_returned']=='1' and int(r['mismatches'])==4),
                       all=ref_distribution(ap))
    shared=sorted({r['guide'] for r in arows if r['released_returned']+r['corrected_returned']=='01'} &
                  {r['guide'] for r in arows if r['released_returned']=='1' and int(r['mismatches'])==4},
                  key=lambda x:int(x[1:]))
    per_active={}; all_t=[]; all_c=[]
    for g in shared:
        ts=[i for i,r in enumerate(arows) if r['guide']==g and r['released_returned']+r['corrected_returned']=='01']
        cs=[i for i,r in enumerate(arows) if r['guide']==g and r['released_returned']=='1' and int(r['mismatches'])==4]
        k=min(len(ts),len(cs)); ts,cs=ts[:k],cs[:k]
        all_t+=ts; all_c+=cs
        per_active[g]=dict(treatment_indices=[int(arows[i]['row_index']) for i in ts],
                           control_indices=[int(arows[i]['row_index']) for i in cs],
                           treatment=ref_distribution(ap[ts]),control=ref_distribution(ap[cs]))
    active_result['exact_guide_mm4_matched']=dict(shared_guides=shared,by_guide=per_active,
                   treatment=ref_distribution(ap[all_t]),control=ref_distribution(ap[all_c]),
                   policy='For each shared guide at four mismatches, first min(n_missed,n_returned) rows by original assay row_index in each arm.')
    compare(result['active'],active_result,'active')
    baseline=read_rows(root/'inputs/ungapped_baseline_total.tsv')
    for r in pg:
        b=next(x for x in baseline if x['guide']==r['guide'])
        values=[int(r[k]) for k in ('released_rna0','corrected_rna0','released_rna1_ungapped','corrected_rna1_ungapped')]
        require(len(set(values))==1 and all(r[k]==b[k] for k in ('released_rna0','corrected_rna0',
                                            'released_rna1_ungapped','corrected_rna1_ungapped')),
                'Ungapped build control differs')
    panel=dict(guides=len(pg),treatment_guides=sum(x['treatment']['rows']>0 for x in perguide.values()),
               released_loci=sum(int(r['released_rna1']) for r in pg),
               corrected_loci=sum(int(r['corrected_rna1']) for r in pg),
               ungapped_released=sum(int(r['released_rna0']) for r in pg),
               ungapped_corrected=sum(int(r['corrected_rna0']) for r in pg),
               deposited_magnitude_baseline=sum(int(r['genomewide_baseline']) for r in baseline))
    compare(result['panel'],panel,'panel')
    require(panel['corrected_loci']-panel['released_loci']==len(raw_t),'Build totals differ from scored arm')
    require(result['cutoff']['exact']=='3/1000' and result['cutoff']['comparison']=='strictly greater','Cutoff changed')
    # Ensure every score in the overlapping assay/treatment rows agrees at its own locus.
    tl={(r['guide'],r['chrom'],r['strand'],r['pam_coord']):(r,float(r['score'])) for r in scored['treatment']}
    overlap=0
    for r in arows:
        if r['released_returned']+r['corrected_returned']=='01':
            tr,p=tl[r['guide'],r['chrom'],r['strand'],r['pam_coord']]
            require(tr['DNA']==r['DNA'] and abs(p-float(r['score']))<2e-12,'Assay/treatment overlap differs')
            overlap+=1
    db.close()
    if not scores_only:
        # Formatting is shared; every supplied numerical value was independently reconstructed above.
        renderer=module('report_formatter',HERE/'render_findings.py')
        layout=json.loads((HERE/'layout/layout_receipt.json').read_text())
        layout_dir=HERE/'layout'
        require(sha(layout_dir/'body_format_source.tex')==layout['format_source_sha256'] and
                sha(layout_dir/'sentence_layout.tex')==layout['tex_sha256'] and
                sha(layout_dir/'sentence_layout.pdf')==layout['pdf_sha256'],'Layout artifact digest differs')
        # Check line count from rendered PDF positions independently of the typesetter's paragraph counter.
        import xml.etree.ElementTree as ET
        bbox=subprocess.run(['pdftotext','-bbox-layout',str(layout_dir/'sentence_layout.pdf'),'-'],
                            capture_output=True,text=True,check=True,timeout=30)
        document=ET.fromstring(bbox.stdout)
        body_lines=[]
        for element in document.iter():
            if element.tag.rsplit('}',1)[-1]=='line' and float(element.attrib['yMin'])<650:
                words=[''.join(w.itertext()) for w in element if w.tag.rsplit('}',1)[-1]=='word']
                if words:
                    body_lines.append(' '.join(words))
        start=next(i for i,line in enumerate(body_lines) if line.startswith('At the published'))
        removal=next(i for i,line in enumerate(body_lines) if line.startswith('The tables and one command'))
        require(removal-start==layout['proposed_lines'] and len(body_lines)-removal==layout['removed_lines'],
                'Independent PDF line count differs')
        require(layout['proposed_lines']<=4 and layout['removed_lines']>=layout['proposed_lines'] and
                layout['build_exit_code']==0,'Printed-line budget fails')
        log=(layout_dir/'sentence_layout.log').read_text()
        require('Overfull \\hbox' not in log,'Overfull layout line')
        sentence=renderer.proposed_sentence(result)
        require(sentence in (layout_dir/'sentence_layout.tex').read_text(),'Sentence differs from layout proof')
        compare((out/'per_guide_scores.tsv').read_text(),renderer.guide_table(result),'paired_guide_table')
        compare((HERE.parent/'FINDINGS.md').read_text(),renderer.render(result,convention,layout),'findings')
        corruption=json.loads((HERE/'corruption_control/receipt.json').read_text())
        original=(root/corruption['input']).read_bytes()
        changed=(HERE/'corruption_control/changed_input'/corruption['input']).read_bytes()
        offsets=[i for i,(a,b) in enumerate(zip(original,changed)) if a!=b]
        require(len(original)==len(changed) and len(offsets)==corruption['changed_bytes']==1 and
                offsets[0]==corruption['byte_offset_0based'],'Corruption byte audit differs')
        require(hashlib.sha256(original).hexdigest()==corruption['original_sha256'] and
                hashlib.sha256(changed).hexdigest()==corruption['changed_sha256'],'Corruption digest differs')
        require(corruption['unmodified_verifier_exit_code']==0 and corruption['corrupted_verifier_exit_code']==1 and
                corruption['rejected_before_scoring'] and corruption['status']=='PASS','Corruption gate failed')
        # Re-execute the rejection, rather than trusting the receipt's exit-code assertion.
        cp=subprocess.run(['/usr/local/bin/python3','-I',str(HERE/'verify.py'),'--integrity-only',
                           '--input-root',str(HERE/'corruption_control/changed_input')],
                          capture_output=True,text=True,timeout=60)
        require(cp.returncode==1 and 'SHA-256 mismatch: '+corruption['input'] in cp.stdout+cp.stderr,
                'Replayed corruption rejection failed')
        print('PASS rendered report, guide table, independent PDF line counts and one-byte corruption replay',flush=True)
    receipt=dict(status='PASS',exit_code=0,independently_replayed_score_rows=sum(len(x) for x in scored.values()),
                  scored_rows_by_arm={k:len(v) for k,v in scored.items()},
                  maximum_absolute_score_difference=replay_max,
                  missed_active_treatment_overlap_rows=overlap,
                  peak_rss_bytes=peak_rss(),wall_seconds=time.monotonic()-begin,
                  findings_numeric_summary_reproduced=not scores_only)
    (out/'verification.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt,sort_keys=True),flush=True)
    print('VERIFIER_EXIT_CODE=0',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root',type=Path,default=HERE)
    parser.add_argument('--integrity-only',action='store_true')
    parser.add_argument('--scores-only',action='store_true',help='Bootstrap replay before rendering the report')
    args=parser.parse_args()
    try:
        verify(args.input_root.resolve(),args.integrity_only,args.scores_only)
    except Exception:
        print('VERIFIER_EXIT_CODE=1',flush=True)
        raise
