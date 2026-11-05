#!/usr/bin/env python3
"""Validate the deposited checkpoint before scoring the literal target TSV pairs."""
import os
import sys
sys.dont_write_bytecode = True
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import argparse
import csv
import hashlib
import json
import time
from datetime import datetime,timezone
from model_numpy import NumpyModel,batches,np
from analysis import metrics,distribution,validation_gate,above

HERE = Path(__file__).resolve().parent
WEIGHT = 'CRISPR_Net_CIRCLE_elevation_SITE_weights.h5'
DATASETS = {
 'CIRCLE_seq_10gRNA_wholeDataset.csv': ('sgRNA_seq','off_seq',False),
 'elevation_6gRNA_wholeDataset.csv': ('crRNA','DNA',False),
 'Kleinstiver_5gRNA_wholeDataset.csv': ('sgRNA_seq','off_seq',True),
 'Listgarten_22gRNA_wholeDataset.csv': ('sgRNA_seq','off_seq',True),
}
EXAMPLES = [
 ('GAGT_CCGAGCAGAAGAAGAATGG','GAGTACCAAGTAGAAGAAAAATTT',0.0),
 ('GTTGCCCCACAGGGCAGTAAAGG','GTGGACACCCCGGGCAGGAAAGG',1.490116e-7),
 ('GGGTGGGGGGAGTTTGCTCCCGG','GTGTGGGGTAAATTTGCTCCCAG',.6078471),
 ('GGGTGGGGGGAGTTTGCTCCAGG','AGGTGGGGTGA_TTTGCTCCAGG',.3451956),
 ('GATGGTAGATGGAGACTCAGNGG','GGTAGGAAATGGAGAGGCAGAGG',5.029583e-5),
 ('GAGTCCGAGCAGAAGAAGAAAGG','GAGTTAGAGCAGAAGAAGAAAGG',.9269089),
]


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return h.hexdigest()


def check_inputs(root):
    if digest(HERE/'pins.json')!='0f755156053bc2b7a6eb67a8b1081d39072c5cee6d6d6e42a64df7127c039657':
        raise ValueError('Pinned digest table changed')
    pins=json.loads((HERE/'pins.json').read_text())
    for rel,want in pins.items():
        if digest(root/rel)!=want: raise ValueError('SHA-256 mismatch: '+rel)
    return pins


def read_data(root,name):
    a,b,upper=DATASETS[name]
    with (root/'inputs'/name).open(newline='') as f:
        rows=list(csv.DictReader(f))
    pairs=[(r[a].upper(),r[b].upper()) if upper else (r[a],r[b]) for r in rows]
    labels=np.array([int(float(r['label'])) for r in rows],dtype=np.int64)
    return pairs,labels


def write_table(path,header,rows):
    with path.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(rows)


def audit_presence(root,targets):
    """Rebuild the TSV's conservative coverage flags from the two complete CSVs."""
    by_spacer={r['spacer']:r['guide_id'] for r in targets}
    wanted={(r['guide_id'],r['off_target_19mer']):set() for r in targets}
    for name in list(DATASETS)[:2]:
        a,b,_=DATASETS[name]
        with (root/'inputs'/name).open(newline='') as f:
            for row in csv.DictReader(f):
                spacer=row[a].replace('_','').replace('-','')[:20]
                if spacer not in by_spacer:continue
                guide=by_spacer[spacer]
                proto=row[b].replace('_','').replace('-','').upper()[:-3]
                kind={19:'gapped19',20:'ungapped20',21:'dnabulge21'}.get(len(proto))
                if kind is None:continue
                for j in range(len(proto)-18):
                    key=(guide,proto[j:j+19])
                    if key in wanted:wanted[key].add(kind)
    kinds={};flags={}
    for row in targets:
        found=wanted[row['guide_id'],row['off_target_19mer']]
        flag='yes' if found else 'no';kind='+'.join(sorted(found)) if found else '-'
        if row['in_published_files']!=flag or row['held_by']!=kind:
            raise ValueError('Target coverage flag differs from source CSVs')
        flags[flag]=flags.get(flag,0)+1;kinds[kind]=kinds.get(kind,0)+1
    return dict(rows=len(targets),by_in_published_files=flags,by_held_by=kinds,
                definition='Strip gaps/padding, remove final three PAM bases, and test contiguous '
                           '19-base protospacer windows for the same 20-base guide spacer.')


def run(root,out):
    root,out=root.resolve(),out.resolve()
    # Even optional reruns stay inside this authorized task directory.
    if not out.is_relative_to(HERE.parent): raise ValueError('Output must stay inside the task directory')
    out.mkdir(parents=True,exist_ok=True)
    pins=check_inputs(root)
    events=[]
    def event(s):
        events.append(dict(event=s,utc=datetime.now(timezone.utc).isoformat()))
        (out/'execution.json').write_text(json.dumps(events,indent=2)+'\n')
        print(s,flush=True)
    event('pinned_inputs_verified')
    model=NumpyModel(root/'inputs'/WEIGHT)
    ep=batches(model,[(a,b) for a,b,_ in EXAMPLES])
    examples=[dict(on_seq=a,off_seq=b,readme_score=v,recomputed_score=float(p),
                   absolute_difference=abs(float(p)-v)) for (a,b,v),p in zip(EXAMPLES,ep)]
    (out/'readme_examples.json').write_text(json.dumps(examples,indent=2)+'\n')
    validation={}
    for name in DATASETS:
        pairs,y=read_data(root,name)
        p=batches(model,pairs)
        validation[name]=metrics(y,p)
        write_table(out/(name+'.scores.tsv'),['data_row_1based','label','score'],
                    ((i,int(label),format(float(score),'.17g'))
                     for i,(label,score) in enumerate(zip(y,p),1)))
        event('validated:'+name)
        print(json.dumps(validation[name],sort_keys=True),flush=True)
    gate=validation_gate(validation)
    validation['gate']=gate
    (out/'validation.json').write_text(json.dumps(validation,indent=2,sort_keys=True)+'\n')
    if not gate['pass_for_scoring']:
        event('validation_failed_no_target_scores');raise RuntimeError('Validation failed')
    event('validation_passed_before_target_scoring')
    with (root/'inputs/unreachable_sequences.tsv').open(newline='') as f:
        targets=list(csv.DictReader(f,delimiter='\t'))
    if len(targets)!=10888: raise ValueError('Unexpected target row count')
    if len({(r['guide_id'],r['off_target_19mer']) for r in targets})!=len(targets):
        raise ValueError('Duplicate target identity')
    if any(len(r['spacer'])!=20 or len(r['off_target_19mer'])!=19 for r in targets):
        raise ValueError('Unexpected sequence length')
    presence=audit_presence(root,targets)
    (out/'presence_audit.json').write_text(json.dumps(presence,sort_keys=True,indent=2)+'\n')
    p=batches(model,[(r['spacer'],r['off_target_19mer']) for r in targets])
    selected=above(p)
    model32=NumpyModel(root/'inputs'/WEIGHT,dtype=np.float32)
    from model_numpy import encode
    pairs=[(r['spacer'],r['off_target_19mer']) for r in targets]
    p32=np.concatenate([model32.predict(encode(pairs[j:j+512]).astype(np.float32))
                        for j in range(0,len(pairs),512)]).astype(np.float64)
    write_table(out/'float32_target_scores.tsv',['data_row_1based','score'],
                ((i,format(float(score),'.17g')) for i,score in enumerate(p32,1)))
    precision_check=dict(rows=len(p),float64_above=int(selected.sum()),
                         float32_above=int(above(p32).sum()),
                         changed_cutoff_decisions=int(np.count_nonzero(selected!=above(p32))),
                         maximum_absolute_score_difference=float(np.max(np.abs(p-p32))))
    (out/'precision_check.json').write_text(json.dumps(precision_check,sort_keys=True,indent=2)+'\n')
    if precision_check['changed_cutoff_decisions']:
        raise RuntimeError('Float32 and float64 threshold decisions differ')
    write_table(out/'target_scores.tsv',
                ['data_row_1based','guide_id','spacer','off_target_19mer',
                 'in_published_files','held_by','score','above_3_over_1000'],
                ((i,r['guide_id'],r['spacer'],r['off_target_19mer'],r['in_published_files'],
                  r['held_by'],format(float(score),'.17g'),int(cut))
                 for i,(r,score,cut) in enumerate(zip(targets,p,selected),1)))
    summary=dict(all=distribution(p),
                 by_guide={g:distribution(p[[r['guide_id']==g for r in targets]])
                           for g in sorted({r['guide_id'] for r in targets})},
                 by_in_published_files={g:distribution(p[[r['in_published_files']==g for r in targets]])
                                       for g in sorted({r['in_published_files'] for r in targets})},
                 input_policy='Literal TSV spacer/off_target_19mer; published encoder separately '
                              'left-pads each to 24. No PAM, gap placement, flanks or alignment added.',
                 cutoff=dict(exact='3/1000',comparison='strictly greater',
                             source='Lin et al. 2020, sensitivity-analysis section before Figure 14',
                             calibration='CIRCLE-seq leave-one-guide-out; not calibrated for this combined-data checkpoint'))
    (out/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    event('target_scoring_completed')
    outputs={str(p.relative_to(out)):digest(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='outputs.sha256.json'}
    (out/'outputs.sha256.json').write_text(json.dumps(outputs,sort_keys=True,indent=2)+'\n')
    return validation,summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root',type=Path,default=HERE)
    parser.add_argument('--out',type=Path,default=HERE/'results')
    a=parser.parse_args()
    run(a.input_root,a.out)
