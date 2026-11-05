#!/usr/bin/env python3
"""Offline replay, with independently coded encoder, convolution and recurrent pass."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='4'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import argparse
import csv
import hashlib
import importlib.util
import json
import subprocess
from fractions import Fraction
import math
import numpy as np
import h5py
from scipy.stats import rankdata

PINS_SHA256='0f755156053bc2b7a6eb67a8b1081d39072c5cee6d6d6e42a64df7127c039657'
WEIGHT_SHA256='a7f6aa381520f5c68fa1f099a6ef3ebc3b8ce846709b97dfde2053f26ca62f80'
TARGET_SHA256='5442bf2b3202c7eae45414769186df9c204c418639a2c92ebd877583d8eed27d'
DATASETS={
 'CIRCLE_seq_10gRNA_wholeDataset.csv':('sgRNA_seq','off_seq',False),
 'elevation_6gRNA_wholeDataset.csv':('crRNA','DNA',False),
 'Kleinstiver_5gRNA_wholeDataset.csv':('sgRNA_seq','off_seq',True),
 'Listgarten_22gRNA_wholeDataset.csv':('sgRNA_seq','off_seq',True),
}


def require(ok,message):
    if not ok: raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def encode_reference(pairs):
    """Independent vectorized formulation of the deposited character encoder."""
    a=np.array([list(x.rjust(24,'-').encode('ascii')) for x,_ in pairs],dtype=np.uint8)
    b=np.array([list(x.rjust(24,'-').encode('ascii')) for _,x in pairs],dtype=np.uint8)
    a=np.where(a==ord('N'),b,a)
    table=np.zeros((256,5),dtype=np.uint8)
    direction=np.zeros(256,dtype=np.int64)
    for i,ch in enumerate('ATGC_'):table[ord(ch),i]=1
    for ch,rank in zip('AGCT_',(5,4,3,2,1)):direction[ord(ch)]=rank
    valid=np.array([ord(c) for c in 'ATGC_-'])
    require(np.isin(a,valid).all() and np.isin(b,valid).all(),'Unsupported nucleotide')
    both=(a!=ord('-')) & (b!=ord('-'))
    return np.concatenate((table[a]|table[b],
                           (both & (direction[a]>direction[b]))[...,None],
                           (both & (direction[a]<direction[b]))[...,None]),axis=-1).astype(np.float64)


class ReferenceModel:
    def __init__(self,path,dtype=np.float64):
        self.w={}
        with h5py.File(path,'r') as f:
            def visit(k,v):
                if isinstance(v,h5py.Dataset):self.w[k]=np.asarray(v,dtype=dtype)
            f.visititems(visit)

    def recurrent(self,x,backward):
        name='backward' if backward else 'forward'
        prefix=f'bidirectional_1/bidirectional_1/{name}_LSTM_out/'
        w,u,b=(self.w[prefix+k+':0'] for k in ('kernel','recurrent_kernel','bias'))
        kernel=np.split(w,4,axis=1)
        recurrent=np.split(u,4,axis=1)
        bias=np.split(b,4)
        if backward:x=x[:,::-1]
        inputs=[x@k+v for k,v in zip(kernel,bias)]
        state=np.zeros((len(x),15),dtype=x.dtype);cell=state.copy();trajectory=[]
        for t in range(24):
            gates=[z[:,t]+state@r for z,r in zip(inputs,recurrent)]
            i,f,g,o=gates
            i=np.minimum(1,np.maximum(0,0.5+i/5))
            f=np.minimum(1,np.maximum(0,0.5+f/5))
            o=np.minimum(1,np.maximum(0,0.5+o/5))
            cell=f*cell+i*np.tanh(g)
            state=o*np.tanh(cell)
            trajectory.append(state)
        result=np.stack(trajectory,axis=1)
        return result[:,::-1] if backward else result

    def predict(self,x):
        features=[x]
        for j,k in enumerate((1,2,3,5),1):
            p=f'conv2d_{j}/conv2d_{j}/'
            padded=np.pad(x,((0,0),((k-1)//2,k//2),(0,0)))
            windows=np.lib.stride_tricks.sliding_window_view(padded,k,axis=1)
            flat=windows.transpose(0,1,3,2).reshape(len(x),24,7*k)
            y=flat@self.w[p+'kernel:0'].reshape(7*k,10)+self.w[p+'bias:0']
            features.append(y.clip(0))
        x=np.concatenate(features,axis=2)
        y=np.concatenate((self.recurrent(x,False),self.recurrent(x,True)),axis=2)
        y=y.reshape(len(x),-1)
        for name in ('dense_1','dense_2'):
            p=f'{name}/{name}/'
            y=(y@self.w[p+'kernel:0']+self.w[p+'bias:0']).clip(0)
        z=(y@self.w['main_output/main_output/kernel:0']+
           self.w['main_output/main_output/bias:0']).ravel()
        # Equivalent logaddexp expression gives an additional sigmoid cross-check.
        return np.exp(-np.logaddexp(0,-z))


def replay(model,pairs,predictions):
    require(len(pairs)==len(predictions),'Prediction row count differs')
    largest=0.
    for start in range(0,len(pairs),512):
        pp=pairs[start:start+512]
        got=model.predict(encode_reference(pp))
        want=predictions[start:start+512]
        largest=max(largest,float(np.max(np.abs(got-want))))
        require(np.allclose(got,want,rtol=2e-11,atol=2e-12),'Independent forward replay differs')
        require(np.array_equal(exact_cut(got),exact_cut(want)),'Cutoff decisions differ')
    return largest


def exact_cut(p):
    return np.array([Fraction.from_float(float(v))>Fraction(3,1000) for v in p])


def compare(actual,expected,path='summary'):
    if isinstance(expected,dict):
        require(isinstance(actual,dict) and set(actual)==set(expected),path+': keys differ')
        for k in expected:compare(actual[k],expected[k],path+'.'+k)
    elif isinstance(expected,list):
        require(len(actual)==len(expected),path+': length differs')
        for j,(a,b) in enumerate(zip(actual,expected)):compare(a,b,path+f'[{j}]')
    elif isinstance(expected,float):
        require(math.isclose(float(actual),expected,rel_tol=2e-12,abs_tol=2e-14),path+': value differs')
    else:require(actual==expected,path+': value differs')


def reference_metrics(y,p):
    P=int(sum(y));N=len(y)-P
    ranks=rankdata(p,method='average')
    rank_sum_twice=int(sum(2*ranks[y==1]))
    auc=Fraction(rank_sum_twice-P*(P+1),2*P*N)
    idx=np.argsort(-p,kind='stable');ys=y[idx];ps=p[idx]
    last=np.r_[np.flatnonzero(np.diff(ps)),len(ps)-1]
    tp=np.cumsum(ys)[last];rec=np.r_[0,tp/P]
    prec=np.r_[1,tp/(last+1)]
    cuts=exact_cut(p)
    return dict(rows=len(y),positive=P,negative=N,auroc=float(auc),auroc_exact=str(auc),
                auprc_trapezoidal=float(np.sum(np.diff(rec)*(prec[1:]+prec[:-1])/2)),
                label_0_mean=float(np.mean(p[y==0])),label_1_mean=float(np.mean(p[y==1])),
                label_0_median=float(np.median(p[y==0])),label_1_median=float(np.median(p[y==1])),
                above_cutoff_label_0=int(cuts[y==0].sum()),above_cutoff_label_1=int(cuts[y==1].sum()),
                first_data_row_1based=1,last_data_row_1based=len(y))


def reference_distribution(scores):
    p=list(map(float,scores));n=len(p)
    if not n:return dict(rows=0)
    rational=list(map(Fraction.from_float,p));threshold=Fraction(3,1000)
    count=sum(v>threshold for v in rational)
    ordered=sorted(p)
    def quantile(q):
        h=(n-1)*q;i=math.floor(h);j=math.ceil(h)
        return ordered[i]+(h-i)*(ordered[j]-ordered[i])
    bounds=['0','0.000001','0.0001','0.001','0.003','0.01','0.1','0.5','1']
    hist=[]
    for j,(a,b) in enumerate(zip(bounds[:-1],bounds[1:])):
        lo,hi=Fraction(a),Fraction(b);inc=j==len(bounds)-2
        hist.append(dict(lower=a,upper=b,upper_inclusive=inc,
                         rows=sum(lo<=v and (v<=hi if inc else v<hi) for v in rational)))
    mean=math.fsum(p)/n
    return dict(rows=n,above_3_over_1000=count,equal_3_over_1000=sum(v==threshold for v in rational),
                fraction_above_exact=str(Fraction(count,n)),percent_above=100*count/n,
                mean=mean,population_sd=math.sqrt(math.fsum((v-mean)**2 for v in p)/n),
                quantiles_linear={str(q):quantile(q) for q in (0,.01,.05,.25,.5,.75,.95,.99,1)},
                histogram=hist,nearest_cutoff_distance=float(min(abs(v-threshold) for v in rational)))


def verify(args):
    require(sha(HERE/'pins.json')==PINS_SHA256,'Pinned digest table changed')
    pins=json.loads((HERE/'pins.json').read_text())
    root=args.input_root.resolve();results=args.results.resolve()
    paths={rel:root/rel for rel in pins}
    if args.target_input:paths['inputs/unreachable_sequences.tsv']=args.target_input
    if args.weight_input:paths['inputs/CRISPR_Net_CIRCLE_elevation_SITE_weights.h5']=args.weight_input
    for rel,want in pins.items():
        require(sha(paths[rel])==want,'SHA-256 mismatch: '+rel)
    require(pins['inputs/unreachable_sequences.tsv']==TARGET_SHA256,'Target pin differs')
    require(pins['inputs/CRISPR_Net_CIRCLE_elevation_SITE_weights.h5']==WEIGHT_SHA256,'Weight pin differs')
    print('PASS: all pinned inputs and sources',flush=True)
    manifest=json.loads((results/'outputs.sha256.json').read_text())
    for rel,want in manifest.items():require(sha(results/rel)==want,'Output digest mismatch: '+rel)
    spec=importlib.util.spec_from_file_location('upstream_encoder',paths['sources/Encoder_sgRNA_off.py'])
    upstream=importlib.util.module_from_spec(spec);spec.loader.exec_module(upstream)
    model=ReferenceModel(paths['inputs/CRISPR_Net_CIRCLE_elevation_SITE_weights.h5'])
    # Check every operation-bearing architecture field, not merely tensor dimensions.
    arch=json.loads(paths['sources/CRISPR_Net_CIRCLE_elevation_SITE_structure.json'].read_text())
    layers={v['name']:v for v in arch['config']['layers']}
    for j,k in enumerate((1,2,3,5),1):
        c=layers[f'conv2d_{j}']['config']
        require(c['kernel_size']==[1,k] and c['padding']=='same' and c['activation']=='linear'
                and c['strides']==[1,1] and c['data_format']=='channels_last','Convolution architecture differs')
        require(layers[f'activation_{j}']['config']['activation']=='relu','Branch activation differs')
    sources=[v[0] for v in layers['concatenate_1']['inbound_nodes'][0]]
    require(sources==['main_input','activation_1','activation_2','activation_3','activation_4'],'Concatenation order differs')
    c=layers['bidirectional_1']['config']['layer']['config']
    require(c['recurrent_activation']=='hard_sigmoid' and c['activation']=='tanh' and
            c['return_sequences'] and c['units']==15 and c['implementation']==1,'LSTM architecture differs')
    require(layers['dense_1']['config']['activation']=='relu' and
            layers['dense_2']['config']['activation']=='relu' and
            layers['main_output']['config']['activation']=='sigmoid','Dense activation differs')
    validation=json.loads((results/'validation.json').read_text())
    max_error=0.;total=0
    for name,(a,b,upper) in DATASETS.items():
        with paths['inputs/'+name].open(newline='') as f:rows=list(csv.DictReader(f))
        pairs=[(r[a].upper(),r[b].upper()) if upper else (r[a],r[b]) for r in rows]
        y=np.array([int(float(r['label'])) for r in rows])
        with (results/(name+'.scores.tsv')).open(newline='') as f:stored=list(csv.DictReader(f,delimiter='\t'))
        require([int(r['data_row_1based']) for r in stored]==list(range(1,len(rows)+1)),name+': row identities differ')
        require(np.array_equal(y,[int(r['label']) for r in stored]),name+': labels differ')
        p=np.array([float(r['score']) for r in stored])
        probes=[pairs[i] for i in np.linspace(0,len(pairs)-1,256,dtype=int)]
        require(np.array_equal(encode_reference(probes),np.array([upstream.Encoder(x,z).on_off_code for x,z in probes])),
                name+': upstream encoder differs')
        error=replay(model,pairs,p);max_error=max(error,max_error);total+=len(p)
        compare(reference_metrics(y,p),validation[name],name)
        print(f'PASS: {name}: {len(rows)} rows; max independent score difference {error:.3g}',flush=True)
    gate=validation['gate']
    require(gate['pass_for_scoring'] and gate['own_label_separation'] and gate['three_published_rounded_benchmarks_match'],'Validation gate not passed')
    for name in list(DATASETS)[:2]:
        require(Fraction(validation[name]['auroc_exact'])>Fraction(9,10) and
                validation[name]['label_1_mean']>validation[name]['label_0_mean'],'Own-label separation failed')
    require(format(validation['Kleinstiver_5gRNA_wholeDataset.csv']['auroc'],'.3f')=='0.995','Five-guide AUROC failed')
    require(format(validation['Listgarten_22gRNA_wholeDataset.csv']['auroc'],'.3f')=='0.995','22-guide AUROC failed')
    require(format(validation['Listgarten_22gRNA_wholeDataset.csv']['auprc_trapezoidal'],'.3f')=='0.317','22-guide PR failed')
    require(gate['kleinstiver_auprc_matches_published']==
            (format(validation['Kleinstiver_5gRNA_wholeDataset.csv']['auprc_trapezoidal'],'.3f')=='0.329'),
            'Five-guide AUPRC discrepancy flag differs')
    with paths['inputs/unreachable_sequences.tsv'].open(newline='') as f:targets=list(csv.DictReader(f,delimiter='\t'))
    with (results/'target_scores.tsv').open(newline='') as f:scored=list(csv.DictReader(f,delimiter='\t'))
    require(len(targets)==len(scored)==10888,'Target row count differs')
    require(len({(r['guide_id'],r['off_target_19mer']) for r in targets})==10888,'Duplicate target identities')
    # Independently rebuild conservative substring presence, not just trust TSV flags.
    target_map={r['spacer']:{} for r in targets}
    for r in targets:target_map[r['spacer']][r['off_target_19mer']]=set()
    for name,(gc,oc,_) in list(DATASETS.items())[:2]:
        with paths['inputs/'+name].open(newline='') as f:
            for r in csv.DictReader(f):
                guide=r[gc].translate(str.maketrans('','','_-'))[:20]
                if guide not in target_map:continue
                proto=r[oc].translate(str.maketrans('','','_-')).upper()[:-3]
                if len(proto) not in (19,20,21):continue
                kind={19:'gapped19',20:'ungapped20',21:'dnabulge21'}[len(proto)]
                for j in range(len(proto)-18):
                    word=proto[j:j+19]
                    if word in target_map[guide]:target_map[guide][word].add(kind)
    flags={};kinds={}
    for r in targets:
        found=target_map[r['spacer']][r['off_target_19mer']]
        flag='yes' if found else 'no';kind='+'.join(sorted(found)) if found else '-'
        require(r['in_published_files']==flag and r['held_by']==kind,'Source-derived coverage differs')
        flags[flag]=flags.get(flag,0)+1;kinds[kind]=kinds.get(kind,0)+1
    presence=json.loads((results/'presence_audit.json').read_text())
    compare(dict(rows=len(targets),by_in_published_files=flags,by_held_by=kinds),
            {k:presence[k] for k in ('rows','by_in_published_files','by_held_by')},'presence audit')
    for i,(a,b) in enumerate(zip(targets,scored),1):
        require(int(b['data_row_1based'])==i and all(b[k]==v for k,v in a.items()),'Target metadata differs')
    pairs=[(r['spacer'],r['off_target_19mer']) for r in targets]
    # All target encodings also compared against the original published encoder.
    require(np.array_equal(encode_reference(pairs),np.array([upstream.Encoder(a,b).on_off_code for a,b in pairs])),
            'Target upstream encoder differs')
    p=np.array([float(r['score']) for r in scored])
    error=replay(model,pairs,p);max_error=max(max_error,error)
    require(np.array_equal(exact_cut(p),[bool(int(r['above_3_over_1000'])) for r in scored]),'Stored decisions differ')
    with (results/'float32_target_scores.tsv').open(newline='') as f:stored32=list(csv.DictReader(f,delimiter='\t'))
    require([int(r['data_row_1based']) for r in stored32]==list(range(1,len(p)+1)),'Float32 row IDs differ')
    p32=np.array([float(r['score']) for r in stored32])
    ref32=ReferenceModel(paths['inputs/CRISPR_Net_CIRCLE_elevation_SITE_weights.h5'],dtype=np.float32)
    # Different legal float32 accumulation orders need a rounding tolerance.
    for start in range(0,len(pairs),512):
        got=ref32.predict(encode_reference(pairs[start:start+512]).astype(np.float32))
        require(np.allclose(got,p32[start:start+512],rtol=2e-5,atol=2e-5),'Float32 replay differs')
        require(np.array_equal(exact_cut(got),exact_cut(p32[start:start+512])),'Float32 decisions differ')
    precision=dict(rows=len(p),float64_above=int(exact_cut(p).sum()),float32_above=int(exact_cut(p32).sum()),
                   changed_cutoff_decisions=int(np.count_nonzero(exact_cut(p)!=exact_cut(p32))),
                   maximum_absolute_score_difference=float(np.max(np.abs(p-p32))))
    compare(precision,json.loads((results/'precision_check.json').read_text()),'precision check')
    require(precision['changed_cutoff_decisions']==0,'Precision changes a threshold decision')
    summary=json.loads((results/'summary.json').read_text())
    compare(reference_distribution(p),summary['all'])
    require(set(summary['by_guide'])=={r['guide_id'] for r in targets},'Guide groups differ')
    require(set(summary['by_in_published_files'])=={r['in_published_files'] for r in targets},'Presence groups differ')
    for field,key in [('by_guide','guide_id'),('by_in_published_files','in_published_files')]:
        for group,want in summary[field].items():
            compare(reference_distribution(p[[r[key]==group for r in targets]]),want,field+'.'+group)
    require(summary['cutoff']['exact']=='3/1000' and summary['cutoff']['comparison']=='strictly greater','Cutoff differs')
    examples=json.loads((results/'readme_examples.json').read_text())
    for e in examples:
        score=float(model.predict(encode_reference([(e['on_seq'],e['off_seq'])]))[0])
        require(math.isclose(score,e['recomputed_score'],rel_tol=2e-11,abs_tol=2e-12),'README replay differs')
        require(math.isclose(abs(score-e['readme_score']),e['absolute_difference'],rel_tol=2e-11,abs_tol=2e-12),'README difference differs')
    events=json.loads((results/'execution.json').read_text())
    event_names=[r['event'] for r in events]
    require(event_names==['pinned_inputs_verified']+['validated:'+n for n in DATASETS]+
            ['validation_passed_before_target_scoring','target_scoring_completed'],'Execution gate ordering differs')
    controls=json.loads((HERE/'corruption_control/receipt.json').read_text())
    require(controls['status']=='PASS' and len(controls['tests'])==2,'Corruption receipt differs')
    for test in controls['tests']:
        name=test['input'];original=paths['inputs/'+name].read_bytes()
        changed_path=HERE/'corruption_control'/(name+'.one_byte_changed')
        changed=changed_path.read_bytes()
        positions=[i for i,(a,b) in enumerate(zip(original,changed)) if a!=b]
        require(len(original)==len(changed)==test['length'] and len(positions)==test['changed_bytes']==1,
                'Corruption is not an equal-length one-byte substitution')
        require(positions==[test['byte_offset_0based']] and
                original[positions[0]]==test['original_byte'] and changed[positions[0]]==test['altered_byte'],
                'Corruption byte receipt differs')
        require(sha(paths['inputs/'+name])==test['original_sha256'] and
                sha(changed_path)==test['altered_sha256'],'Corruption digest receipt differs')
        option='--target-input' if name.endswith('.tsv') else '--weight-input'
        proc=subprocess.run([sys.executable,'-I','-B',str(HERE/'verify.py'),option,str(changed_path)],
                            capture_output=True,text=True)
        require(proc.returncode==test['verifier_exit_status']==1 and
                proc.stderr.strip()==test['verifier_stderr']=='REJECT: SHA-256 mismatch: inputs/'+name,
                'Corrupted input was not rejected as recorded')
    if (HERE.parent/'FINDINGS.md').exists():
        from render_findings import render
        require((HERE.parent/'FINDINGS.md').read_text()==render(),'FINDINGS.md differs from the replayed report')
    print('PASS: one-byte corruption controls and generated findings',flush=True)
    print(f'PASS: all {len(targets)} target scores, every summary number, and all README discrepancies',flush=True)
    print(json.dumps(dict(status='PASS',validation_rows_replayed=total,target_rows_replayed=len(p),
                          largest_absolute_independent_score_difference=max_error,
                          target_above_cutoff=summary['all']['above_3_over_1000'],
                          absent_above_cutoff=summary['by_in_published_files']['no']['above_3_over_1000']),sort_keys=True),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root',type=Path,default=HERE)
    p.add_argument('--results',type=Path,default=HERE/'results')
    p.add_argument('--target-input',type=Path)
    p.add_argument('--weight-input',type=Path)
    args=p.parse_args()
    try:verify(args)
    except (ValueError,OSError,KeyError,AssertionError) as e:
        print('REJECT: '+str(e),file=sys.stderr,flush=True);sys.exit(1)
