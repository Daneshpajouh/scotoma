"""Deterministic summaries; decisions at 3/1000 use exact integer arithmetic."""
from fractions import Fraction
import math
import numpy as np


def above(values, numerator=3, denominator=1000):
    return np.array([v.as_integer_ratio()[0]*denominator >
                     numerator*v.as_integer_ratio()[1] for v in map(float,values)])


def metrics(labels,scores):
    y = np.asarray(labels,dtype=np.int64)
    p = np.asarray(scores,dtype=np.float64)
    if not (np.isfinite(p).all() and np.isin(y,[0,1]).all()):
        raise ValueError('Invalid labels or predictions')
    P,N = int(y.sum()),int(len(y)-y.sum())
    if not P or not N: raise ValueError('Both classes required')
    order = np.argsort(p,kind='stable')
    ss,yy = p[order],y[order]
    ends = np.r_[np.flatnonzero(np.diff(ss))+1,len(ss)]
    starts = np.r_[0,ends[:-1]]
    pos = np.add.reduceat(yy,starts)
    neg = ends-starts-pos
    neg_before = np.cumsum(neg)-neg
    twice_u = int(np.sum(pos*(2*neg_before+neg),dtype=np.int64))
    roc = Fraction(twice_u,2*P*N)
    # Tie-group trapezoidal precision-recall AUC, matching upstream metrics.auc.
    tp = np.cumsum(pos[::-1]); fp = np.cumsum(neg[::-1])
    recall = np.r_[0,tp/P]; precision = np.r_[1,tp/(tp+fp)]
    pr = float(np.trapezoid(precision,recall))
    selected = above(p)
    return dict(rows=len(y),positive=P,negative=N,auroc=float(roc),
                auroc_exact=str(roc),auprc_trapezoidal=pr,
                label_0_mean=float(p[y==0].mean()),label_1_mean=float(p[y==1].mean()),
                label_0_median=float(np.median(p[y==0])),label_1_median=float(np.median(p[y==1])),
                above_cutoff_label_0=int(selected[y==0].sum()),
                above_cutoff_label_1=int(selected[y==1].sum()),
                first_data_row_1based=1,last_data_row_1based=len(y))


def distribution(scores):
    p = np.asarray(scores,dtype=np.float64)
    if len(p)==0: return dict(rows=0)
    if not (np.isfinite(p).all() and ((p>=0)&(p<=1)).all()):
        raise ValueError('Scores must lie in [0,1]')
    cut = above(p)
    boundaries = ['0','0.000001','0.0001','0.001','0.003','0.01','0.1','0.5','1']
    counts = []
    # Bins [low, high), except the final bin includes 1.
    for j,(lo,hi) in enumerate(zip(boundaries[:-1],boundaries[1:])):
        a,b = Fraction(lo),Fraction(hi)
        n = sum(a <= Fraction.from_float(float(v)) and
                (Fraction.from_float(float(v)) <= b if j==len(boundaries)-2
                 else Fraction.from_float(float(v)) < b) for v in p)
        counts.append(dict(lower=lo,upper=hi,upper_inclusive=j==len(boundaries)-2,rows=n))
    assert sum(c['rows'] for c in counts)==len(p)
    quantiles = {str(q):float(np.quantile(p,q,method='linear'))
                 for q in (0,.01,.05,.25,.5,.75,.95,.99,1)}
    return dict(rows=len(p),above_3_over_1000=int(cut.sum()),
                equal_3_over_1000=sum(Fraction.from_float(float(v))==Fraction(3,1000) for v in p),
                fraction_above_exact=str(Fraction(int(cut.sum()),len(p))),
                percent_above=100*int(cut.sum())/len(p),mean=float(p.mean()),
                population_sd=float(p.std()),quantiles_linear=quantiles,histogram=counts,
                nearest_cutoff_distance=float(min(abs(Fraction.from_float(float(v))-Fraction(3,1000))
                                                  for v in p)))


def validation_gate(result):
    """Frozen before target scoring: behavior and published rounded benchmarks."""
    own = all(Fraction(result[k]['auroc_exact']) > Fraction(9,10) and
              result[k]['label_1_mean'] > result[k]['label_0_mean']
              for k in ('CIRCLE_seq_10gRNA_wholeDataset.csv','elevation_6gRNA_wholeDataset.csv'))
    k = result['Kleinstiver_5gRNA_wholeDataset.csv']
    l = result['Listgarten_22gRNA_wholeDataset.csv']
    benchmark = (format(k['auroc'],'.3f')=='0.995' and
                 format(l['auroc'],'.3f')=='0.995' and
                 format(l['auprc_trapezoidal'],'.3f')=='0.317')
    return dict(pass_for_scoring=bool(own and benchmark),own_label_separation=bool(own),
                three_published_rounded_benchmarks_match=bool(benchmark),
                kleinstiver_auprc_matches_published=bool(format(k['auprc_trapezoidal'],'.3f')=='0.329'),
                criteria='Own-data AUROC > 9/10 and positive mean > negative mean in both files; '
                         'both Table 2 AUROCs round to .995 and Listgarten AUPRC rounds to .317. '
                         'Kleinstiver AUPRC and README example discrepancies are retained separately.')
