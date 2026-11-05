#!/usr/local/bin/python3
"""LP support proposals, exact reconstruction, exhaustive integer transform separation."""
import os,sys,json,time,math,subprocess,argparse,hashlib
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'deps'))
import numpy as np
from scipy.optimize import linprog
from checker import Cell
REPO=Path(__file__).resolve().parent.parent

def log(**kw): print(json.dumps(kw),flush=True)

def decode(i,L):
 a=[0]*L
 for k in range(L-1,-1,-1): a[k]=i%3;i//=3
 return tuple(a)

def encode(a):
 v=0
 for x in a:v=v*3+x
 return v

class LightCell(Cell):
 def __init__(self,N,r,m):
  if not all(type(v) is int for v in [N,r,m]) or not (0<=N<=19 and 1<=r<=N+1 and 0<=m<=N):raise ValueError('invalid cell or unsupported N; kernels support N<=19')
  self.N,self.r,self.m=N,r,m;self.L=r-1;self.C=N-self.L
  self._suffix_counts={};self._prefix_signatures={}
  self.library=tuple(((0,)*j+(1,)*(self.L-j)) for j in range(r))
 @lru_cache(None)
 def reflect(self,o):
  i,w=o;a=decode(i,self.L)
  return encode(tuple(1-x if x<2 else 2 for x in reversed(a))),w
 def canonical(self,o):return min(o,self.reflect(o))
 def in_union(self,o):
  i,w=o;return self.dmin_id(i)+w<=self.m
 @lru_cache(None)
 def dmin_id(self,i):
  a=decode(i,self.L);d=sum(x!=1 for x in a);best=d
  for x in a:d+=(x!=0)-(x!=1);best=min(best,d)
  return best
 def raw(self,o):return decode(o[0],self.L),o[1]
 @lru_cache(None)
 def masks(self,i):
  a=decode(i,self.L)
  return tuple(sum(1<<k for k,x in enumerate(a) if x==v) for v in [0,1,2])
 def kernel_id(self,a,b):
  ai,p=a;bi,w=b;az,ao,ax=self.masks(ai);bz,bo,bx=self.masks(bi)
  forced=self.L-((az&bz)|(ao&bo)|(ax&bx)).bit_count();both=(ax&bx).bit_count()
  return self.probability(w,p,forced,both)
 def skernel(self,a,b):
  ref=self.reflect(b)
  return self.kernel_id(a,b) if b==ref else (self.kernel_id(a,b)+self.kernel_id(a,ref))/2
 def size_id(self,o):return self.size(self.raw(o))
 def terms(self,weights):
  return [dict(prefix=''.join('X' if x==2 else str(x) for x in decode(i,self.L)),suffix_weight=w,mass=str(v)) for (i,w),v in weights if v]
 def parse_id(self,terms,packing=False):
  ans={}
  for t in terms:
   if isinstance(t['mass'],float):raise ValueError('float mass')
   a=tuple(2 if x in 'X2' else int(x) for x in t['prefix']);w=t['suffix_weight'];v=F(t['mass'])
   if len(a)!=self.L or any(x not in (0,1,2) for x in a) or type(w)is not int or not 0<=w<=self.C or v<0:raise ValueError('bad term')
   o=(encode(a),w)
   if packing and not self.in_union(o):raise ValueError('packing outside U')
   ans[o]=ans.get(o,F(0))+v
  return [(o,v) for o,v in sorted(ans.items()) if v]
 def expand(self,weights):
  ans={}
  for o,v in weights:
   ref=self.reflect(o)
   for b,u in ([(o,v)] if ref==o else [(o,v/2),(ref,v/2)]):ans[b]=ans.get(b,F(0))+u
  return sorted((o,v) for o,v in ans.items() if v)

def sweep(cell,weights,packing,stem,keep=512):
 if cell.r>=14:
  from screened import sweep as screened_sweep
  return screened_sweep(cell,weights,packing,stem,keep)
 for o,v in weights:
  assert v>=0 and (not packing or cell.in_union(o))
 scale=math.lcm(*(v.denominator*cell.size_id(o) for o,v in weights)) if weights else 1
 inp=HERE/(stem+'.in');out=HERE/(stem+'.out')
 lines=[f'{cell.N} {cell.r} {cell.m} 8 {int(packing)} {keep} {scale} {len(weights)}']
 for (i,w),v in weights:
  val=v*scale/cell.size_id((i,w));assert val.denominator==1
  lines.append(f'{i} {w} {val.numerator}')
 inp.write_text('\n'.join(lines)+'\n')
 subprocess.run([str(HERE/'tensor_general'),str(inp),str(out)],check=True)
 fields={};candidates=[]
 for line in out.read_text().splitlines():
  k,*v=line.split()
  if k=='candidate':candidates.append((F(int(v[2]),scale),(int(v[0]),int(v[1]))))
  else:fields[k]=' '.join(v)
 assert int(fields['errors'])==0 and int(fields['scale'])==scale
 extreme=F(int(fields['extremum']),scale)
 passes=extreme<=1 if packing else extreme>=1
 assert passes==(int(fields['violations'])==0)
 report=dict(passes=passes,extremum=str(extreme),at=list(map(int,fields['arg'].split())),
  total_mass=str(sum((v for _,v in weights),F(0))),orbits_checked=int(fields['checked']),
  total_centre_orbits=int(fields['classes']),violations=int(fields['violations']),
  seconds=fields['seconds'],arithmetic='arbitrary precision integer convolution; Fraction comparisons')
 (HERE/(stem+'.json')).write_text(json.dumps(report,indent=2)+'\n')
 return report,candidates

def exact_solve(A,preferred=None):
 from flint import fmpq,fmpq_mat
 h=len(A);n=len(A[0]);mat=fmpq_mat([[fmpq(v.numerator,v.denominator) for v in row] for row in A]);rhs=fmpq_mat([[1] for _ in range(h)])
 if h==n:
  try:
   sol=mat.solve(rhs,algorithm='dixon');return [F(int(sol[i,0].numerator),int(sol[i,0].denominator)) for i in range(n)]
  except ZeroDivisionError:pass
 aug=fmpq_mat([[fmpq(v.numerator,v.denominator) for v in row]+[fmpq(1)] for row in A])
 rr,rank=aug.rref();ans=[F(0)]*n if preferred is None else list(preferred)
 pivots=[]
 for i in range(h):
  pivot=next((j for j in range(n) if rr[i,j]),None)
  if pivot is None:
   if rr[i,n]:raise ArithmeticError('inconsistent exact system')
   continue
  pivots.append((i,pivot))
 for i,pivot in reversed(pivots):
  val=F(int(rr[i,n].numerator),int(rr[i,n].denominator))
  for j in range(pivot+1,n):
   if rr[i,j]:val-=F(int(rr[i,j].numerator),int(rr[i,j].denominator))*ans[j]
  ans[pivot]=val
 return ans

def reconstruct(cell,active,chosen,result):
 errors=[]
 for cutoff in [1e-8,1e-10,1e-6]:
  cj=np.flatnonzero(result.x>cutoff).tolist();pi=np.flatnonzero(-result.ineqlin.marginals>cutoff).tolist()
  if not cj or not pi:continue
  co=[chosen[j] for j in cj];po=[active[i] for i in pi]
  block=[[cell.skernel(a,b) for b in co] for a in po]
  log(stage='reconstruct',cover=len(co),packing=len(po))
  for pref in [False,True]:
   try:
    cp=[F(float(result.x[j])).limit_denominator(10**9) for j in cj] if pref else None
    pp=[F(float(-result.ineqlin.marginals[i])).limit_denominator(10**9) for i in pi] if pref else None
    cm=exact_solve(block,cp);pm=exact_solve(list(map(list,zip(*block))),pp)
    if any(v<0 for v in cm+pm):raise ArithmeticError('negative exact reconstructed mass')
    if sum(cm,F(0))!=sum(pm,F(0)):raise ArithmeticError('unmatched exact objectives')
    # Independent common-denominator residuals; every equality decision is Fraction.
    Q=math.lcm(*{q.denominator for row in block for q in row})
    integer_block=[[q.numerator*(Q//q.denominator) for q in row] for row in block]
    Dc=math.lcm(*(v.denominator for v in cm));Dp=math.lcm(*(v.denominator for v in pm))
    ic=[v.numerator*(Dc//v.denominator) for v in cm];ip=[v.numerator*(Dp//v.denominator) for v in pm]
    for row in integer_block:
     if F(sum(a*b for a,b in zip(row,ic)),Dc*Q)!=1:raise ArithmeticError('cover residual')
    for col in zip(*integer_block):
     if F(sum(a*b for a,b in zip(col,ip)),Dp*Q)!=1:raise ArithmeticError('packing residual')
    return cell.expand(list(zip(co,cm))),cell.expand(list(zip(po,pm)))
   except ArithmeticError as e:errors.append(str(e))
 raise ArithmeticError(str(errors))

def lift_terms(terms,oldC):
 ans={}
 for t in terms:
  w=t['suffix_weight'];v=F(t['mass'])
  for symbol,wp,part in [('1',w,F(oldC-w,oldC)),('0',w-1,F(w,3*oldC)),('X',w-1,F(2*w,3*oldC))]:
   if part:
    o=(t['prefix']+symbol,wp);ans[o]=ans.get(o,F(0))+v*part
 return [dict(prefix=a,suffix_weight=w,mass=str(v)) for (a,w),v in sorted(ans.items())]

class Numeric:
 def __init__(self,cell):self.cell=cell;self.tabs={}
 def masks(self,orbits):
  ans=[]
  for i,w in orbits:
   a=decode(i,self.cell.L);ans.append([sum((1<<k) for k,x in enumerate(a) if x==j) for j in [0,1,2]]+[w])
  return np.array(ans,dtype=np.uint32)
 def matrix(self,points,centres):
  c=self.cell;cs=self.masks(centres);rs=self.masks([c.reflect(o) for o in centres]);mat=np.empty((len(points),len(centres)))
  for j,o in enumerate(points):
   i,p=o
   if p not in self.tabs:
    self.tabs[p]=np.array([[[float(c.probability(w,p,d,h)) for h in range(c.L+1)] for d in range(c.L+1)] for w in range(c.C+1)])
   a=decode(i,c.L);masks=[sum((1<<k) for k,x in enumerate(a) if x==v) for v in [0,1,2]]
   vals=[]
   for cc in [cs,rs]:
    same=np.bitwise_count((cc[:,0]&masks[0])|(cc[:,1]&masks[1])|(cc[:,2]&masks[2]));both=np.bitwise_count(cc[:,2]&masks[2]);d=c.L-same
    vals.append(self.tabs[p][cc[:,3],d,both])
   mat[j]=(vals[0]+vals[1])/2
  return mat

class WarmLP:
 def __init__(self,numeric,lp_seconds=240.0):
  import highspy
  if lp_seconds<=0:raise ValueError('LP proposal allowance must be positive')
  self.lp_seconds=lp_seconds
  self.highspy=highspy;self.h=highspy.Highs();self.nr=0;self.nc=0;self.numeric=numeric
  for key,value in [('output_flag',False),('solver','ipm' if numeric.cell.r>=15 else 'simplex'),('time_limit',240.0),('dual_feasibility_tolerance',1e-9),('primal_feasibility_tolerance',1e-9),('threads',1)]:self.h.setOptionValue(key,value)
 def solve(self,active,chosen):
  from scipy.sparse import csr_matrix,csc_matrix
  from types import SimpleNamespace
  # Add new columns first, then new rows; preserve the existing simplex basis.
  if len(chosen)>self.nc:
   nc=len(chosen)-self.nc
   if self.nr:
    block=csc_matrix(self.numeric.matrix(active[:self.nr],chosen[self.nc:]))
    self.h.addCols(nc,np.ones(nc),np.zeros(nc),np.full(nc,np.inf),block.nnz,block.indptr.astype(np.int32),block.indices.astype(np.int32),block.data)
   else:
    self.h.addVars(nc,np.zeros(nc),np.full(nc,np.inf));self.h.changeColsCost(nc,np.arange(nc,dtype=np.int32),np.ones(nc))
   self.nc=len(chosen)
  if len(active)>self.nr:
   nr=len(active)-self.nr;block=csr_matrix(self.numeric.matrix(active[self.nr:],chosen))
   self.h.addRows(nr,np.ones(nr),np.full(nr,np.inf),block.nnz,block.indptr.astype(np.int32),block.indices.astype(np.int32),block.data)
   self.nr=len(active)
  self.h.setOptionValue('time_limit',self.h.getRunTime()+self.lp_seconds)
  self.h.run();status=self.h.getModelStatus()
  if status!=self.highspy.HighsModelStatus.kOptimal:raise RuntimeError(str(status))
  solution=self.h.getSolution()
  return SimpleNamespace(x=np.array(solution.col_value),ineqlin=SimpleNamespace(marginals=-np.array(solution.row_dual)),fun=self.h.getObjectiveValue(),success=True)


def run(r,maxrounds=70,resume=False,master=None,lp_seconds=240.0):
 start=time.monotonic();N=19;m=10;c=LightCell(N,r,m);stem=f'cell_{N}_{r}_{m}'
 active=set();chosen=set();seed_weights={} 
 priorpath=HERE/f'cell_{N}_{r-1}_{m}.certificate.json'
 if not priorpath.exists():priorpath=REPO/f'gap_census/cell_{N}_{r-1}_{m}.certificate.json'
 if not priorpath.exists():priorpath=HERE/f'cell_{N}_{r-1}_{m}.bounds.certificate.json'
 if priorpath.exists():
  old=json.loads(priorpath.read_text());oldC=N-r+2
  for key,target in [('cover',chosen),('packing',active)]:
   seed_weights[key]={}
   for o,v in c.parse_id(lift_terms(old[key],oldC),key=='packing'):
    orbit=c.canonical(o);target.add(orbit);seed_weights[key][orbit]=seed_weights[key].get(orbit,F(0))+v
 for j in range(r):chosen.add(c.canonical((encode((0,)*j+(1,)*(c.L-j)),0)))
 # Generate the radius-one proposal pool without scanning the ambient prefix space.
 for a in c.library:
  for k in range(c.L):
   for digit in range(3):
    b=list(a);b[k]=digit;chosen.add(c.canonical((encode(b),0)))
  if c.C:chosen.add(c.canonical((encode(a),1)))
 if r<=14:
  from itertools import product
  binary=product((0,1),repeat=c.L)
 else:
  # At larger r, warm supports plus binary prefixes with <=2 changes seed row generation.
  binary=[]
  for x in [0,1]:
   for u in range(c.L+1):
    for v in range(u,c.L+1):binary.append((x,)*u+(1-x,)*(v-u)+(x,)*(c.L-v))
 for a in binary:
  i=encode(a);d=c.dmin_id(i);w=min(c.C,c.m-d)
  if w>=0:active.add(c.canonical((i,w)))
 if r>=16 and not resume and seed_weights:
  for key,target in [('cover',chosen),('packing',active)]:
   ranked=sorted(seed_weights[key].items(),key=lambda t:(-t[1],t[0]));threshold=sum(seed_weights[key].values(),F(0))*F(3,4);retained=F(0);selected=set()
   for orbit,mass in ranked:
    selected.add(orbit);retained+=mass
    if retained>=threshold:break
   # Retain the inexpensive geometric seeds and the significant transported masses.
   target.intersection_update({o for o in target if o not in seed_weights[key]})
   target.update(selected)
   if key=='cover':
    for j in range(r):target.add(c.canonical((encode((0,)*j+(1,)*(c.L-j)),0)))
   log(stage='significant_warm_seeds',r=r,side=key,selected=len(selected),total_orbits=len(ranked),retained_mass=str(retained))
 if resume:
  old=json.loads((HERE/(stem+'.candidate.json')).read_text())
  active.clear();chosen.clear()
  for j in range(r):chosen.add(c.canonical((encode((0,)*j+(1,)*(c.L-j)),0)))
  for key,target in [('cover',chosen),('packing',active)]:
   for o,v in c.parse_id(old[key],key=='packing'):target.add(c.canonical(o))
 if master:
  recovered=json.loads(Path(master).read_text());assert recovered['r']==r
  active={tuple(o) for o in recovered['active']};chosen={tuple(o) for o in recovered['chosen']}
 active=sorted(active);chosen=sorted(chosen);numeric=Numeric(c);warm=WarmLP(numeric,lp_seconds)
 for it in range(1,maxrounds+1):
  log(stage='LP',r=r,round=it,points=len(active),centres=len(chosen),elapsed=time.monotonic()-start)
  result=warm.solve(active,chosen)
  log(stage='numeric_proposal',r=r,value_display_only=float(result.fun),cover=int(np.count_nonzero(result.x>1e-8)),packing=int(np.count_nonzero(-result.ineqlin.marginals>1e-8)),elapsed=time.monotonic()-start)
  cover,packing=reconstruct(c,active,chosen,result)
  if it==1 and not master:
   supported={c.canonical(o) for o,v in packing}
   drop=[i for i,o in enumerate(active) if o not in supported]
   if len(drop)>len(supported)//2:
    warm.h.deleteRows(len(drop),np.array(drop,dtype=np.int32))
    active[:]=[o for o in active if o in supported];warm.nr=len(active)
    log(stage='removed_unused_seed_rows',r=r,removed=len(drop),retained=len(active))
   essential={c.canonical(o) for o,v in cover}
   for j in range(r):essential.add(c.canonical((encode((0,)*j+(1,)*(c.L-j)),0)))
   remove_columns=[j for j,o in enumerate(chosen) if o not in essential]
   if len(remove_columns)>len(essential)//2:
    warm.h.deleteCols(len(remove_columns),np.array(remove_columns,dtype=np.int32))
    chosen[:]=[o for o in chosen if o in essential];warm.nc=len(chosen)
    log(stage='removed_unused_seed_columns',r=r,removed=len(remove_columns),retained=len(chosen))
  cert=dict(N=N,r=r,m=m,cover=c.terms(cover),packing=c.terms(packing),normalization='total mass, uniform within each orbit',orientation='canonical 0^j 1^(r-1-j) 0^(N-r+1)')
  (HERE/(stem+'.candidate.json')).write_text(json.dumps(cert,indent=2)+'\n')
  reports={};newcounts={}
  for key,weights,ispacking,target in [('cover',cover,False,active),('packing',packing,True,chosen)]:
   report,poss=sweep(c,weights,ispacking,stem+f'.round_{it:02d}.{key}',keep=1024)
   reports[key]=report
   new=[];present=set(target)
   for val,o in poss:
    bad=val>1 if ispacking else val<1
    o=c.canonical(o)
    if bad and o not in present:new.append(o);present.add(o)
    if len(new)>=(512 if r>=16 else 160):break
   target.extend(new);newcounts[key]=len(new)
   log(stage='exact_separation',r=r,round=it,side=key,**report,added=len(new),elapsed=time.monotonic()-start)
  if all(v['passes'] for v in reports.values()):
   v=F(reports['cover']['total_mass']);assert v==F(reports['packing']['total_mass'])
   final=dict(N=N,r=r,m=m,status='VERIFIED-COMPUTATIONALLY',matched_exactly=True,tau_star=str(v),lower_K=math.ceil(v),upper_K=r,integrality_gap_exhibited=(r-1<v<r),**reports)
   if r-1<v<=r:final['K']=r
   (HERE/(stem+'.certificate.json')).write_text(json.dumps(cert,indent=2)+'\n')
   (HERE/(stem+'.verification.json')).write_text(json.dumps(final,indent=2)+'\n')
   log(stage='CERTIFIED',r=r,tau_star=str(v),K=final.get('K'),elapsed=time.monotonic()-start);return final
  (HERE/(stem+'.master.json')).write_text(json.dumps(dict(N=N,r=r,m=m,active=active,chosen=chosen))+'\n')
  if sum(newcounts.values())==0:raise ArithmeticError('separation failed to add exact violations')
 raise RuntimeError('round budget reached')

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('r',type=int);parser.add_argument('--resume',action='store_true');parser.add_argument('--master');parser.add_argument('--lp-seconds',type=float,default=240.0);a=parser.parse_args();run(a.r,resume=a.resume,master=a.master,lp_seconds=a.lp_seconds)
