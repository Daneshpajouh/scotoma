"""Exact small-denominator interval screening followed by exact direct refinement."""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import json,subprocess,time
HERE=Path(__file__).resolve().parent

def decode(i,L):
 a=[0]*L
 for k in range(L-1,-1,-1):a[k]=i%3;i//=3
 return a

def refine(cell,weights,queries,Q,stem):
 D=lcm(*(mass.denominator for _,mass in weights));S=D*Q;H=cell.L+1;C=cell.C
 maxnx={}
 for (i,w),mass in weights:maxnx[w]=max(maxnx.get(w,0),decode(i,cell.L).count(2))
 inp=HERE/(stem+'.refine.in');out=HERE/(stem+'.refine.out')
 with inp.open('w') as f:
  f.write(f'{cell.N} {cell.r} {cell.m} {int(refine.packing)} 8 {S} {len(weights)} {len(queries)}\n')
  for (i,w),mass in weights:
   a=decode(i,cell.L);masks=[sum(1<<k for k,x in enumerate(a) if x==v) for v in [0,1,2]];v=mass*D;assert v.denominator==1
   f.write(' '.join(map(str,masks))+f' {w} {v.numerator}\n')
  vals=[]
  for w in range(C+1):
   for p in range(C+1):
    for d in range(H):
     for h in range(H):
      q=Q*cell.probability(w,p,d,h) if w in maxnx and h<=maxnx[w] and d+h<=cell.L else F(0)
      assert q.denominator==1 and 0<=q<2**64
      vals.append(str(q.numerator))
  f.write(' '.join(vals)+'\n')
  for i,w in queries:f.write(f'{i} {w}\n')
 subprocess.run([str(HERE/'refine'),str(inp),str(out)],check=True)
 fields={};candidates=[]
 for line in out.read_text().splitlines():
  k,*v=line.split()
  if k=='candidate':candidates.append((F(int(v[2]),S),(int(v[0]),int(v[1]))))
  else:fields[k]=' '.join(v)
 assert int(fields['scale'])==S and int(fields['checked'])==len(queries)
 e=F(int(fields['extremum']),S)
 return dict(extremum=str(e),at=list(map(int,fields['arg'].split())),violations=int(fields['violations']),checked=len(queries)),candidates

def sweep(cell,weights,packing,stem,keep=512):
 start=time.monotonic();Q=lcm(*(cell.size_id(o) for o,mass in weights));MOD=(2**31-1)*(2**31-19);D=10**9
 for o,mass in weights:assert mass>=0 and (not packing or cell.in_union(o))
 while True:
  ints=[((o),(mass*D).numerator//(mass*D).denominator if packing else -(-(mass*D).numerator//(mass*D).denominator)) for o,mass in weights]
  bound=Q*sum(v for o,v in ints)
  if bound<MOD:break
  D//=10
  if D<=len(weights):raise ArithmeticError('screen denominator is too small')
 S=D*Q;epsilon=F(len(weights),D);threshold=S-len(weights)*Q if packing else S+len(weights)*Q
 assert 0<threshold<MOD and all(Q%cell.size_id(o)==0 for o,v in ints)
 actualerror=sum((abs(mass-F(v,D)) for (o,mass),(_,v) in zip(weights,ints)),F(0));assert actualerror<=epsilon
 inp=HERE/(stem+'.screen.in');out=HERE/(stem+'.screen.out')
 nz=[(o,v) for o,v in ints if v]
 inp.write_text(f'{cell.N} {cell.r} {cell.m} 8 {int(packing)} {S} {bound} {threshold} 250000 {len(nz)}\n'+''.join(f'{i} {w} {v*(Q//cell.size_id((i,w)))}\n' for (i,w),v in nz))
 subprocess.run([str(HERE/'tensor_screen'),str(inp),str(out)],check=True)
 fields={};allflags=[];priority=[]
 for line in out.read_text().splitlines():
  k,*v=line.split()
  if k=='flag':allflags.append((int(v[0]),int(v[1])))
  elif k=='priority':priority.append((int(v[0]),int(v[1])))
  else:fields[k]=' '.join(v)
 assert int(fields['errors'])==0
 complete=bool(int(fields['complete']));nflag=int(fields['flagged']);refine.packing=packing
 if priority:
  exact,candidates=refine(cell,weights,priority,Q,stem+'.priority')
 else:exact=dict(extremum='0' if packing else '2',at=None,violations=0,checked=0);candidates=[]
 refined=exact['checked']
 if exact['violations']==0 and complete:
  assert len(allflags)==nflag
  if allflags:
   exact,candidates=refine(cell,weights,allflags,Q,stem+'.all_flags');refined+=exact['checked']
  accepted=exact['violations']==0
 else:accepted=False
 if accepted:
  e=F(exact['extremum'])
  # All other packing loads are <1, or all other coverages are >1.
  # The generated LP supports supply a tight constraint in the refined set.
  if nflag and (e>=1 if packing else e<=1):extremum=str(e)
  else:extremum='at most 1' if packing else 'at least 1'
 else:extremum=exact['extremum']
 report=dict(passes=accepted,extremum=extremum,at=exact['at'],total_mass=str(sum((v for o,v in weights),F(0))),orbits_checked=int(fields['checked']),total_centre_orbits=int(fields['classes']),violations=exact['violations'],violation_count_is_lower_bound=not accepted,screen_flagged=nflag,complete_flag_list=complete,direct_evaluations=refined,screen_error_bound=str(epsilon),rounding_denominator=D,crt_absolute_bound=str(bound),crt_modulus=str(MOD),seconds=time.monotonic()-start,arithmetic='exact integer CRT screening; Fraction error bounds; arbitrary precision direct refinement')
 (HERE/(stem+'.json')).write_text(json.dumps(report,indent=2)+'\n')
 return report,candidates[:keep]
