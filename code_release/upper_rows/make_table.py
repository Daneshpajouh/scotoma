"""Exact ten-cell table and compact isometric transport certificates."""
from extend import *
import csv,hashlib
HERE=Path(".")
REPO=Path("..")

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def orbit_counts(N,r,m):
 L=r-1;C=N-L;states={(0,0,0):(1,1)}
 for _ in range(L):
  nxt={}
  for (d,e,b),(count,words) in states.items():
   for dd,de,mult in [(1,-1,1),(0,1,1),(1,0,2)]:
    key=(d+dd,e+de,min(b,e+de));co,wo=nxt.get(key,(0,0));nxt[key]=(co+count,wo+mult*words)
  states=nxt
 hist={}
 for (d,e,b),(co,wo) in states.items():
  a,z=hist.get(d+b,(0,0));hist[d+b]=(a+co,z+wo)
 assert sum(x[0] for x in hist.values())==3**L
 assert sum(x[1] for x in hist.values())==4**L
 points=sum(co*max(0,min(C,m-d)+1) for d,(co,wo) in hist.items())
 uw=sum(wo*sum(math.comb(C,w)*3**w for w in range(max(0,min(C,m-d)+1))) for d,(co,wo) in hist.items())
 return dict(point_orbits=points,centre_orbits=3**L*(C+1),union_words=uw,ambient_words=4**N)

def transport(source_r,target_r,start=0):
 # Source active coordinates map to the varying coordinates of the target subchain.
 lb=source_r-1;Lt=target_r-1
 mapping=list(range(start,start+lb))+[j for j in range(19) if j not in range(start,start+lb)]
 swaps=[j for j in range(start+lb,Lt)]
 for i in range(source_r):
  src=[0]*i+[1]*(lb-i)+[0]*(19-lb);got=[None]*19
  for k,j in enumerate(mapping):got[j]=1-src[k] if j in swaps and src[k]<2 else src[k]
  expected=[0]*(start+i)+[1]*(Lt-start-i)+[0]*(19-Lt)
  assert got==expected
 return dict(source_r=source_r,target_r=target_r,pattern_indices=list(range(start,start+source_r)),source_to_target_coordinate=mapping,swap_0_1_at_target_coordinates=swaps,coordinate_indexing='zero based')

values={1:F(1),2:F(2),3:F(3),4:F(4)};lower_values=dict(values);paths={};exact_cells=set()
for r in range(5,21):
 path=HERE/f'cell_19_{r}_10.certificate.json'
 if not path.exists():path=HERE/f'input_certificates/cell_19_{r}_10.certificate.json'
 if not path.exists():path=REPO/f'gap_census/cell_19_{r}_10.certificate.json'
 if path.exists():
  cert=json.loads(path.read_text());v=sum((F(t['mass']) for t in cert['cover']),F(0));assert v==sum((F(t['mass']) for t in cert['packing']),F(0));values[r]=v;lower_values[r]=v;paths[r]=path;exact_cells.add(r)
 else:
  path=HERE/f'cell_19_{r}_10.bounds.certificate.json'
  if path.exists():
   cert=json.loads(path.read_text());hi=sum((F(t['mass']) for t in cert['cover']),F(0));lo=sum((F(t['mass']) for t in cert['packing']),F(0));assert lo<=hi
   report=json.loads(path.with_name(path.name.replace('.certificate.json','.verification.json')).read_text());assert report['cover']['passes'] and report['packing']['passes']
   values[r]=hi;lower_values[r]=lo;paths[r]=path
# A partition into consecutive sublibraries gives a fractional cover by isometric transport.
upper={0:(F(0),[])}
for r in range(1,21):
 upper[r]=min((upper[r-b][0]+v,upper[r-b][1]+[b]) for b,v in values.items() if b<=r)
rows=[]
for r in range(11,21):
 lower_r=max((b for b in lower_values if b<=r),key=lambda b:(lower_values[b],b));lower=lower_values[lower_r]
 up,partition=upper[r];assert lower<=up<r-1
 spec=dict(N=19,r=r,m=10,status='PROVED',arithmetic='Fraction',packing=dict(base_certificate=str(paths[lower_r]),base_sha256=sha(paths[lower_r]),isometry=transport(lower_r,r),mass=str(lower)),cover=dict(mass=str(up),blocks=[]))
 start=0
 for b in partition:
  block=dict(mass=str(values[b]),isometry=transport(b,r,start))
  if b in paths:block.update(base_certificate=str(paths[b]),base_sha256=sha(paths[b]))
  else:block['base_certificate']='the b canonical library balls, each of mass one'
  spec['cover']['blocks'].append(block);start+=b
 assert start==r
 specpath=HERE/f'cell_19_{r}_10.transport_bounds.json';specpath.write_text(json.dumps(spec,indent=2)+'\n')
 exact=r in exact_cells
 row=dict(N=19,r=r,m=10,status='VERIFIED-COMPUTATIONALLY' if exact or (HERE/f'cell_19_{r}_10.bounds.certificate.json').exists() else 'PROVED',resolution='EXACT_TAU_K_BOUNDED' if exact else 'BOUNDS_ONLY',tau_star=str(values[r]) if exact else '',tau_lower=str(lower),tau_upper=str(up),K_lower=math.ceil(lower),K_upper=r,K_exact='',requested_packing_above_r_minus_1='REFUTED',gap_at_K_equal_r='NOT_CLAIMED',upper_deficit=str(r-up),fractional_cover_below_r_minus_1=True,cover_partition=partition,lower_base_r=lower_r,certificates=(str(paths[r]) if exact else str(specpath)),transport_certificate=str(specpath),**orbit_counts(19,r,10))
 row['local_LP_attempted']=(HERE/f'r{r}.log').exists()
 row['verification_route']='matched full-orbit cover and packing' if exact else 'certified interval and isometric transport'
 row['own_full_orbit_certificate']=exact or (HERE/f'cell_19_{r}_10.bounds.certificate.json').exists()
 ldec=(lower.numerator*10**6)//lower.denominator;udec=-((-up.numerator*10**6)//up.denominator)
 assert F(ldec,10**6)<=lower<=up<=F(udec,10**6)
 row['decimal_bound_lower']=f'{ldec//10**6}.{ldec%10**6:06d}';row['decimal_bound_upper']=f'{udec//10**6}.{udec%10**6:06d}'
 log_path=HERE/f'r{r}.log'
 row['reason']='matching exact certificates; K bounded' if exact else 'proposal LP time limit; certified bounds' if log_path.exists() and 'HighsModelStatus.kTimeLimit' in log_path.read_text() else 'transported bounds; no direct LP run'
 if exact:
  v=values[r];assert lower==up==v;z=(v.numerator*10**6)//v.denominator
  assert F(z,10**6)<v<F(z+1,10**6)
  row.update(decimal_interval_lower=f'{z//10**6}.{z%10**6:06d}',decimal_interval_upper=f'{(z+1)//10**6}.{(z+1)%10**6:06d}',certificate_sha256=sha(paths[r]))
 rows.append(row)
(HERE/'extended_exact_table.json').write_text(json.dumps(rows,indent=2)+'\n')
keys=list(dict.fromkeys(k for row in rows for k in row))
with (HERE/'extended_exact_table.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=keys);wr.writeheader();wr.writerows(rows)
print(json.dumps([dict(r=v['r'],resolution=v['resolution'],K_lower=v['K_lower'],K_upper=v['K_upper'],tau_interval=[v.get('decimal_interval_lower'),v.get('decimal_interval_upper')],cover_partition=v['cover_partition']) for v in rows],indent=2))
