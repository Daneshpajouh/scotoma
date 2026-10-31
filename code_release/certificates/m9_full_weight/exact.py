"""Integer transform certificate driver; all certification arithmetic is exact."""
from pathlib import Path
from fractions import Fraction
from math import comb,lcm,isqrt
import json,sys,subprocess,time
ROOT=Path(__file__).resolve().parent

def prepare(path):
 c=json.loads(Path(path).read_text());n,r,m,D=c['N'],c['r'],c['m'],c['D'];L=r-1;C=n-L
 assert 0<=L<=n and isinstance(D,int) and D>0
 sums=[0]*r;rows=[];Q=1
 for idx,w,v in c['support']:
  assert 0<=idx<3**L and 0<=w<=C and isinstance(v,int) and v>=0
  p=[idx//3**(L-1-k)%3 for k in range(L)]
  ds=[sum(a!=(0 if k<j else 1) for k,a in enumerate(p))+w for j in range(r)]
  assert min(ds)==m and ds.count(m)==1
  sums[ds.index(m)]+=v
  den=2**p.count(2)*comb(C,w)*3**w;Q=lcm(Q,den);rows.append((idx,w,v,den))
 assert sums==[D]*r
 vol=sum(comb(n,k)*3**k for k in range(m+1))
 # Safe absolute bound for every intermediate of either transform.
 assert r*D*Q*(isqrt(vol)+1)*4**n < 2**127
 inp=ROOT/(Path(path).stem+'.in')
 inp.write_text(f'{n} {r} {m} {D} {Q} {len(rows)} 16\n'+''.join(f'{i} {w} {v*(Q//d)}\n' for i,w,v,d in rows))
 return inp

def run(path,direct=False):
 inp=prepare(path);out=inp.with_suffix('.out')
 nclasses=3**(json.loads(Path(path).read_text())['r']-1)*(21-json.loads(Path(path).read_text())['r'])
 kernel='tensor_crt' if nclasses>=100000000 else 'tensor_exact'
 subprocess.run([str(ROOT/kernel),str(inp),str(out)]+(['direct'] if direct else []),check=True)
 ans={p[0]:' '.join(p[1:]) for p in map(str.split,out.read_text().splitlines())}
 ans['badcenters']=[list(map(int,p[1:3])) for p in map(str.split,out.read_text().splitlines()) if p[0]=='badcenter'];ans.pop('badcenter',None)
 ans['file']=str(path);ans['status']='COMPUTED_EXACT';print(json.dumps({k:v for k,v in ans.items() if k!='badcenters'}),flush=True)
 return ans
if __name__=='__main__':
 for p in sys.argv[1:]:
  if p!='--direct':run(p,'--direct' in sys.argv)
