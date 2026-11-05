"""Audit equality certificates, transport maps, and the direction of every K bound."""
import sys,json,math,copy,hashlib
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent

def check_map(mapping):
 source_r=mapping['source_r'];target_r=mapping['target_r'];perm=mapping['source_to_target_coordinate'];swaps=set(mapping['swap_0_1_at_target_coordinates']);indices=mapping['pattern_indices']
 assert sorted(perm)==list(range(19)) and all(0<=j<19 for j in swaps)
 assert len(indices)==source_r and indices==list(range(indices[0],indices[0]+source_r))
 images=set()
 for j in range(source_r):
  src=[0]*j+[1]*(source_r-1-j)+[0]*(20-source_r);dest=[None]*19
  for source,target in enumerate(perm):dest[target]=src[source]
  for k in swaps:dest[k]=1-dest[k]
  images.add(tuple(dest))
 expected={tuple([0]*j+[1]*(target_r-1-j)+[0]*(20-target_r)) for j in indices}
 assert images==expected,'transport does not map the source library to the specified sublibrary'
 return True

def check_transport(row):
 spec=json.loads(Path(row['transport_certificate']).read_text());packing=spec['packing'];check_map(packing['isometry'])
 base=Path(packing['base_certificate']);assert hashlib.sha256(base.read_bytes()).hexdigest()==packing['base_sha256']
 cert=json.loads(base.read_text());assert sum((F(t['mass']) for t in cert['packing']),F(0))==F(packing['mass'])==F(row['tau_lower'])
 cover_mass=F(0);covered=[]
 for block in spec['cover']['blocks']:
  check_map(block['isometry']);covered+=block['isometry']['pattern_indices'];mass=F(block['mass']);cover_mass+=mass
  if 'base_sha256' in block:
   base=Path(block['base_certificate']);assert hashlib.sha256(base.read_bytes()).hexdigest()==block['base_sha256']
   assert sum((F(t['mass']) for t in json.loads(base.read_text())['cover']),F(0))==mass
  else:assert mass==block['isometry']['source_r']
 assert sorted(covered)==list(range(row['r']))
 assert cover_mass==F(spec['cover']['mass'])==F(row['tau_upper'])
 return True

def check_row(row):
 check_transport(row)
 lo=F(row['tau_lower']);hi=F(row['tau_upper']);r=row['r']
 assert lo<=hi<r-1
 assert row['K_lower']==math.ceil(lo) and row['K_upper']==r
 if row.get('K_exact'):
  assert row['K_lower']==row['K_upper']==row['K_exact'], 'K equality has no forcing lower bound'
 if row['tau_star']:
  assert F(row['tau_star'])==lo==hi
  path=Path(row['certificates']);cert=json.loads(path.read_text());vp=sum((F(t['mass']) for t in cert['packing']),F(0));vc=sum((F(t['mass']) for t in cert['cover']),F(0))
  assert vp==vc==lo
  report=json.loads(path.with_name(path.name.replace('.certificate.json','.verification.json')).read_text())
  assert report['cover']['passes'] and report['packing']['passes'] and report['matched_exactly'] and F(report['tau_star'])==lo
  assert report['cover']['orbits_checked']==row['point_orbits'] and report['packing']['orbits_checked']==row['centre_orbits']
  for key in ['cover','packing']:
   if 'screen_flagged' in report[key]:assert report[key]['complete_flag_list'] and not report[key]['violation_count_is_lower_bound']
  assert hashlib.sha256(path.read_bytes()).hexdigest()==row['certificate_sha256']
 return True

if __name__=='__main__':
 rows=json.loads((HERE/'extended_exact_table.json').read_text())
 for row in rows:check_row(row)
 planted=copy.deepcopy(rows[0]);planted['K_exact']=planted['r']
 try:check_row(planted)
 except AssertionError as e:control=dict(name='upper_bound_used_to_claim_K_equal_r',expected='FAIL',observed='FAIL',reason=str(e))
 else:raise AssertionError('wrong-direction K control was incorrectly accepted')
 badmap=copy.deepcopy(json.loads(Path(rows[0]['transport_certificate']).read_text())['packing']['isometry'])
 badmap['source_to_target_coordinate'][0],badmap['source_to_target_coordinate'][1]=badmap['source_to_target_coordinate'][1],badmap['source_to_target_coordinate'][0]
 try:check_map(badmap)
 except AssertionError as e:map_control=dict(name='wrong_transport_coordinate_permutation',expected='FAIL',observed='FAIL',reason=str(e))
 else:raise AssertionError('wrong transport incorrectly accepted')
 report=dict(status='VERIFIED-COMPUTATIONALLY',table_rows=len(rows),rows_passed=len(rows),control=control,transport_control=map_control)
 (HERE/'table_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
