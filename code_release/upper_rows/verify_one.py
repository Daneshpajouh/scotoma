#!/usr/local/bin/python3
"""Verify a final original-orbit certificate without running an LP solver."""
from extend import *

if __name__=='__main__':
 path=Path(sys.argv[1]);cert=json.loads(path.read_text());cell=LightCell(cert['N'],cert['r'],cert['m']);reports={}
 for key,packing in [('cover',False),('packing',True)]:
  weights=cell.parse_id(cert[key],packing);report,_=sweep(cell,weights,packing,path.stem+'.recheck.'+key,keep=0);assert report['passes'];reports[key]=report
 upper=F(reports['cover']['total_mass']);lower=F(reports['packing']['total_mass']);assert upper==lower
 result=dict(status='VERIFIED-COMPUTATIONALLY',N=cell.N,r=cell.r,m=cell.m,tau_star=str(lower),matched_exactly=True,K_lower=math.ceil(lower),K_upper=cell.r,**reports)
 out=HERE/(path.stem+'.recheck.json');out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],r=cell.r,matched_exactly=True,K_lower=result['K_lower'],K_upper=result['K_upper'])))
