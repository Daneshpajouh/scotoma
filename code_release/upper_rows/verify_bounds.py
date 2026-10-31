#!/usr/local/bin/python3
"""Verify a two-sided interval certificate without claiming equality or exact K."""
from extend import *

path=Path(sys.argv[1]);cert=json.loads(path.read_text());cell=LightCell(cert['N'],cert['r'],cert['m']);reports={}
for key,packing in [('cover',False),('packing',True)]:
 weights=cell.parse_id(cert[key],packing);report,_=sweep(cell,weights,packing,path.stem+'.recheck.'+key,keep=0);assert report['passes'];reports[key]=report
lo=F(reports['packing']['total_mass']);hi=F(reports['cover']['total_mass']);assert lo<=hi
result=dict(status='VERIFIED-COMPUTATIONALLY',N=cell.N,r=cell.r,m=cell.m,tau_lower=str(lo),tau_upper=str(hi),K_lower=math.ceil(lo),K_upper=cell.r,**reports)
(HERE/(path.stem+'.recheck.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],r=cell.r,K_lower=math.ceil(lo),K_upper=cell.r,matched_exactly=lo==hi)))
