#!/usr/bin/env python3
"""Make exactly one byte substitution in a deposited input; require rejection."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent


def main():
    out=HERE/'corruption_control'
    out.mkdir(exist_ok=True)
    tests=[]
    for name,option,offset in [
        ('unreachable_sequences.tsv','--target-input',None),
        ('CRISPR_Net_CIRCLE_elevation_SITE_weights.h5','--weight-input',100000),
    ]:
        original=(HERE/'inputs'/name).read_bytes()
        if offset is None:
            # Alter a sequence base while keeping the TSV syntactically valid.
            offset=original.index(b'\n')+1
            offset=original.index(b'\t',offset)+1
        altered=bytearray(original)
        if name.endswith('.tsv'):
            altered[offset]=dict(zip(b'ACGT',b'CGTA'))[original[offset]]
        else:
            altered[offset]=original[offset]^1
        assert len(original)==len(altered)
        assert sum(a!=b for a,b in zip(original,altered))==1
        target=out/(name+'.one_byte_changed')
        target.write_bytes(altered)
        cmd=[sys.executable,'-I','-B',str(HERE/'verify.py'),option,str(target)]
        result=subprocess.run(cmd,capture_output=True,text=True)
        expected='REJECT: SHA-256 mismatch: inputs/'+name
        assert result.returncode==1 and expected in result.stderr,(cmd,result.returncode,result.stderr)
        record=dict(input=name,byte_offset_0based=offset,original_byte=original[offset],
                    altered_byte=altered[offset],changed_bytes=1,length=len(original),
                    original_sha256=hashlib.sha256(original).hexdigest(),
                    altered_sha256=hashlib.sha256(altered).hexdigest(),
                    verifier_exit_status=result.returncode,verifier_stderr=result.stderr.strip(),
                    outcome='REJECTED')
        tests.append(record)
        print(f'PASS: one-byte mutation rejected: {name}',flush=True)
    (out/'receipt.json').write_text(json.dumps(dict(status='PASS',tests=tests),indent=2)+'\n')


if __name__=='__main__':main()
