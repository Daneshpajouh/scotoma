"""Require all six invalid certificates to fail as actual standalone processes."""
from pathlib import Path
import subprocess
import sys
import json
ROOT=Path(__file__).resolve().parent
rows=[]
for control in ['zero-local','empty-polynomial','missing-local','lower-path','zero-tail','wrong-p']:
    r=subprocess.run([sys.executable,'-B',str(ROOT/'verify_w.py'),'--control',control],cwd=ROOT,text=True,capture_output=True)
    assert r.returncode==1,(control,r.returncode,r.stdout,r.stderr)
    rows.append({'control':control,'exit_code':r.returncode,'stderr':r.stderr.strip(),'outcome':'FAILED_AS_REQUIRED'})
    print(control,'exit',r.returncode,r.stderr.strip(),flush=True)
(ROOT/'control_exit_codes.json').write_text(json.dumps(rows,indent=2)+'\n')
