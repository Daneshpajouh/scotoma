#!/usr/bin/env python3
"""Reject a one-byte input mutation through the verifier's actual input gate."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from integrity import check_pins,require,sha
import json
import subprocess


def run():
    pins=check_pins()
    work=HERE/'corruption_control'
    work.mkdir(exist_ok=True)
    badroot=work/'changed_input'
    badroot.mkdir(exist_ok=True)
    target='inputs/corrected_only_alignments.part000.tsv'
    for rel in ['pins.json',*pins]:
        dest=badroot/rel
        dest.parent.mkdir(parents=True,exist_ok=True)
        if rel!=target and not dest.exists():
            dest.symlink_to(HERE/rel)
    raw=bytearray((HERE/target).read_bytes())
    # Change one A/C/G/T at a real DNA alignment column; keep length and TSV syntax.
    first_data=raw.index(b'\n')+1
    dna_start=first_data
    for _ in range(9):
        dna_start=raw.index(b'\t',dna_start)+1
    offset=dna_start
    original=raw[offset]
    raw[offset]=ord('C') if original!=ord('C') else ord('A')
    changed=badroot/target
    changed.write_bytes(raw)
    require(sum(a!=b for a,b in zip(raw,(HERE/target).read_bytes()))==1,'Mutation is not exactly one byte')
    command=['/usr/local/bin/python3','-I',str(HERE/'verify.py'),'--integrity-only']
    good=subprocess.run(command,capture_output=True,text=True,timeout=60)
    bad=subprocess.run(command+['--input-root',str(badroot)],capture_output=True,text=True,timeout=60)
    (work/'unmodified.log').write_text(good.stdout+good.stderr)
    (work/'changed_input.log').write_text(bad.stdout+bad.stderr)
    require(good.returncode==0,'Unmodified control rejected')
    require(bad.returncode!=0 and 'SHA-256 mismatch: '+target in bad.stdout+bad.stderr,
            'One-byte corruption was not rejected at its input seal')
    receipt=dict(status='PASS',input=target,changed_bytes=1,byte_offset_0based=offset,
                  original_byte=original,changed_byte=raw[offset],original_sha256=sha(HERE/target),
                  changed_sha256=sha(changed),unmodified_verifier_exit_code=good.returncode,
                  corrupted_verifier_exit_code=bad.returncode,rejected_before_scoring=True)
    (work/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt,sort_keys=True))


if __name__=='__main__':
    run()
