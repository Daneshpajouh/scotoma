#!/usr/bin/env python3
"""Recheck every completed certificate with the final checker artifact."""
from pathlib import Path
import hashlib
import json
from checker import verify

HERE = Path(__file__).resolve().parent
results = []
files = sorted(HERE.glob('cell_*.certificate.json'))
declared = json.loads((HERE/'exact_gaps.json').read_text())
assert len(declared) == 13, 'the thirteen-row gap inventory is incomplete'
expected = {row['certificate'] for row in declared} | {'cell_5_4_3.certificate.json'}
assert {p.name for p in files} == expected, 'gap certificate inventory is incomplete'
for path in files:
    check = verify(json.loads(path.read_text()))
    old_path = path.with_name(path.name.replace('.certificate.json', '.verification.json'))
    old = json.loads(old_path.read_text())
    current_fields = {k: v for k, v in check.items() if k != 'elapsed_seconds'}
    old_fields = {k: v for k, v in old.items() if k != 'elapsed_seconds'}
    if not check['passes'] or not check['matched_exactly'] or current_fields != old_fields:
        raise AssertionError(('final artifact disagrees', path.name, check, old))
    check['certificate'] = path.name
    check['certificate_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    (HERE/path.name.replace('.certificate.json', '.final_verification.json')).write_text(json.dumps(check, indent=2)+'\n')
    results.append(check)
    print(json.dumps({'certificate': path.name, 'status': 'VERIFIED-COMPUTATIONALLY',
                      'tau_star': check['tau_star'], 'passes': True}), flush=True)
summary = {'status': 'VERIFIED-COMPUTATIONALLY',
           'checker_sha256': hashlib.sha256((HERE/'checker.py').read_bytes()).hexdigest(),
           'completed_matching_pairs': len(results), 'results': results}
(HERE/'final_verification_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
