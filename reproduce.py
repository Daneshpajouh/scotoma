"""Replay deposited evidence in a disposable copy; the deposit stays unchanged."""
from pathlib import Path
import argparse
import hashlib
import os
import json
import resource
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
PY = [sys.executable, '-B']
CXX = os.environ.get('CXX', 'c++')
BOOST = os.environ.get('BOOST_INCLUDE')
if not BOOST and Path('/opt/homebrew/include/boost').is_dir():
    BOOST = '/opt/homebrew/include'
BOOST_FLAGS = ['-I' + BOOST] if BOOST else []
COMMANDS = {
    'm9_full_weight': [('code_release', PY + ['-I', 'certificates/m9_full_weight/verify.py'])],
    'm8_closure': [('code_release', PY + ['-I', 'certificates/m8_closure/verify.py', 'verify', '--workers', '4'])],
    'm9_closure': [('code_release', PY + ['-I', 'certificates/m9_closure/verify.py', 'verify'])],
    'm10_closure': [('code_release', PY + ['-I', 'certificates/m10_closure/verify.py', 'verify'])],
    'm10_r11': [('code_release', PY + ['-I', 'certificates/m10_r11/verify.py', 'verify'])],
    'm10_r11_kernel': [('code_release', PY + ['-I', 'certificates/m10_r11/verify_kernel_certs.py'])],
    'infinite_families': [('code_release', PY + ['-I', 'certificates/infinite_families/verify.py', 'verify'])],
    'thresholds': [('code_release', PY + ['thresholds/verify_all.py', '--full-grid'])],
    'private_sphere': [('code_release', PY + [f'private_sphere/{s}'] + (['--controls'] if s == 'verify_w.py' else []))
                       for s in ['verify_w.py', 'verify_geometry.py', 'verify_companion.py', 'finite_cells.py', 'run_negative_controls.py']],
    'selector': [('code_release', PY + ['selector/proof_checks.py']),
                 ('code_release', PY + ['selector/verify_j2.py', 'verify'])],
    'five_run': [('code_release', PY + ['five_run/verify_k3.py']),
                 ('code_release', PY + ['five_run/independent_literal_check.py'])],
    'cover_families': [('code_release', PY + ['cover_families/verify_j1.py', '--algebra-only'])],
    'monotonicity': [('code_release', PY + ['monotonicity/verify_finite.py', '--table', 'monotonicity/tau_values.csv']),
                     ('code_release', PY + ['monotonicity/audit_map.py']),
                     ('code_release', PY + ['monotonicity/run_controls.py'])],
    'gap_census': [('code_release', PY + ['gap_census/final_verify.py'])],
    'insertion_nofloor': [('code_release', PY + ['-I', 'certificates/insertion_nofloor/verify.py', 'verify'])],
    'pooled': [('code_release', PY + ['-I', 'certificates/pooled_m4_1629/my_private_check.py']),
               ('code_release', PY + ['-I', 'certificates/pooled_m4_1629/space_verify.py', 'data/producers_2026-10-20/bulge_census_full_idx_rna.json'])],
    'scope': [('code_release', PY + ['-I', 'certificates/paper_B_scope/pooled_counterexample.py']),
              ('code_release', PY + ['-I', 'certificates/paper_B_scope/verify_scope.py']),
              ('code_release', PY + ['-I', 'certificates/paper_B_scope/soundness.py'])],
    'cellular': [('code_release/cellular', PY + ['independent_recheck.py'])],
    'multibulge': [('code_release/multibulge', PY + ['independent_counterexample_check.py'])],
    'empirical_gates': [('code_release', PY + [s]) for s in ['verify_minimality.py', 'verify_guide_screen.py',
                         'verify_missed_actives.py', 'make_scheme_lemmas.py', 'make_bulge6_table.py',
                         'verify_distal_bound.py', 'verify_F_independent.py']],
    'r20_m8': [('code_release', ['node', 'certificates/r20_m8/dual20.js']),
               ('code_release', PY + ['-I', 'certificates/r20_m8/independent_check.py'])],
    'floor_rows': [('code_release', PY + [s]) for s in ['verify_floor_certificates.py', 'verify_floor_new_results.py', 'verify_ceiling_reach.py']],
    'ceiling': [('code_release', PY + ['gap_census/checker.py', 'ceiling/' + f]) for f in
                ['cover_19_11_7.converted.json', 'cover_19_12_5.converted.json']],
    'enumerator': [('.', PY + ['-m', 'pytest', '-q', 'tests'])],
    'dials': [('code_release', PY + [s]) for s in ['make_assembly_table.py', 'make_compounding_table.py',
               'make_dials_figure.py']],
    'maize_loci': [('code_release', PY + ['verify_maize_exclusive.py'])],
    'tool_evidence': [('code_release/tool_audit', PY + ['-c',
        'import hashlib,json,pathlib; d=json.loads(pathlib.Path("artifact_sha256.json").read_text()); '
        'assert len(d)==5; assert all(hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()==h for p,h in d.items()); '
        'print("FIVE_DEPOSITED_AUDIT_HASHES_PASS; native pipelines were not rerun")'])],
    'structure': [('code_release', PY + ['verify_path_lemma_wide.py']),
                  ('code_release', PY + ['verify_burst_floor.py', 'data/panel_match.json']),
                  ('code_release', PY + ['verify_transversal_lengths.py'])],
    'validation_budget': [('code_release', PY + ['test_m0_search.py']),
                          ('code_release', PY + ['boundary_exact_test.py']),
                          ('code_release', PY + ['offtarget_budget.py', '--n', '28175', '--k', '1', '10', '25'])],
    'upper_rows': [('code_release/upper_rows', [CXX, '-O3', '-std=c++17', '-pthread', *BOOST_FLAGS,
                     s + '.cpp', '-o', s]) for s in ['tensor_general', 'tensor_screen', 'refine']]
                  + [('code_release/upper_rows', PY + ['verify_one.py', f'cell_19_{r}_10.certificate.json']) for r in range(11,16)]
                  + [('code_release/upper_rows', PY + ['verify_bounds.py', f'cell_19_{r}_10.bounds.certificate.json']) for r in [16,17]]
                  + [('code_release/upper_rows', PY + [s]) for s in ['make_table.py', 'audit_table.py']],
}


def integrity():
    lines = (ROOT / 'SHA256SUMS').read_text().splitlines()
    if len(lines) < 100:
        raise ValueError('incomplete deposit inventory')
    for line in lines:
        expected, name = line.split(maxsplit=1)
        p = ROOT / name
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest() != expected:
            raise ValueError('missing or changed deposit file: ' + name)
    return len(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('package', choices=['integrity', 'all', *COMMANDS])
    args = ap.parse_args()
    print('INTEGRITY_PASS', integrity(), flush=True)
    if args.package == 'integrity':
        return
    env = os.environ.copy()
    for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        env[key] = '1'
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    for key in ['CHANGESEQ_CSV', 'DATASETS_ZIP', 'PANEL_FILE']:
        env.pop(key, None)
    for key in list(env):
        if key.startswith(('CERT_', 'VF_', 'FLOOR_', 'CEILING_', 'COV_')) or key in {
            'PYTHONPATH', 'PYTHONHOME', 'PYTHONSTARTUP', 'SITES129', 'RNA1_IDX', 'MINCERT',
            'BULGE_IDX_DNA', 'MM3_PATH', 'MM4_PATH', 'MM5_PATH', 'MM6_PATH'}:
            env.pop(key, None)
    packages = list(COMMANDS) if args.package == 'all' else [args.package]
    with tempfile.TemporaryDirectory(prefix='scotoma_replay_') as temp:
        work = Path(temp) / 'deposit'
        shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns('__pycache__', '.git'))
        env['MPLCONFIGDIR'] = str(Path(temp) / 'matplotlib')
        for package in packages:
            print('PACKAGE_BEGIN', package, flush=True)
            for folder, command in COMMANDS[package]:
                print('COMMAND', folder, ' '.join(command), flush=True)
                before = resource.getrusage(resource.RUSAGE_CHILDREN)
                started = time.monotonic()
                result = subprocess.run(command, cwd=work / folder, env=env)
                usage = resource.getrusage(resource.RUSAGE_CHILDREN)
                print('REPLAY_RESOURCES ' + json.dumps({'package': package, 'command': command,
                      'exit_code': result.returncode, 'wall_seconds': time.monotonic() - started,
                      'peak_child_rss_bytes': int(usage.ru_maxrss * (1 if sys.platform == 'darwin' else 1024)),
                      'user_seconds': usage.ru_utime - before.ru_utime,
                      'system_seconds': usage.ru_stime - before.ru_stime}), flush=True)
                result.check_returncode()
            print('PACKAGE_PASS', package, flush=True)


if __name__ == '__main__':
    main()
