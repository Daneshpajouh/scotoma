"""Replay a named external-data analysis after verifying the reader's own input."""
from pathlib import Path
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from reproduce import integrity

INPUTS = {
    'changeseq': (878760836, '68f3e1d7eb9e9e18e3ddc859752be8d156378a3e681b0fff329fe679366adc11'),
    'guideseq': (381626484, 'a3398cba2ef452ec37ebebda38a44aa62615c66165a43699691ca210f467bcd2'),
    'newguideseq': (207094458, 'ac9ca55a384a9af5ba37d667766e6723a7c211b5fa1d54875e450db23a13a23b'),
    'cfd-table': (11199, 'eba85c300e484207e5cfe67b8717146cb49e2ee24a4b1b2a4af4960baaa891a5'),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('analysis', choices=INPUTS)
    ap.add_argument('input', type=Path)
    args = ap.parse_args()
    if not args.input.is_file():
        ap.error('external input is missing; see EXTERNAL_INPUTS.md')
    expected_size, expected_hash = INPUTS[args.analysis]
    actual = args.input.resolve()
    if actual.stat().st_size != expected_size:
        raise SystemExit('INPUT_REJECTED: incorrect byte count')
    digest = hashlib.sha256()
    with actual.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    if digest.hexdigest() != expected_hash:
        raise SystemExit('INPUT_REJECTED: incorrect SHA-256')
    print('INTEGRITY_PASS', integrity(), flush=True)
    py = [sys.executable, '-B']
    if args.analysis == 'cfd-table':
        commands = [(py + ['code_release/verify_cfd_gap_table.py', str(actual)], False)]
    else:
        commands = [(py + ['verification/bulge_shares/share_strata.py', str(actual)], False)]
    if args.analysis == 'changeseq':
        commands += [(py + ['code_release/make_sites129.py', '--check', 'code_release/data/sites129.json'], True),
                     (py + ['code_release/verify_missed_activity.py'], True),
                     (py + ['code_release/validate_terminal_law.py'], True),
                     (py + ['code_release/make_rna_position_series.py'], True),
                     (py + ['code_release/read_depth/matched_bin.py'], False),
                     (py + ['verification/cfd_gap_asymmetry/cfd_gap_arms.py', str(actual)], False)]
    with tempfile.TemporaryDirectory(prefix='scotoma_external_') as directory:
        work = Path(directory) / 'deposit'
        shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns('__pycache__', '.git'))
        env = os.environ.copy()
        env['PYTHONDONTWRITEBYTECODE'] = '1'
        env['CHANGESEQ_CSV'] = str(actual)
        for command, stdin in commands:
            print('COMMAND', ' '.join(command), flush=True)
            if stdin:
                with actual.open('rb') as stream:
                    subprocess.run(command, cwd=work, env=env, stdin=stream, check=True)
            else:
                subprocess.run(command, cwd=work, env=env, check=True)
    print('EXTERNAL_ANALYSIS_PASS', args.analysis, flush=True)


if __name__ == '__main__':
    main()
