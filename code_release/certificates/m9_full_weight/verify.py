"""Replay the exact mass-twenty certificate over all ternary centre classes."""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metadata-only', action='store_true',
                        help='check support, regional masses and deposited record without replaying the transform')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    assert 1 <= args.threads <= 16
    for line in (HERE / 'SHA256SUMS').read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == expected, name
    with tempfile.TemporaryDirectory(prefix='m9_full_weight_') as tmp:
        tmp = Path(tmp)
        exact.ROOT = tmp
        input_path = exact.prepare(HERE / 'packing_19_20_9.json')
        lines = input_path.read_text().splitlines()
        header = lines[0].split()
        header[-1] = str(args.threads)
        lines[0] = ' '.join(header)
        input_path.write_text('\n'.join(lines) + '\n')
        deposited = dict(line.split(maxsplit=1) for line in
                         (HERE / 'deposited_verification.txt').read_text().splitlines() if line.strip())
        assert int(deposited['classes']) == 3**19
        assert int(deposited['violations']) == int(deposited['errors']) == 0
        assert Fraction(deposited['nonpattern_max']) <= 1
        record = {'status': 'VERIFIED-COMPUTATIONALLY', 'support_and_regional_masses': 'PASS',
                  'mass': '20', 'centre_classes': 3**19,
                  'transform_replayed': False, 'deposited_record_checked': True}
        if not args.metadata_only:
            binary = tmp / 'tensor_crt'
            subprocess.run([os.environ.get('CXX', 'c++'), '-O3', '-std=c++17', '-pthread',
                            str(HERE / 'tensor_crt.cpp'), '-o', str(binary)], check=True)
            output_path = tmp / 'replay.out'
            subprocess.run([str(binary), str(input_path), str(output_path)], check=True)
            result = dict(line.split(maxsplit=1) for line in output_path.read_text().splitlines())
            assert int(result['classes']) == 3**19
            assert int(result['violations']) == int(result['errors']) == 0
            assert result['nonpattern_max'] == deposited['nonpattern_max']
            record.update(transform_replayed=True, maximum_load='1',
                          nonpattern_max=result['nonpattern_max'])
        print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
