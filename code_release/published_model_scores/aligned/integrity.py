"""Input seals and process limits for the local replay."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '4'
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
import json
import resource

HERE = Path(__file__).resolve().parent
PINS_SHA256 = '4ce48c72ce7beb6264ab9c3f37f298a88f2f59a7f22879b9c33196cfe53484ac'
WEIGHT_SHA256 = 'a7f6aa381520f5c68fa1f099a6ef3ebc3b8ce846709b97dfde2053f26ca62f80'
ENCODER_SHA256 = 'b0bf97d757a75b31607a2fdc76ff8800b9be63a154945a62631ae95158eb616e'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def check_pins(root=HERE):
    require(sha(root/'pins.json') == PINS_SHA256, 'SHA-256 mismatch: pins.json')
    pins = json.loads((root/'pins.json').read_text())
    require(pins['inputs/CRISPR_Net_CIRCLE_elevation_SITE_weights.h5'] == WEIGHT_SHA256,
            'Checkpoint seal differs')
    require(pins['vendor/Encoder_sgRNA_off.py'] == ENCODER_SHA256, 'Encoder seal differs')
    for rel, expected in pins.items():
        require(sha(root/rel) == expected, 'SHA-256 mismatch: '+rel)
    return pins


def peak_rss():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    value = value if sys.platform == 'darwin' else value*1024
    require(value < 2_000_000_000, 'Resident memory limit exceeded')
    return value
