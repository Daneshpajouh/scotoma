#!/bin/sh
# Run the deposit's checks in the order they depend on, so a reviewer does not have to know it.
#
# aligned/verify.py replays the one-byte corruption rejection, and the corrupted copy is BUILT by
# aligned/corruption_control.py rather than shipped: it is written as symlinks into whatever
# directory it runs in, which would dangle for anyone who cloned this repository. So the control
# runs first here. Running verify.py on its own in a fresh checkout stops with a missing
# changed_input tree; that is this order, not a fault.
#
# The checks need numpy and h5py. This picks the first interpreter that has both rather than
# hardcoding a path, because the one named in reproduce.sh is particular to the machine that
# produced the deposit. Override with PYTHON=/path/to/python3.
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

pick() {
    for cand in ${PYTHON:-} python3 /usr/local/bin/python3 /opt/homebrew/bin/python3 python; do
        [ -n "$cand" ] || continue
        command -v "$cand" >/dev/null 2>&1 || continue
        if "$cand" -I -c 'import numpy, h5py' >/dev/null 2>&1; then
            printf '%s' "$cand"; return 0
        fi
    done
    return 1
}

if ! PY=$(pick); then
    echo "No interpreter on PATH has both numpy and h5py, which these checks need." >&2
    echo "Install them, or set PYTHON to an interpreter that has them:" >&2
    echo "  PYTHON=/path/to/python3 $0" >&2
    exit 2
fi
echo "interpreter: $PY ($("$PY" -c 'import sys;print(sys.version.split()[0])'))"

cd "$HERE/aligned"
echo "== input seals and pin chain"
"$PY" -I seal.py --verify
echo "== one-byte corruption control (builds changed_input, then requires its rejection)"
"$PY" -I corruption_control.py
echo "== full replay"
"$PY" -I verify.py
echo "== all three passed"
