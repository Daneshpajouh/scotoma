#!/bin/sh
set -eu
TASK_PACKAGE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$TASK_PACKAGE"
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4 NUMEXPR_NUM_THREADS=4
mkdir -p logs results
nice -n 10 /usr/local/bin/python3 -I produce.py > logs/producer.log 2>&1
nice -n 10 /usr/local/bin/python3 -I verify.py --scores-only > logs/verifier.log 2>&1
nice -n 10 /usr/local/bin/python3 -I corruption_control.py > logs/corruption_control.log 2>&1
nice -n 10 /usr/local/bin/python3 -I check_layout.py > logs/layout.log 2>&1
nice -n 10 /usr/local/bin/python3 -I render_findings.py
if nice -n 10 /usr/local/bin/python3 -I verify.py > logs/final_verifier.log 2>&1; then
    printf 'VERIFIER_EXIT_CODE=0\n'
else
    TASK_VERIFY_EXIT=$?
    printf 'VERIFIER_EXIT_CODE=%s\n' "$TASK_VERIFY_EXIT"
    exit "$TASK_VERIFY_EXIT"
fi
nice -n 10 /usr/local/bin/python3 -I seal.py --write
nice -n 10 /usr/local/bin/python3 -I seal.py --verify
