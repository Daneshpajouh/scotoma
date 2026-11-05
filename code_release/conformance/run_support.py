"""Recorded subprocesses with local output paths and strict expected exit codes."""
import json
import hashlib
import os
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PYTHON = "/usr/local/bin/python3"


def run(label, argv, expected=0, cwd=None, timeout=60):
    evidence = HERE / "evidence"
    evidence.mkdir(exist_ok=True)
    requested_label = label
    attempt = 1
    while (evidence / f"{label}.stdout").exists():
        attempt += 1
        label = f"{requested_label}-r{attempt}"
    cwd = Path(cwd or HERE)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    argv = [str(item) for item in argv]
    source_hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in HERE.glob("*.py")}
    if argv[0].startswith("/Library/TeX/"):
        cache = HERE / "build/tex-cache"
        cache.mkdir(parents=True, exist_ok=True)
        env.update(TEXMFVAR=str(cache / "var"), TEXMFCONFIG=str(cache / "config"))
    try:
        result = subprocess.run(argv, cwd=cwd, env=env, text=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        def as_text(value):
            return value.decode(errors="replace") if isinstance(value, bytes) else value or ""
        result = subprocess.CompletedProcess(argv, 124, as_text(exc.stdout), as_text(exc.stderr) + "\nTimed out\n")
    (evidence / f"{label}.stdout").write_text(result.stdout)
    (evidence / f"{label}.stderr").write_text(result.stderr)
    record = {"label": label, "requested_label": requested_label,
              "argv": argv, "command": shlex.join(argv), "cwd": str(cwd),
              "source_sha256": source_hashes,
              "exit_code": result.returncode, "expected_exit_code": expected,
              "stdout": f"evidence/{label}.stdout", "stderr": f"evidence/{label}.stderr"}
    with (evidence / "runs.jsonl").open("a") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")
    print(f"{label}: exit {result.returncode} (expected {expected})", flush=True)
    if result.returncode != expected:
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"unexpected exit code for {label}")
    return result
