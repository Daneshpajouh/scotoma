#!/usr/bin/env python3
"""Plant a corrupted raw-scan archive; exit non-zero if the verifier accepts it."""
import os, shutil, subprocess, sys, tarfile, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CSVDIR = sys.argv[1] if len(sys.argv) > 1 else None


def run(dst, csvdir):
    cmd = [sys.executable, "-I", "-B", "verify_published_dataset_hole.py"]
    if csvdir:
        cmd.append(csvdir)
    return subprocess.run(cmd, cwd=dst, capture_output=True, text=True)


def empty_archive(dst):
    # replace archive with empty tar.gz
    p = os.path.join(dst, "scan_raw_outputs.tgz")
    os.remove(p)
    with tarfile.open(p, "w:gz") as tf:
        pass


def copy_deposit(dst):
    for name in os.listdir(HERE):
        src = os.path.join(HERE, name)
        if os.path.isfile(src) and not name.endswith(".bak_r7_20261101"):
            shutil.copy(src, os.path.join(dst, name))


print("%-28s %-8s %s" % ("corruption", "verdict", "expected"))
d = tempfile.mkdtemp()
unexpected = 0
try:
    copy_deposit(d)
    empty_archive(d)
    proc = run(d, CSVDIR)
    accepted = proc.returncode == 0
    verdict = "accepted" if accepted else "rejected"
    expected = "rejected"
    bad = accepted
    unexpected += bad
    print("%-28s %-8s %s%s" % ("dataset_empty_scan", verdict, expected,
                               "   <-- UNEXPECTED" if bad else ""))
    if bad:
        sys.stderr.write(proc.stdout[-2000:] + "\n" + proc.stderr[-1000:] + "\n")
finally:
    shutil.rmtree(d)
print("1 corruption, %d unexpected" % unexpected)
sys.exit(1 if unexpected else 0)
