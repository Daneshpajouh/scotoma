# Clean-room runs of this bundle, 2026-10-23

Files copied: only the .py files of this folder (core), or those plus the add-on (addon), into an empty folder that is not the release tree. Every command was run from a second empty folder, with a cleared environment: env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I (Python 3.9.6).

    === lemmas
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/frac_alln_tail_check.py --lemmas)
    real 59.27
    exit: 0
    === centres
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/frac_alln_tail_check.py --centres 6 27 3)
    real 264.58
    exit: 0
    === finite7
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/frac_alln_tail_check.py --finite 7 19)
    real 2.58
    exit: 0
    === finite12
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/frac_alln_tail_check.py --finite 12 29 30 31 32)
    real 232.12
    exit: 0
    === step3_selftest
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/step3_dp.py)
    real 11.60
    exit: 0
    === privsphere_selftest
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/privsphere.py)
    real 0.15
    exit: 0
    === topr_check_dp
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I -c import runpy; g = runpy.run_path('<cleanroom>/rel/topr.py', run_name='imported'); ok, v, i = g['check_dp'](8, 3); print('topr check_dp(8,3):', ok, v, i); raise SystemExit(0 if ok else 1))
    real 0.13
    exit: 0
    === CTRL_checker_alone_finite
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/ctrl_missing/frac_alln_tail_check.py --finite 7 19)
    real 0.08
    exit: 1
    === CTRL_unpatched_topr_import
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I -c import runpy; runpy.run_path('<cleanroom>/ctrl_orig/topr.py', run_name='imported'); print('import ok'))
    real 0.10
    exit: 1
    === topr_run
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/topr_run.py)
    real 164.71
    exit: 0
    === close_r6
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/close_r6.py)
    real 74.38
    exit: 0
    === verify_cert_19_9_9
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/verify_ceiling_certificate.py <cleanroom>/rel/cert_19_9_9.json)
    real 1.28
    exit: 0
    === verify_cert_19_10_7
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/verify_ceiling_certificate.py <cleanroom>/rel/cert_19_10_7.json)
    real 1.21
    exit: 0
    === verify_cert_19_11_6
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/verify_ceiling_certificate.py <cleanroom>/rel/cert_19_11_6.json)
    real 1.01
    exit: 0
    === CTRL_verify_corrupted
    cmd: (cd <cleanroom>/elsewhere && env -i HOME=$HOME PATH=/usr/bin:/bin /usr/bin/python3 -I <cleanroom>/rel/verify_ceiling_certificate.py <cleanroom>/ctrl_bad/cert_19_9_9_bound_plus_one.json)
    real 0.09
    exit: 1
