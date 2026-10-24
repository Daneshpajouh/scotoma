# Reference outputs

Compare your own run with these, line by line. Timings are from one Mac with Python 3.9.6 and will differ on yours.

## From the files in this folder, in a clean-room run

Each of these was made by copying only the four `.py` files into an empty folder and running them from a different empty folder, with a cleared environment and `python3 -I`. `CLEANROOM.md` lists every command and its exit code.

- `tail_and_lemmas.log`: `python3 frac_alln_tail_check.py --lemmas`
- `centres_m6_N27_k3.log`: `python3 frac_alln_tail_check.py --centres 6 27 3`
- `finite_m7_N19.log`: `python3 frac_alln_tail_check.py --finite 7 19`
- `finite_m12_N29-32.log`: `python3 frac_alln_tail_check.py --finite 12 29 30 31 32`
- `step3_dp_selftest.log`: `python3 step3_dp.py`
- `privsphere_selftest.log`: `python3 privsphere.py`

## All 143 finite maxima

`finite_all_143_rows.log` has one line for every (m, n) with m = 6 to 12 and 20 <= n <= n_0 - 1, which is every length the subsection computes. Each line is `step3_dp.scaled_max(N, m, N + 1)`, called by the short wrapper `release_boundary_rows.py` (kept here as a record; to run it, put it beside `step3_dp.py`) and run on 2026-10-23. The wrapper prints M as an exact fraction and a decimal, a maximising centre over {0,1,X}, and `every_r_proved`, which is True exactly when `--finite` prints `every_r=True`. Those runs used `step3_dp.py` before its import line was changed to find `privsphere.py` in its own folder. The computation itself is unchanged. To reprint any line with the checker, use `python3 frac_alln_tail_check.py --finite m N`.
