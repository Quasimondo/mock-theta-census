# Scripts for addendum v14 (P J ~ Phi on the sector -90 deg < arg t < -12.44 deg; scan of straight contours)

Python 3 with mpmath; numpy and scipy for the scan. Run from inside this folder. `run_log.txt` is the output of the commands listed in it.

- `u1_checks.py [a b c d e]` — checks of the proof: (a) the potential U and Lemma 2, (b) Euler-Maclaurin with remainder,
  (c) the uniform bound of Lemma 3, (d), (e) the remainder of the series against the kernel integral.
  Parts a, b: seconds. Part c: about 5 minutes. Part d: about 20 minutes (most of it the direction -80 degrees). Part e: 3 minutes.
- `u2_scan.py` — leading-order scan: for which directions a straight contour descends (20 seconds).
- `borel.py`, `phi_exact.json` — from v13: the kernel integral P J(x) and the exact coefficients of Phi.
