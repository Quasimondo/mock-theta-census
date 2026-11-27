# Scripts for addendum v15 (the descent inequality; P J ~ Phi on both sides of the real x-axis)

Python 3 with mpmath; numpy and scipy for w0, w1. Run from inside this folder. `run_log.txt` is the output of the commands listed in it.

- `w0_concavity.py` — on how many admissible straight contours the leading exponent is concave (seconds).
- `w1_descent.py` — Proposition 1: the exact formula for the exponent, the bounds on X and Y, the margin kappa for delta = 10 degrees (one minute).
- `w2_checks.py [a b c d]` — checks for Theorem 2 with contour angle 10 degrees: (a) inversion relation for |q| > 1, (b) Euler-Maclaurin
  on the contour on both sides, (c) the uniform bound, (d) remainder of the series against the kernel integral.
  Parts a, b: seconds. Part c: about 5 minutes. Part d: about 15 minutes.
- `borel.py`, `phi_exact.json` — from v13: the kernel integral P J(x) on both sides of the real axis, and the exact coefficients of Phi.
