# Scripts for addendum v16 (the wedge around the real x-axis)

Python 3 with mpmath, numpy, scipy. Run from inside this folder. `run_log.txt` is the output of `python3 z1_checks.py a b c d e`.

- `fad.py` — Faddeev's integral for log Phi_b in the variable w = 2 pi b y: mpmath version (in the strip), numpy version (any w, via the
  functional equation in the direction i/b and the inversion relation); the product formulas for comparison.
- `z1_checks.py [a b c d e]` — (a) integral = products, functional equations at real x; (b) expansion in the strip with remainder;
  (c) the leading-term bound up to the strip edge; (d) the uniform bound on the contour; (e) the kernel integral at real x against the
  Borel sum and the series.  Parts a, b: about one minute (mpmath). Parts c, d, e: about one minute (double precision).
- `borel.py`, `phi_exact.json` — from v13.
