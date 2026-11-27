# Scripts of addendum v23 (bound for the kernel integral at real x; the volume rule at even roots of unity)

Python 3 with numpy and mpmath. Run from inside this folder: `python3 p1_unitary.py > run_log.txt` (about 3 minutes).

- `p1_unitary.py [a b c]`
  - a: |Φ_b(y)| = 1 for real b, y: Re log Φ_b at 200 random real points (Faddeev's integral, `fad.py`).
  - b: |P·J(a/b)| at x = 1/X computed over the real line, against the integral over the tilted line L and against the bound 1/(2√2·sin πa).
  - c: Theorem 15: the finite sum F(z; e^{2πi/X}) divided by the saddle-point expression, X = 400 … 3200.
- `n1_zroots.py` (from v22: integrand and quadrature), `fad.py` (from v16), `run_log.txt`.
