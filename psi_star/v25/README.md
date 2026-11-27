# Scripts of addendum v25 (the formulas at roots of unity for X = 1, 2, 3)

Python 3 with numpy and mpmath. Run from inside this folder: `python3 r1_small.py > run_log.txt` (about 25 seconds).

- `r1_small.py`
  - a: Theorems 9 and 11 at x = 1/X, X = 1 … 5, for ψ₁* and ψ₀*: the exact value against the three integrals over L (double precision, Faddeev's integral).
  - b: Theorem 14 at X = 2 for five values of a.
  - c: the radial limits of ψ₁*, ψ₀* at q → 1 and q → e^{2πi/3}.
- `h1_roots.py` (from v20), `k1_odd.py` (from v21), `n1_zroots.py` (from v22), `fad.py` (from v16), `run_log.txt`.
