# Scripts of addendum v20 (ψ₁*, ψ₀* at even roots of unity)

Python 3 with numpy and mpmath. Run from inside this folder: `python3 h1_roots.py > run_log.txt` (about one minute).

- `h1_roots.py [a b c]`
  - a: Theorem 9. The finite sums at q₀ = e^{2πi/X} (mpmath, exact) against P·J + P·I(+) − P·I(−) at x = 1/X, X = 4 … 60, in double precision with Φ_b from Faddeev's integral. `PI_np` integrates over the line L (used for J, and for the decaying ψ₀* integral when X ≤ 8); `PI_w` integrates the state integrals over the descent contours of v19 (over L they cancel by many orders of magnitude).
  - b: the radial limits of Lemmas R1, R2 at X = 8.
  - c: Theorem 10: exact finite sums divided by 3^{−1/4}(X/2)^{1/2}e^{Cl₂X/2π}e^{2πi/3 − iπX/6} (times e^{2πi/3} for ψ₀*), X = 100 … 3200.
- `fad.py` (from v16): Faddeev's integral in the variable w.
- `run_log.txt`.
