# Scripts of addendum v21 (ψ₁*, ψ₀* at odd roots of unity)

Python 3 with numpy and mpmath. Run from inside this folder: `python3 k1_odd.py > run_log.txt` (about two minutes).

- `k1_odd.py [a b c]`
  - a: Theorem 11. A = Σ_{n<X} (−1)ⁿ q₀^{n(n+1)/2} (−q₀;q₀)ₙ (mpmath, exact) against P·[3J − I(+) − I(−)] at x = 1/X, X odd from 5 to 61, and the same for ψ₀*; integrals in double precision with `h1_roots.PI_np` (line L) and `h1_roots.PI_w` (descent contours of v19).
  - b: Lemma R2′ at X = 5 along x = 1/X + iε: W(p), p·W₃(p), F_a, F_b, ψ·P_± and the identities (V), (VI) near p = −1 (up to 300 digits), and ψ₁*(q) against A/3.
  - c: Theorem 12: A/3 against ⅓·3^{−1/4}(X/2)^{1/2}e^{Cl₂X/2π}e^{2πi/3 − iπX/6}, X = 101 … 3201.
- `h1_roots.py` (from v20), `fad.py` (from v16).
- `run_log.txt`.
