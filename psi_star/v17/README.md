# Scripts of addendum v17 (factorial constants; P·J is the Borel sum of Φ)

Python 3 with numpy, scipy (imported in g1_shift.py), mpmath. Run from inside this folder. `run_all.sh` runs everything and writes `run_log.txt` (about 6 minutes).

- `g1_shift.py`
  - (1) the shift identity P·J = Σ_{n<N} (−1)ⁿ τₙ + (−1)^N e^{t/8} ∫ H·D_N ds for N = 1, 2, 5, at real x and on both sides, ψ₁* and ψ₀*;
  - (2) |P·J − Σ_{n<N} (−1)ⁿ τₙ| as a function of N at real x;
  - (3) the elementary constants of the proof (c₃, L, λ, B₀).
  - Functions: `integral(t, x, N, case)` (Gauss–Legendre on the contour with δ = 10°, Φ_b from `fad.logPhi_np`; valid for |arg x| < 90°), `tau(n, t, case)`, `tx(T, theta_degrees)`.
- `g2_checks.py [a b c d e]`
  - a: the bound |D_N| ≤ exp(−ν N + L c₁ |s|) on the contour (Corollary G) and sup |h′| (Lemma G);
  - b: |τₙ| ≤ (e^c/2) λⁿ on |t| = c/n (Lemma T);
  - c: aₖ = Σ_n (−1)ⁿ [tᵏ] τₙ in exact rational arithmetic against `phi_exact.json`; size of aₖ against A₀ k! B₀ᵏ;
  - d: Gevrey ratios (|P·J − Σ_{k<K} aₖ tᵏ| / (K! |t|ᴷ))^{1/K};
  - e: P·J against the Borel sums in the directions −120°, −100°, −80°.
- `fad.py` (from v16): Faddeev's integral in the variable w.
- `borel.py`, `phi_exact.json` (from v13): Borel–Padé–Laplace, kernel integral with the product formulas, exact coefficients a₀ … a₃₀₀.
- `run_all.sh`, `run_log.txt`.

Precision: everything at real x and the fast integrator are double precision (about 10⁻¹³). Parts d and e use mpmath with 30 digits for the directions −200° and 5°, where the fast integrator does not apply.
