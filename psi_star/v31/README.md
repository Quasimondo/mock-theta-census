# Scripts of addendum v31 (the cocycle Ω_γ continues across (−½, 0))

Python 3 with mpmath and numpy. Run from inside this folder: `python3 c1_continuation.py > run_log.txt` (about 9 minutes).

- `c1_continuation.py [a b c d e]`
  - a: Proposition 1, 𝐉̃(−a; x)·𝐒ᶜ(a, x) = Δ(a, x)·𝐉(a; x) for Im x > 0, with 𝐒ᶜ(a, x) = conj 𝐒(−ā, −x̄) by quadrature at 30 digits (6 minutes).
  - b: the dual shift 𝐒(a + 1, x) = 𝐂·𝐒(a, x), for Im x > 0, real x and Im x < 0 (double precision).
  - c: Theorem 22 (i), 𝐑(a′, x′) = 𝐑(a, x) for Im x < 0, 𝐑 = 𝐒ᶜ·𝐒⁻¹, at six points; first row and block separately.
  - d: Ω⁽²⁾ = 𝐒ᶜ(a, x)⁻¹·𝐒ᶜ(a′, x′) against the series above the real axis and against Ω⁽¹⁾ = 𝐒(a, x)⁻¹·𝐒(a′, x′) below it; the value of Ω_γ on the interval (−½, 0) from state integrals at positive moduli, approached from both sides; a circle around −0.25.
  - e: the algebra of Lemma U at random Laurent polynomials, and an explicit 𝐕 ≠ 1 for |q| < 1 (which shows that |q| > 1 is needed).
- In double precision the matrix 𝐒(a, x) is computed by `Sany`: Re a is first reduced to [0, 1) with the dual shift, and for Im x ≥ 0 the line of integration is steeper (δ₀ up to 0.75) than in v30, because the integrands at the reflected variable −ā are large.
- `frames.py`, `npint.py`, `lowint.py`, `fad.py`, `common.py`, `a1_formula.py`: from v30, unchanged.
