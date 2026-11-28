# Scripts of addendum v32 (the cocycle Ω_γ at x = 0: one formula on the whole domain)

Python 3 with mpmath and numpy. Run from inside this folder: `python3 d1_zero.py > run_log.txt` (about 12 minutes; part a at 30 digits takes most of it).

- `d1_zero.py [a b c d]`
  - a: Theorem 23 for Im x > 0, max|𝐒(a, x+1)⁻¹·𝐒ᶜ(a′, x′−1) − 𝐉⁻¹Λ⁻¹𝐉′|, by quadrature at 30 digits (two points) and in double precision (four points, two of them with |x| ≤ 0.07).
  - b: the constant of the proof, i·δ(a, x+1)/δ(a′, x′−1) = λ(a, x).
  - c: the new formula Ω⁽³⁾ against Ω⁽¹⁾ = 𝐒(a,x)⁻¹𝐒(a′,x′) for Im x < 0 and on the positive real axis, and against Ω⁽²⁾ = 𝐒ᶜ(a,x)⁻¹𝐒ᶜ(a′,x′) on the interval (−½, 0).
  - d: the value at x = 0, and circles around 0 of radius 0.05, 0.1, 0.2: Ω⁽³⁾ against the series on the upper half, against the old formulas on the lower half, the mean value, and the negative Fourier modes. The old formulas need state integrals at the small moduli x and x′ and lose accuracy as |x| → 0 (at |x| = 0.05 they are useless in double precision); Ω⁽³⁾ only needs moduli near 1.
- `smat.py`: 𝐒(a, x) in double precision (the function `Sany` of `psi_star/v31/c1_continuation.py`), 𝐒ᶜ, and the three formulas Ω⁽¹⁾, Ω⁽²⁾, Ω⁽³⁾.
- `frames.py`, `npint.py`, `lowint.py`, `fad.py`, `common.py`, `a1_formula.py`: from v30, unchanged.
