# ψ₁* and the figure-eight knot (`fig8.py`)

Comparison of ψ₁*(e^{2πi/X}) (exact values, addendum-v3 §1) with the Kashaev invariant
J_N = Σ_{n<N} |(q;q)ₙ|², q = e^{2πi/N}, of the figure-eight knot 4₁ (2026-11-26).

- Control: J_N / (3^{−1/4} N^{3/2} e^{Vol·N/2π}) → 1 (1.0056, 1.0028, 1.0014, 1.0007 at N = 100 … 800).
- Same exponential at N = X/2: |ψ₁*(e^{2πi/X})| / J_{X/2} · X = 1.9908, 1.9954, 1.9977, 1.9989 at X = 200 … 1,600,
  tending to exactly 2. So ψ₁* carries the 4₁ volume and the 4₁ torsion 3^{−1/4}, with weight ½ instead of 3/2.
- Different perturbative series: in the common variable the first correction is 0.17641i·t for 4₁ (Zagier's
  11/(72√−3)·h with h = −2t) and (1/4 + i/(18√3))·t for ψ₁*; the imaginary parts differ by a factor of ≈ 5.5.

Reading: ψ₁* is not J_{4₁}(q²) in disguise, but it belongs to the same arithmetic family (regular ideal
tetrahedron, field ℚ(√−3), torsion 3^{−1/4}). Whether it is one of the weight-½ q-series that Garoufalidis–Zagier
attach to 4₁ ("Knots and their related q-series") is not yet checked.

## First-principles derivation of the constants (`phib_asym.py`, `saddle_B.py`)

- `phib_asym.py`: checks log Φ_b(z/2πb) ~ Li₂(−e^z)/ħ + Σ B₂ₙ(½)ħ^{2n−1}/(2n)!·Li_{2−2n}(−e^z), ħ = 2πib², against
  the product formula (errors scale as ħ³, ħ⁵).
- `saddle_B.py`: stationary phase through first order for the four state integrals of addendum-v7, with the exact
  prefactor. It reproduces the saddles (z = 2πi/3, 4πi/3), |c| = (4/3)^{1/4} and all four first-order coefficients
  of addendum-v6 §2 (∓1/24 ∓ i/(18√3)) to 10⁻³², i.e. exactly; nothing fitted. The overall phase agrees up to a
  factor i common to all four, from a square-root branch in the prefactor that was not pinned down.

## The Garoufalidis–Zagier recipe applied to ψ₁* (`gz_recipe.py`, `check_term.py`, `check_term2.py`)

Garoufalidis–Zagier, *Knots, perturbative series and quantum modularity* (SIGMA 20 (2024) 055, arXiv:2111.06645),
§6.3 define the perturbative Φ-series of a (generalized) Nahm sum by writing Pochhammer symbols with their expansion
(6.8), replacing the sum by an integral and expanding formally at the stationary points. ψ₁* = Σ(−1)ⁿq^{n(n+1)/2}(−q;q)ₙ
is not of their form (6.6) (it has (−q;q)ₙ = (q²;q²)ₙ/(q;q)ₙ), but it is one of the "even more general Nahm sums" they
mention, and the recipe applies verbatim:

- `check_term2.py`: for the factor (−q^{n+1};q)∞ the expansion (6.8) is used with Bernoulli argument 0 (checked
  numerically; arguments ±1 fail).
- `check_term.py`: the expansion of log T(n) to first order agrees with the exact summand to O(h²) (10⁻⁹ … 10⁻¹⁰).
- `gz_recipe.py`: the stationary points are x = 2πi/3, 4πi/3 with values −S, −S̄ (S = π²/3 − i·Cl₂(π/3)); at the one
  with value −S̄ the first correction is **1/4 + i/(18√3)**, exactly the growing-sector constant measured in
  addendum-v3 §5.3 and rechecked in `psi_star/v3/check_law.py`; the other gives its conjugate, the decaying sector of
  addendum-v5 §5. So ψ₁*'s perturbative series are the Garoufalidis–Zagier Φ-series of a rank-1 generalized Nahm sum.
