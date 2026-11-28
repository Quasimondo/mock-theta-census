# Scripts of addendum v30 (the 3 × 3 matrix of state integrals; the q-frame under Γ₀(2); the relations (R3), (R4))

Python 3 with mpmath and numpy. Run from inside this folder. `run_log.txt` is the output of all seven scripts in the order below (about 40 minutes; `m3_modular.py` takes 25 of them).

- `r1_letter.py [N]` (§1): exact integer series to order N = 300. The relations (A0), (A1) of Ramanujan's last letter in the form used, (R1)–(R4) of v24, and four further relations of the two fifth-order families.
- `m1_matrix.py` (§2): Theorem 20, 𝐉̃·𝐒 = Δ·𝐉, with 𝐒 by quadrature at 30 digits; det 𝐒, det 𝐉, P elementary; Lemma W and the equation of H.
- `m2_extend.py` (§3, §6): the matrix 𝐒 for real x and for Im x < 0 in double precision. Against the series for Im x > 0; Fourier modes on circles around points of the positive real axis; for Im x < 0 the equations L·𝒥 = 1, L·𝒾 = 0 and I(u; x) = I(u; 1/x).
- `m4_block.py` (§5, §6): the 2 × 2 block of Ω_γ on both half-planes by series (for Im x < 0 Theorem C at 1/x). Fourier modes on circles around 0.5, −0.25, −0.2, −0.7, and the jump across the negative real axis at seven points.
- `m3_modular.py` (§4): Theorem 21 (T), (γ), the formula for λ, and the relation E₂² = −1, which expresses F(1/z) by F(z), F(zq).
- `m5_jump.py` (§5): the full Ω_γ just above the negative real axis (series) and just below (integrals). Part (a) compares the integrals of `lowint.py` with those of `npint.py` (printed as `nan` where Re x < 0, since `npint.py` needs Re x > 0) and with the series of the block.
- `m6_circle.py` (§5): Fourier modes of the full Ω_γ on circles around −0.2 and −0.75.

Modules:

- `frames.py`: the q-frame 𝐉, the dual frame 𝐉̃, the factor Δ, P in elementary form, and 𝐒 by quadrature (mpmath). The line of integration is e^{−iδ₀}·ℝ with δ₀ = min(0.4, (π/2 − |arg b|)/2).
- `npint.py`: I(u), J(u) and 𝐒 for Re x > 0 (also real x and Im x < 0), numpy, Φ_b from Faddeev's integral in `fad.py`.
- `lowint.py`: the same for Im x < 0 (also Re x < 0), Φ_b from the product formula of Φ_{1/b}.
- `fad.py`: from v27, with one correction. Faddeev's integral is taken on the line Im σ = min(1, Re(1/x))/2 with step 0.16 times that height; before it was Im σ = ½ with step 0.08, which is wrong for Re(1/x) ≤ ½ (|x| large). The earlier addenda used it for 0 < x ≤ 1 only.
- `common.py`, `a1_formula.py`: from v18 (`common.py` with more subdivision points when the line is longer than 10).
