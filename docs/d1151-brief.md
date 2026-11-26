# Brief: identify R = (D − θ₄·f)/4

For a fresh session. The task is self-contained, and everything needed is in this brief and in the public repository https://github.com/Quasimondo/mock-theta-census (start with `docs/d1151.md`). The goal is a closed form, or at least a structural identification, of one q-series R(q). An analytic derivation is especially welcome, because a blind computational search has already been tried and failed (section 4).

Notation: (a;q)ₙ = ∏_{k=0}^{n−1}(1 − aqᵏ), (q)ₙ = (q;q)ₙ, (−q)ₙ = (−q;q)ₙ = ∏_{k=1}^{n}(1 + qᵏ).

## 1. The objects

- **D:** the one mock-theta-type candidate of a census of 6,327 double sums.

  D(q) = Σ_{m,n≥0} q^{m² − mn + n²} / ((−q)ₘ (−q)ₙ)

  The exponent is the A₂ quadratic form (as in the Andrews–Schilling–Warnaar A₂ Rogers–Ramanujan identities, which use (q)ₙ denominators). Without the cross term the sum is f₀(q)², where f₀(q) = Σ q^{n²}/(−q)ₙ is Ramanujan's fifth-order mock theta function. So D is an A₂ deformation of f₀².
- **θ₄(q)** = Σ_{n∈ℤ} (−1)ⁿ q^{n²} = (q)∞/(−q)∞.
- **f(q)** = Σ q^{n²}/(−q)ₙ²: Ramanujan's third-order mock theta function.
- **R(q)** := (D(q) − θ₄(q)·f(q)) / 4.

Coefficients of D, q⁰ … q⁶⁰:
`1, 3, -4, 7, -7, 7, -11, 19, -22, 22, -30, 42, -45, 49, -65, 80, -87, 95, -111, 138, -157, 164, -196, 239, -257, 280, -324, 374, -416, 447, -510, 596, -649, 695, -786, 903, -979, 1049, -1187, 1335, -1454, 1567, -1724, 1947, -2121, 2264, -2512, 2802, -3025, 3247, -3574, 3945, -4267, 4583, -5007, 5524, -5959, 6368, -6950, 7624, -8203`

Coefficients of R, q⁰ … q¹⁰⁰:
`0, 1, 0, 0, 0, -1, 1, -1, 1, 0, 1, -2, 2, -1, 1, -2, 2, -2, 2, -4, 4, -2, 3, -5, 4, -2, 4, -6, 6, -5, 5, -7, 8, -5, 5, -9, 8, -7, 9, -11, 11, -10, 10, -13, 14, -11, 13, -17, 15, -13, 16, -19, 20, -18, 19, -24, 25, -21, 22, -28, 28, -25, 28, -32, 32, -32, 32, -38, 41, -36, 40, -48, 46, -43, 47, -53, 54, -52, 55, -63, 65, -59, 63, -74, 74, -71, 78, -85, 87, -85, 87, -98, 102, -97, 103, -116, 116, -112, 118, -129, 134`

**Caution when computing D yourself:** because of the cross term, the terms with small exponent in row m lie near n ≈ m/2, not on the axes. Sum over all (m, n) with m² − mn + n² < N; a box bounded by the axes alone misses terms. One of our own scripts had exactly this bug.

## 2. What is known (with status)

| Fact | Status |
|---|---|
| Every coefficient of D − θ₄f is divisible by 4 (so R has integer coefficients) | verified to q²⁵⁹, not proved |
| D/(θ₄f) → 1 along rays to every even-order root tested (1/2, 1/4, 1/6, 1/8, 3/8, 1/10, 1/12, 5/12): 22 digits, change between two t values 10⁻¹⁰⁹ | numerical |
| Growth fingerprint of D at 19 roots: exact law with κ = 1/24, weight ½ at all even-order roots; bounded at odd-order roots. Identical to f's fingerprint. | numerical (classifier) |
| Fingerprint of R: bounded at the 7 odd-order roots. At all 12 even-order roots: rational exponent κ = 1/60, weight ½, but no exact law (power-series corrections). This is the behaviour of θ₄·R₀, θ₄·R₁ below. | numerical (classifier) |
| At q = 1: D(e^{−t}) = Σ aₖtᵏ with k!·aₖ integers (4, −12, 212, −9588, 813380, …, computed exactly to k = 40). Large order: aₖ ≈ K(−1)ᵏΓ(k+½)A^{−k−½}(1 + s₁A/(k−½) + …) with A = π²/60, πK² = 8(5+√5)/25, s₁ = (35+6√5)/300. Since θ₄f is exponentially small at q = 1, this is also R's expansion there. | exact coefficients; constants numerical (10–12 digits, PSLQ) |
| f₀ and f₁ (fifth order) diverge at q = 1 with the same Borel singularity π²/60 | numerical, 12 digits |

## 3. The model case that worked (for orientation)

For the sign-twisted fifth-order functions ψ₁\* = Σ (−1)ⁿ q^{n(n+1)/2}(−q)ₙ and ψ₀\* = Σ (−1)ⁿ q^{(n+1)(n+2)/2}(−q)ₙ, these identities were found and proved (`docs/addendum-v2.md`):

  2ψ₀\* = f₀\* − θ₄·R₀,   2ψ₁\* = f₁\* + θ₄·R₁

with f₀\* = Σ (−1)ⁿq^{n²}/(−q)ₙ, f₁\* = Σ (−1)ⁿq^{n²+n}/(−q)ₙ, R₀ = Σ (−1)ⁿq^{n²}/(q)ₙ and R₁ = Σ (−1)ⁿq^{n²+n}/(q)ₙ. R₀, R₁ are Ramanujan's entire Rogers–Ramanujan function A_q(a) = Σ (−a)ⁿq^{n²}/(q)ₙ at a = 1 and a = q. The proof idea:
1. Extend the sum to n ∈ ℤ.
2. Write (−q)ₙ = (−q)∞/(−qⁿ⁺¹)∞, and use the partial fraction 1/(−xq)∞ = (1/(q)∞) Σⱼ (−1)ʲ q^{j(j+1)/2} / ((q)ⱼ (1 + xq^{j+1})).
3. At the special value z = −1 an Appell-type sum becomes modular and a theta factor vanishes.

Near q = −1, θ₄R₁ behaves exactly like R does at the even-order roots: κ = 1/60 with power-series corrections. That is the reason to expect R to be θ₄ times a Rogers–Ramanujan-type series, or a combination of such terms. D = θ₄·f + 4R would then mirror the ψ\* identities with the roles distributed differently.

## 4. Already ruled out (exact linear algebra on coefficients; do not repeat)

- **Single sums times eta quotients:** R is not, up to a factor ±1, ±2, ±4 and a monomial, a single Eulerian sum Σ (±1)ⁿ q^{(An²+Bn)/2}/(factor)ₙ times an eta quotient ∏ (q^a;q^a)∞^{e_a} (a ∈ {1, 2, 3, 4, 5, 6, 10}, e_a ∈ {−1, 0, 1}). The sums cover A ≤ 8, with factor types (q;q), (q²;q²), (q;q²), (−q;q), (−q;q²), (q⁴;q⁴), (−q;q) in the numerator, 1/(−q;q)² and 1/((q;q)(−q;q)).
- **θ₄ times one of those sums:** no match; also not with θ₄(q²), θ₃, θ₄², (q;q)∞, (q²;q²)∞ or θ₄/(−q;q)∞ as the multiplier.
- **Infinite product:** R/(q·θ₄) is not an infinite product ∏(1 − qⁿ)^{−cₙ} with periodic exponents.
- **Linear combinations with q-power multipliers:** R is not a rational linear combination, to q²⁵⁰, of these 20 blocks, each also multiplied by q^j for j = −2, …, 3: θ₄R₀, θ₄R₁, f₀\*, f₁\*, f₀, f₁, ψ₀\*, ψ₁\*, ψ₀, ψ₁, f, θ₄, 1, θ₄G, θ₄H, θ₄f(−q), f₀(−q), f₁(−q), f₀\*(−q), f₁\*(−q). Here G, H are the Rogers–Ramanujan functions.
- **Pairwise products:** R is not a linear combination of all pairwise products of f₀, f₁, f₀\*, f₁\*, ψ₀\*, ψ₁\*, ψ₀, ψ₁, R₀, R₁, G, H, θ₄ and f (to q²⁵⁰).
- **Along the ray to −1:** R is not θ₄·(αR₀ + βR₁) with constant α, β (relative misfit about 10⁻⁵).

The same searches do reproduce (I), (II) and f₀(q) − 2ψ₀(q) − 2f₀(−q) + θ₄(q)G(q) = 0, so the method works. R simply lies outside that space.

## 5. Suggested routes

1. **Analytic (preferred).** Write D = Σₘ q^{m²}/(−q)ₘ · Φ(q^{−m}), where Φ(z) = Σₙ zⁿ q^{n²}/(−q)ₙ is a two-variable f₀. Try:
   - the partial-fraction / bilateral method of section 3 on Φ, or a Bailey-pair argument;
   - known two-variable expansions of Φ (Appell–Lerch or Hecke-type, as for f₀ in Hickerson–Mortenson, arXiv:1208.1421);
   - an A₂ Bailey lemma (Andrews–Schilling–Warnaar 1999) adapted to (−q) denominators.

   The aim is a derivation that explains why θ₄·f appears and what the rest is.
2. **q-difference equation.** Introduce D(x, y) = Σ x^m y^n q^{m²−mn+n²}/((−q)ₘ(−q)ₙ), find the functional equations, and see whether θ₄f + 4R separates naturally.
3. **Asymptotic identification.** Near q = −1, R ≈ θ₄·e^{π²/(60t)}·C·(1 + σ₁t + σ₂t² + …). Compute σ₁, σ₂ to high precision and compare with the same data for candidate Rogers–Ramanujan-type series. (This route may be pursued in parallel elsewhere, so the results can be cross-checked.)

## 6. What to deliver

- **An identity, or a precise negative result.** For an identity: the formula, a proof or proof sketch, and an exact coefficient check to at least q²⁵⁰ against the R coefficients above (recompute them independently). For a negative result: what space was searched.
- **Clear labels:** proved / verified exactly to q^N / numerical / guess, for every claim. Please keep these distinct; the repository is public, and the earlier documents follow this convention.
- **Literature:** if you recognize D, R or the decomposition from the literature, say where. That would be the most valuable outcome of all.
