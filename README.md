# mock-theta-census

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24443194.svg)](https://doi.org/10.5281/zenodo.24443194)

A complete numerical census of short Ramanujan-style q-series, classified as **modular**, **mock theta**, **false/quantum type** or **neither** by their radial asymptotics at roots of unity. Along the way it found one pair that does not appear in the literature we checked:

  **ψ₁\*(q) = Σₙ≥0 (−1)ⁿ q^{n(n+1)/2} (−q;q)ₙ**, **ψ₀\*(q) = Σₙ≥0 (−1)ⁿ q^{(n+1)(n+2)/2} (−q;q)ₙ**

These are Ramanujan's fifth-order mock theta functions ψ₁, ψ₀ with a sign twist (−1)ⁿ. They are bounded at every root of unity and terminate at even-order roots, but they are neither modular nor mock theta functions. What they are is now largely worked out (items 3–7 below):

- **Identities.** Each is half of a sign-twisted fifth-order function plus or minus θ₄ times Ramanujan's entire Rogers–Ramanujan function.
- **Quantum modularity.** Near q = 1 they obey a transformation law of the kind known for the Kashaev invariant of a knot: exponential growth governed by the volume of a regular ideal tetrahedron, a new dual sum W at the transformed root of unity, and three sectors. The forward law is now a theorem, inside the unit disk (v19) and at the roots of unity e^{2πi/X} (v20, v21).
- **Exact formulas, proved.** ψ₁\*, ψ₀\* and their companions f₁\*, f₀\* equal the Borel sum of their perturbative series plus state integrals of Faddeev's quantum dilogarithm times the dual series. The proofs (addenda v12, v14–v17, v19–v21) were written by one AI model and checked numerically by another; they have not been refereed.
- **A bridge.** ψ₁\* = F(−1) and Ramanujan's ψ₁ = F(1) are two members of one family F(z) = Σ zⁿq^{n(n+1)/2}(−q;q)ₙ with one exact formula for all z. Its growth is the volume of an ideal tetrahedron with shape W² + W = 1/z. At z = 1 the tetrahedron is flat, and the formula becomes, term by term, the classical transformation law of Ramanujan's fifth-order mock theta functions: its kernel integral is a Mordell integral (v24).

As far as we could find, neither the functions, the identities, the dual sums nor the laws appear in the literature. **We'd be glad to hear if they do.**

Two other findings: the same sign twist gives bounded companions for several other mock theta functions (fifth-order F₀, F₁, eighth-order T₀; for the third-order ψ it gives the known −σ\*/2 of Andrews–Dyson–Hickerson). And a census of 6,327 double sums found one mock-theta-like double sum, **D(q) = Σ_{m,n≥0} q^{m²−mn+n²} / ((−q;q)ₘ (−q;q)ₙ)**, now explained by a proved identity **D = θ₄·f + 4R**, with R = Σ_{d>k≥0} (−1)^{d+k+1} q^{d²}/(−q;q)ₖ a sum of theta tails.

**Interactive pages:** [quasimondo.github.io/mock-theta-census](https://quasimondo.github.io/mock-theta-census/)
- [Explorer](https://quasimondo.github.io/mock-theta-census/explorer-v4.html): mock theta functions on the unit disk, with a radial probe toward any root of unity, including ψ₁\* and ψ₀\* and the double sum D
- [Fingerprints](https://quasimondo.github.io/mock-theta-census/fingerprints-v3.html): the classifier on 51 known functions, both censuses, the sign twists, and the double sum D
- [Formula sheet](docs/formulas.md): every formula in LaTeX, Unicode, plain text and Mathematica

**Status labels.** Every item or bullet under Main results carries one of four labels: **proved** (a complete argument is written in the cited addendum), **proved in part** (the label says which part), **numerical** (verified, usually to 15–50 digits) or **heuristic** (a derivation or fit not verified to that standard). "Exact to qᴺ" is a series check that accompanies a proof, not a status. From addendum v3 on, the addenda were written by Claude Fable 5.2, and each one was rechecked independently by Claude Opus 5.5 with separate code (the `check_*.py` scripts and the note at the top of each addendum). The analytic estimates in the proofs were read but not checked line by line, and no mathematician has refereed them yet. Computations used Arb ball arithmetic (python-flint) and mpmath at 30–300 digits, and exact integer arithmetic for series identities.

A one-page summary of the state of the ψ\* work, with a status table, is in [`docs/psi-star-handoff-v21.md`](docs/psi-star-handoff-v21.md).

**For referees:** the theorems on ψ₁\*, ψ₀\* and their proofs are collected in one document, in logical order with one statement of the results, a map of the proof and a list of the places most in need of checking: [`docs/psi-star-theorems-v1.md`](docs/psi-star-theorems-v1.md).

## Main results

1. **A validated classifier (numerical).** Along q = e^{−t}ζ, fit log|g| = c/t + α log t + β + γt (+ exponentially small terms) at 19 roots of unity ζ, and test whether κ = c·m²/4π² is rational and whether the fit is exact.
   - On a labeled catalog of 51 known functions, it classifies every one correctly except a sum of two modular forms with different q-offsets, which is a documented limitation.
   - It recovers known normalizations (q^{−1/24}f, q^{−1/60}f₀, q^{−1/168}F₀ of order 7) and dilogarithm values of non-modular Nahm sums to 11–12 digits.
2. **A census (numerical).** All 9,936 formulas Σ ±q^{(An²+Bn)/2} × (one or two q-product factors from a fixed list) are classified. They reduce to 8,869 distinct functions: **25 mock, 28 modular, 1,357 bounded at all 19 roots**, and the rest rejected.
   - 23 of the 25 mock hits are known functions, rescalings, or specializations of universal mock theta functions. This includes the elementary identity **H(x) = (1−x²)·g₂(x,q) + 1 − x** (proved by telescoping), which identifies two hits as g₂(∓ω, q).
   - The other two hits are ½ψ + ½ψ\*, i.e. mock only in Ramanujan's original sense.
3. **Identities for ψ\* (proved, `docs/addendum-v2.md`; exact to q³⁹⁹):**
   - (I) 2ψ₀\* = f₀\* − θ₄·R₀,  (II) 2ψ₁\* = f₁\* + θ₄·R₁

   with f₀\* = Σ(−1)ⁿq^{n²}/(−q;q)ₙ, f₁\* = Σ(−1)ⁿq^{n²+n}/(−q;q)ₙ, θ₄ = Σ(−1)ⁿq^{n²}, R₀ = Σ(−1)ⁿq^{n²}/(q;q)ₙ, R₁ = Σ(−1)ⁿq^{n²+n}/(q;q)ₙ (R₀, R₁ are Ramanujan's entire Rogers–Ramanujan function at a = 1, q). These are the z = −1 analogues of the classical bilateral-series identities for f₀ + 2ψ₀ and f₁ + 2ψ₁; the method of relating fifth-order functions to modular bilateral series is that of Bajpai, Kimport, Liang, Ma and Ricci (*Bilateral series and Ramanujan's radial limits*, Proc. AMS 2015).
4. **Values and asymptotics at roots of unity.**
   - **Radial limits (proved, `docs/addendum-v3.md` §1):** at even order m the series terminates; at odd order m the terms repeat with a factor −2 every m steps, and the limit is one third of the first m terms.
   - **Resurgence (numerical, 15–30 digits):** at q = −1, ψ₁\*(−e^{−t}) = Σ cₖ tᵏ/k! with integer cₖ = 1, −1, 27, −1237, …, Borel singularity at −π²/60 and Stokes data in ℚ(√5), matching R₀, R₁ exactly (`docs/addendum-v2.md`); at q = i the singularity is at −π²/240; at q → 1 at π²/3 ∓ i·Cl₂(π/3), with closed-form Stokes constants (`docs/addendum-v3.md` §3).
5. **Quantum modularity of ψ\*: the three laws** (`docs/addendum-v3.md`–`v5.md`; found numerically and rechecked independently; the forward law is now proved, inside the disk in v19 and at the roots e^{2πi/X} in v20, v21, see item 6).
   - **Forward law (proved in part: the leading term, v20, v21; the correction s̄₁ is numerical):** As X → ∞, ψ₁\*(e^{2πi/X}) = C̄·t^{−1/2}·e^{−S̄/t}·W(e^{−πiX})·(1 + s̄₁t + …) + O(1), with t = −2πi/X, S̄ = π²/3 + i·Cl₂(π/3), C̄ = √(π/√3)·e^{5πi/12}, s̄₁ = 1/4 + i√3/54, and a **new dual sum W(ξ) = Σ(−1)ʲ ξ^{j²+j}(ξ;ξ²)ⱼ**. The growth rate is Cl₂(π/3)/2π: the volume of the regular ideal tetrahedron, half that of the figure-eight knot complement. This explains why earlier fits under Γ₀(2) found a weight of ≈ 0.508 instead of ½.
   - **The W family (proved, v4):** W = G(−q) has exactly one companion W₃ = G(−q³) in G(z) = Σ zʲq^{j²}(q;q²)ⱼ, with bilateral identities (V), (VI) of the same type as (I), (II) and Ramanujan's ψ(q) = Σq^{n(n+1)/2} as prefactor.
   - **Converse and mirror laws (numerical; v4, v5):** (W, W₃) goes back to (ψ₀\*, ψ₁\*) through a 2×2 matrix with golden-ratio entries and Rogers–Ramanujan phases; near q = −1 the pair goes to itself under x ↦ (x−1)/(2x−1) ∈ Γ₀(2) with weight ½.
   - **Three sectors (numerical, v5; the forward law is proved in v19):** inside the unit disk (along horocycles) all three laws hold to about 10⁻²⁹, and the decaying sector of the forward law carries q̃·W₃(q̃). So the sectors (trivial, growing, decaying) of ψ\* pair with the basis (1, W, q̃W₃).
   - **The same pattern for F₀\*, F₁\* (numerical; v8, v9):** their duals are 𝒥(r) = Σ r^{2m²+m}/(−r²;r²)ₘ and 𝒦(r) = 2Σ r^{(n+1)²}(−r²;r²)ₙ, unsigned fifth-order-type sums at a half shift, with a new bilateral identity (VII) 𝒥 + 𝒦 = ψ(r)(r²;r⁴)∞·Σ r^{m(m+1)/2}/(r²;r²)ₘ (exact to r³²⁹, not proved), and 𝒥 goes back to (q^{−1/120}(F₀\* − 1), q^{71/120}F₁\*), the exponents of Ramanujan's unsigned pair. Both families are closed two-way dualities.
6. **Exact formulas and their proof** (`docs/addendum-v6.md`–`v17.md`, `v19.md`–`v21.md`).
   - **One function behind both sides (proved, v6; an instance of the rank-one Nahm theory, see Related work):** R₀, R₁, the bilateral parts of the W family, and X, Y of (III), (IV) are all values of the Nahm-type sum H(z;q) = Σ zᵐq^{m²/4}/(q;q)ₘ, with Casoratian identities X·N₀ + N₁·Y = 1 and E₊′E₋ + E₋′E₊ = 2.
   - **State integrals (proved, v7 + v12 §4.1):** the bilinear combinations 𝓑₁…𝓑₄ that carry the growing and decaying sectors are one-dimensional state integrals ∫Φ_b(y)·e^{−πiy²/2 + 2πuy}dy of Faddeev's quantum dilogarithm, b² = x. Their stationary-phase constants reproduce the measured ones exactly (`psi_star/knot/saddle_B.py`).
   - **The exact formula (proved, v11 + v12):** ψ₁\*(q) = P·J + (κ/2)·[𝓑₂·p·W₃(p) − 𝓑₄·W(p)], where P·J = P·∫ e^{−πiy² + πy/b}/(Φ_b(y)(1 + e^{2πy/b}))dy is a state integral with a kernel, p the dual nome; likewise for ψ₀\*, f₁\*, f₀\*. The kernel is of the same kind as the tanh kernel in the trivial-connection state integral of the figure-eight knot (Garoufalidis–Gu–Mariño–Wheeler, arXiv:2111.04763); their proof does not carry over, and v12 uses a geometric series of kernel-free integrals and a three-term recursion whose generating sums are W, W₃.
   - **Borel summability (proved, v13–v17, using the Nevanlinna–Sokal theorem):** the perturbative series Φ has a closed form (301 exact coefficients); P·J is analytic and has Φ as its asymptotic expansion, with factorial constants, on the sector −210° ≤ arg t ≤ 10° (which includes |q| = 1, where ψ₁\* does not exist); hence **P·J is the Borel sum of Φ** in the directions −120° … −80°. So for |q| < 1: ψ₁\* = (Borel sum of its perturbative series) + (state integrals) × (dual series), as a theorem.
   - **The forward law as a theorem (proved, v19):** steepest-descent contours with closed-form heights give 𝓑₄, 𝓑₃ = (4/3)^{1/4}e^{∓iπ/12}e^{−i·Cl₂/t}(1 + O(t)) and 𝓑₂, 𝓑₁ = (4/3)^{1/4}e^{±iπ/12}e^{+i·Cl₂/t}(1 + O(t)), uniformly for 0 ≤ arg x ≤ 11° including tangential approach. Hence, for 0 < arg x ≤ 10°, ψ₁\* = 𝒮[Φ₁] + (π/√3)^{1/2}t^{−1/2}[e^{5πi/12}e^{−S̄/t}W(p)(1 + O(t)) + e^{−5πi/12}e^{−S/t}p·W₃(p)(1 + O(t))], and the same for f₁\*, ψ₀\*, f₀\*: the law of item 5 with its constants, a Borel sum, exact dual functions and error terms. Checked against ψ₁\* summed directly up to |ψ₁\*| ≈ 10²⁷, where the first correction comes out as s̄₁ = 1/4 + i√3/54 (`psi_star/v19/check_thm68.py`). The converse and mirror laws are still numerical.
   - **At roots of unity e^{2πi/X} (proved, v20, v21, v25):** for even X ≥ 2 (X = 2: v25), ψ₁\*(e^{2πi/X}) = P·[J + I(1/(2b)) − I(−1/(2b))] exactly, with Faddeev's Φ_b at real b (the dual series become 1 there). Hence ψ₁\*(e^{2πi/X}) = 3^{−1/4}·(X/2)^{1/2}·e^{Cl₂(π/3)·X/2π}·e^{2πi/3 − iπX/6}·(1 + O(1/X)), and the same times e^{2πi/3} for ψ₀\*: the law as measured on the unit circle is a theorem for even X. For odd X ≥ 1 (v21; X = 1, 3: v25) the radial limit, one third of a period by the odd-root rule, equals P·[J − ⅓I(1/(2b)) − ⅓I(−1/(2b))], because W → ⅓ and p·W₃ → −⅓ at p = −1, and its asymptotics is one third of the same expression. So the forward law ψ₁\*(e^{2πi/X}) ≈ C̄·t^{−1/2}·e^{−S̄/t}·W(e^{−πiX}) is a theorem for every X ≥ 4. Checked against the exact sums up to X = 3201; the next term is −2πi·s̄₁/X for both parities (`psi_star/v20/check_thm10.py`, `psi_star/v21/check_thm12.py`).
7. **The z-family: a bridge to Ramanujan's mock theta functions** (`docs/addendum-v18.md`, `docs/addendum-v22.md`, `docs/addendum-v23.md`, `docs/addendum-v24.md`, `census3/zfamily/`, `psi_star/knot/`).
   - **Exact formula for every z (proved, from v12):** F(z) = i·P·e^{−πia}·[J(a/b) − I((a−1)/b)·𝒲(T) + p·T·I(a/b)·𝒲(p²T)], z = e^{2πia}, with T = e^{πi(2a−1)/x}: the dual side carries a Jacobi variable transformed as a → a/x. ψ₁\* is the case a = ½.
   - **Tetrahedron volumes (proved in part: at the even roots of unity e^{2πi/X} for every z ≠ 1 on the unit circle, constant and phase included, v23 Theorem 15; elsewhere numerical, 10⁻²⁶ … 10⁻³⁶):** the growth of F(z) at roots of unity is Vol(z)/2π, where Vol(z) is the volume of the ideal tetrahedron with shape W, W² + W = 1/z; the dominant saddle is the root with Im W > 0. Both state integrals are proved to have the predicted asymptotics for every z: the decaying one in v19 (Theorem 7), the growing one in v22 (Theorem 13, at real x and in a small sector), on a contour whose heights are closed-form (the ray Im w = 2π, the imaginary axis, the edge of the strip, the vertical line through the saddle, the real axis). At even roots of unity F(z; e^{2πi/X}) = i·P·e^{−πia}·[J(a/b) − I((a−1)/b) + z^X·I(a/b)] exactly (v22, Theorem 14; checked to 10⁻⁹ … 10⁻¹³). The kernel integral is bounded at real x by an elementary argument: for real b Faddeev's dilogarithm has modulus 1 on the real line, so |P·J(a/b)| ≤ 1/(2√2·sin πa) (v23, Lemma J). Together these give F(z; e^{2πi/X}) = (1 + W)·W^{1/2}·(π/(t(1 + 2W)))^{1/2}·e^{(V(w₀) + π²/12)/t}·(1 + O(1/X)) for even X → ∞ (v23, Theorem 15). Odd X for z ≠ −1 is open. For example z = −1, i, ω, −ω give 1.014942, 0.696701, 0.865129, 0.486618.
   - **z = 1 and z = q (proved in part: the laws with the kernel integral, given classical Rogers–Ramanujan identities; the Mordell integral, given in addition the known law of χ₀, χ₁ and two series identities checked to q³⁰⁰; everything verified to 10⁻³⁰; `docs/addendum-v24.md`):** the tetrahedron flattens (Vol ≈ 2π·ln φ·|a|) and the formula of the family becomes, term by term, the classical transformation law of Ramanujan's fifth-order functions ψ₁, f₁ (z = 1) and ψ₀, f₀ (z = q). The state integrals are pure powers times the Rogers–Ramanujan S-matrix entries, for example ∫ e^{−πiy²}/Φ_b(y) dy = e^{−πi/4}·(2/√5)·sin(2π/5)·e^{−πi(b² + b⁻²)/30}; the dual series are φ₀, φ₁, F₀, F₁; and the kernel integral is a Mordell integral, 2i·q^{11/60}·P·J(1/b) = (60t/π)^{1/2}·L(2/5, 10t) with L as in Gordon–McIntosh. So the trivial sector P·J(a/b) of the family interpolates between a Mordell integral (z = 1) and the Borel sum of item 6 (z = −1). A proof of the Mordell identity from the integral itself, which would reprove the fifth-order laws, is open.
   - **Figure-eight knot (proved in part: the ratio 2/X, given the known asymptotics of J; the rest numerical):** ψ₁\*(e^{2πi/X}) has the exponential growth and the torsion 3^{−1/4} of the Kashaev invariant of 4₁ at N = X/2 (|ψ₁\*|/J_{X/2} → 2/X), but weight ½ instead of 3/2 and a different perturbative series, so it is not that invariant in disguise. With v20 the ratio |ψ₁\*|/J_{X/2} → 2/X is proved for even X (given the known asymptotics of J); at odd X the radial limit is one third of that size. The Garoufalidis–Zagier recipe for perturbative series of (generalized) Nahm sums, applied to ψ₁\*, reproduces the measured constants exactly: ψ₁\* behaves like a rank-1, one-tetrahedron example of their framework.
8. **Sign twists: a pattern beyond ψ\* (numerical; identities (III), (IV) proved).** Inserting (−1)ⁿ into the other catalog mock theta functions gives more companions that are bounded at all 19 roots and are not a known function at ±q: the 5th-order F₀, F₁ and the 8th-order T₀, and the family F(z) at z = i, −ω (now covered by item 7). The twist of the 3rd-order ψ is also bounded, but it is a known function: −σ\*(q)/2, the Andrews–Dyson–Hickerson companion of σ(q). Three more twists turn out to be false theta functions. None of the dense bounded twists becomes sparse after multiplying by an eta quotient, so none of them is a single Hecke-type double sum over such a product (`census3/hecke/`). F₀\*, F₁\* follow the ψ\* pattern (proved, `docs/addendum-v3.md` §4): (III) F₀\* + Σₙ≥₁ q^{n²}(q;q²)ₙ = θ₄(q²)(−q;q)∞·X(q²) and (IV) qF₁\* − Σₙ≥₀ q^{(n+1)²}(q;q²)ₙ = θ₄(q²)(−q;q)∞·Y(q²), with X = Σ(−1)ⁿq^{n²}/(q;q)₂ₙ, Y = Σ(−1)ⁿq^{(n+1)²}/(q;q)₂ₙ₊₁. Details and the full table: [`docs/twists.md`](docs/twists.md).
9. **The double sum d1151, and an identity for it.** The rank-2 census left one candidate,

     D(q) = Σ_{m,n≥0} q^{m²−mn+n²} / ((−q;q)ₘ (−q;q)ₙ).

   - **Observed (numerical):** the growth fingerprint of Ramanujan's third-order f(q) (κ = 1/24 at all 12 even-order roots), and at q = 1 an exact expansion whose large-order behaviour carries golden-ratio constants (π²/60, π·K² = 8(5+√5)/25, s₁ = (35+6√5)/300; numerical).
   - **The identity (proved):** D = θ₄(q)·f(q) + 4R(q) with R = Σ_{d>k≥0} (−1)^{d+k+1} q^{d²}/(−q;q)ₖ. Each diagonal m − n = d of D is ±f plus a finite sum of 1/(−q)ₖ (a three-line telescoping lemma), and the signs are those of θ₄. The same argument gives an identity for the whole family Σ q^{m²−mn+n²}/((zq)ₘ(q/z)ₙ) in terms of the rank generating function, a theta series and tail sums. So D is not a new mock theta function in the strict sense; whether the general-z identity is known is open.

   Details: [`docs/d1151.md`](docs/d1151.md) (how it was found), [`docs/d1151-identity.md`](docs/d1151-identity.md) (the identity and its proof). Census results and scripts: `census2/`.
10. **One-parameter families at roots of unity (numerical; the threshold rule is heuristic)** (`census3/`). F(z), g₂(x,q) and g₃(x,q) at all z, x of order ≤ 12 (136 functions).
    - g₂, g₃ behave like mock theta functions everywhere, with κ = (1−4a)²/16 and (1−6a)²/24 at q = 1 for x = e^{2πia} below a = 1/4, 1/6.
    - F(z) is bounded at even-order roots. At odd order m it grows exactly when z^m is closer than a₀ = 0.107366 of a turn to 1, so it is bounded everywhere only for z = −1, ±i, −ω, −ω̄. a₀ is where a dilogarithm saddle-point value changes sign; its predicted κ matches the fits to 4–5 digits (heuristic).
    - Negative result: ψ₀\*, ψ₁\* have no q → 1/q partner (like σ ↔ σ\*) among 1,110 testable bounded census functions (`census3/qm/`).

## Open questions

- **Is f\* a Nahm-type or indefinite-theta object?** R is a rank-one Nahm sum (see Related work); whether f₁\*, f₀\* (or F₀\*, F₁\*) can be written as Nahm sums, generalized Nahm sums or indefinite theta series decides how much of the structure is genuinely new.
- **Is any of this known?** ψ₀\*, ψ₁\*, the family F(z), the dual sums W, W₃, 𝒥, 𝒦, the identities (I)–(VII), the laws of items 5–7, the general-z identity for Σ q^{m²−mn+n²}/((zq)ₘ(q/z)ₙ) (`docs/d1151-identity.md` §3).
- **An independent reader** for the proofs in addenda v12, v14–v17 and v19–v21, and for the formulas of v18. The consolidated document `docs/psi-star-theorems-v1.md` is meant as the starting point (its §0.5 lists where to look first).
- **From exact formulas to asymptotic laws:** the tetrahedron-volume rule of item 7 at odd roots of unity for z ≠ −1 (even ones are done, v23); radial limits at the other roots of unity e^{2πih/k}, h ≠ 1 (e^{2πi/X} is done, v20, v21); the converse and mirror laws of item 5; Borel summability in directions outside −120° … −80°, and the Stokes jump at arg t = −17.1°.
- **The rest of the structure:** Theorem A of v12 for |q| > 1; the trivial sector for general z; a proof from the integral itself that the kernel integral at z = 1 is the Mordell integral (v24 proves it from the known law of χ₀, χ₁); whether 𝒲(T; p) is an Appell–Lerch sum; the full 3×3 matrix (a matrix-valued quantum modular form with a Jacobi variable in the sense of Garoufalidis–Zagier).
- **Smaller items:** the root e^{2πi/3} (x = ⅓ lies outside the range of the bounds used in v20, v21) and the correction −2πi·s̄₁/X at roots of unity as a theorem; a proof of (VII); W, W₃ and 𝒥 in their bounded directions; identities for the twist T₀\*; a proof of the F(z) threshold rule of item 10.

## Related work

The asymptotic and exact-formula part of this project uses the framework of quantum modularity for knots and Nahm sums of Garoufalidis, Zagier and collaborators. Two references describe the closest known results:

- S. Garoufalidis, D. Zagier, *Knots, perturbative series and quantum modularity*, SIGMA 20 (2024), 055 (arXiv:2111.06645). Their recipe for perturbative series of Nahm sums reproduces our measured constants for ψ₁\* (`psi_star/knot/gz_recipe.py`).
- C. Wheeler, *Modular q-difference equations and quantum invariants of hyperbolic three-manifolds*, PhD thesis, Bonn 2023. It treats rank-one Nahm sums f_{A,m,n} = Σ e(nk) q^{Ak²/2+km}/(q;q)ₖ in general. Their asymptotics come from the Nahm equation 1 − X = e(n)X^A (§4.7). They are quantum modular (Theorems 5 and 6, §8.7). The proof factorises a state integral ∫ e(Az²/τ + uz/τ)/Φ(z;τ) dz into q- and q̃-series (§9.1, §10.1, worked out for A = 4). Nahm sums satisfy quadratic relations (Proposition 7).

**What this means for our results.** R₀ and R₁ in identities (I) and (II) are the rank-one Nahm sums f_{2,0,½} and f_{2,1,½}. Their Nahm equation 1 − X + X² = 0 has the roots e^{±iπ/3} of the regular ideal tetrahedron, so the volume Cl₂(π/3) in the laws of ψ\* comes from them. The sum H(z;q) of item 6 is the dual member with A = ½. The kernel-free state integrals I(u) of v7 and v12 (Theorem C) are of the type Wheeler factorises, and the Casoratian identities of v6 are quadratic relations of the kind in his Proposition 7. Those parts are instances of known methods, here for A = 2 and A = ½ with a phase, which are not written out there.

We found no treatment of the following:
- the sign-twisted fifth-order functions ψ₀\*, ψ₁\*, f₀\*, f₁\*, F₀\*, F₁\*;
- the identities (I)–(VII);
- the dual sums W, W₃, 𝒥, 𝒦, which are not Nahm sums;
- the state integral with the kernel 1/(1 + e^{2πy/b}) and its Borel summability;
- the forward laws as theorems, inside the disk and at the roots of unity e^{2πi/X};
- the z-family and its limit to the transformation laws of Ramanujan's fifth-order functions ψ₀, ψ₁, f₀, f₁ (laws of the type Watson found for third order; for fifth order due to Zwegers and Gordon–McIntosh).

In short: ψ\* = ½(f\* ± θ₄·R) splits into a rank-one Nahm part R, which is covered by the known theory, and a sign-twisted fifth-order part f\*. The f\* part is not a Nahm sum but is quantum modular as well, with its own non-Nahm duals and a kernel state integral. The new material is the f\* side, the single exact formula that joins both parts, and the bridge to Ramanujan's mock theta functions.

Also checked, without overlap: Folsom–Ono–Rhoades and later work on radial limits of mock theta functions, Lovejoy's quantum q-series identities (which contain the unsigned ψ₁, not ψ₁\*), Folsom–Metacarpa on (anti)quantum q-series identities (third order only), Hikami's mock and false theta functions as quantum invariants, and Lovejoy–Osburn's mock theta double sums. Pointers to anything we missed are very welcome.

## Layout

| Folder | Contents |
|---|---|
| `census/` | Evaluator and labeled catalog (`catalog.py`), per-root probe (`probe.py`), generator and dedupe (`search.py`), staged search (`stage.py`, `filt.py`, `dclass.py`), q→1 dilogarithm pre-filter (`dilog.py`), identification tools (`ident.py`, `shadow.py`, `companions.py`, …), plus all result data (`cands.json`, `labels.json`, `stageC_summary.json`, `*.jsonl`) |
| `census2/` | Rank-2 (double-sum) census: generator (`gen2.py`), evaluators (`sum2_fast.py`, `sum2.py`), staged runner (`run2.py`), report (`report2.py`); results of the first run in `census2/results/`, the analysis of d1151 in `census2/d1151/`; sign-twist scripts are in `census/twist.py` and `census/twist_checks/` |
| `psi_star/` | Quantum-modularity tests with Kontsevich–Zagier control (`qmf*.py`); exact identities, expansion at q = −1 and large-order analysis (`s1`–`s11`) |
| `census3/` | One-parameter families F(z), g₂(x,q), g₃(x,q) at all roots of unity of order ≤ 12 (`zfam.py`, `zreport.py`, `zextra.py`, results and `README.md`); `zfamily/`: growth rates of F(z) as tetrahedron volumes; `hecke/`: sparsity test of the bounded twists against eta quotients and an f_{a,b,c} matcher; `qm/`: radial limits of the 1,357 bounded census functions and the q → 1/q duality search |
| `docs/` | `d1151.md` (the double sum d1151, how it was found), `d1151-identity.md` (its identity and proof), `d1151-brief.md` (the brief that led to it); `twists.md` (sign twists of the catalog functions); GitHub Pages site (`index.html`, `explorer-v4.html`, `fingerprints-v3.html`); `formulas.md` (every formula in LaTeX, Unicode, plain text and Mathematica, with equivalent forms and coefficients), `handoff-v1.md` (full account of method and results), `addendum-v2.md` (identities, proofs, resurgence analysis), `addendum-v3.md` (odd-root rule, identities (III), (IV), resurgence at q = i and q → 1, the law with W), `addendum-v4.md` (the W family, identities (V), (VI), converse and mirror laws), `addendum-v5.md` (the laws inside the disk, the third sector), `addendum-v6.md` (one function H behind both sides, Casoratians), `addendum-v7.md` (the bilinear combinations as state integrals), `addendum-v8.md` (the law for F₀\*, F₁\*, identity (VII)), `addendum-v9.md` (𝒥 back to F₀\*, F₁\*; status table), `addendum-v10.md` (the one-sided law as an exact formula with a Borel sum), `addendum-v11.md` (the trivial sector as a state integral with a kernel), `addendum-v12.md` (proof of the exact formula), `addendum-v13.md` (closed form of Φ, the kernel integral as Borel sum on both sides of the real axis), `addendum-v14.md`, `addendum-v15.md`, `addendum-v16.md` (P·J ~ Φ proved on a 220° sector), `addendum-v17.md` (factorial constants; P·J is the Borel sum of Φ), `addendum-v18.md` (the z-family: exact formula, tetrahedron volumes, the limit z → 1), `addendum-v19.md` (asymptotics of the state integrals; the forward laws as theorems), `addendum-v20.md` (exact formula and asymptotics at even roots of unity), `addendum-v21.md` (the same at odd roots of unity), `addendum-v22.md` (the growing state integral for every z; F(z) at even roots of unity), `addendum-v23.md` (bound for the kernel integral at real x; the volume rule at even roots of unity as a theorem), `addendum-v24.md` (z = 1 and z = q: state integrals as S-matrix entries, the fifth-order laws, the Mordell integral), `addendum-v25.md` (the formulas at roots of unity for X = 1, 2, 3), `psi-star-handoff-v21.md` (status of v3–v21 for a fresh session), `psi-star-theorems-v1.md` (the theorems on ψ₁\*, ψ₀\* and their proofs in one place, for referees), `psi-star-brief.md`, `psi-star-brief-v2.md` (the briefs that led to them), `psi-star-stress-test-brief-v1.md` (brief for an adversarial numerical test of every inequality in the proofs; not yet carried out) |
| `psi_star/theorems/` | Source of `docs/psi-star-theorems-v1.md`: the hand-written overview (`overview.md`) and the script that appends the proofs from the addenda (`assemble.py`) |
| `psi_star/formulas_check.py` | Checks every equivalent form in `docs/formulas.md` to 50 digits |
| `psi_star/s12_odd_limits.py` | Radial limits of ψ₁\*, ψ₀\* at odd-order roots vs. the closed form (first period)/3 |
| `psi_star/v25/` | Scripts of `docs/addendum-v25.md` (the root-of-unity formulas at X = 1, 2, 3) |
| `psi_star/v24/` | Scripts of `docs/addendum-v24.md` (series identities to q³⁰⁰, the laws at z = 1 and z = q, the Mordell integrals) |
| `psi_star/v23/` | Scripts of `docs/addendum-v23.md` (unitarity of Φ_b on the real line, bound for the kernel integral, Theorem 15 against the finite sums) |
| `psi_star/v22/` | Scripts of `docs/addendum-v22.md` (contour tests for the growing integral, F(z) at even roots of unity) |
| `psi_star/v21/` | Scripts of `docs/addendum-v21.md` (odd roots of unity), plus the independent check `check_thm12.py` |
| `psi_star/v20/` | Scripts of `docs/addendum-v20.md` (even roots of unity), plus the independent check `check_thm10.py` |
| `psi_star/v19/` | Scripts of `docs/addendum-v19.md` (contours and asymptotics of the state integrals), plus the independent check `check_thm68.py` |
| `psi_star/v18/` | Scripts of `docs/addendum-v18.md`, plus independent checks (`check_rr.py`, `check_saddle.py`, `check_thmZ.py`) |
| `psi_star/v17/` | Scripts of `docs/addendum-v17.md`, plus the independent check `check_v17.py` |
| `psi_star/v16/` | Scripts of `docs/addendum-v16.md` (`fad.py`: Faddeev's integral), plus the independent check `check_faddeev.py` |
| `psi_star/v15/` | Scripts of `docs/addendum-v15.md`, plus independent checks (`check_prop1.py`, `check_thm2.py`) |
| `psi_star/v14/` | Scripts of `docs/addendum-v14.md`, plus the independent check `check_thm1.py` |
| `psi_star/v13/` | Scripts of `docs/addendum-v13.md` (`phi_exact.json`: 301 exact coefficients), plus the independent check `check_borel.py` |
| `psi_star/v12/` | Scripts of `docs/addendum-v12.md` (`run_all.sh`, `t1`–`t5`, `common.py`), plus the independent check `check_thmAC.py` |
| `psi_star/v11/` | Script of `docs/addendum-v11.md`, plus the independent check `check_kernel.py` (with `tfun.py`) |
| `psi_star/v10/` | Scripts of `docs/addendum-v10.md`, plus independent checks (`phi_psi1.py`, `check_exact.py`) |
| `psi_star/knot/` | Comparison of ψ₁\* with the Kashaev invariant of the figure-eight knot (`fig8.py`); the Garoufalidis–Zagier recipe applied to ψ₁\* (`gz_recipe.py`); stationary phase of the four state integrals (`saddle_B.py`); `README.md` |
| `psi_star/v9/` | Scripts of `docs/addendum-v9.md`, plus independent checks (`phiJ.py`, `check_Jloop.py`) |
| `psi_star/v8/` | Scripts of `docs/addendum-v8.md`, plus independent checks (`check_VII.py`, `check_F0_law.py`) |
| `psi_star/v7/` | Script of `docs/addendum-v7.md` (`state_integral.py`), plus the independent check `check_state_integral.py` |
| `psi_star/v6/` | Scripts of `docs/addendum-v6.md`, plus independent checks (`check_casoratian.py`, `check_B1.py`) |
| `psi_star/v5/` | Scripts of `docs/addendum-v5.md`, plus independent checks (`check_M.py`, `check_psi_mod.py`, `check_converse_disk.py`) |
| `psi_star/v4/` | Scripts of `docs/addendum-v4.md`, plus independent checks (`check_V_VI.py`, `check_laws.py`, `conv_trend.py`, `check_W_stokes.py`) |
| `psi_star/v3/` | Scripts of `docs/addendum-v3.md` (`run_log.txt` is their output), plus the independent checks `check_III_IV.py` and `check_law.py` |

## Reproducing

```bash
pip install -r requirements.txt
cd census
python3 catalog.py                       # identity checks of the catalog formulas
python3 probe.py f theta3 "Rogers-Ramanujan G" "Nahm A=4,B=0" | python3 fmt.py
cd ../psi_star
python3 s1_identities.py 400             # identities (I), (II) exactly to q^400
python3 s2_expand.py 120 && python3 s6_largeorder.py && python3 s7_id.py
python3 qmf_scale.py                     # transformation test with Kontsevich control
cd v3
python3 check_III_IV.py && python3 check_law.py   # identities (III), (IV) and the law with W
cd ../v12 && python3 check_thmAC.py                # the kernel identity (Theorems A, C)
cd ../v18 && python3 check_saddle.py              # tetrahedron-volume rule with its constant
cd ../v19 && python3 check_thm68.py quick         # the forward law (Theorem 8) against the q-series
cd ../v20 && python3 check_thm10.py               # the law at even roots of unity (Theorem 10)
cd ../v21 && python3 check_thm12.py               # the law at odd roots of unity (Theorem 12)
cd ../v22 && python3 check_thm13.py               # the growing saddle contour for every z (Theorem 13)
cd ../v23 && python3 check_lemmaJ.py              # the kernel bound and the volume rule (Lemma J, Theorem 15)
python3 p1_unitary.py c                          # the volume rule at even roots of unity for six z (Theorem 15)
cd ../v24 && python3 s1_series.py && python3 s2_law.py c   # z = 1: series identities; kernel integral = Mordell integral
python3 check_thm17.py                           # independent check: Proposition M, Theorem 17, the classical laws
cd ../theorems && python3 assemble.py --check     # docs/psi-star-theorems-v1.md is current with the addenda
```

Each `psi_star/vN/` folder has the scripts of `docs/addendum-vN.md` with a `README.md` and `run_log.txt`, plus the independent `check_*.py` scripts.

The full census takes a few hours on 2 cores. `stage.py` is resumable, and `census/labels.json` holds the final label for every candidate, so nothing needs re-running to inspect results. Scripts expect to run from their own folder.

## How this was made

This was an exploratory project by Mario Klingemann, an artist and creative technologist, not a professional mathematician. The computations and much of the code were developed in conversation with AI assistants (Anthropic's Claude Opus 5.5 and Claude Fable 5.2). The addenda from v3 on, including the proofs, were written by Claude Fable 5.2 and rechecked by Claude Opus 5.5 with independent code. Every claim in the documents is backed by a script in this repository. Corrections are very welcome.

## How to cite

Mario Klingemann, *mock-theta-census: a numerical census of Ramanujan-style q-series and the sign-twisted fifth-order functions ψ₀\*, ψ₁\**, Zenodo, 2026. https://doi.org/10.5281/zenodo.24443194 (all versions; each release also has its own version DOI on Zenodo)

GitHub's "Cite this repository" button (from `CITATION.cff`) gives the same reference in BibTeX and APA.

## License

Code: MIT (see `LICENSE`). Data and documents: CC BY 4.0.
