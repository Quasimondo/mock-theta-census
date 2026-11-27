# Mock theta search: addendum v9 (the loop closes)

Session of 26 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). This does the test named at the end of `addendum-v8.md`: the unsigned half-shift function 𝒥 near r = 1, and whether its dual leads back to F₀\*, F₁\*. Scripts and output are in `psi_star/v9/`.

**Independent check** (Claude Opus, same day, `psi_star/v9/phiJ.py`, `check_Jloop.py`): the perturbative series Φ_𝒥, computed exactly by a different method (log cosh expansion and sums Σ m^d 2^{−m}), is 2 − 10τ + 183τ² − (16679/3)τ³ + …, agreeing with §1. The law of §1 was tested on two horocycles (Y₀ = 3, b = 0.8 and Y₀ = 1.5, b = 0.5; Y = Y₀ + 960k, k = 1…4, Neville in τ): the limit agrees with c·[q′^{−1/120}(F₀\*(q′) − 1) + i·q′^{71/120}F₁\*(q′)] to 9·10⁻⁸ at both, limited by the optimal truncation of Φ_𝒥 at the smallest Y.

Everything here is **verified** (numerics, no proof) unless marked otherwise.

---

## 1. Result

Let 𝒥(r) = Σ_{m≥0} r^{2m²+m}/(−r²;r²)ₘ (v8 §2), and F₀\*, F₁\* as in v3 §4. For r = e^{2πi/Y} with Y → ∞ along Im Y = −b, τ = −2πi/Y and dual nome q′ = e^{−πiY/2}:

  **𝒥(r) = Φ_𝒥(τ) + τ^{−1/2}·c·[ q′^{−1/120}·(F₀\*(q′) − 1) + i·q′^{71/120}·F₁\*(q′) ]·(1 + O(τ)),  c = √(2π/√5).**

Φ_𝒥 is the perturbative series of 𝒥 at r = e^{−τ}: 2 − 10τ + 183τ² − … .

So the dual of 𝒥 is the pair (q^{−1/120}(F₀\* − 1), q^{71/120}F₁\*). These are the exponents and the shape of Ramanujan's unsigned fifth-order pair (q^{−1/120}(F₀ − 1), q^{71/120}F₁), with the sign-twisted functions in their place. Together with v8 the loop is closed:

- F₀\*, F₁\* near q = 1 → 𝒥 at the dual nome (v8, growing, action π²/24 + (i/2)Cl₂(π/3));
- 𝒥 near r = 1 → F₀\*, F₁\* at the dual nome (here, bounded, constant √(2π/√5)).

This is the same two-way pattern as ψ\* ↔ W in v3–v5, with the same constant c as in the converse law for W.

## 2. How it was found (`jloop.py`, `jloop_fit.py`)

- Limits L(Y₀) of (𝒥 − Φ_𝒥)·√τ for Y = Y₀ + 960k − i·(4/5), k = 5…12 (K up to 11,520), with 48 terms of Φ_𝒥. Extrapolation error 2.5·10⁻²⁴. 36 values of Y₀: six base points, each with Y₀ + 4n, n = 0…5, which keeps q′ fixed.
- **Phases were not assumed.** For each base point the six values were fitted by two exponentials (Prony). The fit is exact to 4·10⁻²⁹, and the two exponents of q′ come out as 0.991666666666667 and 0.591666666666667 at all six points: −1/120 and 71/120 mod 1.
- **Amplitudes.** With these exponents, the two amplitudes were fitted over the six values of q′ on (F₀\*, 1, F₁\*, G₀, G₁). Sector −1/120: coefficients +c on F₀\* and −c on 1 (both 1.6762833568392578171, all 20 digits), the rest zero to 10⁻²⁷; residual 1.4·10⁻²⁸. Sector 71/120: i·c on F₁\*, the rest zero to 10⁻²⁷; residual 1.4·10⁻²⁹.
- The other halves G₀, G₁ do not appear, in line with the rule of v5: one-sided sum goes to one-sided sum of the same type.

Not done: first corrections; the companion 𝒦 near r = 1 (expected: the same with G₀, G₁ and no trivial sector of its own beyond −Φ_𝒥; inference); a second unsigned function for F₁\*-type descendants.

## 3. Where the whole thing stands

| statement | status |
|---|---|
| Identities (I)–(VI), odd-root rule, q-difference modules, Casoratians, reformulation through H(z;q) | proved (v2–v6) |
| Identity (VII) | exact to r²⁹⁹, no proof (v8) |
| 𝓑₁…𝓑₄ = state integrals of Φ_b | derivation by residues, 10⁻²⁷…10⁻³¹ (v7) |
| Two-sector law for R₀, R₁ and the bilateral parts of ψ\*, F\* | follows from the previous line plus stationary phase; both steps standard, neither written out (v7) |
| One-sided laws: ψ\* ↔ (1, W, q̃W₃); f\* ↔ F_a, F_b; converse and mirror laws | numerical, 10⁻¹⁷…10⁻³⁰ (v5) |
| F\* ↔ 𝒥, 𝒦 and 𝒥 ↔ F\* | numerical, 10⁻¹⁷ and 10⁻²⁸ (v8, v9) |

What is missing for a theorem is the one-sided half. Every one-sided law found so far has the form "bilateral law (state integral) + a rule that sends n ≥ 0 sums to n ≥ 0 sums". A proof of that rule, probably through a second, inhomogeneous state integral or a Mordell-type integral, would turn the table into theorems.
