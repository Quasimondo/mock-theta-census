# ψ₁\*, ψ₀\*: the theorems and their proofs in one place (v1)

Consolidated from addenda v2–v4, v12, v14–v17, v19–v21, v25, v26 of this repository (27 Nov 2026). Written by Claude Fable 5.2 with Mario Klingemann.

**Status.** Every proof below was written by one AI model. The numerical statements were rechecked by a second model with separate code (`psi_star/v*/check_*.py`). A numerical stress test of the displayed inequalities and constants by the second model (`docs/psi-star-stress-test-v1.md`) found no violation; its findings are answered in addendum v26, included below. The analytic estimates have not been checked line by line by anyone, and no mathematician has refereed them. "Proved" in this document means: a complete argument is written here. §0.5 lists where I think a referee should look first.

**What this document is.** §0 states the results once, in one notation, and gives the map of the proof. Parts 1–7 and the appendix are the proofs, taken unchanged from the addenda and put in logical order (headings demoted and prefixed with their source, e.g. "[v21 §2]", nothing else edited; references such as "§4" inside a section mean that section of the same addendum). The addenda keep their theorem numbers, so "Theorem 6" means the same here and there. Part 1 (from v12) uses a few symbols differently from the rest; §0.6 has the concordance. Numerical checks are not repeated here; each addendum has its table and its scripts.

---

## 0. Results and map

### 0.1 Objects and notation

- **The functions.** For |q| < 1:
  ψ₁\* = Σ_{n≥0} (−1)ⁿ q^{n(n+1)/2} (−q;q)ₙ, ψ₀\* = Σ (−1)ⁿ q^{(n+1)(n+2)/2} (−q;q)ₙ, f₁\* = Σ (−1)ⁿ q^{n²+n}/(−q;q)ₙ, f₀\* = Σ (−1)ⁿ q^{n²}/(−q;q)ₙ.
- **Variables.** x with Re x > 0; b = √x, β = arg b; q = e^{2πix} = e^{−t}, t = −2πix, θ = arg t = 2β − 90°; dual nome p = e^{−πi/x} = e^{−2π²/t}. |q| < 1 means Im x > 0, that is −90° < θ < 90°. Fractional powers are defined through x: q^{α} = e^{2πixα}, p^{α} = e^{−πiα/x} (so p^{n²/2} in Theorem C is e^{−πin²/(2x)}).
- **Faddeev's dilogarithm.** Φ_b(y) = (−q^{1/2}e^{2πby}; q)∞/(−p·e^{2πy/b}; p²)∞ for Im x > 0; for all Re b > 0 it is given by Faddeev's integral (Part 2, Lemma Φ), which also defines it for real x.
- **Integrals.** L = e^{−iδ}·ℝ with 0 < δ < 90° − |β| (δ = 10° from Part 3 on), and for u ∈ ℂ
  I(u) = ∫_L e^{−πiy² + 2πuy}/Φ_b(y) dy,  J(u) = ∫_L e^{−πiy² + 2πuy}/(Φ_b(y)·(1 + e^{2πy/b})) dy.
  Prefactor P = 2^{−1/2}·e^{π²/(12t) + t/6}/b, which equals (−q;q)∞/(b·q^{1/8}·(p;p²)∞) for Im x > 0.
- **Rescaled variables (Parts 3–7).** y = b·s, w = 2π·b·y = i·t·s. The contour is s = r·e^{−iγ}, γ = β + δ, on which w = ±ρ·e^{iα}, α = β − δ.
- **Dual series.** W(p) = Σ (−1)ʲ p^{j²+j} (p;p²)ⱼ, W₃(p) = Σ (−1)ʲ p^{j²+3j} (p;p²)ⱼ, F_a(p) = Σ_{n≥1} p^{2n²−n}/(p;p²)ₙ, F_b(p) = Σ_{n≥1} p^{2n²−3n+1}/(p;p²)ₙ.
- **Perturbative series.** τₙ(t) = q^{n²+n}/(2·(−q;q)ₙ) is analytic at t = 0, and Φ₁(t) = Σ_k aₖtᵏ with aₖ = Σ_{n≥0} (−1)ⁿ·[tᵏ]τₙ (absolutely convergent; a₀ = 1/3, a₁ = 1/9, a₂ = 5/162). Φ₀ is defined in the same way from q^{n²}/(2·(−q;q)ₙ). 𝒮_{θ₀}[Φ](t) = (1/t)·∫₀^{∞·e^{iθ₀}} e^{−ξ/t}·B(ξ) dξ with B(ξ) = Σ aₖξᵏ/k! is the Borel sum in the direction θ₀.
- **Constants.** Cl₂ = Cl₂(π/3) = 1.0149416…; S = π²/3 − i·Cl₂, S̄ its conjugate; A = (π/√3)^{1/2}.

### 0.2 The theorems

**Theorem I (exact formulas; Part 1; v12).** For Im x > 0:

- ψ₁\*(q) = P·[ J(1/(2b)) + I(1/(2b))·p·W₃(p) − I(−1/(2b))·W(p) ]
- f₁\*(q) = 2P·[ J(1/(2b)) − I(1/(2b))·F_b(p) + I(−1/(2b))·F_a(p) ]
- ψ₀\*(q) = q^{1/2}·P·[ J(b + 1/(2b)) + I(b + 1/(2b))·p·W₃(p) − I(b − 1/(2b))·W(p) ]
- f₀\*(q) = 2q^{1/2}·P·[ J(b + 1/(2b)) − I(b + 1/(2b))·F_b(p) + I(b − 1/(2b))·F_a(p) ]

Moreover I(u) = c·e^{−πiu²/2}·D(u) with an explicit double series D that is bilinear in q-series and p-series (Theorem C), so P·I(±1/(2b)) = (κ/2)·𝓑₂,₄ in the notation of v6/v7. Theorems A, B, C of Part 1 give the same for every u.

**Theorem II (the kernel integral is a Borel sum; Parts 2–5; v15–v17).** There are constants A′, t₀, R > 0 such that:

- P·J(1/(2b)) is analytic in t on −210° < θ < 10°, 0 < |t| < t₀ (this includes real x, |q| = 1), and for all K ≥ 0 and all such t, also on the closed sector: |P·J(1/(2b)) − Σ_{k<K} aₖtᵏ| ≤ A′·K!·(240·|t|)ᴷ (Theorem 4).
- For every θ₀ ∈ [−120°, −80°] and Re(e^{iθ₀}/t) > 1/R: P·J(1/(2b)) = 𝒮_{θ₀}[Φ₁](t) (Theorem 5, using the Nevanlinna–Sokal theorem).

The same holds for q^{1/2}·P·J(b + 1/(2b)) and Φ₀.

**Theorem III (the state integrals; Part 6; v19 Theorem 6).** As x → 0 with 0 ≤ arg x ≤ 11°, uniformly:

- P·I(−1/(2b)) = −A·e^{5πi/12}·t^{−1/2}·e^{−S̄/t}·(1 + O(t)) (growing),
- P·I(1/(2b)) = +A·e^{−5πi/12}·t^{−1/2}·e^{−S/t}·(1 + O(t)) (decaying; this one for 0 ≤ arg x ≤ 60°),
- q^{1/2}·P·I(b ± 1/(2b)) = e^{∓2πi/3}·P·I(±1/(2b))·(1 + O(t)).

Equivalently 𝓑₄, 𝓑₃ = (4/3)^{1/4}·e^{∓iπ/12}·e^{−i·Cl₂/t}·(1 + O(t)) and 𝓑₂, 𝓑₁ = (4/3)^{1/4}·e^{±iπ/12}·e^{+i·Cl₂/t}·(1 + O(t)).

**Theorem IV (the forward law inside the disk; Part 6; v19 Theorem 8).** For 0 < arg x ≤ 10° and |t| < R, however x tends to 0 in this sector:

  ψ₁\* = 𝒮[Φ₁](t) + A·t^{−1/2}·[ e^{5πi/12}·e^{−S̄/t}·W(p)·(1 + O(t)) + e^{−5πi/12}·e^{−S/t}·p·W₃(p)·(1 + O(t)) ],

and the corresponding formulas for f₁\*, ψ₀\*, f₀\* (Part 6, §6.1). Here 𝒮 is the Borel sum in the direction of t.

**Theorem V (roots of unity; Part 7; v20, v21).** Let q₀ = e^{2πi/X}, x = 1/X, and let the integrals be taken with Faddeev's Φ_b at real b.

- X ≥ 2 even (X = 2 by Theorem 18 of v25): the finite sum ψ₁\*(q₀) = Σ_{n<X/2} (−1)ⁿ q₀^{n(n+1)/2} (−q₀;q₀)ₙ equals P·[J(1/(2b)) + I(1/(2b)) − I(−1/(2b))] (Theorem 9).
- X ≥ 1 odd (X = 1, 3 by Theorem 18 of v25): lim_{r→1⁻} ψ₁\*(r·q₀) = A_X/3 = P·[J(1/(2b)) − ⅓·I(1/(2b)) − ⅓·I(−1/(2b))], with A_X the sum of the first X terms (Theorem 11; the first equality is the odd-root rule of v3).
- As X → ∞: the radial value is w_X·3^{−1/4}·(X/2)^{1/2}·e^{Cl₂·X/2π}·e^{2πi/3 − iπX/6}·(1 + O(1/X)) with w_X = 1 for even X and ⅓ for odd X (Theorems 10, 12).
- The same for ψ₀\* with the extra factor e^{2πi/3} in the asymptotics.

### 0.3 Map of the proof

1. **Part 1 → Theorem I.** Residues at the kernel poles on one side give the one-sided q-series; on the moved contour the kernel is a geometric series of kernel-free integrals I(u ∓ m/b) (Theorem A). These obey a three-term recursion whose generating sums are W, W₃ or F_a, F_b (Theorem B). I(u) is evaluated by two Gaussian linearisations and the Fourier transform of Φ_b (Theorem C).
2. **Part 2: properties of Φ_b.** Faddeev's integral equals both product formulas (Lemma Φ). In the strip |Im w| < π it has an expansion uniform in x up to the real axis (Lemma S) and a leading-term bound up to the edge (Lemma F). Outside the strip the functional equation brings the argument back at negligible cost on the contours used (Lemma D, Lemma B).
3. **Part 3: descent.** Along the straight contour with δ = 10° the leading exponent of the kernel integrand decreases linearly, with margin κ ≥ 0.185, for −210° ≤ θ ≤ 10° (Proposition 1).
4. **Part 4: uniform bound.** |H| ≤ C·e^{−κ|s|} for the integrand H of P·J on the whole contour, uniformly on the sector: away from the real x-axis by Euler–Maclaurin for the products (v14 Lemma 1, v15 Lemmas 1′–3′), near it by Part 2 (v16 §4). Analyticity of P·J on the sector follows (the argument is written out in [v26 §1]).
5. **Part 5 → Theorem II.** Shifting the contour across N kernel poles gives P·J = Σ_{n<N} (−1)ⁿτₙ + (moved integral), and the functional equation turns the moved integrand into H·D_N with |D_N| ≤ e^{−0.544·N + 0.149·|s|}. With N = ⌊0.04/|t|⌋ and Cauchy estimates for τₙ this gives the factorial bounds; Nevanlinna–Sokal gives the Borel sum.
6. **Part 6 → Theorems III, IV.** Explicit descent contours for the kernel-free integrals with closed-form heights, Lemma B, and Laplace's method. Theorem IV is Theorem I + II + III.
7. **Part 7 → Theorem V.** Radial limits: of the q-series (Lemma R1; odd-root rule of v3, Appendix), of W, p·W₃ at p = ±1 (Lemmas R2, R2′, the latter through identities (V), (VI) of v4, Appendix), continuity of the integrals at real x (Lemma R3; for x₀ ≤ 1 in [v25 §1]). Then Theorem III at real x.

Only the uniform bound of Part 4 and Part 2 are used from the "asymptotic expansion" addenda v14–v16; the expansions proved there (Theorems 1–3) are superseded by Theorem 4 and are needed only for the remark aₖ = coefficients of v13.

### 0.4 What is quoted from the literature

- The Nevanlinna–Sokal theorem (A. D. Sokal, J. Math. Phys. 21 (1980) 261–263). Its statement was taken from later papers that quote it, not from the original. Watson's original theorem would also suffice (the sector has opening 220°).
- Classical facts used without proof: the Jacobi triple product and the transformation of η and of the theta functions; Euler's and the q-binomial series; the Fourier integral ∫ e^{−iuσ}/(2σ·sinh πσ) dσ = −log(1 + e^{u}); the product forms of the Rogers–Ramanujan functions; Stirling-type bounds; "Laplace's method" in the form of Lemma L of Part 6, whose proof is only sketched.
- For the corollary on the figure-eight knot only: the asymptotics of the Kashaev invariant of 4₁ (quoted from memory).

### 0.5 Where a referee should look first

A numerical stress test of every displayed inequality, identity and uniform constant (`docs/psi-star-stress-test-v1.md`, Claude Opus 5.5, 27 Nov 2026) found no violation. Its findings were places where a statement was incomplete; they are answered in addendum v26, which is included below (in Parts 4 and 6 and in the last part). The list, with what the stress test and v26 changed:

1. **Part 1, §4** (Theorem C): two applications of Fubini with complex Gaussians, and the normalisation of the Fourier transform of Φ_b, which was not compared with the literature. (Tested numerically; the closing lemma of v12 §4.1 was not tested.)
2. **Part 3** (Proposition 1): the exact formulas for ψ/ρ and the bounds for X, Y. Everything after it rests on κ ≥ 0.185. (The measured rate is ≥ 0.97; the constant 0.185 is sharp at θ = 10° by construction.)
3. **Part 4**: the Euler–Maclaurin remainder (v14 Lemma 1) for complex step, and the constants "uniform on the sector" C, t₀. (Now measured: C = 0.60 for |t| ≤ 0.4; see [v26 §5]. The proof of analyticity is restated in [v26 §1].)
4. **Part 2, Lemma F**: the contour shift to Im σ = ½ and the bound for r₁ there; Lemma D's margin. (The numerical remark after Lemma F is corrected in [v26 §4].)
5. **Part 5**: the justification of the contour shift (Lemma E) at real x; the margin 0.175 − 0.149 between the decay of H and the growth of D_N. (Measured: about 0.97 against 0.06.)
6. **Part 6**: Lemma L is the standard argument, not written in full; it holds for |t| ≤ t₁, and for the integrals with the factor e^{w} the descent contours need |t| < π/3 and |t| < π at real x ([v26 §2]). The statements "the integrand tends to 0 at infinity in the region swept" are brief. Theorem 7 now has ε(a) = 60° ([v26 §3]).
7. **Part 7**: Lemma R3 is brief; its range is extended to x₀ ≤ 1 in [v25 §1], which removes the edge case X = 4 and the restriction X ≥ 4, X ≥ 5.
8. **Appendix**: the odd-root rule uses identity (II) of v2 and standard growth bounds for θ₄ and R₁.

Not tested by the stress test: the closing lemma of v12 §4.1, Lemma L and Lemma E as separate statements, Theorem 8 as a whole (only its ingredients), and the quoted theorems.

### 0.6 Concordance

| here | source | remark |
|---|---|---|
| Part 1 | `addendum-v12.md` §1–§5 | Theorems A, B, C, Corollary |
| Part 2 | `addendum-v16.md` §1–§3; `addendum-v19.md` §2 | Lemmas Φ, S, F, D; Lemmas B, L |
| Part 3 | `addendum-v15.md` §1 | Proposition 1 |
| Part 4 | `addendum-v14.md` §2.1–§2.3; `addendum-v15.md` §2 up to Lemma 3′; `addendum-v16.md` §4 | uniform bound, analyticity |
| Part 5 | `addendum-v17.md` §1–§5 | Theorems 4, 5 |
| Part 6 | `addendum-v19.md` §1, §3–§6 | Theorems 6, 7, 8 |
| Part 7 | `addendum-v20.md` §1–§3; `addendum-v21.md` §1–§3 | Theorems 9–12 |
| Appendix | `addendum-v2.md` §2; `addendum-v3.md` §1; `addendum-v4.md` §1 | identities (I), (II); odd-root rule; identities (V), (VI) |

References inside the proofs such as "v16 §4" or "Lemma F of v16" point to the sources in this table.

**Symbols that differ in Part 1 (v12).**

| Part 1 | rest of the document |
|---|---|
| ξ = y/b | s |
| v = e^{2πby} | e^{w} |
| w = e^{2πy/b} (and, in §4, a Fourier variable) | w′ = e^{2πs} in Parts 4–5; elsewhere w = 2π·b·y |
| t = p·e^{2πiu/b} (dual variable of Theorem B) | T in v18; elsewhere t = −2πix |
| c = 2^{−1/2}e^{−iπ/4}q^{−1/24}p^{1/12} | c, C also denote generic constants |
| κ (constant of v10 in the Corollary) | κ ≥ 0.185 is the descent margin of Part 3 |
| Ψ(−1), q·Ψ(−q), Ψ̌(−1), Ψ̌(−1/q) | ψ₁\*, ψ₀\*, f₁\*, f₀\* |

In Parts 4–5 the letters A, B, C, K are also used for constants; Theorems A, B, C always refer to Part 1.

---
