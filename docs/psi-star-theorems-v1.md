# ψ₁\*, ψ₀\*: the theorems and their proofs in one place (v1)

Consolidated from addenda v2–v4, v12, v14–v17, v19–v21, v25, v26, v27 of this repository (27 Nov 2026). Written by Claude Fable 5.2 with Mario Klingemann.

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

1. **Part 1, §4** (Theorem C): two applications of Fubini with complex Gaussians, and the normalisation of the Fourier transform of Φ_b, which was not compared with the literature. (Tested numerically in both rounds; in the closing lemma of v12 §4.1 the constant c is the exact leading constant, and its remark on the factor is reworded in [v27 §4].)
2. **Part 3** (Proposition 1): the exact formulas for ψ/ρ and the bounds for X, Y. Everything after it rests on κ ≥ 0.185. (The measured rate is ≥ 0.97; the constant 0.185 is sharp at θ = 10° by construction.)
3. **Part 4**: the Euler–Maclaurin remainder (v14 Lemma 1) for complex step, and the constants "uniform on the sector" C, t₀. (Now measured: C = 0.60 for |t| ≤ 0.4; see [v26 §5]. The proof of analyticity is restated in [v26 §1].)
4. **Part 2, Lemma F**: the contour shift to Im σ = ½ and the bound for r₁ there; Lemma D's margin. (The numerical remark after Lemma F is corrected in [v26 §4].)
5. **Part 5**: the justification of the contour shift (Lemma E) at real x; the margin 0.175 − 0.149 between the decay of H and the growth of D_N. (Measured: about 0.97 against 0.06.)
6. **Part 6**: Lemma L is the standard argument, not written in full; it holds for |t| ≤ t₁, and for the integrals with the factor e^{w} the descent contours need |t| < π/3 and |t| < π at real x ([v26 §2]). The statements "the integrand tends to 0 at infinity in the region swept" are brief. Theorem 7 now has ε(a) = 60° ([v26 §3]).
7. **Part 7**: Lemma R3 is brief; its range is extended to x₀ ≤ 1 in [v25 §1], which removes the edge case X = 4 and the restriction X ≥ 4, X ≥ 5.
8. **Appendix**: the odd-root rule uses identity (II) of v2 and standard growth bounds for θ₄ and R₁.

A second round (`docs/psi-star-stress-test-v2.md`) tested what the first had left: the closing lemma of v12 §4.1, Lemma L on a model, Theorem 8 as a whole for ψ₁\* (uniform up to |p| = 0.77), and the corrections of v26; again no violation. Still not tested: Lemma E as a separate statement, Theorem 8 for f₁\*, ψ₀\*, f₀\*, and the quoted theorems.

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

## Part 1. The kernel identity (Theorem I)

Local notation: see §0.6. Statements in §1, proofs in §2–§5.

*From `docs/addendum-v12.md`, §1–§5; text unchanged.*

### [v12 §1] Notation and statements

x in the upper half-plane, b = √x with β = arg b ∈ (0, π/2); q = e^{2πix}, p = e^{−πi/x}, q̃ = p²; c_b = i(b + 1/b)/2. For a point y write ξ = y/b, v = e^{2πby} = e^{2πxξ}, w = e^{2πy/b} = e^{2πξ}. As in v7,

  Φ_b(y) = (−q^{1/2}v; q)∞ / (−p·w; p²)∞.

Fix δ with 0 < δ < π/2 − β and let L = e^{−iδ}·ℝ. For u ∈ ℂ put

  g_u(y) = e^{−πi·y² + 2π·u·y} / Φ_b(y),  **I(u) = ∫_L g_u(y) dy,  J(u) = ∫_L g_u(y) / (1 + e^{2πy/b}) dy.**

The integral of v11 is J(1/(2b)). Further

- Ψ(z) = Σ_{n≥0} q^{n(n+1)/2}·zⁿ·(−q;q)ₙ, Ψ̌(z) = Σ_{n≥0} q^{n²+n}·zⁿ/(−q;q)ₙ, so ψ₁\* = Ψ(−1), ψ₀\* = q·Ψ(−q), f₁\* = Ψ̌(−1), f₀\* = Ψ̌(−1/q);
- 𝒲(t) = Σ_{k≥0} (−t)ᵏ·p^{k²+k}·(p;p²)ₖ, so W(p) = 𝒲(1), W₃(p) = 𝒲(p²);
- 𝒢₀(t) = Σ_{k≥0} t^{−k}·p^{2k²+k}/(p;p²)ₖ, 𝒢₁(t) = Σ_{k≥0} t^{−k}·p^{2k²+3k}/(p;p²)ₖ₊₁, so 𝒢₀(1) = F_b(p) − F_a(p), 𝒢₁(1) = F_a(p)/p;
- D(u) = Σ_{k,n≥0} (−1)^{kn+n}·e^{πiu(bk + n/b)}·q^{k²/4}·p^{n²/2+n} / ((q;q)ₖ·(p²;p²)ₙ);
- c = 2^{−1/2}·e^{−iπ/4}·q^{−1/24}·q̃^{1/24};
- P = (−q;q)∞/(b·q^{1/8}·(p;p²)∞) (prefactor of v11), κ as in v10, 𝓑₁…𝓑₄ as in v6.

**Theorem A (kernel integral = one-sided sum + dual sum of kernel-free integrals).** For every u ∈ ℂ, with both sums over m absolutely convergent,

  (A↑) J(u) = −i·b·q^{1/8}·e^{πiub}·((p;p²)∞/(−q;q)∞)·Ψ(e^{2πiub}) − Σ_{m≥1} (−1)ᵐ·I(u − m/b),

  (A↓) J(u) = +i·b·q^{1/8}·e^{−πiub}·((p;p²)∞/(2(−q;q)∞))·Ψ̌(e^{−2πiub}) + Σ_{m≥0} (−1)ᵐ·I(u + m/b).

**Theorem B (the dual sums).** For every u ∈ ℂ:

  (B1) I(u) = p^{−1}·e^{2πiu/b}·[ I(u + 1/b) + p^{−1}·I(u + 2/b) ], and likewise I(u) = q^{1/2}·e^{2πiub}·[ I(u + b) + q^{1/2}·I(u + 2b) ].

  (B2) With t = p·e^{2πiu/b}:

  - Σ_{m≥1} (−1)ᵐ·I(u − m/b) = I(u)·[𝒲(t/p²) − 1] − (t/p)·I(u + 1/b)·𝒲(t),
  - Σ_{m≥0} (−1)ᵐ·I(u + m/b) = I(u)·𝒢₀(t) − I(u + 1/b)·𝒢₁(t).

  (B3) If t = 1, that is u ∈ 1/(2b) + b·ℤ:

  - Σ_{m≥1} (−1)ᵐ·I(u − m/b) = I(u)·p·W₃(p) − I(u − 1/b)·W(p),
  - Σ_{m≥0} (−1)ᵐ·I(u + m/b) = I(u)·F_b(p) − I(u − 1/b)·F_a(p).

**Theorem C (the kernel-free integral).** For every u ∈ ℂ,

  I(u) = c·e^{−πi·u²/2}·D(u).

In particular I(1/(2b)) = c·p^{1/8}·𝓑₂, I(−1/(2b)) = c·p^{1/8}·𝓑₄, I(b + 1/(2b)) = −i·c·q^{−1/4}·p^{1/8}·𝓑₁, I(b − 1/(2b)) = +i·c·q^{−1/4}·p^{1/8}·𝓑₃, and c·p^{1/8} = κ/(2P).

**Corollary (the four exact formulas).**

- ψ₁\*(q) = P·J(1/(2b)) + (κ/2)·[ 𝓑₂·p·W₃(p) − 𝓑₄·W(p) ]   (the identity of v11)
- f₁\*(q) = 2P·J(1/(2b)) − κ·[ 𝓑₂·F_b(p) − 𝓑₄·F_a(p) ]
- ψ₀\*(q) = P·q^{1/2}·J(b + 1/(2b)) − i·q^{1/4}·(κ/2)·[ 𝓑₁·p·W₃(p) + 𝓑₃·W(p) ]
- f₀\*(q) = 2P·q^{1/2}·J(b + 1/(2b)) + i·q^{1/4}·κ·[ 𝓑₁·F_b(p) + 𝓑₃·F_a(p) ]

All of §1 is **proved** (§2–§5) and **verified** (§6).

---

### [v12 §2] Proof of Theorem A

**Bounds for 1/Φ_b.** Two forms are used.

- (2.1) 1/Φ_b(y) = (−p·w; p²)∞/(−q^{1/2}v; q)∞. Where |w| ≤ R and |q^{1/2}v| ≤ r < 1, it is bounded by (−|p|R; |p|²)∞/(r; |q|)∞.
- (2.2) The inversion relation Φ_b(y)·Φ_b(−y) = Φ_b(0)²·e^{πi·y²}, Φ_b(0)² = q^{1/24}·q̃^{−1/24} (classical: Jacobi triple product and the theta transformation; checked in §6), gives 1/Φ_b(y) = e^{−πi·y²}·Φ_b(−y)/Φ_b(0)², and Φ_b(−y) = (−q^{1/2}/v; q)∞/(−p/w; p²)∞ is bounded where |1/v| ≤ R and |p/w| ≤ r < 1.

So g_u = O(|e^{−πiy² + 2πuy}|) in regions of the first kind, which include the left sector {Re(by) ≤ 0, Re(y/b) ≤ 0} (arg y from π/2 + β to 3π/2 − β), and g_u = O(|e^{−2πiy² + 2πuy}|) in regions of the second kind, which include the right sector (arg y from −π/2 + β to π/2 − β). Since |e^{−πiy²}| = exp(π·|y|²·sin(2·arg y)), both ends of L lie in directions of Gaussian decay, and I(u), J(u) converge for all u.

**Poles.** The kernel has simple poles at yₙ = i·b·(n + ½), n ∈ ℤ, that is ξ = i(n + ½), with residue −b/(2π). 1/Φ_b has poles at −c_b − i·b·k − i·n/b (k, n ≥ 0), all below L, with Re ξ ≤ −Im x/(2|x|²). The kernel poles with n ≥ 0 are above L, those with n < 0 below.

**Upper route.** Choose ε > 0 with e^{2π|Re x|ε}·|q|^{1/2} < 1 and let L′ be the path that comes down the line ξ = ε + iτ (τ from +∞ to 0) and continues along b·ε + e^{−iδ}·s (s ≥ 0). On the strip 0 ≤ Re ξ ≤ ε, Im ξ ≥ 0 we have |w| ≤ e^{2πε} and |q^{1/2}v| ≤ e^{2π|Re x|ε}|q|^{1/2} < 1, so (2.1) applies there and in the left sector. Close the region between the left ray of L and L′ by the segment Im ξ = N (N a positive integer). On that segment w = e^{2π·Re ξ} > 0, so |1/(1 + w)| ≤ 1, and arg y stays in [π/2 + β − o(1), π − δ], so |g_u| ≤ C·e^{−c·N²}. Between the two right rays there are no poles and Gaussian decay. Hence

  J(u) = 2πi·Σ_{n≥0} Res_{yₙ} + ∫_{L′} g_u/(1 + w) dy.

At yₙ: w = −1, v = q^{n+½}, so 1/Φ_b = (p;p²)∞·(−q;q)ₙ/(−q;q)∞, and e^{−πiyₙ² + 2πuyₙ} = q^{(n+½)²/2}·e^{2πiub(n+½)}. The residue sum is the first term of (A↑).

On L′ we have Re ξ ≥ ε, so 1/(1 + w) = Σ_{m≥1} (−1)^{m−1}·w^{−m} with Σ|w|^{−m} ≤ 1/(e^{2πε} − 1). Since ∫_{L′}|g_u||dy| < ∞, dominated convergence gives ∫_{L′} g_u/(1 + w) = Σ_{m≥1} (−1)^{m−1}·∫_{L′} g_u·w^{−m}, absolutely convergent. Finally g_u·w^{−m} = g_{u−m/b} has no poles between L′ and L, and the same bounds as above let the contour go back: ∫_{L′} g_{u−m/b} = I(u − m/b). This is (A↑). ∎

**Lower route.** Mirror image. Choose ε < Im x/(2|x|²), so that the strip −ε ≤ Re ξ ≤ 0, Im ξ ≤ 0 contains no pole of 1/Φ_b and |p/w| ≤ |p|e^{2πε} < 1, |1/v| ≤ e^{2π|Re x|ε} there; (2.2) applies in the strip and in the right sector. L″ runs along −b·ε + e^{i(π−δ)}·s (from s = ∞) and then down the line ξ = −ε − iτ. On the closing segments Im ξ = −N again w > 0. The contour is moved clockwise, so J(u) = −2πi·Σ_{n<0} Res + ∫_{L″}. At y = −i·b·(j + ½): v = q^{−j−½}, (−q^{−j}; q)∞ = 2·q^{−j(j+1)/2}·(−q;q)ⱼ·(−q;q)∞, which gives the terms q^{j²+j}/(−q;q)ⱼ of Ψ̌ and the factor ½. On L″, |w| ≤ e^{−2πε} and 1/(1 + w) = Σ_{m≥0} (−1)ᵐwᵐ. This is (A↓). ∎

Remark. Subtracting the two routes shows that neither one-sided sum alone is the integral: (A↑) − (A↓) says residue sum(n ≥ 0) + residue sum(n < 0) = Σ_{m∈ℤ} (−1)ᵐ·I(u + m/b), a bilateral dual sum.

---

### [v12 §3] Proof of Theorem B

**(B1).** From Φ_b(y + i/b) = Φ_b(y)/(1 + p^{−1}w) one gets g_u(y + i/b) = p^{−1}·e^{2πiu/b}·w·(1 + p^{−1}w)·g_u(y). There is no pole between L and L + i/b and the decay is Gaussian in the strip, so ∫_L g_u(y + i/b) dy = I(u); and w·g_u = g_{u+1/b}. The second relation follows in the same way from Φ_b(y + ib) = Φ_b(y)/(1 + q^{1/2}v). ∎

**(B2).** Put aₘ = I(u − m/b). (B1) at u − m/b reads aₘ = t·p^{2m−2}·aₘ₋₁ + t·p^{2m−3}·aₘ₋₂ for all m ∈ ℤ. By Theorem A, Σ_{m≥1}|aₘ| < ∞, so Ã(z) = a₀ + Σ_{m≥1} aₘzᵐ is continuous on |z| ≤ 1, and the recursion gives

  Ã(z) = t·z·(1 + pz)·Ã(p²z) + a₀ + t·p^{−1}·z·a₋₁.

Iterating (the remainder t^K·z^K·p^{K(K−1)}·(−pz;p²)_K·Ã(p^{2K}z) tends to 0, because Ã is bounded on the disk):

  Ã(z) = Σ_{k≥0} tᵏzᵏp^{k(k−1)}·(−pz; p²)ₖ·[ a₀ + t·p^{2k−1}·z·a₋₁ ].

At z = −1 this is a₀·𝒲(t/p²) − (t/p)·a₋₁·𝒲(t), which is the first formula. For the second put bⱼ = I(u + j/b) = a₋ⱼ, with Σ_{j≥0}|bⱼ| < ∞ again by Theorem A; the recursion solved for the lowest index is bⱼ = −p·bⱼ₋₁ + t^{−1}p^{2j−1}·bⱼ₋₂, so B(z) = Σ_{j≥0} bⱼzʲ satisfies B(z) = b₀ + b₁z/(1 + pz) + t^{−1}p³z²·B(p²z)/(1 + pz), and iteration at z = −1 gives b₀·𝒢₀(t) − b₁·𝒢₁(t). ∎

**(B3).** For t = 1: 𝒲(p^{−2}) = Σ (−1)ᵏp^{k²−k}(p;p²)ₖ = 1 − W + p·W₃ (shift k → k + 1 and use (p;p²)ₖ₊₁ = (p;p²)ₖ(1 − p^{2k+1})), and a₁ = a₀ + p^{−1}a₋₁. So a₀·(p·W₃ − W) − p^{−1}a₋₁·W = a₀·p·W₃ − a₁·W. In the same way 𝒢₀(1) = F_b − F_a, 𝒢₁(1) = F_a/p, and b₀(F_b − F_a) − p^{−1}b₁F_a = I(u)·F_b − I(u − 1/b)·F_a. ∎

So W, W₃, F_a, F_b are the four generating sums of one three-term recursion, read in its two directions.

---

### [v12 §4] Proof of Theorem C

The integral I(u) cannot be closed in either half-plane (see the remark after the closing lemma). The integral of v7 can, and I(u) is a Fourier transform away from it.

#### [v12 §4.1] Closing lemma (this also completes v7)

Let C be the contour of v7 with both angles equal to δ: from ∞·e^{i(π−δ)} through 0 to ∞·e^{iδ}. Consider ∫_C Φ_b(y)·e^{−Aπi·y² + 2πλy} dy with A = ½ (any λ), or with A = 0 and 2πλ = 2πi·w where Im(w·e^{−iδ}) < 0 and Re(w·b) > 0.

**Claim.** The integral equals 2πi times the sum of the residues at all poles c_b + i·b·k + i·n/b (k, n ≥ 0), summed first over k, then over n.

Proof. In the variable ξ the poles are at Re ξ = (n + ½)·Im x/|x|², Im ξ = k + ½ + (n + ½)·Re x/|x|². Close C by the vertical line Re ξ = σ_N := N·Im x/|x|², which runs between the pole columns n = N − 1 and n = N, and by a horizontal line Im ξ = T₀ + K (K integer, T₀ fixed so that the line meets no pole with n < N).

- On Re ξ = σ_N, |w| = |p|^{−2N}, and |(−p·w; p²)∞| ≥ Π_j |1 − |p|^{2j+1−2N}| = |p|^{−N²}·(|p|;|p|²)_N·(|p|;|p|²)∞. So |1/(−p·w; p²)∞| ≤ C·|p|^{N²} = C·exp(−π·σ_N²·|x|²/Im x), uniformly on the whole line.
- For all y, |(−q^{1/2}v; q)∞| ≤ C·exp(π·(σ·Re x − τ·Im x)₊²/Im x), ξ = σ + iτ (triple product: (−a·r^{1/2}; r)∞ ≤ C_r·exp((log a)²/(2·log(1/r)))).
- |e^{−Aπiy²}| = exp(Aπ·Im(xξ²)) = exp(Aπ·[Im x·(σ² − τ²) + 2στ·Re x]).

Adding the exponents on the vertical line, with τ\* = σ_N·Re x/Im x: for τ ≥ τ\* the total is at most −π(1 − A)·σ_N²·|x|²/Im x − πA·Im x·(τ − τ\*)²; for τ ≤ τ\* (down to the point where the line meets C) it is −π(1 − A)·Im(y²). There arg y ∈ [δ, π/2 − β], and Im(y²) = Im(xξ²) is concave in τ, so its minimum is at an end of the segment: the total is at most −c·σ_N² with c = π(1 − A)·min(|x|²/Im x, |x|·sin 2δ/cos²(δ − β)). (The second value belongs to the point on C; `t5_bounds.py` confirms that the supremum on the line is this bound times a factor between 1 and 1.5.) The linear factor costs e^{O(σ_N + |τ|)} for A = ½, and for A = 0 it supplies the decay e^{−2πτ·Re(wb)} as τ → ∞. So the integral over the vertical line tends to 0 as N → ∞ whenever A < 1. On the horizontal line with N fixed, 1/(−p·w;p²)∞ is bounded by periodicity in ξ → ξ + i, |v| is bounded, and the integrand is O(e^{−cK²}) for A = ½, O(e^{−cK}) for A = 0. Let K → ∞, then N → ∞. ∎

Consequences:

- **v7 is proved**: I_½(v) := ∫_C Φ_b(s)·e^{−πi·s²/2 + 2πvs} ds = −i·b·((q)∞/(q̃)∞)·e^{−πi·c_b²/2 + 2πv·c_b}·Σ_{k,n} (−1)^{kn}·Zᵏq^{k²/4}/(q)ₖ·Z̃ⁿq̃^{n²/4}/(q̃)ₙ, Z = i·e^{2πivb}·q^{1/4}, Z̃ = −i·e^{2πiv/b}·q̃^{1/4}, for all v. (Residues as computed in v7.)
- **Fourier transform of Φ_b**: F(w) := ∫_C Φ_b(y)·e^{2πiwy} dy = −i·b·((q)∞/(q̃)∞)·e^{2πi·c_b·w}/Φ_b(−w − c_b). The residue sum is two q-binomial series. Both sides are analytic in w on the half-plane H = {Im(w·e^{−iδ}) < 0}, where the integral converges, so the formula holds on all of H.
- **Fourier transform of 1/Φ_b**: G(w) := ∫_L e^{2πiwy}/Φ_b(y) dy = (i/b)·((q̃)∞/(q)∞)·e^{−2πi·c_b·w}·Φ_b(w + c_b) on H. Same proof, closing downward with the lines Re(by) = N·Im x, on which |1/(−q^{1/2}v; q)∞| ≤ C·|q|^{N²/2}; the roles of the two products are exchanged.

These are the classical Fourier transform formulas for the quantum dilogarithm (Faddeev–Kashaev–Volkov). I did not compare normalisations with the literature; the proof above is included so that complex b is covered.

Remark. The estimate needs A < 1 for Φ_b·e^{−Aπiy²}. (As far as I remember, without having re-read it, this is the condition B > A of Garoufalidis–Kashaev for B = 1.) I(u) is, by (2.2), an integral of Φ_b·e^{−2πiy²}, so A = 2: it cannot be closed, which is why v11 found that "closing naively" fails and why step (C) needs a different tool.

#### [v12 §4.2] Two Gaussian linearisations

Let Σ_η = {e^{iδ}(ρ − iη) : ρ ∈ ℝ} with η > 0 large, and K(λ) = ∫ Φ_b(z)·e^{πi·z² + 2πλz} dz over any line of direction e^{iδ} below all poles of Φ_b (left end in the left sector at angle π + δ, right end at angle δ; both are directions of Gaussian decay).

**(X)** For all y, e^{−πiy²} = e^{−iπ/4}·∫_{Σ_η} e^{πiσ² + 2πiσy} dσ. Insert this into I(u). For σ ∈ Σ_η and y = r·e^{−iδ} ∈ L, |e^{2πi(σ − iu)y}| = e^{2πr(η + Re(u·e^{−iδ}))}, independent of ρ; |e^{πiσ²}| ≤ C·e^{−π·sin 2δ·ρ² + 2πη|ρ|}; |1/Φ_b| is bounded on the left ray and O(e^{−π·sin 2δ·r²}) on the right ray. For η > −Re(u·e^{−iδ}) the double integral converges absolutely, so by Fubini

  I(u) = e^{−iπ/4}·∫_{Σ_η} e^{πiσ²}·G(σ − iu) dσ = e^{−iπ/4}·(i/b)·((q̃)∞/(q)∞)·e^{πi(iu − c_b)² + 2πi·c_b²}·K(−u − 2i·c_b)

(substitute z = σ − iu + c_b).

**(Y)** For all s, e^{−πis²/2} = √2·e^{−iπ/4}·∫_{Σ_η} e^{2πiσ² + 2πiσs} dσ. Insert this into I_½(v). On the left ray of C the factor |e^{2πi(σ − iv)s}| decays, independently of ρ; on the right ray s = r·e^{iδ} the three Gaussian exponents add up to −π·sin 2δ·[ρ² + (|ρ| − r)²] plus linear terms. Again Fubini applies, and with F from §4.1 and the inversion relation,

  I_½(v) = √2·e^{−iπ/4}·(−i·b)·((q)∞/(q̃)∞)·Φ_b(0)^{−2}·e^{2πi(iv − c_b)² − 2πi·c_b²}·K(−2v − i·c_b).

**Combining.** The two K's coincide for v = u/2 + i·c_b/2 (the two lines are parallel, below all poles, with decay between them). Dividing, and inserting the formula for I_½:

  I(u) = (i/(√2·b))·((q̃)∞/(q)∞)·Φ_b(0)²·e^{πi·c_b²}·e^{−πi·u²/2}·D(u),

because for v = u/2 + i·c_b/2 one has Z = e^{πiub}, Z̃ = −p·e^{πiu/b}, which turns the double sum of v7 into D(u). With (q)∞/(q̃)∞ = q^{−1/24}q̃^{1/24}/√(−ix), √(−ix) = e^{−iπ/4}b, Φ_b(0)² = q^{1/24}q̃^{−1/24} and e^{πi·c_b²} = −i·q^{−1/8}q̃^{1/8}, the constant is c. ∎

**Special values.** At u = ±1/(2b), b ± 1/(2b) the parity split of D(u) in k gives 𝓑₂, 𝓑₄, 𝓑₁, 𝓑₃ exactly as in v7 (the map v = u/2 + i·c_b/2 sends these u to u₂, u₄, u₁, u₃ of v7). And c·p^{1/8} = κ/(2P): from the definitions, κ/(2P) = (b/(4i·q^{1/8}))·θ₄(q)·(p;p²)∞²/((−q;q)∞²·ψ(p)), and by the eta transformations θ₄(q)/(−q;q)∞² · (p;p²)∞²/ψ(p) = (q)∞/(−q)∞³ · (p;p²)∞³/(p²;p²)∞ = 2^{3/2}·q^{1/12}·p^{5/24}/√(−ix). Hence κ/(2P) = 2^{−1/2}·e^{−iπ/4}·q^{−1/24}·p^{5/24} = c·p^{1/8}. So κ of v10 is exactly the constant the proof produces.

---

### [v12 §5] The corollary, and consistency

- **ψ₁\*.** (A↑) at u = 1/(2b): e^{πiub} = i, e^{2πiub} = −1, so the residue term is ψ₁\*/P. By (B3) and Theorem C the dual sum is c·p^{1/8}·[𝓑₂·p·W₃ − 𝓑₄·W] = (κ/(2P))·[…]. This is the identity of v11. ∎
- **f₁\*.** (A↓) at the same u: residue term f₁\*/(2P), dual sum (κ/(2P))·[𝓑₂F_b − 𝓑₄F_a].
- **ψ₀\*, f₀\*.** u = b + 1/(2b): e^{2πiub} = −q, e^{πiub} = i·q^{1/2}, t = 1 again; I(u), I(u − 1/b) are ∓i·c·q^{−1/4}p^{1/8} times 𝓑₁, 𝓑₃.
- **Consistency.** Subtracting the first two formulas: 2ψ₁\* − f₁\* = κ·[𝓑₂·(p·W₃ + F_b) − 𝓑₄·(W + F_a)] = κ·ψ(p)·[𝓑₂·P₋(p) − 𝓑₄·P₊(p)] by (V), (VI). By the inversion formula of v6 (i·q^{1/4}N₁ = ½(𝓑₂E₋ − 𝓑₄E₊)) the right side is θ₄R₁. So the contour argument reproduces identity (II) of v2; the same with (I) for ψ₀\*, f₀\*. Nothing new, but it ties the two closings, (V), (VI), v2 and v6 together.
- **tanh form.** Since 1/(1 + w) = ½(1 − tanh(πy/b)), J = ½I − ½Z with Z(u) = ∫_L tanh(πy/b)·g_u dy, and the first formula becomes ψ₁\* = −(P/2)·Z(1/(2b)) + (κ/2)·[𝓑₂·(p·W₃ + ½) − 𝓑₄·W].

---

## Part 2. Faddeev's dilogarithm

Lemma Φ and the strip lemmas S, F, D (v16), then the bounds in the form used for the state integrals (Lemma B) and Laplace's method (Lemma L) from v19.

*From `docs/addendum-v16.md`, §1–§3; text unchanged.*

### [v16 §1] Φ_b by Faddeev's integral

For Re b > 0 and |Im z| < Re(b + 1/b)/2 let

  Φ_b(z) := exp( ∫_{ℝ+i0} e^{−2izξ} / (4·sinh(bξ)·sinh(ξ/b)) · dξ/ξ ).

**Lemma Φ.**

- (i) For Im b² > 0 this is the product formula of v12. For Im b² < 0 it is Φ_{1/b}, the function used in v13 and v15.
- (ii) It is analytic in b on Re b > 0. Hence the two functional equations Φ_b(z − i·b^{±1}/2) = (1 + e^{2π·b^{±1}·z})·Φ_b(z + i·b^{±1}/2) and the inversion relation Φ_b(z)·Φ_b(−z) = Φ_b(0)²·e^{πiz²}, Φ_b(0)² = e^{−t/24 + π²/(6t)}, hold for all Re b > 0, and Φ_b continues meromorphically in z.

Proof. (i) Let Ψ denote the exponential of the integral and let Im b² > 0.

- The exponents of Ψ(z − ib/2) and Ψ(z + ib/2) differ by −∫_{ℝ+i0} e^{−2izξ}/(2ξ·sinh(ξ/b)) dξ, which is log(1 + e^{2πbz}) by the classical Fourier integral ∫_{ℝ+i0} e^{−iuσ}/(2σ·sinh πσ) dσ = −log(1 + e^{u}), |Im u| < π. (Both sides of the latter vanish as Re u → −∞ and have the same derivative, by ∫ e^{−iuσ}/sinh(πσ) dσ = −i·tanh(u/2) as a principal value.) The integral is symmetric under b → 1/b, so the same holds with 1/b.
- The product formula satisfies the same two equations (v12 §3), and it has no zero or pole in the strip |Im z| < Re(b + 1/b)/2.
- So R = Ψ/(product) is analytic and non-vanishing in the strip and invariant under z → z + ib and z → z + i/b wherever both points are in the strip. The strip is higher than Re b, so R extends to an entire function with these two periods. They are independent over ℝ because b² is not real. An entire doubly periodic function is constant, and R → 1 as Re z → −∞ in the strip. So R = 1.

For Im b² < 0 apply this to 1/b. (ii) The integral is analytic in b, both sides of each identity are analytic in b, and they agree for Im b² > 0. ∎

**The variable w.** As before w = 2π·b·y = i·t·s for y = b·s. For Re x > 0 (substituting ξ = π·b·σ and turning the contour back to the real line, which crosses no pole):

  Λ(w) := log Φ_b(w/(2πb)) = ∫_{ℝ+i0} e^{−iwσ} / (4σ·sinh(πσ)·sinh(πxσ)) dσ,  |Im w| < π·(1 + Re x).

In this variable the two functional equations shift w by 2πi·x and by 2πi; the second one reads

  Φ_b(y + i/b) = Φ_b(y)/(1 + p^{−1}·w′),  w′ = e^{2πs} = e^{w/x},  p = e^{−πi/x} = e^{−2π²/t},

and the inversion relation reads Λ(w) + Λ(−w) = −t/24 + π²/(6t) + w²/(2t).

---

### [v16 §2] Two uniform lemmas in the strip

Throughout: |arg x| ≤ 60° and |x| ≤ ¼. Constants do not depend on x in this set; in particular they do not depend on whether x is real.

**Lemma S (expansion in the strip).** Let K ≥ 2 and 0 < ε₁ < π. For |Im w| ≤ π − ε₁ and any Re w,

  −Λ(w) = Σ_{j=0}^{K−1} B_{2j}(½)·t^{2j−1}·Li_{2−2j}(−e^{w})/(2j)! + O(|t|^{2K−1}),

with Li₂(−e^{w}) the principal value in the strip.

Proof. Write 1/sinh u = Σ_{k<K} a_k·u^{2k−1} + r_K(u), a_k = 2(1 − 2^{2k−1})·B_{2k}/(2k)!. On the rays arg u = arg x and arg x + 180° one has |r_K(u)| ≤ C_K·|u|^{2K−1}: for |u| ≤ 1 by the Laurent series, for |u| ≥ 1 because 1/sinh u is bounded there and each term is at most C·|u|^{2K−3}. Insert u = πxσ:

- the terms give a_k·(πx)^{2k−1}·I_k(w), I_k(w) = ∫_{ℝ+i0} e^{−iwσ}·σ^{2k−2}/(4·sinh πσ) dσ = (i/2)·(−1)^{k−1}·Li_{2−2k}(−e^{w}) (residues at σ = i·n for Re w < 0, then analytic continuation in the strip); with πx = i·t/2 this is the stated series;
- the remainder is ∫_ℝ e^{−iwσ}·r_K(πxσ)/(4σ·sinh πσ) dσ, regular at σ = 0 for K ≥ 2, and bounded by C_K·|x|^{2K−1}·∫ |σ|^{2K−2}·e^{(π−ε₁)|σ|}/(4·|sinh πσ|) dσ. On the real σ-line |e^{−iwσ}| = e^{σ·Im w} does not depend on Re w. ∎

**Lemma F (leading term up to the edge).** For Re w ≤ 0 and |Im w| ≤ π,

  | Λ(w) + Li₂(−e^{w})/t | ≤ C₁·e^{Re w/2}.

Proof. Λ(w) + Li₂(−e^{w})/t = ∫ e^{−iwσ}·r₁(πxσ)/(4σ·sinh πσ) dσ with r₁(u) = 1/sinh u − 1/u, first over ℝ + i0 (for |Im w| < π, and by continuity at |Im w| = π, where both sides still converge). Move the contour to Im σ = ½: no pole is crossed, since r₁ is regular at 0, the zeros of sinh πσ start at σ = i, and the poles of r₁(πxσ) are at σ = i·n/x with Im(i/x) ≥ 2. On the new line, with σ = a + i/2:

- |e^{−iwσ}| = e^{Re w/2 + a·Im w} ≤ e^{Re w/2}·e^{π|a|}, and |sinh πσ| = cosh πa, so e^{π|a|}/cosh(πa) ≤ 2;
- |r₁(u)| ≤ C·min(|u|, 1/|u|) for u = πxσ on the line (u stays at distance ≥ 1 from the poles i·π·n, and |Re u| grows like |u|/2);
- hence the integrand is at most e^{Re w/2}·C·min(|x|, 1/(|x|·|σ|²)), and the integral of the minimum over the line is bounded by a constant. ∎

Numerically the supremum of the left side times e^{−Re w/2} is 0.347 (attained at w = ±iπ), the same for real x and for the neighbouring directions, and independent of |t| (§6).

---

### [v16 §3] The contour

Take the contour of v15 with δ = 10°, for |β| ≤ 30° (−150° ≤ θ ≤ −30°): w = ±ρ·e^{iα}, α = β − 10° ∈ [−40°, 20°], γ = β + 10° ∈ [−20°, 40°].

**Lemma D (outside the strip).** Let w be on the left half of the contour and N the integer with |Im w − 2πN| ≤ π. Then

  Λ(w) = Λ(w − 2πi·N) + E,  |E| ≤ C·e^{−c/|t|},  E = 0 if N = 0,

with absolute constants c, C > 0. The same holds for −w when w is on the right half.

Proof. Iterating the functional equation in the direction i/b:

- N > 0: Φ_b(y) = Φ_b(y − iN/b)/Π_{j<N} (1 + p^{2j+1}·w′),
- N = −M < 0: Φ_b(y) = Φ_b(y + iM/b)·Π_{j<M} (1 + p^{−(2j+1)}·w′).

On the left half |w′| = e^{2π·Re s} = e^{−2πρ·cos γ/|t|}. Further |p| = e^{−2π²·sin 2β/|t|} and 2|N| − 1 ≤ |Im w|/π = ρ·|sin α|/π, so every factor satisfies

  |p^{±(2j+1)}·w′| ≤ exp( −(2πρ/|t|)·(cos γ − |sin α·sin 2β|) ).

For |β| ≤ 30° the bracket is at least 0.383, and N ≠ 0 forces ρ ≥ π/|sin α| ≥ 4.88. So each factor is 1 + O(e^{−11/|t|}) and there are at most ρ of them, each with its own factor e^{−2.4ρ/|t|}. On the right half replace w by −w and w′ by 1/w′. ∎

For real x this is easy to see directly: |p| = 1, and the factors are 1 + (a number of modulus e^{−2π|Re s|}).

*From `docs/addendum-v19.md`, §2; text unchanged.*

### [v19 §2] Bounds for Φ_b (from v16)

**Lemma B.** For 0 ≤ arg x ≤ 60°, |x| ≤ ¼, with constants independent of x:

- (i) In the closed strip |Im w| ≤ π, for every Re w: 1/|Φ_b| ≤ C·|exp(Li₂(−e^{w})/t)|.
- (ii) For |Im w| ≤ π − ε₁: 1/Φ_b = exp(Li₂(−e^{w})/t)·(1 + O(|t|)), uniformly.
- (iii) On the ray w = iπ − ρ·e^{−iδ₁}, ρ ≥ 0, with δ₁ + arg x ≤ 90°: 1/|Φ_b| ≤ C·2^{N}·|exp(Li₂(−e^{w})/t)|, N ≤ 1 + ρ/(2π).

Proof. (i) For Re w ≤ 0 this is Lemma F of v16. For Re w ≥ 0 the inversion relation gives −Λ(w) = Λ(−w) + t/24 − π²/(6t) − w²/(2t), Lemma F at −w gives Λ(−w) = −Li₂(−e^{−w})/t + O(1), and Li₂(−e^{w}) + Li₂(−e^{−w}) = −π²/6 − w²/2 in the strip. (ii) Lemma S of v16 with K = 2: the term after the leading one is t·B₂(½)·Li₀(−e^{w})/2, and Li₀(−e^{w}) = −e^{w}/(1 + e^{w}) is bounded at distance ≥ ε₁ from ±iπ. (iii) The functional equation in the direction i/b, N times (N the integer with |Im w − 2πN| ≤ π): 1/Φ_b(y) = Π_{j<N} (1 + p^{2j+1}·w′)/Φ_b(y − iN/b), w′ = e^{w/x}. Here |p| ≤ 1 and |p·w′| = |e^{(w − iπ)/x}| = exp(−(ρ/|x|)·cos(δ₁ + arg x)) ≤ 1, so each factor is at most 2; the shifted point is in the strip with the same Re w ≤ 0 and the same e^{w}; apply (i). ∎

**Lemma L (Laplace's method; standard).** Let Γ be a contour of finitely many smooth pieces, w\* an interior point of a straight piece with direction d, |Im w\*| < π, V′(w\*) = 0, and suppose for all τ in the range considered:

- (L1) H_τ(w) − H_τ(w\*) ≤ −c·min(|w − w\*|², |w − w\*|) on Γ (V continued along Γ);
- (L2) the integrand E satisfies |E| ≤ C·e^{C|w|}·|e^{V/t}| on Γ and E = e^{V/t}·(1 + O(|t|)) for |w − w\*| ≤ r₁.

Then ∫_Γ E dw = (2πt/(−V″(w\*)))^{1/2}·e^{V(w\*)/t}·(1 + O(|t|)), the square root being the one whose quotient by d has positive real part; the O is uniform.

Proof. Outside |w − w\*| ≤ |t|^{1/3} the integrand is O(e^{−c′|t|^{−1/3}}) times e^{cos ε·H_τ(w\*)/|t|}, by (L1), (L2). Inside, with w = w\* + s·d: V = V\* + V″s²d²/2 + V‴s³d³/6 + O(s⁴); (L1) gives Re(−V″d²/t) ≥ c/|t|; expand exp((V‴s³d³/6 + O(s⁴))/t) = 1 + V‴s³d³/(6t) + O(s⁴/|t| + s⁶/|t|²); the odd term integrates to 0 and the Gaussian moments of the rest are O(|t|). ∎

---

## Part 3. Descent of the leading exponent

The contour angle is δ = 10° from here on.

*From `docs/addendum-v15.md`, §1; text unchanged.*

### [v15 §1] The descent inequality

#### [v15 §1.1] Setting

t = |t|·e^{iθ} with −270° < θ < 90°, β = (θ + 90°)/2 = arg b. A contour angle δ is admissible if 0 < δ < 90° − |β|. Put α = β − δ and γ = β + δ; both lie in (−90°, 90°). The contour is s = r·e^{−iγ} (r real), and

  w := i·t·s = |t|·r·e^{iα},  ρ := |w|.

So w runs over the line e^{iα}·ℝ; the right half of the contour (r > 0) has Re w > 0. Let

  G(w) = −w²/2 + Li₂(−e^{w}) + π²/12 for Re w ≤ 0,  G(w) = −w² − Li₂(−e^{−w}) − π²/12 for Re w ≥ 0

(the two expressions agree for |Im w| < π; G is the function U of v14 continued to complex argument). The leading exponent of the integrand of P·J is ψ(ρ)/|t| with

  ψ(ρ) = Re( e^{−iθ}·G(w) ) − π·ρ·cos γ.

The first term comes from the Gaussian and the quantum dilogarithm, the second from the kernel, |1/(2cosh πs)| ≈ e^{−π|Re s|}.

#### [v15 §1.2] Proposition 1

Let ℓ = X + i·Y be

  ℓ = ∫₀¹ log(1 + e^{−λw}) dλ on the right half,  ℓ = ∫₀¹ log(1 + e^{λw}) dλ on the left half

(principal logarithm; the exponentials have modulus ≤ 1). Then:

- **(i) Exact formula.**
  - right half: ψ(ρ)/ρ = −(π − Y)·cos γ − X·sin γ − ρ·sin 2δ,
  - left half: ψ(ρ)/ρ = −(π + Y)·cos γ + X·sin γ − (ρ/2)·sin 2δ.
- **(ii) Bounds.** |Y| ≤ π/2 and λ_α ≤ X ≤ ln 2, where λ_α = ln(1 − e^{−(π/2)·cot|α|}) for α ≠ 0 and λ₀ = 0.
- **(iii) Descent.** ψ(ρ) ≤ −κ·ρ − c·ρ²·sin 2δ with κ = (π/2)·cos γ − max(ln 2, −λ_α)·|sin γ| and c = 1 on the right half, c = ½ on the left half.
- **(iv) One contour angle for a wide sector.** For δ = 10° and −60° ≤ β ≤ 50°, that is −210° ≤ θ ≤ 10°: κ ≥ 0.185.

**Proof.**

(i) G′(u) = −u − log(1 + e^{u}) for Re u ≤ 0 and G′(u) = −2u − log(1 + e^{−u}) for Re u ≥ 0, and G(0) = 0. Integrating along the segment from 0 to w: G(w) = −w²/2 − w·ℓ on the left half, G(w) = −w² − w·ℓ on the right half. Now e^{−iθ} = i·e^{−2iβ}, so e^{−iθ}·w = ±i·ρ·e^{−iγ} and e^{−iθ}·w² = i·ρ²·e^{−2iδ}, with Re(i·e^{−2iδ}) = sin 2δ and Re(i·e^{−iγ}·(X + iY)) = X·sin γ − Y·cos γ. Inserting gives the two formulas.

(ii) For |z| ≤ 1 the point 1 + z lies in the closed disk of radius 1 around 1, so |arg(1 + z)| ≤ π/2 and log|1 + z| ≤ ln 2; averages inherit both. For the lower bound: z = e^{∓λw} has modulus e^{−λρ·cos α} and argument of size λρ·|sin α|. If the real part of z is ≥ 0, then |1 + z| ≥ 1. Otherwise λρ·|sin α| > π/2, hence |z| < e^{−(π/2)·cot|α|} and |1 + z| ≥ 1 − |z| > e^{λ_α}. So log|1 + z| ≥ λ_α at every point, and X ≥ λ_α.

(iii) cos γ > 0 and π ∓ Y ≥ π/2 bound the first term by −(π/2)·cos γ; |X| ≤ max(ln 2, −λ_α) bounds the second.

(iv) Here γ ∈ [−50°, 60°] and α ∈ [−70°, 40°]. For γ ≥ 0 we have |α| ≤ 40° < 66.2°, where −λ_α ≤ ln 2, so κ ≥ (π/2)·cos 60° − ln 2·sin 60° = 0.185. For γ < 0: −λ_α ≤ −λ_{70°} = 0.831 and κ ≥ (π/2)·cos 50° − 0.831·sin 50° = 0.373. ∎

Remarks.

- On the real contour of v14 (α = 0) the formula of (i) is U(τ)/τ rewritten, and (iii) gives descent for γ < 66.2°. v14 got γ < 77.56° there by concavity, and the scan of v14 §4 says 77.56° is the truth for all admissible contours. The bounds of (ii) are not sharp together, and they do not need to be: one good contour per direction is enough.
- Concavity of ψ in ρ holds on 2101 of 2760 admissible grid contours, on all of them only for −85° ≤ θ ≤ −60°, and on none for θ ≤ −200° or θ ≥ 55° (`w0_concavity.py`). So it cannot be the general argument.

---

## Part 4. The uniform bound for the kernel integrand

Needed later: Lemma 3′ (away from the real x-axis) and the first bullet of the proof of Theorem 3 (near it), which together give |H| ≤ C·e^{−κ|s|} on the closed sector, and the Corollary (analyticity). The expansions themselves (Theorems 1–3) are superseded by Theorem 4 of Part 5.

*From `docs/addendum-v14.md`, §2.1–§2.3 (contour, Euler–Maclaurin lemma, second form of the integrand); text unchanged.*

#### [v14 §2.1] The contour and the integrand

In the sector, 0 < β < 45°, so δ = β is an admissible angle for the line L = e^{−iδ}·ℝ of v12 (the condition there is 0 < δ < 90° − β). J does not depend on the admissible δ: turning L keeps both ends in directions of Gaussian decay and crosses no pole.

On this L write y = b·s and s = −i·τ/t. Then τ = i·t·s runs over the real line, e^{2πby} = e^{τ}, and (v13 §1)

  P·J(1/(2b)) = e^{t/8} · ∫ H ds,  H = e^{−τ²/(2t)} · R(τ) · Π(τ) / (2·cos(πτ/t)),

  R(τ) = (−q;q)∞/(−q^{1/2}e^{τ}; q)∞,  Π(τ) = (−p·w; p²)∞/(p;p²)∞,  w = e^{2πs} = e^{−2πiτ/t}.

Two facts about the geometry, both from θ < 0: Re s = τ·|sin θ|/|t| and Im s = −τ·cos θ/|t|. For the ψ₀\* integral the integrand has the extra factor e^{τ} and the prefactor is e^{t/8 − t/2}.

#### [v14 §2.2] Lemma 1 (Euler–Maclaurin with remainder)

Let u ∈ ℂ with 1 + u·e^{−ty} in the right half-plane or in the disk |z − 1| ≤ 1 − d for all y ≥ 0 (some d > 0), and let F(y) = log(1 + u·e^{−ty}) with the principal logarithm. Then for a ∈ (0, 1] and K ≥ 1

  Σ_{n≥0} F(n + a) = −Li₂(−u)/t + Σ_{k=1}^{K} B_k(a)·(−t)^{k−1}·Li_{2−k}(−u)/k! + ρ_K,  |ρ_K| ≤ C_K·|t|^{K−1}·|u|·(1 + |u|)^{K−1}/(cos θ·d^K).

Proof. Euler–Maclaurin with offset a: the sum is ∫₀^∞ F − Σ_{k≤K} B_k(a)·F^{(k−1)}(0)/k! plus a remainder bounded by (sup|B̃_K|/K!)·∫₀^∞ |F^{(K)}|. Here ∫₀^∞ F = −Li₂(−u)/t (substitute u′ = u·e^{−ty}; the path stays in the region where log(1 + u′) is the principal branch), and F^{(k)}(y) = −(−t)ᵏ·Li_{1−k}(−u·e^{−ty}). For k ≥ 1, Li_{1−k}(−v) is v times a polynomial of degree k − 1 in v, divided by (1 + v)ᵏ. So |F^{(K)}(y)| ≤ C_K·|t|ᴷ·|u|e^{−y·Re t}·(1 + |u|)^{K−1}/d^K, and ∫₀^∞ e^{−y·Re t} dy = 1/(|t|·cos θ). ∎

**The hypothesis holds for u = e^{τ} with τ ≤ τ₁ := (π/4)·cot|θ|, with d = d₀ := 1 − e^{−(π/4)·cot|θ|}.** The point u·e^{−ty} has modulus e^{τ − φ·cot|θ|} and argument φ = y·|Im t| ≥ 0. For φ ≤ 90° its real part is ≥ 0. For φ > 90° its modulus is below e^{τ₁ − (π/2)cot|θ|} = 1 − d₀. In the sector, cot|θ| ≥ tan ε, so d₀ is bounded below.

**Consequences.** With a = ½, u = e^{τ} for the denominator of R and a = 1, u = 1 for the numerator, for τ ≤ τ₁ and every K:

  (2.1) log R(τ) − τ²/(2t) = U(τ)/t + Σ_{k=1}^{K} (−t)^{k−1}·cₖ(τ)/k! + O(|t|ᴷ),

  **U(τ) = −τ²/2 + Li₂(−e^{τ}) + π²/12**,  cₖ(τ) = B_k(1)·Li_{2−k}(−1) − B_k(½)·Li_{2−k}(−e^{τ}),

uniformly in τ ≤ τ₁ and in the sector. (Apply the lemma with K + 1; c₁ = −½·ln 2.)

#### [v14 §2.3] Lemma 2 (the other form of the integrand)

For all τ,

  e^{−τ²/(2t)}·R(τ)·Π(τ) = 2^{−1/2}·e^{−π²/(12t) + t/12} · e^{−τ²/t} · (−q^{1/2}e^{−τ}; q)∞/(−p/w; p²)∞.

Proof. R·Π = (−q;q)∞/((p;p²)∞·Φ_b(bs)). Insert 1/Φ_b(y) = e^{−πiy²}·Φ_b(−y)/Φ_b(0)² (v12, 2.2), e^{−πi·x·s²} = e^{−τ²/(2t)}, (−q;q)∞/(p;p²)∞ = 2^{−1/2}e^{π²/(12t) + t/24} (v6) and Φ_b(0)^{−2} = e^{t/24 − π²/(6t)}. ∎

By Lemma 1 with u = e^{−τ} (τ ≥ 0, so again d ≥ d₀), the logarithm of the right side is [−τ² − Li₂(−e^{−τ}) − π²/12]/t + O(1), and this bracket is U(τ) again, by the reflection formula Li₂(−e^{τ}) + Li₂(−e^{−τ}) = −π²/6 − τ²/2.

*From `docs/addendum-v15.md`, §2 up to Lemma 3′; text unchanged.*

### [v15 §2] Theorem 2

#### [v15 §2.1] Statement

For Im x < 0 (|q| > 1) let Φ_b := Φ_{1/b}, with the right side given by the product formula of v12 (legitimate, because Im(1/b²) > 0). Explicitly, with t₁ = t for Re t > 0 and t₁ = −t for Re t < 0, q₁ = e^{−t₁}, p₁ = e^{−2π²/t₁}, w′ = e^{2πs}:

  1/Φ_b(b·s) = [ (−p₁·w′; p₁²)∞ / (−q₁^{1/2}·e^{w}; q₁)∞ ]^{±1},  + for Im x > 0, − for Im x < 0.

This is the continuation used in v13. All identities of v12 for Φ_b that were stated for Im x > 0 hold for Im x < 0 by applying them to 1/b; in particular Φ_b(y)·Φ_b(−y) = Φ_b(0)²·e^{πiy²} with Φ_b(0)² = e^{−t/24 + π²/(6t)}. Let C₀ = 2^{−1/2}·e^{π²/(12t) + t/24} and

  P·J := e^{t/8} · ∫ H ds,  H = C₀·e^{−w²/(2t)} / (Φ_b(b·s)·2·cosh πs),

over the contour of §1. For Im x > 0 this is P·J(1/(2b)) of v12 (v13 §1, with C₀ = (−q;q)∞/(p;p²)∞ by v6). For Im x < 0 it is the definition. The ψ₀\* version has the extra factor e^{w} and the prefactor e^{t/8 − t/2}.

**Theorem 2.** Fix ε > 0 and K ≥ 0. There are C and t₀ > 0 such that for all t with −210° ≤ θ ≤ 10°, |θ + 90°| ≥ ε and |t| ≤ t₀:

  | P·J − Σ_{k<K} aₖ·tᵏ | ≤ C·|t|ᴷ,

and the same for the ψ₀\* version with its coefficients.

In terms of x: ε ≤ arg x ≤ 100° and −120° ≤ arg x ≤ −ε.

#### [v15 §2.2] Proof

The contour is that of §1 with δ = 10° (admissible because |β| ≤ 60°; for Im x > 0 the integral does not depend on the admissible δ, v14 §2.1). The proof follows v14 §2 step by step; I state each lemma in the form needed now and give only what changes.

**Lemma 1′ (Euler–Maclaurin).** Let Re t₁ > 0, ζ₀ ∈ ℂ and ζ(y) = ζ₀ − t₁·y. Suppose that for all y ≥ 0 the point ζ(y) has distance ≥ η from every i·π·(2m + 1) and Re ζ(y) ≤ η/2. Then, with F(y) = log(1 + e^{ζ(y)}) (principal), a ∈ (0, 1], K ≥ 1:

  Σ_{n≥0} F(n + a) = −Li₂(−e^{ζ₀})/t₁ + Σ_{k=1}^{K} B_k(a)·(−t₁)^{k−1}·Li_{2−k}(−e^{ζ₀})/k! + ρ_K,  |ρ_K| ≤ C(K, η)·|t₁|^{K−1}/cos(arg t₁).

Proof as in v14, Lemma 1. The hypotheses give: 1 + e^{ζ} never lies on (−∞, 0] (that would need Im ζ ≡ π mod 2π with 0 ≤ Re ζ ≤ η/2, a point within η/2 of the lattice), |1 + e^{ζ}| ≥ c(η) > 0 (periodicity and compactness), and |e^{ζ}| ≤ e^{η/2}. These are all that the proof uses. ∎

**Lemma A (geometry of the paths).** Let η = π·min(cos α, |cos θ|). For w on the left half of the contour, the path w − t₁·y satisfies the hypotheses of Lemma 1′ with this η; for w on the right half, so does −w − t₁·y; and for any w with |w| ≤ η/4, the path w − t₁·y satisfies them with 3η/4.

Proof. Both −e^{iα} and −t₁ point into the open left half-plane, so the paths from points ∓ρe^{iα} lie in the closed convex cone spanned by these two directions, which meets the imaginary axis only at 0. A point i·π·(2m + 1) has distance at least π·sin(angle to the nearer edge) from the cone, and the sines of the angles between the edges and the imaginary axis are cos α and |cos θ|. A path starting within η/4 of the cone stays within η/4 of it. ∎

In the sector of the theorem, cos α ≥ cos 70° and |cos θ| ≥ sin ε, so η is bounded below.

**Lemma 2′ (two forms and their expansions).**

- Left form, used on the left half and for |w| ≤ η/4:

  (2.1′) log( C₀·e^{−w²/(2t)}/Φ_b(b·s) ) = G(w)/t − ½·ln 2 + t/24 − Σ_{k=2}^{K} B_k(½)·(−t)^{k−1}·Li_{2−k}(−e^{w})/k! + log D + O(|t|ᴷ),

  with D = (−p₁·w′; p₁²)∞^{±1}.
- Right form, used on the right half: C₀·e^{−w²/(2t)}/Φ_b(b·s) = C₀·Φ_b(0)^{−2}·e^{−w²/t}·Φ_b(−b·s) by the inversion relation, and its logarithm is G(w)/t + O(1) + log D′, D′ = (−p₁/w′; p₁²)∞^{∓1}.

Proof. Apply Lemma 1′ with a = ½ to the q₁-product, using Lemma A. For Im x < 0 the q₁-product stands in the numerator and t₁ = −t. This changes the sign of the logarithm and replaces (−t)^{k−1} by t^{k−1}; the two changes cancel for even k (and in the leading term), and B_k(½) = 0 for odd k. So (2.1′) has the same form on both sides of the real axis. (For Im x > 0 it is (2.1) of v14 with the numerator (−q;q)∞ replaced by its exact value C₀·(p;p²)∞.) ∎

**The dual nome.** |p₁| = e^{−2π²·|cos θ|/|t|}. On the left half |w′| ≤ 1, on the right half |1/w′| ≤ 1, and for |s| ≤ |t|^{−1/2} one has |w′| ≤ e^{2π|t|^{−1/2}}. So D = 1 + O(e^{−c/|t|}) wherever the left form is used, and D′ = 1 + O(|p₁|) on the right half.

**Lemma 3′ (uniform bound).** |H| ≤ C·e^{ψ(ρ)/|t|} ≤ C·e^{−κ·|s|} on the whole contour, with κ ≥ 0.185.

Proof. Lemma 2′ with K = 1 gives |C₀·e^{−w²/(2t)}/Φ_b| ≤ C·exp(Re(G(w)/t)), and Re(G/t) = Re(e^{−iθ}G)/|t|. On rays with |γ| ≤ 60°, 1/|2cosh πs| ≤ C·e^{−π|Re s|} and π·|Re s| = π·ρ·cos γ/|t|. Proposition 1 (iii), (iv) and ρ/|t| = |s| finish the proof. ∎

*From `docs/addendum-v16.md`, §4; text unchanged.*

### [v16 §4] Theorem 3

**Theorem 3.** Let P·J = e^{t/8}·∫ H ds with H = C₀·e^{−w²/(2t)}/(Φ_b(b·s)·2·cosh πs), C₀ = 2^{−1/2}·e^{π²/(12t) + t/24}, over the contour of §3, and Φ_b as in §1. For every K there are C and t₀ > 0 such that for −150° ≤ θ ≤ −30° and |t| ≤ t₀

  | P·J − Σ_{k<K} aₖ·tᵏ | ≤ C·|t|ᴷ,

and the same for the ψ₀\* version (extra factor e^{w}, prefactor e^{t/8 − t/2}) with its coefficients.

By Lemma Φ (i) this P·J is the one of Theorem 2 of v15 for Im x ≠ 0, and for Im x > 0 the one of v12.

Proof. The structure is that of v14 §2 and v15 §2.2. Only the two inputs about Φ_b change.

- **Uniform bound (replaces Lemma 3′).** Left half: by Lemma D and Lemma F (at w − 2πi·N, which has the same real part and the same e^{w}), −Λ(w) = Li₂(−e^{w})/t + O(1). So log(C₀·e^{−w²/(2t)}/Φ_b) = G(w)/t + O(1) with G(w) = −w²/2 + Li₂(−e^{w}) + π²/12. Right half: by the inversion relation C₀·e^{−w²/(2t)}/Φ_b(b·s) = 2^{−1/2}·e^{−π²/(12t) + t/12}·e^{−w²/t}·Φ_b(−b·s), and by Lemmas D and F at −w, Λ(−w) = −Li₂(−e^{−w})/t + O(1); the result is again G(w)/t + O(1), now with G(w) = −w² − Li₂(−e^{−w}) − π²/12. With the window bound for 1/cosh and Proposition 1 of v15 (κ ≥ 0.185 for δ = 10° and |β| ≤ 30°): |H| ≤ C·e^{ψ(ρ)/|t|} ≤ C·e^{−κ|s|} on the whole contour.
- **Expansion in the window (replaces (2.1′)).** For |w| ≤ |t|^{1/2}, Lemma S with ε₁ = π/2 gives log(C₀·e^{−w²/(2t)}/Φ_b) = G(w)/t − ½·ln 2 + t/24 + Σ_{j=1}^{K−1} B_{2j}(½)·t^{2j−1}·Li_{2−2j}(−e^{w})/(2j)! + O(|t|^{2K−1}). This is (2.1′) of v15 term by term, without a dual factor.
- The rest (Lemma 4′, splitting at |s| = |t|^{−1/2}, the moments, the identification of the coefficients) is unchanged. All constants are uniform in |β| ≤ 30°. ∎

**Corollary.** P·J is analytic in t on −210° < θ < 10°, and for every K the bound |P·J − Σ_{k<K} aₖtᵏ| ≤ C_K·|t|ᴷ holds uniformly on the closed sector −210° ≤ θ ≤ 10°, |t| ≤ t₀.

Proof. Theorem 2 of v15 with ε = 60° covers |θ + 90°| ≥ 60°, Theorem 3 the rest. Analyticity: the integrand is analytic in b by Lemma Φ, the contour angle is fixed, and the bound |H| ≤ C·e^{−κ|s|} is locally uniform. ∎

In terms of x the sector is −120° ≤ arg x ≤ 100°. For real x > 0 this is a statement about |q| = 1, where neither ψ₁\* nor f₁\* exists.

*From `docs/addendum-v26.md`, §1 (the proof of analyticity, restated); text unchanged.*

### [v26 §1] Analyticity of P·J (the Corollary of v16 §4)

**Corollary (the statement of v16, with the range of |t| written out).** P·J is analytic in t on −210° < θ < 10°, 0 < |t| < t₀, where t₀ is the constant of the uniform bound.

Proof. With y = b·s the contour s = r·e^{−iγ}, γ = β + 10°, is y = |b|·r·e^{−10°·i}: the line L = e^{−10°·i}·ℝ of v12, for every t. Since w = 2πby and t = −2πib², e^{−w²/(2t)} = e^{−πiy²}, and

  P·J = e^{t/8}·(C₀/b)·∫_L e^{−πiy²}/(Φ_b(y)·2cosh(πy/b)) dy,  C₀ = 2^{−1/2}·e^{π²/(12t) + t/24},  b = (it/2π)^{1/2}.

- *The integrand is analytic in b for fixed y ∈ L.* Φ_b(y) is analytic in (b, y) for Re b > 0 (Lemma Φ). Its zeros y = −i·[(m + ½)·b + (n + ½)/b] have arguments in [−90° − |β|, −90° + |β|], the poles of the kernel y = i·b·(n + ½) have arguments ±90° + β. For |β| ≤ 60° neither set meets L, whose points have arguments −10° and 170°.
- *Domination.* The bound |H| ≤ C·e^{−κ|s|} (Lemma 3′ of v15 and v16 §4), uniform on the closed sector for |t| ≤ t₀, reads |integrand| ≤ C·e^{−κ|y|/|b|} in the variable y, and |b| is bounded above and below on compact subsets of the punctured sector.
- An integral over a fixed contour of functions analytic in a parameter, dominated locally uniformly, is analytic (Fubini and Morera). ∎

No contour other than L is used, so the questions raised in F3 (descent for δ ≠ 10°, independence of the contour for Im x ≤ 0) do not arise. The stress test's S3.7 (κ > 0 for δ between 5° and 15°) shows that the other route would work as well.

---

## Part 5. Factorial bounds and Borel summability (Theorem II)

The coefficients called ãₖ in the proof of Theorem 4 are the aₖ of §0.1, by definition; the last step of that proof (ãₖ = aₖ) only identifies them with the coefficients computed in v13.

*From `docs/addendum-v17.md`, §1–§5; text unchanged.*

### [v17 §1] The shift identity

Notation as in v16 §4: t in the sector 𝕊 = {−210° ≤ θ ≤ 10°, 0 < |t| ≤ t₀}, θ = arg t, β = (θ + 90°)/2, contour Γ: s = r·e^{−iγ}, r ∈ ℝ, with γ = β + 10°, so w = i·t·s = ±ρ·e^{iα}, α = β − 10°, ρ = |t|·|s|. On 𝕊: α ∈ [−70°, 40°], γ ∈ [−50°, 60°].

  H(s) = C₀·e^{ts²/2}/(Φ_b(b·s)·2·cosh πs),  C₀ = 2^{−1/2}·e^{π²/(12t) + t/24},  P·J = e^{t/8}·∫_Γ H ds.

Let

  τₙ(t) = q^{n²+n}/(2·(−q;q)ₙ),  D_N(w) = Π_{j=1}^{N} g(ζⱼ),  g(ζ) = ζ²/(1 + ζ),  ζⱼ = e^{−ωⱼ},  ωⱼ = w + t·(j − ½),

with q = e^{−t}. Written out, D_N(w) = q^{N²}·e^{−2Nw}/(−q^{1/2}e^{−w}; q)_N.

**Proposition 2 (shift identity).** For t ∈ 𝕊 and every integer N ≥ 0 with (N − ½)·|t|·cos γ < π·cos α (in particular whenever N·|t| ≤ 1):

  P·J = Σ_{n<N} (−1)ⁿ·τₙ(t) + R_N(t),  R_N(t) = (−1)^N·e^{t/8}·∫_Γ H(s)·D_N(w) ds.

For the ψ₀\* version (extra factor e^{w} in H, prefactor e^{t/8 − t/2}): τₙ = q^{n²}/(2·(−q;q)ₙ) and D_N gets the extra factor q^{−N}.

Proof.

- **Transformation of H.** For s → s − iN: cosh(π(s − iN)) = (−1)^N·cosh πs; e^{t(s−iN)²/2} = e^{ts²/2}·e^{−Nw}·q^{N²/2}; and by the functional equation in the direction b (v16 Lemma Φ (ii)), iterated, Φ_b(y − iNb) = Φ_b(y)·Π_{j=1}^{N} (1 + q^{1/2−j}·e^{w}) = Φ_b(y)·q^{−N²/2}·e^{Nw}·(−q^{1/2}e^{−w}; q)_N. Together: H(s − iN) = (−1)^N·H(s)·D_N(w). For ψ₀\*: e^{w} → e^{w}·q^{−N}.
- **Poles between Γ and Γ − iN.** The zeros of cosh πs at s = −i(n + ½), n = 0 … N − 1. The zeros of Φ_b(b·s) are at s = −i(k + ½) − (2l + 1)·π/t, k, l ≥ 0 (the lattice −c_b − i·b·k − i·l/b; from the products for Im x ≠ 0, and for real x from the functional equations, since Φ_b has no zero in the strip of v16 §1). Such a point lies on Γ − iλ only if (λ − k − ½)·|t|·cos γ = (2l + 1)·π·cos α, which is excluded for 0 ≤ λ ≤ N by the hypothesis. On 𝕊, π·cos α ≥ π·cos 70° = 1.07.
- **Residues.** At s = −i/2: e^{ts²/2} = e^{−t/8}, the derivative of 2·cosh πs is −2πi, and Φ_b(−ib/2) = 2·C₀. (Functional equation at z = 0: Φ_b(−ib/2) = 2·Φ_b(ib/2); inversion: Φ_b(ib/2)·Φ_b(−ib/2) = Φ_b(0)²·e^{t/8} = e^{π²/(6t) + t/12}; so Φ_b(−ib/2)² = 4·C₀², and the sign is fixed by the product formula and continuity in b.) Hence e^{t/8}·(−2πi)·Res_{s=−i/2} H = ½ = τ₀. By the transformation, the residue at s = −i(n + ½) is (−1)ⁿ times this times Dₙ(t/2), and at w = t/2 one has ζⱼ = qʲ, so Dₙ(t/2) = q^{n²+n}/(−q;q)ₙ. For ψ₀\* the extra factor at the pole is e^{w}·e^{−t/2} = q^{−n}.
- **The ends (Lemma E).** For fixed t ∈ 𝕊 and λ₀ > 0, sup_{0≤λ≤λ₀} |H(s − iλ)| → 0 as s → ∞ on Γ. Proof: put s̃ = s − iλ, w̃ = i·t·s̃ = w + tλ. On the left end Re w̃ → −∞ and |e^{2πs̃}| = e^{−2π|s|cos γ} → 0, and Φ_b(b·s̃) → 1 uniformly in λ: for Im x > 0 and Im x < 0 from the two product formulas (both arguments e^{w̃} and e^{2πs̃} tend to 0); for real x from Lemma F of v16 after bringing w̃ back into the strip with the functional equation in the direction i/b, which costs at most O(|s|) factors 1 + O(e^{−2π|s|cos γ}) (as in Lemma D of v16, where |p| = 1). So |H(s̃)| ≤ C·|e^{ts̃²/2}| = C·exp(−½·|t|·|s|²·sin 20° + O(|s|)), because arg(t·s²) = θ − 2γ = −90° − 20°. On the right end use the inversion relation: 1/Φ_b(y) = Φ_b(−y)·e^{−πiy²}/Φ_b(0)², e^{−πiy²} = e^{ts̃²/2}, so |H(s̃)| ≤ C·|e^{ts̃²}|·|Φ_b(−b·s̃)| and Φ_b(−b·s̃) → 1 as before. ∎
- **Cauchy's theorem** on the parallelogram with corners ±r·e^{−iγ}, ±r·e^{−iγ} − iN, r → ∞: the short sides vanish by Lemma E; the contour is clockwise, so ∫_Γ H ds = ∫_Γ H(s − iN) ds + Σ_{n<N} (−2πi)·Res_{s=−i(n+½)} H. ∎

Numerically (§6): the two sides agree to 10⁻¹³ for N = 1, 2, 5 at real x and on both sides.

---

### [v17 §2] The size of D_N

**Lemma G.** Let c₃ = 2·e^{−1/2}·sin(½) = 0.5816 and L = 2 + 1/c₃ = 3.7195. If every point of the segment [0, ω] has distance ≥ 1 from all points i·π·(2m + 1), then

  |g(e^{−ω})| ≤ ½·e^{L·|ω|}.

Proof. h(ω) = log(2·g(e^{−ω})) = −2ω − log((1 + e^{−ω})/2) is analytic along the segment, h(0) = 0, and h′(ω) = −2 + 1/(1 + e^{ω}). For ω′ on the segment choose m with |Im ω′ − π(2m + 1)| ≤ π and put η = ω′ − i·π·(2m + 1) = a + i·b. Then 1 + e^{ω′} = 1 − e^{η}, |η| ≥ 1, |b| ≤ π, and |e^{η} − 1|² = 4·e^{a}·(sinh²(a/2) + sin²(b/2)).

- a ≤ −1: |1 − e^{η}| ≥ 1 − e^{−1} = 0.632. a ≥ 1: |e^{η} − 1| ≥ e − 1.
- |a| < 1 and 1 ≤ |b| ≤ π: sin²(b/2) ≥ sin²(½), so |e^{η} − 1|² ≥ 4·e^{−1}·sin²(½).
- |a| < 1 and |b| < 1: sin(|b|/2) ≥ |b|·sin(½) (concavity on [0, ½]) and sinh²(a/2) ≥ a²/4 ≥ a²·sin²(½); the bracket is at least sin²(½)·(a² + b²) ≥ sin²(½).

So |1 + e^{ω′}| ≥ c₃ and |h′| ≤ L on the segment; hence |h(ω)| ≤ L·|ω|. ∎

(The true infimum of |1 − e^{η}| on that set is 1 − 1/e, attained at η = −1; numerically the supremum of |h′| near our contours is 2.53, at α = −70°.)

**Corollary G.** Let c₁ = 0.04, t ∈ 𝕊 and N·|t| ≤ c₁. Then for w on Γ

  |D_N(w)| ≤ exp( −ν·N + L·c₁·|s| ),  ν = ln 2 − L·c₁ = 0.5444,  L·c₁ = 0.1488.

Proof. |ωⱼ − w| = |t|·(j − ½) ≤ c₁, so the segment [0, ωⱼ] stays within c₁ of the line e^{iα}ℝ. The points i·π·(2m + 1) have distance ≥ π·cos α ≥ 1.0745 from that line. So Lemma G applies to each ωⱼ: |D_N| ≤ 2^{−N}·exp(L·Σⱼ|ωⱼ|) ≤ 2^{−N}·exp(L·N·(|w| + c₁)), and N·|w| = N·|t|·|s| ≤ c₁·|s|. ∎

**Corollary R (N residues approximate P·J exponentially well; proved).** Take t₀ ≤ 0.01. For t ∈ 𝕊 and N ≤ c₁/|t|:

  | P·J − Σ_{n<N} (−1)ⁿ·τₙ(t) | ≤ C₂·e^{−ν·N}.

Proof. By the uniform bound (v15 Lemma 3′ with ε = 60° for |θ + 90°| ≥ 60°; v16 §4 for |θ + 90°| ≤ 60°) |H| ≤ C_U·e^{−κ|s|} on Γ with κ = 0.185; for ψ₀\*, |e^{w}| ≤ e^{|t||s|}, so the exponent is κ₀ = κ − t₀ ≥ 0.175 in both cases. With Proposition 2 and Corollary G:

  |R_N| ≤ e^{t₀/2 + c₁}·C_U·e^{−νN}·∫_Γ e^{−(κ₀ − L·c₁)|s|} |ds| = e^{t₀/2 + c₁}·C_U·e^{−νN}·2/(κ₀ − L·c₁),

where e^{t₀/2} covers the prefactor e^{t/8} (or e^{t/8 − t/2}) and e^{c₁} the extra factor q^{−N} of the ψ₀\* version. ∎

Numerically the decay is like 2^{−N}, and it continues far beyond N = 0.04/|t| (§6).

---

### [v17 §3] The residues as functions of t

**Lemma T.** Let c = 0.3, m_c = 1 + e^{−c}·cos c = 1.7077, λ = e^{c}/m_c = 0.7904. Then τₙ is analytic for |t| < π/n and

  |τₙ(t)| ≤ (e^{c}/2)·λⁿ for |t| ≤ c/n,  hence |τ_{n,k}| ≤ (e^{c}/2)·λⁿ·(n/c)ᵏ for the Taylor coefficients τ_{n,k} = [tᵏ] τₙ.

The same holds for the ψ₀\* version.

Proof. For 1 ≤ j ≤ n and |t| ≤ c/n the number u = j·t has |u| ≤ c < π/2, so |1 + e^{−u}| ≥ Re(1 + e^{−u}) = 1 + e^{−Re u}·cos(Im u) ≥ m_c. And |q^{n²+n}| ≤ e^{|t|(n²+n)} ≤ e^{c(n+1)} (for ψ₀\*: e^{cn}). Cauchy's estimate on the circle |t| = c/n gives the coefficients. ∎

Consequently ãₖ := Σ_{n≥0} (−1)ⁿ·τ_{n,k} converges absolutely, and with Σ_{n≥0} nᵏ·λⁿ ≤ k!/(1 − λ)^{k+1} (because nᵏ ≤ (n+1)…(n+k)):

  |ãₖ| ≤ A₀·k!·B₀ᵏ,  A₀ = e^{c}/(2(1 − λ)) = 3.22,  B₀ = 1/(c·(1 − λ)) = 15.91.

---

### [v17 §4] Theorem 4

**Theorem 4.** Let t₀ ≤ 0.01 be such that the uniform bound of v15/v16 holds on 𝕊. There is A such that for all t ∈ 𝕊 and all K ≥ 0

  | P·J − Σ_{k<K} aₖ·tᵏ | ≤ A·K!·(B·|t|)ᴷ,  B = 240,

for the ψ₁\* and the ψ₀\* version, with the coefficients aₖ of Φ (v13). Moreover aₖ = Σ_{n≥0} (−1)ⁿ·[tᵏ] τₙ.

Proof. Fix t ∈ 𝕊 and put N = ⌊c₁/|t|⌋, so that c₁/2 ≤ c₁ − |t| < N·|t| ≤ c₁. By Proposition 2,

  P·J − Σ_{k<K} ãₖtᵏ = R_N + Σ_{n<N} (−1)ⁿ·[τₙ(t) − Σ_{k<K} τ_{n,k}tᵏ] − Σ_{n≥N} (−1)ⁿ·Σ_{k<K} τ_{n,k}tᵏ.

- **(I)** By Corollary R, |R_N| ≤ C₂·e^{ν}·e^{−ν·c₁/|t|} ≤ C₂·e^{ν}·K!·(|t|/(ν·c₁))ᴷ, and 1/(ν·c₁) = 45.9. (e^{−a/|t|} ≤ K!·(|t|/a)ᴷ for every K.)
- **(II)** For n < N: n·|t|/c ≤ c₁/c < ½, so the Taylor series converges and by Lemma T the remainder is at most (e^{c}/2)·λⁿ·Σ_{k≥K} (n|t|/c)ᵏ ≤ e^{c}·λⁿ·(n|t|/c)ᴷ. Summed over n: at most (e^{c}/(1 − λ))·K!·(B₀·|t|)ᴷ.
- **(III)** For n ≥ N put u = n·|t|/c ≥ u₀ := c₁/(2c) = 1/15. Then Σ_{k<K} uᵏ ≤ uᴷ·Σ_{m=1}^{K} u₀^{−m} ≤ 1.08·(u/u₀)ᴷ = 1.08·(2n|t|/c₁)ᴷ. With Lemma T and the sum over n: at most (0.73/(1 − λ))·K!·(B₁·|t|)ᴷ, B₁ = 2/(c₁·(1 − λ)) = 238.5.

So |P·J − Σ_{k<K} ãₖtᵏ| ≤ A·K!·(240·|t|)ᴷ. In particular Σ ãₖtᵏ is an asymptotic expansion of P·J on 𝕊; by the Corollary of v16 so is Σ aₖtᵏ; asymptotic expansions are unique, so ãₖ = aₖ. ∎

Remarks.

- ãₖ = aₖ is also what v13 §2 observed directly (Φ is the termwise expansion of the residue series), and it is checked exactly for k < 13 in §6.
- The constant 240 is far from the truth. Numerically the ratios approach 1/|S| = 0.2905, the distance of the nearest Borel singularity (§6). Watson's theorem needs only some finite B.
- A depends on C_U and t₀, which v15/v16 did not make explicit.

---

### [v17 §5] Theorem 5: Borel summability

#### [v17 §5.1] Statement and proof

**Theorem (Nevanlinna; Sokal 1980; quoted).** Let f be analytic in the disc C_R = {z : Re(1/z) > 1/R} and suppose f(z) = Σ_{k<K} cₖzᵏ + R_K(z) with |R_K(z)| ≤ A·σᴷ·K!·|z|ᴷ for all K ≥ 0 and all z ∈ C_R. Then B(ζ) = Σ cₖζᵏ/k! converges for |ζ| < 1/σ, continues analytically to the strip {ζ : dist(ζ, ℝ₊) < 1/σ}, is bounded there by a constant times e^{|ζ|/R} (on every narrower strip), and f(z) = (1/z)·∫₀^∞ e^{−ζ/z}·B(ζ) dζ for z ∈ C_R.

Reference: A. D. Sokal, *An improvement of Watson's theorem on Borel summability*, J. Math. Phys. 21 (1980) 261–263, https://doi.org/10.1063/1.524408. I took the statement from quotations of it in later papers (for example arXiv:2511.04509 §4.2, arXiv:1011.6490 §4), not from the original. Watson's original version (a sector of opening greater than 180°) would also do here.

**Theorem 5.** Let B(ξ) = Σ aₖξᵏ/k! be the Borel transform of Φ (ψ₁\* or ψ₀\* version), and R = t₀. Then:

- (i) B is analytic in the set {ξ : dist(ξ, Σ₀) < 1/240}, Σ₀ = {−120° ≤ arg ξ ≤ −80°} ∪ {0}, and |B(ξ)| ≤ C·e^{|ξ|/R} on {dist(ξ, Σ₀) ≤ 1/241}.
- (ii) For every θ₀ ∈ [−120°, −80°] and every t with Re(e^{iθ₀}/t) > 1/R:

  P·J(t) = 𝒮_{θ₀}[Φ](t) := (1/t)·∫₀^{∞·e^{iθ₀}} e^{−ξ/t}·B(ξ) dξ.

Proof. Fix θ₀ and put f(z) = P·J(e^{iθ₀}·z). For z ∈ C_R: |z| < R and |arg z| < 90°, so t = e^{iθ₀}z lies in the open sector −210° < arg t < 10°, |t| < t₀, where P·J is analytic (Corollary of v16) and Theorem 4 holds. So the hypotheses hold with cₖ = aₖ·e^{ikθ₀} and σ = 240. The Borel transform of f is B(e^{iθ₀}ζ); substituting ξ = e^{iθ₀}ζ gives (ii), and the union of the strips over θ₀ gives (i). ∎

The union of the discs in (ii) contains every t with −210° < arg t < 10° and |t| < R·cos(dist(arg t, [−120°, −80°])).

#### [v17 §5.2] What this says about the earlier formulas

- **|q| < 1 (proved).** For Im x > 0, −90° < arg t < 10° and t in a disc of (ii), the Corollary of v12 becomes

  ψ₁\*(q) = 𝒮_{θ₀}[Φ₁](t) + (κ/2)·[𝓑₂·p·W₃(p) − 𝓑₄·W(p)],  f₁\*(q) = 2·𝒮_{θ₀}[Φ₁](t) − κ·[𝓑₂·F_b(p) − 𝓑₄·F_a(p)],

  ψ₀\*(q) = 𝒮_{θ₀}[Φ₀](t) − i·q^{1/4}·(κ/2)·[𝓑₁·p·W₃(p) + 𝓑₃·W(p)],  f₀\*(q) = 2·𝒮_{θ₀}[Φ₀](t) + i·q^{1/4}·κ·[𝓑₁·F_b(p) + 𝓑₃·F_a(p)],

  with any θ₀ ∈ [−120°, −80°] within 90° of arg t, and Φ₁, Φ₀ the perturbative series of ψ₁\*, ψ₀\* at q → 1. This is the three-sector formula of v10 with 𝒯 = 𝒮_{θ₀}[Φ₁].
- **The caveat.** v10 defined 𝒯 as the Borel sum *along the direction of t* (for −90° < arg t < −17.1°). For −90° < arg t ≤ −80° that is literally what is proved. For −80° < arg t < −17.1° the two definitions agree if B has no singularity in −80° < arg ξ ≤ arg t and is exponentially bounded there. Numerically that is the case (the nearest singularities are S and its conjugate, at arg ∓17.1°; v13). It is not proved: Theorem 5 only reaches the directions −120° … −80°.
- **|q| = 1 and |q| > 1 (proved for P·J).** P·J is the Borel sum of Φ also for real x and for Im x < 0 (−210° < arg t ≤ −90°). What is missing on that side is the other half: Theorem A of v12 for |q| > 1.

---

## Part 6. The state integrals and the forward law (Theorems III, IV)

Lemmas B and L are in Part 2.

*From `docs/addendum-v19.md`, §1; text unchanged.*

### [v19 §1] Setting

x → 0 with 0 ≤ arg x ≤ 60°, |x| ≤ ¼; t = −2πix, so arg t = −90° + ε with ε = arg x, τ = tan ε. Φ_b is Faddeev's integral (v16 §1); for real x > 0 the state integrals are defined with it. As in v18, with w = 2π·b·y and u = a′/b:

  I(a′/b) = (1/(2πb))·∫ e^{Q(w)/t}/Φ_b dw,  Q(w) = −w²/2 − 2πi·a′·w,  V(w) = Q(w) + Li₂(−e^{w}).

**Height.** |e^{V/t}| = exp(cos ε·H_τ(w)/|t|), H_τ = −Im V + τ·Re V. (Re(V/t) = Re(i·e^{−iε}·V)/|t|.)

V′(w) = −w − 2πi·a′ − log(1 + e^{w}), V″(w) = −(1 + 2e^{w})/(1 + e^{w}). Li₂(−e^{w}) is principal in the strip |Im w| ≤ π; for Re w < 0 it is analytic and 2πi-periodic.

*From `docs/addendum-v19.md`, §3–§6; text unchanged.*

### [v19 §3] Statement for the integrals

**Theorem 6.** Uniformly for |x| → 0:

- (decaying, 0 ≤ arg x ≤ 60°) I(1/(2b)) = (t/2π)^{1/2}·3^{−1/4}·e^{+iπ/12}·b^{−1}·exp((−5π²/12 + i·Cl₂)/t)·(1 + O(t)), and I(b + 1/(2b)) = e^{−2πi/3}·(the same)·(1 + O(t));
- (growing, 0 ≤ arg x ≤ 11°) I(−1/(2b)) = (t/2π)^{1/2}·3^{−1/4}·e^{−iπ/12}·b^{−1}·exp((−5π²/12 − i·Cl₂)/t)·(1 + O(t)), and I(b − 1/(2b)) = e^{+2πi/3}·(the same)·(1 + O(t)),

with principal roots. For Im x > 0 this is item 1 of the summary, by Theorem C of v12: 𝓑₂,₄ = I(±1/(2b))/(c·p^{1/8}), 𝓑₁,₃ = ±i·q^{1/4}·I(b ± 1/(2b))/(c·p^{1/8}), c·p^{1/8} = 2^{−1/2}·e^{−iπ/4}·e^{t/24}·e^{−5π²/(12t)}, and t^{1/2}/b = (2π)^{1/2}·e^{−iπ/4}.

**Theorem 7.** Let 0 < a < 1, W the root of W² + W = e^{−2πia} with Im W < 0, w₁ = Log W. There is ε(a) > 0 such that for 0 ≤ arg x ≤ ε(a), uniformly,

  I(a/b) = (1/(2πb))·(2πt·(1 + W)/(1 + 2W))^{1/2}·e^{V(w₁)/t}·(1 + O(t)),  −Im V(w₁) = −Vol(z) < 0,

with the principal root. For a = ½ one can take ε = 60°.

---

### [v19 §4] The decaying integrals (proof of Theorem 7 and of the first half of Theorem 6)

a′ = a ∈ (0, 1). Write η₁ = Im w₁ = arg W ∈ (−π, 0), r₁ = Re w₁.

- **Critical point.** W(1 + W) = e^{−2πia} and both arguments lie in (−π, 0), so arg W + arg(1 + W) = −2πa exactly: w₁ + 2πia + Log(1 + W) = 0.
- **Contour.** The horizontal line w = r + iη₁. It lies in the strip, at distance π − |η₁| from its edges.
- **Height at τ = 0.** d/dr Im V(r + iη₁) = Im V′ = −η₁ − 2πa − arg(1 + e^{r+iη₁}) = arg(1 + W) − arg(1 + e^{r−r₁}·W). For Im W < 0 the argument of 1 + λW decreases strictly from 0 to arg W as λ runs from 0 to ∞. So the derivative is negative for r < r₁, positive for r > r₁, has a simple zero at r₁, and tends to arg(1 + W) < 0 and to arg(1 + W) − arg W > 0 at the two ends. Hence Im V − Im V(w₁) ≥ c·min((r − r₁)², |r − r₁|): (L1) for τ = 0.
- **τ > 0.** Re V(r + iη₁) is smooth, has zero derivative at r₁ (critical point) and tends to −∞ quadratically at both ends, so Re V − Re V(w₁) ≤ M·min((r − r₁)², |r − r₁|) for some M. (L1) follows for τ ≤ c/(2M). For a = ½ (W = e^{−2πi/3}): d/dr Re V = −r − ½·ln(1 − e^{r} + e^{2r}), positive for r < 0 and negative for r > 0, so M = 0 and (L1) holds for every τ ≥ 0.
- **(L2).** Lemma B (i) on the line and (ii) near w₁; e^{Q/t} is exact. For I(b + 1/(2b)) the integrand has the extra factor e^{w}, which is W at the saddle and at most e^{|w|} elsewhere.
- **Moving the contour** from L of v12 (through 0, slope −δ) to the line: the zeros of Φ_b lie at w = t(k + ½) − iπ(2l + 1), k, l ≥ 0, below L and below the line Im w = −π; the region between L and the horizontal line contains none of them, and the integrand tends to 0 at infinity there (for Im x > 0 by the Gaussian factor, v12 §2; for real x by Lemma B and its mirror image: on the left the size is e^{r(η + 2πa)/|t|} with η ≥ η₁ > −2πa, on the right e^{r(2η + 2πa)/|t|} with η ≤ η₁ < −πa).
- **Result.** Lemma L with d = 1 and V″(w₁) = −(1 + 2W)/(1 + W). For a = ½: (1 + W)/(1 + 2W) = e^{iπ/6}/√3 and V(w₁) = −5π²/12 + i·Cl₂. ∎

The two valley conditions used above, −2πa < η₁ < −πa, are automatic: the first is arg(1 + W) < 0, the second arg W < arg(1 + W).

---

### [v19 §5] The growing integrals at z = −1 (proof of the second half of Theorem 6)

a′ = −½: V(w) = −w²/2 + iπw + Li₂(−e^{w}), saddle w₀ = 2πi/3 (W = e^{2πi/3}), V(w₀) = −5π²/12 − i·Cl₂, H_τ(w₀) = Cl₂ − τ·5π²/12 ≥ 0.19 for τ ≤ 0.2.

**Contour Γ** = Γ₁ ∪ Γ₂ ∪ Γ₃:

- Γ₁: w = iπ − ρ·e^{−iπ/6}, ρ from ∞ to 0 (in the left half-plane, above the strip);
- Γ₂: w = iη, η from π to 0;
- Γ₃: w = r, r from 0 to ∞.

**Heights, exactly.**

- Γ₃: Li₂(−e^{r}) is real and ≤ −π²/12, so Im V = π·r and Re V ≤ −π²/12 − r²/2. Hence H_τ ≤ −π·r.
- Γ₂: Im V(iη) = −Cl₂(π − η) and Re V(iη) = ¾η² − πη − π²/12 (Re Li₂ on the unit circle is a quadratic). So H_τ(iη) = g(η) := Cl₂(π − η) + τ·(¾η² − πη − π²/12), and g′(η) = f(η) − τ·(π − 3η/2) with f(η) = ln(2cos(η/2)). f is concave with f(0) = ln 2, f(2π/3) = 0, f′(2π/3) = −√3/2. A concave function lies above its chord and below its tangent: f(η) ≥ (ln 2/π)·(π − 3η/2) on [0, 2π/3] and f(η) ≤ −(1/√3)·(3η/2 − π) on [2π/3, π). So for τ ≤ 0.2 < ln 2/π = 0.2206: g′ ≥ 0.0206·(π − 3η/2) below 2π/3 and g′ ≤ −0.377·(3η/2 − π) above, hence g(η) − g(2π/3) ≤ −0.015·(η − 2π/3)².
- Γ₁: Im Q = ρ²·sin(π/6)·cos(π/6) = 0.433·ρ² and Re Q = −π²/2 − ρ²/4. Further −e^{w} = e^{−ρcos(π/6)}·e^{iρ/2} has modulus < 1 and, for ρ ≤ 2π, argument in [0, π], where Im Li₂ ≥ 0; for ρ ≥ 2π, |Li₂(−e^{w})| ≤ (π²/6)·e^{−2π·cos(π/6)} < 0.01. And |Re Li₂| ≤ π²/6. So Im V ≥ 0.433·ρ² − 0.01, Re V ≤ −π²/3, and H_τ ≤ 0.01 − 0.433·ρ².

Together with H_τ(w₀) ≥ 0.19 this is (L1) on all of Γ for 0 ≤ τ ≤ 0.2, that is 0 ≤ arg x ≤ 11.3°.

**(L2).** Γ₁: Lemma B (iii) with δ₁ = 30°. Γ₂, Γ₃: Lemma B (i); near w₀, (ii) with ε₁ = π/6. For I(b − 1/(2b)) the extra factor e^{w} is e^{2πi/3} at the saddle, has modulus ≤ 1 on Γ₁, Γ₂ and e^{r} on Γ₃.

**Moving the contour** from L to Γ: on the left the region between the two rays contains only poles of Φ_b (zeros of the integrand); on the right the wedge between L and the real axis lies above L and contains no zero of Φ_b; decay at infinity as in §4.

**Result.** Lemma L on Γ₂, d = −i, V″(w₀) = −√3·e^{iπ/6}: the root is −i·(2πt/V″(w₀))^{1/2} = (2πt)^{1/2}·3^{−1/4}·e^{−iπ/12} (principal t^{1/2}). ∎

Remark. Along the imaginary axis the two saddles are the minimum and the maximum of one Clausen function: Im V(iη) = −Cl₂(π − η) on (−π, π). The same axis is the ascent path of the decaying integral, which is why that one uses the horizontal line.

---

### [v19 §6] Consequences

#### [v19 §6.1] Theorem 8: the forward laws (proved)

By v12, ψ₁\* = P·J(1/(2b)) + P·I(1/(2b))·p·W₃(p) − P·I(−1/(2b))·W(p) (since (κ/2)·𝓑 = P·I), with P = 2^{−1/2}·e^{π²/(12t) + t/6}/b. With Theorem 6 and 1/x = −2πi/t:

  P·I(1/(2b)) = +A·e^{−5πi/12}·t^{−1/2}·e^{−S/t}·(1 + O(t)),  P·I(−1/(2b)) = −A·e^{+5πi/12}·t^{−1/2}·e^{−S̄/t}·(1 + O(t)),

with A = (π/√3)^{1/2}, S = π²/3 − i·Cl₂ (checked numerically, §7). For 0 < arg x ≤ 10° the direction of t lies in (−90°, −80°], where P·J = 𝒮[Φ] by Theorem 5 of v17. Hence, for |t| < R:

- ψ₁\* = 𝒮[Φ₁] + A·t^{−1/2}·[ e^{5πi/12}·e^{−S̄/t}·W(p)·(1 + O(t)) + e^{−5πi/12}·e^{−S/t}·p·W₃(p)·(1 + O(t)) ],
- f₁\* = 2𝒮[Φ₁] − 2A·t^{−1/2}·[ e^{5πi/12}·e^{−S̄/t}·F_a(p)·(1 + O(t)) + e^{−5πi/12}·e^{−S/t}·F_b(p)·(1 + O(t)) ],
- ψ₀\* = 𝒮[Φ₀] − A·t^{−1/2}·[ e^{iπ/12}·e^{−S̄/t}·W(p)·(1 + O(t)) + e^{−iπ/12}·e^{−S/t}·p·W₃(p)·(1 + O(t)) ],
- f₀\* = 2𝒮[Φ₀] + 2A·t^{−1/2}·[ e^{iπ/12}·e^{−S̄/t}·F_a(p)·(1 + O(t)) + e^{−iπ/12}·e^{−S/t}·F_b(p)·(1 + O(t)) ],

The first exponential grows, the second decays; each O(t) is a factor of an exact coefficient and does not depend on p. The constants S̄ = π²/3 + i·Cl₂ and A·e^{5πi/12} are those of v5. The exponent S of the decaying term is the Borel singularity of Φ found in v13.

(For ψ₀\*: q^{1/2}·P·I(b ± 1/(2b)) = e^{∓2πi/3}·P·I(±1/(2b))·(1 + O(t)), and e^{∓2πi/3}·e^{∓5πi/12} = −e^{∓iπ/12}.)

#### [v19 §6.2] What remains exact and what is asymptotic

The formulas of v12 are exact. Theorem 8 replaces the three coefficient functions by their leading behaviour. On a horocycle (p fixed) it gives the two-term law of v5 with relative error O(t) in each sector; radially (p → 0) it gives ψ₁\* = 𝒮 + A·e^{5πi/12}·t^{−1/2}·e^{−S̄/t}·(1 + O(t)).

#### [v19 §6.3] The z-family

- Decaying sector: Theorem 7, all 0 < a < 1.
- Growing sector for a ≠ ½: not proved. The vertical line through the saddle is again the local descent path (d/dη Im V = −r₀ − ln|1 + e^{r₀+iη}| is monotone in η ∈ (0, π)), and the real axis works on the right (Im V = 2π|a′|·r). The left part needs a case analysis in a, because the left valley starts at Im w > 2π|a′|, which is above the strip for |a′| > ½, and the saddle is not on the imaginary axis.
- A statement at roots of unity (the volume rule as measured) needs in addition the radial limit of the formula of v18 and a bound for J(a/b) at real x for a ≠ ½.

*From `docs/addendum-v26.md`, §2–§3 (Lemma L completed, the e^{w}-versions, Theorem 7 with ε(a) = 60°); text unchanged.*

### [v26 §2] Lemma L and the integrals with the factor e^{w}

**Lemma L (v19 §2), completed.** Add to the hypotheses: every piece of Γ is a segment, a ray or a line, and |t| ≤ t₁ := c·cos ε₁/(2C), where ε ≤ ε₁ is the range of arg x, c is the constant of (L1) and C the constant in the factor e^{C|w|} of (L2). If (L2) holds without that factor there is no condition on |t| beyond the one of the bounds used.

Why. By (L1), (L2), on Γ: |E| ≤ C·e^{C|w\*|}·exp(cos ε·H_τ(w\*)/|t|)·exp(C·|w − w\*| − cos ε·c·min(|w − w\*|², |w − w\*|)/|t|). For |t| ≤ t₁ and |w − w\*| ≥ 1 the last exponent is at most −(c·cos ε/(2|t|))·|w − w\*|, which gives convergence and the tail estimate of the proof. Without the condition the integral over Γ need not converge.

**The e^{w}-versions at real x.** For I(b ± 1/(2b)) the integrand has the extra factor e^{w}, and (L2) holds with e^{|w|} (C = 1). At real x (t = −i|t|) the sizes on the right are:

- *Decaying integral (v19 §4), line Im w = η₁.* Between L and the line, and on the line, the size is e^{r·(1 + (2η + 2πa)/|t|)} with η ≤ η₁. So the integrand tends to 0 on the right, and the line is an admissible contour, exactly when |t| < −(2η₁ + 2πa). For a = ½ (η₁ = −2π/3) this is **|t| < π/3**, that is X ≥ 7 at x = 1/X.
- *Growing integral (v19 §5), contour Γ.* On Γ₃ and in the wedge between L and the real axis the size is at most e^{r·(1 − π/|t|)}. The condition is **|t| < π**, that is X ≥ 3.
- On the left the factor e^{w} only helps. For Im x > 0 the Gaussian factor gives decay for every t (v12 §2), and there is no condition.

v19 §4 proved the decay at infinity for the integrals without e^{w} and for the others only stated (L2). The condition is a missing hypothesis in "Moving the contour", as F2 says.

**Consequences.** None for the theorems: Theorem 6 is a statement for |x| → 0; Theorems 9, 11, 18 are stated with the integrals over L, which exist for every real x in their range. For numerical work: at small X the integrals I(b ± 1/(2b)) must be taken over L (the scripts of v20, v21 do this for X ≤ 8).

---

### [v26 §3] Theorem 7 with an explicit sector

Notation of v19 §4: 0 < a < 1, W the root of W² + W = e^{−2πia} with Im W < 0, w₁ = r₁ + iη₁ = Log W, V(w) = −w²/2 − 2πia·w + Li₂(−e^{w}).

**Lemma V.** On the line w = r + iη₁, Re V(w) < Re V(w₁) for r ≠ r₁.

Proof. Let R(r) = Re V(r + iη₁). Since V′(w) = −w − 2πia − log(1 + e^{w}),

  R′(r) = −r − ln|1 + e^{r + iη₁}| = −½·ln g(e^{r}),  g(λ) = λ²·|1 + λe^{iη₁}|² = λ⁴ + 2cλ³ + λ²,  c = cos η₁.

g(|W|) = |W(1 + W)|² = 1. The claim is that g < 1 on (0, |W|) and g > 1 on (|W|, ∞); then R increases up to r₁ and decreases after it. Now g(0) = 0, g(∞) = ∞ and g′(λ) = 2λ·(2λ² + 3cλ + 1).

- If c ≥ −2√2/3, the quadratic 2λ² + 3cλ + 1 has no sign change on λ > 0 (its discriminant is 9c² − 8 ≤ 0, or c ≥ 0). So g increases strictly and takes the value 1 once.
- If −1 < c < −2√2/3, g increases on (0, λ₋), decreases on (λ₋, λ₊) and increases on (λ₊, ∞), with λ± = (3|c| ± (9c² − 8)^{1/2})/4 and ½ < λ₋ < 2^{−1/2}. At the local maximum, g(λ₋) = λ₋²·(1 − 2|c|λ₋ + λ₋²) ≤ ½·(1 − (4√2/3)·½ + ¼) < 0.16, because 1 − (4√2/3)λ + λ² decreases on [½, 2^{−1/2}]. So g < 1 on (0, λ₊], and on [λ₊, ∞) it increases strictly to ∞: again the value 1 is taken once. ∎

**Theorem 7′.** Theorem 7 of v19 holds with ε(a) = 60° for every 0 < a < 1: uniformly for 0 ≤ arg x ≤ 60°, |x| → 0,

  I(a/b) = (1/(2πb))·(2πt·(1 + W)/(1 + 2W))^{1/2}·e^{V(w₁)/t}·(1 + O(t)).

Proof. v19 §4 proves (L1) on the line for τ = 0: Im V − Im V(w₁) ≥ c·min((r − r₁)², |r − r₁|). By Lemma V, for every τ ≥ 0: H_τ(w) − H_τ(w₁) = −(Im V − Im V(w₁)) + τ·(Re V − Re V(w₁)) ≤ −c·min((r − r₁)², |r − r₁|). So (L1) holds for all τ ≥ 0 with the constant of τ = 0 (in the words of v19: M = 0). The rest of the proof of v19 §4 is for 0 ≤ arg x ≤ 60°, the range of Lemma B. ∎

Remarks.

- R is concave on the whole line exactly in the first case of the proof, cos η₁ ≥ −2√2/3, which is a ≤ 0.807674. For larger a it is not concave but still has its maximum at w₁.
- The constants depend on a and degenerate as a → 0 and a → 1 (the stress test: c₀ = 0.0048 at a = 0.01; at a → 1 the line approaches the edge of the strip).
- This removes the item "Theorem 7 gives no explicit ε(a)" from the list for referees. Theorem 13 of v22 (the growing integral for general a) keeps its own unspecified ε(a).

---

## Part 7. Roots of unity (Theorem V)



*From `docs/addendum-v20.md`, §1–§3 (even X); text unchanged.*

### [v20 §1] Three lemmas

Throughout X = 2M is even, q₀ = e^{2πi/X}, so q₀^M = −1.

**Lemma R1 (radial limit of the q-series).** As ρ → 1⁻: ψ₁\*(ρq₀) → Σ_{n<M} (−1)ⁿ q₀^{n(n+1)/2} (−q₀;q₀)ₙ, and likewise for ψ₀\*.

Proof. Let q = ρq₀, δ = 1 − ρ, X·δ ≤ 1. It suffices to show Σ_{n≥M} |q|^{n(n+1)/2}·|(−q;q)ₙ| → 0; the terms with n < M are continuous.

- *Blocks.* For m ≥ 0 let B_m = Π_{j=1}^{X} (1 + q^{mX+j}). Put c = ρ^{(m+1)X} and c_j = ρ^{mX+j}, so 0 ≤ c_j − c ≤ X·δ. Since Π_{j=1}^{X} (1 + c·q₀ʲ) = 1 − c^X (X even), expanding the product gives |B_m − (1 − c^X)| ≤ (2 + Xδ)^X − 2^X ≤ X²·3^{X−1}·δ. And 0 ≤ 1 − c^X ≤ min(1, (m + 1)·X²·δ). So with C = X²·(1 + 3^{X−1}): |B_m| ≤ min(1 + Cδ, (m + 1)·Cδ).
- *M ≤ n < X.* (−q;q)ₙ contains 1 + q^M = 1 − ρ^M ≤ Mδ; the other factors are at most 2. These finitely many terms tend to 0.
- *n ≥ X.* Write n = BX + r, 0 ≤ r < X, B ≥ 1. Then |(−q;q)ₙ| ≤ 2^X·Π_{m<B} |B_m|. Let B₁ = ⌊1/(2Cδ)⌋. For B ≤ B₁ every factor (m + 1)Cδ is at most ½, so the product is at most min(Cδ, 2^{−B}); the sum over B ≤ B₁ tends to 0. For B > B₁ the product is at most 2^{−B₁}·e^{Cδ·B}, and |q|^{n²/2} ≤ e^{−δB²X²/2}; the sum over B of e^{CδB − δB²X²/2} is O(δ^{−1/2}), so this part is O(δ^{−1/2}·2^{−1/(2Cδ)}) → 0. Each B stands for X values of n. ∎

**Lemma R2 (the dual series at p = 1).** As p → 1 with |p| < 1 and |1 − p| ≤ (3/2)·(1 − |p|): W(p) → 1 and W₃(p) → 1.

Proof. W(p) − 1 = Σ_{j≥1} (−1)ʲ p^{j²+j} (p;p²)ⱼ. Let δ = 1 − |p|, K = 3/2. |1 − pᵏ| ≤ k·Kδ and |1 − pᵏ| ≤ 1 + |p|ᵏ. Let j₁ = ⌊1/(2Kδ)⌋.

- j ≤ j₁: |(p;p²)ⱼ| ≤ (2j − 1)!!·(Kδ)ʲ ≤ min(Kδ, √2·e^{−j}) (each factor (2i + 1)Kδ ≤ 1, and (2j − 1)!! ≤ √2·(2j/e)ʲ). The sum over these j tends to 0.
- j > j₁: |(p;p²)ⱼ| ≤ √2·e^{−j₁}·Π_{i≥j₁} (1 + |p|^{2i+1}) ≤ √2·e^{−j₁}·exp(|p|^{2j₁+1}/(1 − |p|²)) ≤ C·exp(−(1/δ)·[1/(2K) − e^{−1/K + 2δ}/(2 − δ)]). For K = 3/2 the bracket tends to 1/3 − 0.257 = 0.077 > 0. With Σ_j |p|^{j²} = O(δ^{−1/2}) this part tends to 0.

W₃ has the extra factors p^{2j} of modulus < 1. ∎

**Lemma R3 (continuity at real x).** Let 0 < x₀ ≤ ¼. The integrals I(u), J(u) of v12 over L = e^{−10°·i}·ℝ, with Faddeev's Φ_b, for u ∈ {±1/(2b), b ± 1/(2b)}, are continuous in x on {Im x ≥ 0} near x₀.

Proof. Φ_b(y) is analytic in (b, y) off its poles and zeros, and the zeros stay below L. By Lemmas F and D of v16 and the inversion relation, on L one has |1/Φ_b| ≤ C on the left half and |1/Φ_b| ≤ C·|e^{−w²/(2t)}| on the right half, with C uniform near x₀ (|t| is bounded away from 0 there, so the terms Li₂/t are bounded). The Gaussian factor e^{−w²/(2t)} of the integrand has modulus exp(−ρ²·sin 20°/(2|t|)) on L, and the remaining factors are e^{O(ρ)}. Dominated convergence. (For X = 4 the point x₀ = ¼ is on the edge of the range |x| ≤ ¼ of v16; the condition actually used there is Im(i/x) ≥ 2, which holds in a neighbourhood of ¼.) ∎

**The path.** For x = 1/X + iε, ε → 0⁺: q = q₀·e^{−2πε} is radial, and p = e^{−πi/x} has |p| = exp(−πεX²/(1 + ε²X²)) and arg p = πε²X³/(1 + ε²X²) (mod 2π; X even). So p → 1 and |1 − p|/(1 − |p|) → 1: Lemma R2 applies.

---

### [v20 §2] Theorem 9

**Theorem 9.** For even X ≥ 4 the two identities of the summary hold.

Proof. For Im x > 0, v12 gives (with (κ/2)·𝓑₂ = P·I(1/(2b)), (κ/2)·𝓑₄ = P·I(−1/(2b)), −i·q^{1/4}·(κ/2)·𝓑₁ = q^{1/2}·P·I(b + 1/(2b)), i·q^{1/4}·(κ/2)·𝓑₃ = q^{1/2}·P·I(b − 1/(2b)), from Theorem C there)

  ψ₁\*(q) = P·[ J(1/(2b)) + I(1/(2b))·p·W₃(p) − I(−1/(2b))·W(p) ],

  ψ₀\*(q) = q^{1/2}·P·[ J(b + 1/(2b)) + I(b + 1/(2b))·p·W₃(p) − I(b − 1/(2b))·W(p) ].

The integrals with the product formula for Φ_b and any admissible δ are those with Faddeev's integral and δ = 10° (v16 Lemma Φ; v14 §2.1), and P = 2^{−1/2}·e^{π²/(12t) + t/6}/b. Let x = 1/X + iε, ε → 0⁺, and use Lemmas R1, R2, R3. ∎

Remarks.

- f₁\*, f₀\* have no such limit: their terms have (−q;q)ₙ in the denominator, and on the other side F_a(p), F_b(p) have poles at p = 1.
- The formula holds for every even X ≥ 4, not only asymptotically: at X = 4, ψ₁\*(i) = 2 − i is the sum of P·J = 0.276 − 0.146i, P·I(1/(2b)) = −0.539 − 0.092i and −P·I(−1/(2b)) = 2.263 − 0.763i.

---

### [v20 §3] Theorem 10

**Theorem 10.** As X → ∞ through even integers, the two asymptotic formulas of the summary hold.

Proof. In Theorem 9, by Theorem 6 of v19 at real x (arg x = 0 is included there):

- −P·I(−1/(2b)) = A·e^{5πi/12}·t^{−1/2}·e^{−S̄/t}·(1 + O(t)), A = (π/√3)^{1/2}, S̄ = π²/3 + i·Cl₂ (v19 §6.1);
- P·I(1/(2b)) = O(|t|^{−1/2}·e^{−Cl₂/|t|});
- |P·J(1/(2b))| ≤ C (v16, Theorem 3 with K = 0; in fact P·J = 1/3 + O(1/X)).

At t = −2πi/X: t^{−1/2} = (X/2π)^{1/2}·e^{iπ/4} and e^{−S̄/t} = e^{−i·S̄·X/2π} = e^{Cl₂·X/2π}·e^{−iπX/6}. So the leading term is (π/√3)^{1/2}·(X/2π)^{1/2}·e^{5πi/12 + iπ/4}·e^{Cl₂X/2π − iπX/6} = 3^{−1/4}·(X/2)^{1/2}·e^{2πi/3}·e^{Cl₂X/2π − iπX/6}, and the other two terms are smaller by e^{−Cl₂X/2π}. For ψ₀\*: −e^{πi/X}·P·I(b − 1/(2b)) = e^{2πi/3}·(−P·I(−1/(2b)))·(1 + O(t)) by Theorem 6. ∎

**Corollary (figure-eight).** The Kashaev invariant of 4₁ is ⟨4₁⟩_N = Σ_{k<N} Π_{j≤k} |1 − e^{2πij/N}|², and ⟨4₁⟩_N = 3^{−1/4}·N^{3/2}·e^{N·Vol(4₁)/2π}·(1 + O(1/N)) with Vol(4₁) = 2·Cl₂(π/3) (**quoted** from memory: Andersen–Hansen 2006; I did not look it up). With Theorem 10 at N = X/2:

  |ψ₁\*(e^{2πi/X})| / ⟨4₁⟩_{X/2} = (2/X)·(1 + O(1/X)).

Both have the growth rate Cl₂(π/3)/2π per unit of X and the torsion factor 3^{−1/4}; they differ by one power of N. (Verified directly in v18 §2.5: the ratio times X/2 is 0.99887, 0.99943, 0.99972 at X = 802, 1602, 3202.)

*From `docs/addendum-v21.md`, §1–§3 (odd X); text unchanged.*

### [v21 §1] The dual series at p = −1

For Im x > 0 let p = e^{−πi/x}. With x = 1/X + iε, X odd, ε → 0⁺ (v20 §1): p = −ρ·e^{iφ} with ρ = exp(−πεX²/(1 + ε²X²)) and φ = πε²X³/(1 + ε²X²). Put δ = 1 − ρ. Then δ → 0 and 0 < φ ≤ δ² for small ε (φ/δ² → 1/(πX)).

Recall (v4, proved): (V) W + F_a = ψ·P₊, (VI) p·W₃ + F_b = ψ·P₋, with F_a = Σ_{n≥1} p^{2n²−n}/(p;p²)ₙ, F_b = Σ_{n≥1} p^{2n²−3n+1}/(p;p²)ₙ, P₊ = Σ p^{2k²+k}/(p²;p²)ₖ, P₋ = Σ p^{2k²−k}/(p²;p²)ₖ, ψ(p) = Σ_{n≥0} p^{n(n+1)/2}.

**Lemma R2′.** Along this path: F_a(p) → −⅓, F_b(p) → ⅓, ψ(p)·P_±(p) → 0. Hence W(p) → ⅓ and p·W₃(p) → −⅓.

Proof.

- *Values at p = −1.* (p;p²)ₙ = 2ⁿ, p^{2n²−n} = (−1)ⁿ, p^{2n²−3n+1} = (−1)^{n+1}. So F_a(−1) = Σ_{n≥1} (−½)ⁿ = −⅓ and F_b(−1) = −Σ_{n≥1} (−½)ⁿ = ⅓.
- *Lower bound for the Pochhammer symbol.* For odd k: 1 − pᵏ = 1 + ρᵏ·e^{ikφ}. If k ≤ 1/δ then kφ ≤ δ and ρᵏ ≥ ¼, so |1 − pᵏ| ≥ 1 + ρᵏ·cos(kφ) ≥ 9/8. If k > 1/δ then |1 − pᵏ| ≥ 1 − ρᵏ with ρᵏ ≤ e^{−1}, and Π_{k odd > 1/δ} (1 − ρᵏ) ≥ exp(−1.59·Σ_{k odd > 1/δ} ρᵏ) ≥ exp(−0.3/δ) for small δ. Let n₁ be the number of odd k ≤ 1/δ, so n₁ ≥ 1/(2δ) − 1. Then |(p;p²)ₙ| ≥ (9/8)ⁿ for n ≤ n₁ and |(p;p²)ₙ| ≥ e^{−0.3/δ} for n > n₁.
- *F_a, F_b.* For n ≤ n₁ the terms are bounded by (8/9)ⁿ and converge to the terms at p = −1: dominated convergence. For n > n₁: ρ^{2n²−n} ≤ e^{−n(1−δ)} and ρ^{(2n−1)(n−1)} ≤ e^{−(n−1)(1−δ)} (because n ≥ 1/(2δ), so 2n − 1 ≥ 1/δ − 1), so the rest is at most C·e^{0.3/δ}·e^{−(1−δ)/(2δ)} → 0.
- *P_±.* |1 − p^{2j}| ≥ 1 − Q^{j} with Q = ρ², so |P₊| ≤ Σ Q^{k²}/(Q;Q)ₖ = G(Q) and |P₋| ≤ Σ Q^{k²−k}/(Q;Q)ₖ = G(Q) + H(Q) ≤ 2G(Q) (Rogers–Ramanujan functions). From the product, log G(Q) = Σ_r (Q^r + Q^{4r})/(r·(1 − Q^{5r})) ≤ 2·log(1/(1 − Q)) + π²/(15λ) with λ = −log Q ≥ 2δ. So |P_±| ≤ 2·δ^{−2}·e^{π²/(30δ)}.
- *ψ.* With p = e^{2πiτ}, τ = ½ + s, s = (φ + iℓ)/(2π), ℓ = −log ρ: ψ(p) = ½·e^{−πiτ/4}·θ₂(τ), θ₂(τ) = (−iτ)^{−1/2}·θ₄(−1/τ), and −1/τ = −2 + σ with σ = 4s/(1 + 2s). θ₄ has period 2, and θ₄(σ) = (−iσ)^{−1/2}·θ₂(−1/σ) with |θ₂(u)| ≤ 4·e^{−π·Im u/4} for Im u ≥ 1. Here Im(−1/σ) = Im s/(4|s|²) = (π/2)·ℓ/(φ² + ℓ²) ≥ (π/2)·(1 − 2δ)/δ. So |ψ(p)| ≤ C·δ^{−1/2}·e^{−π²/(8δ)}.
- Together |ψ·P_±| ≤ C·δ^{−5/2}·exp(−(π²/8 − π²/30)/δ) → 0. ∎

---

### [v21 §2] Theorem 11

**Theorem 11.** For odd X ≥ 5 the identity of the summary holds, and likewise

  lim_{r→1⁻} ψ₀\*(r·q₀) = A₀/3 = e^{πi/X}·P·[ J(b + 1/(2b)) − ⅓·I(b + 1/(2b)) − ⅓·I(b − 1/(2b)) ].

Proof. For x = 1/X + iε the formulas of v12 (in the form of v20 §2) hold. As ε → 0⁺: q = q₀·e^{−2πε} is radial, so the left sides tend to A/3, A₀/3 (v3 §1, Theorem; q₀ is a primitive X-th root of unity, X odd). On the right, the integrals are continuous at x = 1/X (v20, Lemma R3; 1/X ≤ ¼), and W(p) → ⅓, p·W₃(p) → −⅓ (Lemma R2′). ∎

Equivalently A = P·[3·J(1/(2b)) − I(1/(2b)) − I(−1/(2b))]: at odd order the full period sum A is three times the kernel integral minus the two state integrals, while at even order the half-period sum is the kernel integral plus one state integral minus the other (v20).

---

### [v21 §3] Theorem 12

**Theorem 12.** As X → ∞ through odd integers, the asymptotic formulas of the summary hold for the radial limits A/3 and A₀/3.

Proof. As for Theorem 10 of v20: −⅓·P·I(−1/(2b)) = ⅓·A′·e^{5πi/12}·t^{−1/2}·e^{−S̄/t}·(1 + O(t)) with A′ = (π/√3)^{1/2} (v19, Theorem 6 at real x), the other two terms are O(1) and O(|t|^{−1/2}e^{−Cl₂/|t|}), and at t = −2πi/X the leading term is ⅓·3^{−1/4}·(X/2)^{1/2}·e^{2πi/3}·e^{Cl₂X/2π − iπX/6}. For ψ₀\* the growing integral carries the extra factor e^{2πi/3}. ∎

*From `docs/addendum-v25.md`, §1 (X = 1, 2, 3); text unchanged.*

### [v25 §1] Continuity at real x up to x = 1

**Lemma R3′.** Let 0 < x₀ ≤ 1 and let u(x) be continuous and bounded near x₀. The integrals I(u), J(u) of v12 over L = e^{−10°·i}·ℝ, with Faddeev's Φ_b, are continuous in x on {Im x ≥ 0} near x₀.

Proof. The proof of Lemma R3 (v20 §1) needs, for x in a neighbourhood N of x₀ in {Im x ≥ 0}, a constant C with |1/Φ_b| ≤ C on the left half of L and |1/Φ_b| ≤ C·|e^{−w²/(2t)}| on the right half. It took this from Lemmas F and D of v16, which were stated for |x| ≤ ¼. Their proofs use the size of x as follows.

- *Lemma F* (|Λ(w) + Li₂(−e^{w})/t| ≤ C₁·e^{Re w/2} for Re w ≤ 0, |Im w| ≤ π). The contour of Faddeev's integral is moved from ℝ + i0 to Im σ = ½. (i) No pole is crossed if the poles σ = i·n/x of r₁(πxσ) lie above the line: Im(i/x) = Re(1/x) > ½. (ii) On the line, u = πxσ must stay away from the poles iπn (n ≠ 0) of 1/sinh u: for real x the distance is at least π·(1 − x/2). For x ∈ N both hold, with Re(1/x) ≥ 0.9 and distance ≥ 1.5 if N is small. The other steps (|r₁(u)| ≤ C·min(|u|, 1/|u|) on the line, and the integral of min(|x|, 1/(|x|·|σ|²))) do not use the size of x. So Lemma F holds on N with a constant C₁(N).
- *Lemma D* (outside the strip). The factors of the products are 1 + η with |η| ≤ exp(−(2πρ/|t|)·0.383), there are at most ρ of them, and ρ ≥ 4.88. For |t| ≤ 2π·1.01 this gives |η| ≤ e^{−0.37ρ} ≤ 0.17 and |E| ≤ Σ |log(1 + η)| ≤ 1.2·ρ·e^{−0.37ρ} ≤ 1.2. So E is bounded (no longer exponentially small, which is not needed here).

With these two bounds and the inversion relation the proof of Lemma R3 applies word for word; the factor e^{2πuy} is e^{O(ρ)} for bounded u. ∎

**Theorem 18.** (a) Theorem 9 of v20 holds for X = 2. (b) Theorem 11 of v21 holds for X = 1 and X = 3. (c) Theorem 14 of v22 holds for X = 2.

Proof. In the three proofs the restriction on X entered only through Lemma R3 (x₀ = 1/X ≤ ¼); Lemma R3′ removes it. The other ingredients have no restriction: Lemma R1 of v20 is stated for X = 2M with any M ≥ 1; Lemmas R2 (v20) and R2′ (v21) are statements about p → ±1 along the path x = 1/X + iε, for any X; the odd-root rule of v3 is proved for every odd m, and its proof covers m = 1 (ζ = 1: (−1;1)ₙ = 2ⁿ, the terms of f₁\* satisfy term_{n+1} = −½·term_n, and θ₄(e^{−t}) → 0). The values: for X = 1, A = A₀ = 1; for X = 2 the sums stop at n = 0 because 1 + q₀ = 0; for X = 3, A = 1 − ω(1 + ω) + (1 + ω)(1 + ω²) = 3 and A₀ = ω − (1 + ω) + 1 = 0. P = 2^{−1/2}·e^{π²/(12t) + t/6}/b at t = −2πi is 2^{−1/2}·e^{−7πi/24}. ∎

Remarks.

- b = 1 is the self-dual point of Φ_b. There the three integrals are taken with Φ₁, and the identity 1 = P·[3·J(½) − I(½) − I(−½)] says that the divergent series ψ₁\*(1) = Σ (−2)ⁿ, summed as 1/(1 + 2) = ⅓ by the radial limit, is given by the formula of v12 at its most symmetric point.
- At X = 3 the radial limit of ψ₀\* vanishes, so there the kernel integral is exactly one third of the sum of the two state integrals.
- Lemma R3′ also removes the remark on the edge case X = 4 in v20.

---

## Appendix. Identities and the odd-root rule

Used in Part 7: identity (II) (and (I)) in the odd-root rule, and identities (V), (VI) in Lemma R2′.

*From `docs/addendum-v2.md`, §2; text unchanged.*

### [v2 §2] New identities for ψ₀\*, ψ₁\*

Notation:
- f₀\*(q) = Σ (−1)ⁿ q^{n²}/(−q;q)ₙ, f₁\*(q) = Σ (−1)ⁿ q^{n²+n}/(−q;q)ₙ (sign-twisted fifth-order f₀, f₁)
- R₀(q) = Σ (−1)ⁿ q^{n²}/(q;q)ₙ, R₁(q) = Σ (−1)ⁿ q^{n²+n}/(q;q)ₙ (Rogers–Ramanujan series Σ zⁿq^{n²}/(q)ₙ at z = −1 and z = −q)
- θ₄(q) = Σ (−1)ⁿ q^{n²} = (q;q)∞/(−q;q)∞

**Identities (proved; verified to q³⁰⁰ in `s1_identities.py`):**

  (I) 2ψ₀\*(q) = f₀\*(q) − θ₄(q)·R₀(q)
  (II) 2ψ₁\*(q) = f₁\*(q) + θ₄(q)·R₁(q)

**Proof sketch.** Extend the sum to all n ∈ ℤ using (−q)₋ₘ = q^{m(m−1)/2}/(2(−q)ₘ₋₁). The negative half of Σ(−1)ⁿq^{n(n+1)/2}(−q)ₙ is −½f₁\*, so the bilateral sum is B = ψ₁\* − ½f₁\*. Write (−q)ₙ = (−q)∞/(−qⁿ⁺¹)∞ and use the partial fraction 1/(−xq)∞ = (1/(q)∞) Σⱼ (−1)ʲ q^{j(j+1)/2} / ((q)ⱼ(1 + xq^{j+1})), valid for all x. With A(w) = Σₙ wⁿ q^{n(n−1)/2}/(1+qⁿ) and Θ(w) = Σₙ wⁿ q^{n(n−1)/2}, the general bilateral sum is

  Σₙ zⁿ q^{n(n+1)/2}(−q)ₙ = ((−q)∞/(q)∞) · { A(z)·z⁻¹ Σⱼ q^{j(j+1)} z⁻ʲ/(q)ⱼ + Θ(z)·(a Hecke-type double sum) }.

At z = −1 two things happen: Θ(−1) = 0, and A(−1) is an Appell sum at a point where its theta function vanishes, so it is modular: A(−1) = −½θ₄². This gives B = ½θ₄R₁, which is (II). The case z = −q gives (I). ∎

**Context.** For the unsigned functions the bilateral series are known to be modular (e.g. Σ_{n∈ℤ} q^{n²}/(−q)ₙ = f₀ + 2ψ₀; see "Bilateral series and Ramanujan's radial limits", arXiv:2202.12141, and Mc Laughlin, arXiv:1906.11997). (I) and (II) are the sign-twisted analogues: the bilateral series is no longer modular, but equals θ₄ times a Rogers–Ramanujan series at z = −1. At z = +1 the same formula has the modular Rogers–Ramanujan H(q) in that slot.

**Double-sum forms (verified to q³⁰⁰, `s11_checks.py`):**

  (−q)∞·ψ₁\* = Σ_{k even, n≥0} (−1)ⁿ q^{k(k+1)/2 + n² + n + kn}/(q)ₖ
  (−q)∞·ψ₀\* = Σ_{k odd, n≥0} (−1)ⁿ q^{k(k+1)/2 + n² + kn}/(q)ₖ

Equivalently (−q)∞·ψ₁\* = Σ_{j≥m≥0} (−1)^{m+j} q^{m²+j²+j}/(q)₂ₘ: a Nahm-type sum weighted by tails of a false theta function. This explains why no eta-quotient multiplier makes ψ\* sparse (v1 §7, fact 4): the Bailey pair with βₙ = (−1)ⁿ has no short αₙ (`s5` check, not included).

*From `docs/addendum-v3.md`, §1; text unchanged.*

### [v3 §1] The odd-root rule (question 3): proved

**Theorem.** Let ζ be a primitive m-th root of unity, m odd. Then

  lim_{r→1⁻} ψ₁\*(rζ) = ⅓·A, A = Σ_{n=0}^{m−1} (−1)ⁿ ζ^{n(n+1)/2} (−ζ;ζ)ₙ,
  lim_{r→1⁻} ψ₀\*(rζ) = ⅓·A₀, A₀ = Σ_{n=0}^{m−1} (−1)ⁿ ζ^{(n+1)(n+2)/2} (−ζ;ζ)ₙ.

**Proof.** It uses (II): 2ψ₁\* = f₁\* + θ₄R₁ (and (I) for ψ₀\*), valid for |q| < 1.

*(a) A finite identity.* For odd m, ∏_{k=0}^{m−1}(1+ζᵏ) = 2, so (−ζ;ζ)_{m−1} = 1 and

  (−ζ;ζ)_{m−1−n} = 1/∏_{j=1}^{n}(1+ζ^{−j}) = ζ^{n(n+1)/2}/(−ζ;ζ)ₙ.

Also (m−1−n)(m−n)/2 ≡ n(n+1)/2 mod m and (−1)^{m−1−n} = (−1)ⁿ. Replacing n by m−1−n in A gives

  A = B := Σ_{n<m} (−1)ⁿ ζ^{n²+n}/(−ζ;ζ)ₙ, and in the same way A₀ = B₀ := Σ_{n<m} (−1)ⁿ ζ^{n²}/(−ζ;ζ)ₙ.

*(b) f₁\* at the root.* At q = ζ the terms of f₁\* = Σ(−1)ⁿq^{n²+n}/(−q;q)ₙ satisfy term_{n+m} = −½·term_n, because (−ζ;ζ)_{n+m} = 2(−ζ;ζ)ₙ. So the series converges there, to B/(1+½) = ⅔B.

*(c) Radial continuity of f₁\*.* For q = rζ group the factors of (−q;q)ₙ in blocks of m. Within a block rᵏ varies by a factor 1+O(1−r), and for equal radii s the block product is exactly ∏ₖ(1+sζᵏ) = 1+s^m. Hence |(−q;q)_{mJ}| ≥ ∏_{j≤J}(1+r^{jm²})·(1−c(1−r))^J. If r^{Jm²} ≥ ½ this is ≥ 1.4^J; otherwise r^{n²} ≤ 2^{−J}. Either way the n-th term is bounded by C·θ^{n/m} with θ < 1, uniformly for r near 1. By dominated convergence f₁\*(rζ) → ⅔B.

*(d) θ₄R₁ → 0.* The same block comparison for (q;q)ₙ, where the block product is 1−s^m, gives |(q;q)_{mJ}| ≥ (Q;Q)_J(1−ct)^{mJ} with Q = e^{−m²t}. Then

  |R₁(ζe^{−t})| ≤ C Σ_J Q^{J²−ηJ}/(Q;Q)_J ≤ C′ (−Q^{1/2};Q)∞ ≤ C′ exp((1+ε)/(m²t)),

using J²−ηJ ≥ J²/2 − η²/2 and Euler's Σ Q^{J²/2}/(Q;Q)_J = (−Q^{1/2};Q)∞. On the other hand θ₄ = (q;q)∞²/(q²;q²)∞ vanishes at odd-order roots: |θ₄(ζe^{−t})| ≤ C t^{−1/2} exp(−π²/(4m²t)) (standard, from the transformation of η). Since 1 < π²/4, the product tends to 0.

By (II), ψ₁\*(rζ) → ½·⅔B = ⅓A. The proof for ψ₀\* is the same with (I), f₀\*, R₀. ∎

So Abel summation with ratio −2 gives the right value because the series f₁\*, which has the same value at the root by the reversal (a), genuinely converges there with ratio −½, and ⅔·½ = ⅓.

**Corollary (proved).** ψ₁\*(ζe^{−t}) − ½f₁\*(ζe^{−t}) = O(exp(−(π²/4−1−ε)/(m²t))). The perturbative series of ψ₁\* at an odd root is therefore that of ½f₁\*, which can be computed term by term (used in §3).

**Checks** (`q3_oddroot.py`):
- A = B and A₀ = B₀ at all primitive roots of odd order ≤ 25: largest difference 7.6·10⁻³⁷.
- Identity (II) at complex q near seven odd roots: residual ≤ 4·10⁻⁵⁵.
- f₁\* → ⅔A and ψ₁\* → ⅓A with error proportional to t.
- **Verified, not proved:** the actual size of θ₄R₁ is e^{−π²/(3m²t)} (m²t·log|θ₄R₁| = −3.13 … −3.18 at m²t = 1/32, tending to −π²/3 = −3.29), i.e. R₁ itself decays like e^{−π²/(12m²t)} at odd roots.

*From `docs/addendum-v4.md`, §1; text unchanged.*

### [v4 §1] The family G(z), the companion W₃, identities (V) and (VI)

Define, for |q| < 1,

  G(z) = Σ_{j≥0} zʲ q^{j²} (q;q²)ⱼ,  W = G(−q),  W₃ = G(−q³).

**q-difference equation (proved, checked to q³¹⁹).** From (q;q²)ⱼ = (q;q²)ⱼ₋₁(1 − q^{2j−1}):

  G(z) = 1 + zq·G(zq²) − zq²·G(zq⁴).

So every G(−q^{2k+1}), k ∈ ℤ, is a combination of 1, W, W₃ with coefficients in ℤ[q, q⁻¹]. Examples: G(−q⁵) = (W − 1 + q²W₃)/q³ and G(−1/q) = 1 − W + qW₃. The W family is therefore the three-dimensional module ⟨1, W, W₃⟩, and **W₃ is the only companion**.

The same holds on the ψ\* side: F(z) = 1 + zq·F(zq) + zq²·F(zq²), so ψ₁\* = 1 − ψ₀\* − q²F(−q²), and every qʲF(−qʲ) lies in ⟨1, ψ₁\*, ψ₀\*⟩ (`f_relation.py`).

**Theorem (proved, checked to q³¹⁹; `w_identities.py`).** Let

  P₊(q) = Σ_{k≥0} q^{2k²+k}/(q²;q²)ₖ,  P₋(q) = Σ_{k≥0} q^{2k²−k}/(q²;q²)ₖ.

Then

  **(V)**  W(q) + Σ_{n≥1} q^{2n²−n}/(q;q²)ₙ = ψ(q)·P₊(q),

  **(VI)**  q·W₃(q) + Σ_{n≥1} q^{2n²−3n+1}/(q;q²)ₙ = ψ(q)·P₋(q).

P± are the two-variable Rogers–Ramanujan function Σ zᵏQ^{k²}/(Q;Q)ₖ at Q = q², z = q^{±1}. They are not functions of q² alone, and they are not modular (see §3: their asymptotic series do not terminate).

*Proof.* Extend the sum for W to negative j using (q;q²)₋ₙ = (−1)ⁿq^{n²}/(q;q²)ₙ. The term j = −n of Σ(−1)ʲq^{j²+j}(q;q²)ⱼ becomes q^{2n²−n}/(q;q²)ₙ, so the left side of (V) is the bilateral sum Σ_{j∈ℤ}(−1)ʲq^{j²+j}(q;q²)ⱼ. Write (q;q²)ⱼ = (q;q²)∞/(q^{2j+1};q²)∞ and expand, with Q = q²,

  1/(yQ;Q)∞ = (1/(Q;Q)∞)·Σ_{k≥0} (−1)ᵏQ^{k(k+1)/2} / ((Q;Q)ₖ(1 − yQ^{k+1})),  y = q^{2j−1}.

This leaves the inner sums Aₖ = Σ_{j∈ℤ}(−1)ʲq^{j²+j}/(1 − q^{2j+2k+1}). Shifting j by k and using x^{−k}/(1−x) = 1/(1−x) + x⁻¹ + … + x^{−k} with x = q^{2j+1}, every extra term is a multiple of Σ_{j∈ℤ}(−1)ʲq^{j²+j} = 0 (the terms j and −1−j cancel). Hence Aₖ = (−1)ᵏq^{k²}·A₀. Finally A₀ = Σ_{j∈ℤ}(−1)ʲq^{j(j+1)}/(1 − q^{2j+1}) = (q²;q²)∞²/(q;q²)∞² = ψ(q)², by the partial-fraction expansion of the reciprocal theta function. Collecting factors gives (V). For (VI) the same steps with q^{j²+3j} = q^{j²+j}·x/q give Bₖ = (−1)ᵏq^{k²−2k−1}A₀. ∎

Compare (I), (II) of addendum-v2: there the prefactor is θ₄ and the sums are sign-twisted Rogers–Ramanujan sums; here the prefactor is ψ and the sums carry a half-integer shift. The one-sided halves Σ_{n≥1} q^{2n²−n}/(q;q²)ₙ and Σ_{n≥1} q^{2n²−3n+1}/(q;q²)ₙ play the role of f₁\*, f₀\*. I did not find them among the classical mock theta functions (short look only).

---

## Corrections and constants after the stress test

From addendum v26, which answers `docs/psi-star-stress-test-v1.md` (its §1 is in Part 4, its §2–§3 in Part 6), and from addendum v27, which answers the second round, `docs/psi-star-stress-test-v2.md` (its other sections concern the z-family, which is not part of this document).

*From `docs/addendum-v26.md`, §4–§5; text unchanged.*

### [v26 §4] Smaller corrections

1. **Lemma F of v16, the numerical remark.** It said that the supremum of |Λ(w) + Li₂(−e^{w})/t|·e^{−Re w/2} is 0.347, "the same for real x and for the neighbouring directions, and independent of |t|". Correct is: 0.347 = ½·ln 2 for real x and for |arg x| ≤ 10°; on the whole domain |arg x| ≤ 60°, |x| ≤ ¼ the supremum measured by the stress test is 0.511 (at |x| = ¼, arg x = ±60°, w = −0.441 ± iπ), and it depends on |t| and arg x. The lemma itself (existence of C₁) is not affected; C₁ = 0.52 works numerically.
2. **Notation in Theorem C (v12).** All fractional powers of q and p are defined through x: q^{α} = e^{2πixα}, p^{α} = e^{−πiα/x}. In D(u) the factor p^{n²/2} is e^{−πi·n²/(2x)}; with a principal power of the number p the odd-n terms can change sign.
3. **Lemma B (iii) of v19.** The factor 2^{N} is what the proof gives; numerically the factor 2 suffices (stress test, S2.8).
4. **(L1) on the contour Γ of v19 §5.** Proved for τ ≤ 0.2 (arg x ≤ 11.3°). Numerically it holds for τ < 0.3085 (arg x < 17.15°) (stress test, S6.3); the proof loses the difference in the chord bound for ln(2cos(η/2)).
5. **Lemma R2 of v20.** The bracket [1/(2K) − e^{−1/K + 2δ}/(2 − δ)] with K = 3/2 is positive for δ < 0.104, not for all δ; the proof uses it for δ → 0 only. (The argument gives e^{−1/K + δ}, for which the range is δ < 0.172.)
6. **Constants without slack.** Lemma A (equality at ρ = 0, θ = −30°), Lemma D (0.383022 at β = −30°) and Proposition 1 (iv) (κ = 0.185115 at θ = 10°) are sharp at an edge of their domains. The statements are correct as written.

---

### [v26 §5] Constants measured by the stress test

The proofs of Parts 4 and 5 use constants "uniform on the sector" without values. The stress test measured them (double precision; details in its table):

| constant | where | measured |
|---|---|---|
| C_U in \|H\| ≤ C_U·e^{−0.185\|s\|} | Lemma 3′, v16 §4 | 0.60, for \|t\| from 0.4 down to 0.001; so t₀ = 0.4 works |
| C_ψ in \|H\| ≤ C_ψ·e^{ψ/\|t\|} | the same | 0.89 |
| true descent rate inf(−ψ/ρ) | Proposition 1 (proves 0.185) | ≥ 0.97 |
| sharp constant in Lemma G / Corollary G | v17 §2 (proves L = 3.72) | 1.49 |
| C₂ of Corollary R | v17 §2 | 48 |
| A in Theorem 4 with B = 240 | v17 §4 | 0.37 |
| C₁ of Lemma F | v16 §2 | 0.52 |
| constants of Lemma B (i), (ii) | v19 §2 | e^{0.40}; 0.090 (ε₁ = π/6), 0.043 (ε₁ = π/2) |
| O(t) constants of Theorem 6 | v19 §3 | ≤ 0.27 |

So the margin "0.175 against 0.149" between the decay of H and the growth of D_N in v17 is an effect of the bounds; the measured rates are about 0.97 against 0.06.

*From `docs/addendum-v27.md`, §4 (the remark in the closing lemma of v12); text unchanged.*

### [v27 §4] The closing lemma of v12 §4.1: the remark on the factor

The sentence in the proof, "(The second value belongs to the point on C; `t5_bounds.py` confirms that the supremum on the line is this bound times a factor between 1 and 1.5.)", is to be read as follows, and is replaced by:

"(The second value belongs to the point on C. Numerically, for the product of the two q-factors and the Gaussian, without the linear factor, the supremum of the modulus on the line is e^{−c·σ_N²} times a factor that stays bounded in N: between 1.0 and 1.5 in the two cases of `t5_bounds.py`. In terms of the exponent E with the linear factor included: sup E/(−c·σ_N²) → 1 as N → ∞, from below (0.94 at N = 2, 0.98 to 1.00 at N = 8, 0.9998 at N = 32), and sup E + c·σ_N² ≤ 0.141·(σ_N + 1) in the samples of the stress test. So c is the exact leading constant, and the difference is the lower-order term e^{O(σ_N)} that the proof allows.)"

The proof is not affected.

---
