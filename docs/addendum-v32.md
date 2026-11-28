# Mock theta search: addendum v32 (the cocycle Ω_γ at x = 0: one formula on the whole domain)

Session of 28 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). Task: the gap in Theorem 22 of v31 that the independent check found (Claude Opus, note at the top of v31): the two functions glued there both exclude x = 0, so the proof gave holomorphy on ℂ \ ((−∞, −½] ∪ {0}) only. This addendum closes the gap. Scripts and output are in `psi_star/v32/`.

Labels: **proved** (a complete argument is written here or in an earlier addendum; nobody else has checked it), **verified** (numerics).

Notation of v30, v31: 𝐒(a, x) the matrix of state integrals, holomorphic for x ∈ ℂ \ (−∞, 0]; 𝐒ᶜ(a, x) = conj 𝐒(−ā, −x̄), holomorphic for x ∈ ℂ \ [0, ∞); γ = (1 0; 2 1), x′ = x/(1 + 2x), a′ = a/(1 + 2x); Ω_γ is defined for Im x > 0 by 𝐉(a′; x′) = Λ(a, x)·𝐉(a; x)·Ω_γ(a, x), Λ = diag(1, λ, −λ).

---

## 0. Summary

1. **The objection is correct.** v31 proved Ω⁽¹⁾ = Ω⁽²⁾ on the lower half-plane and glued them. The union of the two domains is ℂ \ ((−∞, −½] ∪ {0}). The point x = 0 was not covered, and the statement "holomorphic on ℂ \ (−∞, −½]" was not proved there.
2. **Theorem 23 (proved).** For Im x > 0,

     Ω_γ(a, x) = 𝐒(a, x + 1)⁻¹ · 𝐒ᶜ(a′, x′ − 1),

   and the right side is holomorphic in (a, x) for all x ∈ ℂ \ (−∞, −½]. So Ω_γ continues to the whole domain cx + d ∉ (−∞, 0], the point 0 included, by one formula. At x = 0:

     Ω_γ(a, 0) = 𝐒(a, 1)⁻¹ · conj 𝐒(−ā, 1),

   with state integrals at the modulus 1 (b = 1).
3. **How.** γ = T·S·T²·S·T in SL₂(ℤ). Read from the right: T keeps the q-frame, S takes it to the dual frame, T² keeps the dual frame up to a sign, S takes it back, T keeps the q-frame. The first S is applied at x + 1 (Theorem 20 there), the second at x′ − 1 (Proposition 1 of v31 there). For x near 0 these two moduli are near 1 and near −1, far from the cuts of 𝐒 and 𝐒ᶜ.
4. **Consequences (proved).** Theorem 22 holds as stated. The formulas Ω⁽¹⁾ and Ω⁽²⁾ of v31 are restrictions of the same function, so their equality on the lower half-plane (Theorem 22 (i), proved there with Lemma U) now has a second proof that does not use |q| > 1. Lemma U and Proposition U of v31 are unaffected.
5. **A numerical point.** Near x = 0 the old formulas need state integrals at small moduli and lose all accuracy in double precision; the new one does not. This explains why the independent check could not decide the question numerically.

---

## 1. Theorem 23

**Theorem 23 (proved).** Let Im x > 0, a ∈ ℂ. Then Ω_γ(a, x) = 𝐒(a, x + 1)⁻¹·𝐒ᶜ(a′, x′ − 1). The right side is holomorphic in (a, x) for x ∈ ℂ \ (−∞, −½].

Proof. Note x′ − 1 = −(x + 1)/(2x + 1), which lies in the upper half-plane when x does.

- (1) Theorem 21 (T): 𝐉(a; x + 1) = C_T·𝐉(a; x) and 𝐉(a′; x′) = C_T·𝐉(a′; x′ − 1).
- (2) Theorem 20 at (a, x + 1): 𝐉̃(a; x + 1)·𝐒(a, x + 1) = Δ(a, x + 1)·𝐉(a; x + 1).
- (3) Proposition 1 of v31 at (a′, x′ − 1): 𝐉̃(−a′; x′ − 1)·𝐒ᶜ(a′, x′ − 1) = Δ(a′, x′ − 1)·𝐉(a′; x′ − 1).
- (4) The two dual frames: −1/(x′ − 1) = (2x + 1)/(x + 1) = −1/(x + 1) + 2 and −a′/(x′ − 1) = a/(x + 1). So p and ζ̃ are the same for (−a′, x′ − 1) and for (a, x + 1), and p^{1/2} changes sign. As in the proof of Theorem 21, 𝐉̃(−a′; x′ − 1) = E·𝐉̃(a; x + 1), E = diag(1, 1, −1).
- Together: 𝐉(a′; x′) = C_T·Δ(a′, x′−1)⁻¹·E·Δ(a, x+1)·C_T · 𝐉(a; x) · 𝐒(a, x+1)⁻¹·𝐒ᶜ(a′, x′−1).
- The left factor. With d = δ(a, x + 1)/δ(a′, x′ − 1) it is C_T·diag(1, d, −d)·C_T. The lower block of C_T is c = [ α, β ; β, α ] with α = (1 + i)/2, β = (1 − i)/2, and c·diag(d, −d)·c = d·(α² − β²)·diag(1, −1) = i·d·diag(1, −1). So the factor is diag(1, i·d, −i·d).
- i·d = λ(a, x). The parts of d = δ(a, X₁)/δ(a′, X₂), X₁ = x + 1, X₂ = x′ − 1:
  - √X₂/√X₁ = i·(1 + 2x)^{−1/2} (principal values; X₂ = X₁·w with w = −1/(2x + 1) in the upper half-plane, and √w = i/√(2x + 1));
  - X₁ − X₂ = 2xx′ + 2, which gives e^{−5πi·xx′/6}·e^{−5πi/6};
  - 1/X₁ − 1/X₂ = 2, which gives e^{−πi/12};
  - e^{−πi(a − a′)}, and a²/X₁ − a′²/X₂ = 2aa′, which gives e^{−πi·aa′}.

  So d = i·e^{−11πi/12}·(1 + 2x)^{−1/2}·e^{−5πi·xx′/6}·e^{−πi(a − a′)}·e^{−πi·aa′}, and i·d = e^{πi/12}·(the same) = λ.
- Hence 𝐉(a′; x′) = Λ·𝐉(a; x)·[𝐒(a, x+1)⁻¹·𝐒ᶜ(a′, x′−1)]. This is the equation that defines Ω_γ, and Λ·𝐉 is invertible.
- The domain. 𝐒(a, x + 1) is holomorphic and invertible for x + 1 ∉ (−∞, 0], that is x ∉ (−∞, −1]. 𝐒ᶜ(a′, x′ − 1) is holomorphic for x ≠ −½ and x′ − 1 ∉ [0, ∞); now −(x + 1)/(2x + 1) ≥ 0 exactly for x ∈ [−1, −½). Together: x ∉ (−∞, −½]. ∎

**Corollary (proved).**

- Ω_γ(a, ·) continues to a holomorphic function on ℂ \ (−∞, −½], jointly holomorphic in (a, x); Theorem 22 of v31 holds as stated there, and Ω_γ(a, 0) = 𝐒(a, 1)⁻¹·conj 𝐒(−ā, 1).
- On ℂ \ (−∞, 0] this function is Ω⁽¹⁾ = 𝐒(a, x)⁻¹·𝐒(a′, x′), and on ℂ \ ((−∞, −½] ∪ [0, ∞)) it is Ω⁽²⁾ = 𝐒ᶜ(a, x)⁻¹·𝐒ᶜ(a′, x′).
- In particular Ω⁽¹⁾ = Ω⁽²⁾ for Im x < 0, which is Theorem 22 (i) of v31.

Proof. All three formulas are holomorphic on connected domains that contain the upper half-plane, and they agree there (Theorem 21, Proposition 2 of v31, Theorem 23). ∎

Remarks.

- **Two proofs of Theorem 22 (i).** The one in v31 uses the dual shift and Lemma U (for |q| > 1 a Laurent matrix commuting with the q-shift is the identity). The one here uses only that three holomorphic functions agree on the upper half-plane. They are independent, and the numerical checks of v31 §6 and of the independent check support the common statement.
- **Why v31 missed 0.** Ω⁽¹⁾ applies S at x and x′, Ω⁽²⁾ applies S⁻¹ at x and x′. Both are singular at x = 0, the fixed point of γ. The word T·S·T²·S·T moves the two applications of S to x + 1 and x′ − 1.
- **Other elements.** For γⁿ, n ≥ 2, the analogous short word does not exist (S·T^{2n}·S = −γ^{−n}, and T·γ^{−n}·T = −γ only for n = 1). The domain of Ω_g for general g ∈ Γ₀(2) stays open (item 2 of v31 §7).

---

## 2. Numerical checks (verified; `psi_star/v32/d1_zero.py`, `run_log.txt`)

Ω⁽³⁾ denotes 𝐒(a, x+1)⁻¹·𝐒ᶜ(a′, x′−1).

| check | result |
|---|---|
| Theorem 23, Im x > 0: max\|Ω⁽³⁾ − 𝐉⁻¹Λ⁻¹𝐉′\|, quadrature at 30 digits, x = 0.3 + 0.2i and −0.15 + 0.1i | 6·10⁻²⁹, 6·10⁻³⁰ |
| the same in double precision, x = 0.05 + 0.05i, −0.05 + 0.03i, 0.4 + 0.5i, −0.3 + 0.2i | ≤ 8·10⁻¹⁴ |
| i·δ(a, x+1)/δ(a′, x′−1) = λ(a, x), three points | ≤ 8·10⁻³¹ |
| Ω⁽³⁾ against Ω⁽¹⁾ for Im x < 0, six points (Re x from −0.7 to 0.5) | ≤ 4·10⁻¹¹ |
| Ω⁽³⁾ against Ω⁽²⁾ on the interval, x = −0.1, −0.2, −0.3, −0.4 | ≤ 5·10⁻¹³ |
| Ω⁽³⁾ against Ω⁽¹⁾ on the positive real axis, x = 0.2, 0.5, 1 | ≤ 4·10⁻¹³ |
| circle \|x\| = 0.05 (16 points): Ω⁽³⁾ against the series on the upper half; mean value against Ω⁽³⁾(0); negative Fourier modes | 2·10⁻¹⁴; 2·10⁻¹³; ≤ 6·10⁻¹³ |
| circle \|x\| = 0.2 (32 points): the same three; and Ω⁽³⁾ against the old formulas on the lower half | 3·10⁻¹⁴; 2·10⁻¹¹; ≤ 2·10⁻¹⁰; 9·10⁻¹⁰ |
| the old formulas on the lower half of \|x\| = 0.1 and \|x\| = 0.05 | 2·10⁻⁵; useless (error of order 100) |

The value at 0 for a = 0.3: the first row of Ω_γ(0.3, 0) is (1, 0.4612731569 − 0.3395605587i, 0.0197242516 + 0.4838433601i).

The last row of the table is a statement about the numerics of Ω⁽¹⁾ and Ω⁽²⁾, not about the functions: at small |x| the integrals at the moduli x and x′ are exponentially large and cancel in the quotient. Ω⁽³⁾ agrees with the series on the upper half of each circle to 10⁻¹⁴ and satisfies the mean value property.

---

## 3. Open

1. The domain of Ω_g for all g ∈ Γ₀(2) (unchanged).
2. A proof that Ω_γ does not continue across (−∞, −½) (unchanged).
3. The Taylor expansion of Ω_γ at x = 0, which the formula of Theorem 23 makes computable, and its relation to the asymptotic series of the earlier addenda. Not looked at.

---

## 4. Scripts (`psi_star/v32/`)

- `d1_zero.py [a b c d]`: the checks of §2 (12 minutes; part a at 30 digits takes most of it).
- `smat.py`: 𝐒, 𝐒ᶜ in double precision and the three formulas for Ω_γ.
- `frames.py`, `npint.py`, `lowint.py`, `fad.py`, `common.py`, `a1_formula.py`: from v30, unchanged. `run_log.txt`, `README.md`.
