# Mock theta search: addendum v21 (ψ₁\* and ψ₀\* at odd roots of unity)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). Item 1 of `psi-star-addendum-v20.md` §5 for the roots e^{2πi/X} with X odd. Scripts and output are in `psi_star/v21/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v21/check_thm12.py`, own code): Theorem 12 against A/3 computed from one period of the series at X = 201 … 3201 (up to 430 digits): Richardson extrapolation of X·(ratio − 1) gives 0.2015332627 − 1.570796327i for ψ₁\* (= −2πi·s̄₁ to all 10 digits) and 0.2015332627 + 0.5235987757i for ψ₀\* (= −2πi·(−1/12 + i√3/54)), the same values as at even X. Lemma R2′ on the path x = 1/5 + iε, with W and W₃ summed directly at 150 digits: |W − ⅓| = 0.0089, 0.0026, 0.00088 and |p·W₃ + ⅓| = 0.019, 0.0053, 0.0018 for 1 − |p| = 0.076, 0.023, 0.0078, falling like 1 − |p|. Theorem 11 (the integrals at real b) was not recomputed; Lemma R2′ was read.

Labels: **proved** (a complete argument is written here or in the addenda cited; nobody else has checked it), **verified** (numerics).

---

## 0. Summary

1. **Theorem 11 (proved; §1–§2).** Let X ≥ 5 be odd, q₀ = e^{2πi/X}, x = 1/X, b = X^{−1/2}, t = −2πi/X, and A = Σ_{n<X} (−1)ⁿ·q₀^{n(n+1)/2}·(−q₀;q₀)ₙ (one period of the series). Then

   lim_{r→1⁻} ψ₁\*(r·q₀) = A/3 = P·[ J(1/(2b)) − ⅓·I(1/(2b)) − ⅓·I(−1/(2b)) ],

   and the same for ψ₀\* with A₀ (exponent (n+1)(n+2)/2), the integrals at u = b ± 1/(2b) and the factor e^{πi/X}. The first equality is the odd-root rule of addendum v3 §1.
2. **One statement for all X.** With v20: the radial value of ψ₁\* at e^{2πi/X} is P·[J + I(1/(2b))·w₃ − I(−1/(2b))·w], where (w, w₃) are the radial values of (W, p·W₃) at the dual point p₀ = (−1)^X: (1, 1) for even X and (⅓, −⅓) for odd X.
3. **Theorem 12 (proved; §3).** As X → ∞ through odd integers, the radial limit is

   ⅓·3^{−1/4}·(X/2)^{1/2}·e^{Cl₂(π/3)·X/2π}·e^{2πi/3 − iπX/6}·(1 + O(1/X)),

   and the same times e^{2πi/3} for ψ₀\*. So the forward law on the unit circle, ψ₁\*(e^{2πi/X}) ≈ C̄·t^{−1/2}·e^{−S̄/t}·W(e^{−πiX}), is a theorem for every X, with W(1) = 1 and W(−1) = ⅓.
4. **The new lemma (R2′).** As p → −1 along the path that x = 1/X + iε induces: W(p) → ⅓ and p·W₃(p) → −⅓. Proof through the identities (V), (VI) of v4: F_a, F_b converge at p = −1 as geometric series with ratio −½ and are continuous there, and ψ(p)·P_±(p) → 0 because the theta function ψ vanishes like e^{−π²/(8δ)} while P_± grow at most like e^{π²/(30δ)}. This also proves the radial limit of W at ξ = −1 that v4 §2 had verified.
5. **Numerically (verified):** A = P·[3J − I(+) − I(−)] to 10⁻⁹ … 10⁻¹³ for X = 5 … 61, both functions; in Theorem 12, X·(ratio − 1) tends to −2πi·(1/4 + i√3/54) = 0.2015 − 1.5708i for ψ₁\*, the first correction s̄₁ of the forward law.
6. **Not covered:** other roots of unity e^{2πih/k} with h ≠ 1; X = 3; z ≠ −1.

---

## 1. The dual series at p = −1

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

## 2. Theorem 11

**Theorem 11.** For odd X ≥ 5 the identity of the summary holds, and likewise

  lim_{r→1⁻} ψ₀\*(r·q₀) = A₀/3 = e^{πi/X}·P·[ J(b + 1/(2b)) − ⅓·I(b + 1/(2b)) − ⅓·I(b − 1/(2b)) ].

Proof. For x = 1/X + iε the formulas of v12 (in the form of v20 §2) hold. As ε → 0⁺: q = q₀·e^{−2πε} is radial, so the left sides tend to A/3, A₀/3 (v3 §1, Theorem; q₀ is a primitive X-th root of unity, X odd). On the right, the integrals are continuous at x = 1/X (v20, Lemma R3; 1/X ≤ ¼), and W(p) → ⅓, p·W₃(p) → −⅓ (Lemma R2′). ∎

Equivalently A = P·[3·J(1/(2b)) − I(1/(2b)) − I(−1/(2b))]: at odd order the full period sum A is three times the kernel integral minus the two state integrals, while at even order the half-period sum is the kernel integral plus one state integral minus the other (v20).

---

## 3. Theorem 12

**Theorem 12.** As X → ∞ through odd integers, the asymptotic formulas of the summary hold for the radial limits A/3 and A₀/3.

Proof. As for Theorem 10 of v20: −⅓·P·I(−1/(2b)) = ⅓·A′·e^{5πi/12}·t^{−1/2}·e^{−S̄/t}·(1 + O(t)) with A′ = (π/√3)^{1/2} (v19, Theorem 6 at real x), the other two terms are O(1) and O(|t|^{−1/2}e^{−Cl₂/|t|}), and at t = −2πi/X the leading term is ⅓·3^{−1/4}·(X/2)^{1/2}·e^{2πi/3}·e^{Cl₂X/2π − iπX/6}. For ψ₀\* the growing integral carries the extra factor e^{2πi/3}. ∎

---

## 4. Numerical checks (verified; `k1_odd.py`, `run_log.txt`)

| check | result |
|---|---|
| Theorem 11, ψ₁\*: A against P·[3J − I(+) − I(−)] at x = 1/X, X = 5, 7, 9, 15, 25, 41, 61 | relative difference 7·10⁻¹⁰, 4·10⁻¹⁰, 2·10⁻¹⁰, 6·10⁻¹², 8·10⁻¹³, 4·10⁻¹², 5·10⁻¹² |
| Theorem 11, ψ₀\*, the same X | 7·10⁻¹¹, 4·10⁻⁹, 2·10⁻¹³, 4·10⁻¹³, 8·10⁻¹³, 3·10⁻¹², 5·10⁻¹² |
| Lemma R2′ at X = 5, ε = 10⁻², 10⁻³, 3·10⁻⁴, 10⁻⁴ (δ = 0.54, 0.076, 0.023, 0.0078) | \|W − ⅓\| = 0.38, 0.0089, 0.0026, 0.00088; \|p·W₃ + ⅓\| = 0.096, 0.019, 0.0053, 0.0018; \|ψ·P₊\| = 0.40, 6·10⁻⁹, 3·10⁻³⁰, 2·10⁻⁹⁰; φ/δ² → 0.064 = 1/(πX) |
| identities (V), (VI) at these points | to working precision (10⁻⁴¹ … 10⁻²⁹⁶) |
| odd-root rule of v3 on the same path: \|ψ₁\*(q) − A/3\| | 0.93, 0.022, 0.0064, 0.0021 |
| Theorem 12: X·(exact/asymptotic − 1), X = 101 … 3201 | ψ₁\*: 0.1970 − 1.5874i … 0.2014 − 1.5713i (limit −2πi·s̄₁ = 0.2015 − 1.5708i with s̄₁ = 1/4 + i√3/54); ψ₀\*: 0.2081 + 0.5381i … 0.2017 + 0.5240i |

The integrals at real x are computed as in v20 (double precision, Faddeev's integral, descent contours for the state integrals). The series W, W₃ near p = −1 have terms of size e^{0.3/δ} with alternating signs; part b uses up to 300 digits.

---

## 5. Open

1. Roots e^{2πih/k} with h ≠ 1: x → h/k is not small, so the asymptotics of v15–v19 do not apply directly. The exact formula of v12 still holds for Im x > 0.
2. X = 3 (x = ⅓ is outside the range |x| ≤ ¼ in which the bounds of v16 were proved; the formula presumably holds, not checked).
3. The second-order term −2πi·s̄₁/X as a theorem (it is the first correction of 𝓑₄, which the method of v19 gives but which is not written out).
4. The remaining items of v19 §8 and v20 §5.

---

## 6. Scripts (`psi_star/v21/`)

- `k1_odd.py [a b c]`: a: Theorem 11; b: Lemma R2′ and the identities on the path; c: Theorem 12 (two minutes in all).
- `h1_roots.py` (from v20: the integrators), `fad.py` (from v16), `run_log.txt`, `README.md`.
