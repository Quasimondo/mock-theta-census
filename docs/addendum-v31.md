# Mock theta search: addendum v31 (the cocycle Ω_γ continues across (−½, 0): proof)

Session of 28 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). Task: open item 1 of v30. There Ω_γ(a, x) = 𝐒(a, x)⁻¹·𝐒(a′, x′) was proved holomorphic on ℂ \ (−∞, 0], and seen numerically to continue across (−½, 0). This addendum proves the continuation. Scripts and output are in `psi_star/v31/`.

Labels: **proved** (a complete argument is written here or in an earlier addendum; nobody else has checked it), **verified** (numerics), **from memory**.

Notation of v30: 𝐒(a, x) is the matrix of the six state integrals, holomorphic in (a, x) for x ∈ ℂ \ (−∞, 0] (Proposition E); 𝐉, 𝐉̃, Δ are the q-frame, the dual frame and the Jacobi factor; γ = (1 0; 2 1), x′ = x/(1 + 2x), a′ = a/(1 + 2x); p = e^{−πi/x}, T = e^{πi(2a−1)/x}, z = e^{2πia}, q = e^{2πix}.

---

## 0. Summary

1. **Theorem 22 (proved).** For every a ∈ ℂ the function Ω_γ(a, ·) continues holomorphically to ℂ \ (−∞, −½]. On the interval −½ < x < 0 it is

     Ω_γ(a, x) = conj[ 𝐒(−ā, −x)⁻¹ · 𝐒(−ā/(1 + 2x), −x/(1 + 2x)) ],

   with state integrals at the two positive moduli −x and −x/(1 + 2x). So the domain is the one of Garoufalidis–Zagier's definition, cx + d ∉ (−∞, 0] (from memory), and the q-frame 𝐉 is a matrix-valued holomorphic quantum Jacobi form on Γ₀(2) in that sense, for the generators T and γ.
2. **How.** There is a second matrix of state integrals, the reflection 𝐒ᶜ(a, x) = conj 𝐒(−ā, −x̄), holomorphic for x ∉ [0, ∞). It gives a second formula Ω_γ = 𝐒ᶜ(a, x)⁻¹·𝐒ᶜ(a′, x′) in the upper half-plane, and this one continues across (−½, 0) into the lower half-plane. The point is that the two continuations agree below the real axis. That follows from a uniqueness statement which holds only for |q| > 1: a matrix 𝐕(z) of Laurent series that commutes with the q-shift of the inhomogeneous equation is the identity (Lemma U). It rests on the fact, found in v30 §6, that the equation L·y = c has no Laurent series solution when |q| > 1.
3. **By-product (proved).** For Im x < 0 the kernel integral is characterised by its two difference equations: 𝒥(·; x) is the only entire function y with L·y = 1 and y(a) − y(a + 1) = 𝒾₀(a) (Proposition U). For Im x > 0 this is false, because the bilateral sum B(z) can be added.
4. **Not proved.** That Ω_γ does not continue across (−∞, −½) is numerical only (v30 §5).

---

## 1. The reflected matrix

**Definition.** 𝐒ᶜ(a, x) := conj 𝐒(−ā, −x̄). It is holomorphic in (a, x) for −x̄ ∈ ℂ \ (−∞, 0], that is for x ∈ ℂ \ [0, ∞).

**Proposition 1 (proved).** For Im x > 0 and every a:

  𝐉̃(−a; x)·𝐒ᶜ(a, x) = Δ(a, x)·𝐉(a; x).

Proof. Take the complex conjugate of Theorem 20 at the point (−ā, −x̄), which also lies in the upper half-plane.

- The entries of 𝐉 are series in q^{1/4} and ζ^{±1} with real coefficients, and conj q^{1/4}(−x̄) = q^{1/4}(x), conj ζ(−ā) = ζ(a). So conj 𝐉(−ā; −x̄) = 𝐉(a; x).
- The entries of 𝐉̃ are series in p^{1/2} and ζ̃^{±1} with real coefficients, and conj p^{1/2}(−x̄) = p^{1/2}(x), conj ζ̃(−ā; −x̄) = 1/ζ̃(a; x). So conj 𝐉̃(−ā; −x̄) = 𝐉̃(−a; x).
- conj δ(−ā, −x̄) = e^{−πi/4}·q^{−5/24}p^{1/24}e^{−πia − πia²/(2x)}/(2·√(−x)), and √(−x) = −i·√x for Im x > 0. So conj δ(−ā, −x̄) = δ(a, x). ∎

**Proposition 2 (proved).** Put Ω⁽²⁾(a, x) := 𝐒ᶜ(a, x)⁻¹·𝐒ᶜ(a′, x′).

- For Im x > 0, Ω⁽²⁾ = Ω_γ.
- Ω⁽²⁾ is holomorphic in (a, x) for x ∈ ℂ \ ((−∞, −½] ∪ [0, ∞)).

Proof. As in the proof of Theorem 21: p, ζ̃ are the same at (a′, x′) and (a, x) and p^{1/2} changes sign, so 𝐉̃(−a′; x′) = E·𝐉̃(−a; x) with E = diag(1, 1, −1). Proposition 1 at both points gives 𝐉(a′; x′) = Δ(a′, x′)⁻¹·E·Δ(a, x)·𝐉(a; x)·Ω⁽²⁾ = Λ·𝐉(a; x)·Ω⁽²⁾, the equation that defines Ω_γ, and 𝐉 is invertible. For the domain: x′ ∈ [0, ∞) exactly for x ∈ [0, ∞) ∪ (−∞, −½), and x′ = ∞ at x = −½. ∎

So Ω_γ, a function on the upper half-plane, has two continuations: Ω⁽¹⁾ = 𝐒(a, x)⁻¹·𝐒(a′, x′) on ℂ \ (−∞, 0] (through the positive real axis, Theorem 21) and Ω⁽²⁾ (through (−½, 0)). Both are defined on the lower half-plane. They glue to one function on ℂ \ (−∞, −½] if and only if they agree there.

---

## 2. The two shifts

**The q-shift (v30).** 𝐒(a + x, x) = 𝐒(a, x)·𝐀(z), 𝐀 = [ 1, 0, −1/(zq²) ; 0, 0, 1/(zq²) ; 0, 1, −1/q ], on the whole cut plane. The entries of 𝐀 are rational in z, q with real coefficients, and conj z(−ā) = z(a), conj q(−x̄) = q(x). Hence also 𝐒ᶜ(a + x, x) = 𝐒ᶜ(a, x)·𝐀(z).

**The dual shift (proved).** 𝐒(a + 1, x) = 𝐂(a, x)·𝐒(a, x) on the cut plane, with

  𝐂 = [ 1, −1, 0 ; 0, p, −p/T ; 0, −1, 0 ].

Proof. From the definitions: 𝒥(a + 1) = 𝒥(a) − 𝒾₀(a) (the kernel identity J(u) + J(u + 1/b) = I(u)), 𝒾₁(a + 1) = −𝒾₀(a), and 𝒾₀(a + 1) = p·𝒾₀(a) − (p/T)·𝒾₁(a) (the recursion I((a−1)/b) = T·[I(a/b) + I((a+1)/b)/p] of v29 §3). T does not change under a → a + x, so the same relations hold in both columns. They hold for Im x > 0 and continue in x. ∎

Consequently 𝐒ᶜ(a + 1, x) = 𝐂ᶜ(a, x)·𝐒ᶜ(a, x) with 𝐂ᶜ(a, x) = conj 𝐂(−ā − 1, −x̄)⁻¹. Its entries are rational in conj p(−x̄) = p and conj T(−ā − 1, −x̄) = p⁴/T. So **𝐂 and 𝐂ᶜ depend on (a, x) only through p and T**, and p, T are the same at (a′, x′) and at (a, x) (because 1/x′ = 1/x + 2 and a′/x′ = a/x).

---

## 3. Uniqueness for |q| > 1

**Lemma U0 (proved; v30 §6 for c = 1).** Let |q| > 1 and c ∈ ℂ. If y = Σ_{n∈ℤ} cₙzⁿ converges on some annulus and y(z) − zq·y(zq) − zq²·y(zq²) = c, then c = 0 and y = 0.

Proof. The equation is cₙ = qⁿ(1 + qⁿ)·cₙ₋₁ + c·δ_{n,0}. If cₙ ≠ 0 for some n ≥ 0, then all later coefficients are non-zero and grow like |q|^{n²}, so the series diverges for every z ≠ 0. Otherwise the equation at n = 0 gives c₋₁ = −c/2, and c₋ₘ₋₁ = c₋ₘ·qᵐ/(1 + q⁻ᵐ) for m ≥ 1. If c ≠ 0 these grow like |q|^{m²/2}, and the series diverges for every z. If c = 0 all coefficients vanish. ∎

**Lemma U (proved).** Let |q| > 1. Let 𝐕(z) = [ 1, v ; 0, V ] (v a row of length 2, V a 2 × 2 matrix) have entries holomorphic on ℂ \ {0}, and let 𝐀(z)·𝐕(zq) = 𝐕(z)·𝐀(z). Then 𝐕 = 1.

Proof. Write 𝐀 = [ 1, α ; 0, A ] with α = (0, −1/(zq²)) and A = [ 0, 1/(zq²) ; 1, −1/q ]. The equation splits into

  A(z)·V(zq) = V(z)·A(z),  v(zq) + α·V(zq) = α + v(z)·A(z).

- *The block.* Let V = [ g, u ; r, s ]. The four entries of A·V(zq) = V·A are
  - u(z) = r(zq)/(zq²), so r(z) = zq·u(z/q);
  - s(z) = g(zq) − r(zq)/q = g(zq) − zq·u(z);
  - g(z)/(zq²) − u(z)/q = s(zq)/(zq²), which becomes g(z) − g(zq²) = zq·u(z) − zq²·u(zq);
  - r(z)/(zq²) − s(z)/q = u(zq) − s(zq)/q, which becomes u(z/q) − q·u(zq) = g(zq) − g(zq²) − zq·u(z) + zq²·u(zq).

  With g = Σ gₙzⁿ, u = Σ uₙzⁿ the last two read gₙ(1 − q^{2n}) = q(1 − qⁿ)·uₙ₋₁ and uₙ(q⁻ⁿ − q^{n+1}) = gₙqⁿ(1 − qⁿ) − q(1 − qⁿ)·uₙ₋₁.
  - n = 0: the second gives u₀(1 − q) = 0, so u₀ = 0.
  - n ≠ 0: the first gives gₙ(1 + qⁿ) = q·uₙ₋₁; inserted in the second, uₙ = −q^{n+1}(1 − qⁿ)·uₙ₋₁/((1 + qⁿ)(1 − q^{2n+1})).
  - So uₙ = 0 for n ≥ 0. For n = −m ≤ −1 the recursion runs downwards, u₋ₘ₋₁ = −u₋ₘ·(1 + q⁻ᵐ)(1 − q^{1−2m})·q^{m−1}/(1 − q⁻ᵐ), with no vanishing factor; if u₋₁ ≠ 0 the coefficients grow like |q|^{m²/2} and u diverges for every z. So u = 0.
  - Then gₙ = 0 for n ≠ 0, r = 0, s = g₀: V = g₀·1 with a constant g₀.
- *The row.* With V = g₀ the second equation is v(zq) = v(z)·A(z) + (1 − g₀)·α(z), that is v₁(zq) = v₂(z) and v₂(zq) = v₁(z)/(zq²) − v₂(z)/q − (1 − g₀)/(zq²). Eliminating v₂: v₁(z) − zq·v₁(zq) − zq²·v₁(zq²) = 1 − g₀. By Lemma U0, g₀ = 1 and v₁ = 0; then v₂ = 0. ∎

---

## 4. Theorem 22

**Theorem 22 (proved).**

- (i) For Im x < 0 and every a: 𝐑(a′, x′) = 𝐑(a, x), where 𝐑(a, x) := 𝐒ᶜ(a, x)·𝐒(a, x)⁻¹.
- (ii) Hence Ω⁽¹⁾ = Ω⁽²⁾ on the lower half-plane, and Ω_γ(a, ·) continues to a holomorphic function on ℂ \ (−∞, −½], jointly holomorphic in (a, x). On −½ < x < 0 it is given by the formula of the summary.

Proof. (i) Fix x with Im x < 0; then Im x′ < 0 too. All matrices below have the form [ 1, ∗ ; 0, ∗ ], are entire in a, and are invertible (the determinant of 𝐒 is elementary and has no zero).

- *𝐑 has period x in a.* 𝐑(a + x, x) = 𝐒ᶜ(a)·𝐀·𝐀⁻¹·𝐒(a)⁻¹ = 𝐑(a, x). The same at x′: 𝐑(α + x′, x′) = 𝐑(α, x′).
- Put 𝐕(a) := 𝐒ᶜ(a, x)⁻¹·𝐑(a′, x′)·𝐒(a, x).
- *𝐕(a + 1) = 𝐕(a).* (a + 1)′ = a′ + 1 − 2x′. By the period x′ and the dual shift at x′, 𝐑((a+1)′, x′) = 𝐑(a′ + 1, x′) = 𝐂ᶜ(a′, x′)·𝐑(a′, x′)·𝐂(a′, x′)⁻¹. By §2, 𝐂(a′, x′) = 𝐂(a, x) and 𝐂ᶜ(a′, x′) = 𝐂ᶜ(a, x). With 𝐒(a + 1) = 𝐂·𝐒(a) and 𝐒ᶜ(a + 1) = 𝐂ᶜ·𝐒ᶜ(a) the four factors cancel.
- *𝐕(a + x) = 𝐀(z)⁻¹·𝐕(a)·𝐀(z).* (a + x)′ = a′ + x′, so 𝐑((a + x)′, x′) = 𝐑(a′, x′), and 𝐒, 𝐒ᶜ are multiplied by 𝐀(z) on the right.
- So 𝐕 is an entire function of a with period 1, that is a matrix of Laurent series in z converging on ℂ \ {0}, and 𝐀(z)·𝐕(zq) = 𝐕(z)·𝐀(z). Here |q| > 1. By Lemma U, 𝐕 = 1, which is 𝐑(a′, x′) = 𝐒ᶜ(a, x)·𝐒(a, x)⁻¹.

(ii) 𝐒ᶜ(a, x)⁻¹·𝐒ᶜ(a′, x′) = 𝐒(a, x)⁻¹·𝐑(a, x)⁻¹·𝐑(a′, x′)·𝐒(a′, x′) = 𝐒(a, x)⁻¹·𝐒(a′, x′) for Im x < 0. Ω⁽¹⁾ is holomorphic on ℂ \ (−∞, 0], Ω⁽²⁾ on ℂ \ ((−∞, −½] ∪ [0, ∞)) (Proposition 2), they agree on both half-planes, and the union of the two domains is ℂ \ (−∞, −½]. For real −½ < x < 0 the points −x̄ = −x and −x̄′ = −x/(1 + 2x) are positive, and Ω⁽²⁾ is the stated formula. ∎

Remarks.

- **Where |q| > 1 is used.** Only in Lemma U. For Im x > 0 the equality 𝐑(a′, x′) = 𝐑(a, x) also holds, for a different reason: there 𝐑 = 𝐉̃(−a; x)⁻¹·𝐉̃(a; x) is a function of ζ̃ and p^{1/2}, and the sign of p^{1/2} cancels. Lemma U itself fails for |q| < 1: the bilateral sum B(z) of v28 is a Laurent solution of L·y = 0, and for example 𝐕 = [ 1, −F(z), −F(zq) ; 0, 2, 0 ; 0, 0, 2 ] satisfies 𝐀(z)·𝐕(zq) = 𝐕(z)·𝐀(z).
- **The inverse generator.** Ω_{γ⁻¹}(a, x) = Ω_γ(γ⁻¹(a, x))⁻¹ is holomorphic for x/(1 − 2x) ∉ (−∞, −½], that is on ℂ \ [½, ∞), again the domain cx + d ∉ (−∞, 0]. The two are mirror images: Ω_γ(a, x) = conj Ω_{γ⁻¹}(−ā, −x̄) on ℂ \ (−∞, −½] (this is Ω⁽¹⁾ = Ω⁽²⁾ written with the definition of 𝐒ᶜ).
- **Other elements of Γ₀(2).** Ω_g is a product of the matrices Ω_γ^{±1} at points of the orbit (Theorem 21), so it is holomorphic where all those points lie in the domains above. That this is always the full domain cx + d ∉ (−∞, 0] was not checked.
- **The cut (−∞, −½].** Numerically Ω_γ jumps there (v30 §5). A proof is not attempted.

---

## 5. A characterisation of the kernel integral below the real axis

**Proposition U (proved).** Let Im x < 0. If y is an entire function of a with

  y(a) − zq·y(a + x) − zq²·y(a + 2x) = 1  and  y(a) − y(a + 1) = 𝒾₀(a),

then y = 𝒥(·; x).

Proof. 𝒥 satisfies both (v30, Proposition Q (i)). The difference d = y − 𝒥 is entire, has period 1, and satisfies L·d = 0. As a Laurent series in z it vanishes by Lemma U0. ∎

This answers, in one sense, item 4 of v30 §8 ("a closed form of the trivial row for Im x < 0 other than the integral"): there is no series, but the two difference equations determine the function. In the upper half-plane they do not: 𝒥 + κ·B(z) with the bilateral sum B and any constant κ satisfies both.

---

## 6. Numerical checks (verified; `psi_star/v31/c1_continuation.py`, `run_log.txt`)

| check | result |
|---|---|
| Proposition 1: max\|𝐉̃(−a; x)·𝐒ᶜ(a, x) − Δ·𝐉(a; x)\|, 𝐒ᶜ by quadrature at 30 digits, three points (two with Re x < 0) | 4·10⁻²⁹, 3·10⁻²⁷, 9·10⁻²⁴ |
| the dual shift 𝐒(a + 1) = 𝐂·𝐒(a): x = 0.3 + 0.2i, x = 0.6 (real), x = 0.5 − 0.15i, x = −0.2 − 0.3i | 3·10⁻¹⁵, 2·10⁻¹⁴, 5·10⁻¹², 6·10⁻¹⁰ |
| Theorem 22 (i): max\|𝐑(a′, x′) − 𝐑(a, x)\| at six points with Im x < 0, first row (size 1 to 17) | ≤ 2·10⁻¹¹ |
| the same, 2 × 2 block, relative to its size (4 to 4·10⁶) | ≤ 3·10⁻¹² |
| Ω⁽²⁾ against the series 𝐉⁻¹Λ⁻¹𝐉′ above the axis, x = −0.25 + 0.05i, −0.15 + 0.02i, −0.4 + 0.1i | ≤ 2·10⁻¹³ |
| Ω⁽²⁾ against Ω⁽¹⁾ below the axis, five points, Im x = −0.03 … −0.4 | ≤ 9·10⁻¹² |
| the same on the lower half of the circle \|x + 0.25\| = 0.1 (8 points, down to Im x = −0.02) | ≤ 4·10⁻⁹ |
| on the interval, x₀ = −0.1, −0.2, −0.3: Ω_γ(x₀) from the formula of Theorem 22; distance of the series at x₀ + iε, ε = 0.02, 0.01, 0.005 | halves with ε at each point (0.104, 0.054, 0.027 at x₀ = −0.1) |
| the same from below, Ω⁽¹⁾(x₀ − iε), ε = 0.04, 0.02 | 0.25, 0.12 at x₀ = −0.1; 0.23, 0.11 at −0.2; 0.26, 0.15 at −0.3 |
| Lemma U: the reduction of A·V(zq) = V·A to the two scalar equations, at random Laurent polynomials g, u; the example 𝐕 ≠ 1 for \|q\| < 1 | 4·10⁻¹⁵; 3·10⁻¹⁵ |

Example: Ω_γ(0.3, −0.2) has first row (1, 0.4789227 − 0.2993521i, −0.3277870 + 0.2317086i).

A numerical remark. The integrals at the reflected Jacobi variable −ā have large integrands (|e^{2πuy}| with |u| up to 4), which costs digits in double precision. The script first reduces Re a to [0, 1) with the dual shift and uses a steeper line (δ₀ up to 0.75).

---

## 7. Open

1. A proof that Ω_γ does not continue across (−∞, −½).
2. The domain of Ω_g for all g ∈ Γ₀(2).
3. Unchanged from v30: the formula for general second argument y of 𝔉(z, y; q) (SL₂(ℤ) instead of Γ₀(2)); whether the entries of 𝐒⁻¹ are single state integrals; the comparison of the relations (A0), (A1) with the printed text.

---

## 8. Scripts (`psi_star/v31/`)

- `c1_continuation.py [a b c d e]`: the checks of §6 (9 minutes; part a at 30 digits takes 6 of them).
- `frames.py`, `npint.py`, `lowint.py`, `fad.py`, `common.py`, `a1_formula.py`: from v30, unchanged. `run_log.txt`, `README.md`.
