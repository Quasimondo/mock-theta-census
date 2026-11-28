# Mock theta search: addendum v29 (Casoratians of the trivial series with the homogeneous solutions: Appell–Lerch type objects on both sides)

Session of 28 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). Question 1 of `docs/related-ggmw.md` §3 (Claude Opus 5.5): in Garoufalidis–Gu–Mariño–Wheeler the Appell–Lerch type sums are the Casoratians of the trivial and the non-trivial series. What is the analogue here, and does it explain the products H(·;p²)·𝒲 in Theorem C of v12? Scripts and output are in `psi_star/v29/`.

Labels: **proved** (a complete argument is written here or in the addenda cited; nobody else has checked it), **verified** (numerics), **interpretation**.

Notation: L·y(z) = y(z) − zq·y(zq) − zq²·y(zq²); F(z) = Σ_{n≥0} zⁿq^{n(n+1)/2}(−q;q)ₙ (L·F = 1); G(z) = F̌(1/z)/(2z) (L·G = −1); B = F + G = Σ_{n∈ℤ} zⁿq^{n(n+1)/2}(−q;q)ₙ (L·B = 0). I(u), J(u), P, T = e^{πi(2a−1)/x}, 𝒲(T) = 𝒲(T; p) as in v12, v18.

---

## 0. Summary

1. **The Casoratian on the q-side (proved; §1).** D(z) = B(z)·F(zq) − B(zq)·F(z) equals G(z)·F(zq) − G(zq)·F(z), and

   D(z) = −Σ_{k≥0} (−1)ᵏ·zᵏ·q^{k(k+3)/2}·B(zq^{k+1}) = −Σ_{n∈ℤ} (−1)ⁿ·q^{n(n+3)/2}·sₙ·zⁿ,  sₙ = Σ_{j≤n} (−1)ʲ·(−q;q)ⱼ.

   **It is of Appell–Lerch type, one level up**: a classical Appell–Lerch sum is Σ_{k≥0} xᵏ·ϑ(zqᵏ), a sum over the shifts of a theta function, the solution of a first-order equation; D is the same sum over the shifts of B, the solution of the second-order equation. It is not a partial theta function. Its Laurent coefficients are the partial sums of the divergent series Σ_{j∈ℤ} (−1)ʲ(−q;q)ⱼ, whose negative half sums to −½ times Rogers' false theta function Σ q^{n(3n+1)/2}(1 − q^{2n+1}).
2. **Two identities (proved; §1).** At z = −1 and z = 1: f₀\*ψ₁\* − f₁\*ψ₀\* = −2q·Σ_{n∈ℤ} q^{n(n+3)/2}·sₙ and f₁ψ₀ − f₀ψ₁ = −2q·Σ_{n∈ℤ} (−1)ⁿ·q^{n(n+3)/2}·sₙ (the second for Ramanujan's functions).
3. **The inverse matrix (proved; §2).** The first row of the inverse of the 3 × 3 matrix is (1, D_{h₂}/C, −D_{h₁}/C), Casoratians over a theta function, as in GGMW's Theorem 4.
4. **The dual side: yes, it explains Theorem C (proved; §3).** The state-integral part of Theorem Z, K(a) = −I((a−1)/b)·𝒲(T) + pT·I(a/b)·𝒲(p²T), is a Casoratian in the shift a → a + 1 of the dual trivial series 𝒲 with the dual homogeneous solutions I((a+m)/b). The products H(·;p²)·𝒲 are its two cross terms. It satisfies the first-order equation K(a) + K(a+1) = −I(a/b) and is the sum over shifts K(a) = Σ_{j≥1} (−1)ʲ·I((a−j)/b). So 𝒲 is the dual *trivial* series (it solves a second-order inhomogeneous equation), and the Appell–Lerch type object is K.
5. **Kernel integral and Casoratian are a pair (proved; §3).** J(a/b) + J((a+1)/b) = +I(a/b) and K(a) + K(a+1) = −I(a/b): −J and K solve the same first-order equation, and Theorem Z says that their difference from each other is the function e^{πia}·F(z)/(iP), which only changes sign under a → a + 1. At z = 1 this is Zwegers' pair: K becomes fifth-order mock theta functions (Appell–Lerch sums) with pure-power coefficients, J the Mordell integral (v24).
6. **A correction to v28 §3.** B(z) is not "a theta function times a rank-one Nahm sum" for general z. It is a combination of the two Nahm-type solutions with elliptic coefficients; it is a single product at z = −1 (verified: exactly there, not at z = 1 or at generic z).

---

## 1. The Casoratian on the q-side

**Proposition C1.** Let h solve L·h = 0 and D_h(z) = h(z)·F(zq) − h(zq)·F(z). Then D_h(z) + zq²·D_h(zq) = −h(zq).

Proof. zq²·F(zq²) = F(z) − 1 − zq·F(zq) and zq²·h(zq²) = h(z) − zq·h(zq). Insert both in zq²·D_h(zq) = zq²·h(zq)F(zq²) − zq²·h(zq²)F(zq); the terms with zq cancel and −D_h(z) − h(zq) remains. ∎

**Proposition C2.** For h = B, write D = D_B. For z ≠ 0 and |q| < 1:

- (i) D(z) = G(z)·F(zq) − G(zq)·F(z);
- (ii) D(z) = −Σ_{n∈ℤ} (−1)ⁿ·q^{n(n+3)/2}·sₙ·zⁿ with sₙ = Σ_{j≤n} (−1)ʲ·(−q;q)ⱼ, where (−q;q)₋ₘ = q^{m(m−1)/2}/(2(−q;q)ₘ₋₁);
- (iii) D(z) = −Σ_{k≥0} (−1)ᵏ·zᵏ·q^{k(k+3)/2}·B(zq^{k+1});
- (iv) for n < 0: sₙ = −½·Σ_{k≥−n−1} (−1)ᵏ·q^{k(k+1)/2}/(−q;q)ₖ, a tail of σ₋ := Σ_{k≥0} (−1)ᵏq^{k(k+1)/2}/(−q;q)ₖ = Σ_{n≥0} q^{n(3n+1)/2}(1 − q^{2n+1}); in particular s₋₁ = −½σ₋. For n → +∞, sₙ does not converge: s₂ₘ and s₂ₘ₊₁ tend to two limits that differ by (−q;q)∞.

Proof. (i): B = F + G, and the terms F(z)F(zq) cancel. (ii): D is a convergent Laurent series on z ≠ 0 (a product of a series in 1/z and a series in z). Put D = Σ dₙzⁿ and B(zq) = Σ bₙqⁿzⁿ with bₙ = q^{n(n+1)/2}(−q;q)ₙ for all n ∈ ℤ. Proposition C1 gives dₙ = −q^{n+1}·dₙ₋₁ − bₙqⁿ. With dₙ = (−1)ⁿq^{n(n+3)/2}·eₙ this becomes eₙ − eₙ₋₁ = −(−1)ⁿ(−q;q)ₙ, because n(n+1)/2 + n − n(n+3)/2 = 0. So eₙ = c − sₙ with a constant c = lim_{n→−∞} eₙ (the sum defining sₙ converges at −∞, where (−q;q)ⱼ decays like q^{j²/2}). c = 0: for N < 0 the coefficient d_N of (i) is Σ_{n≥0} g_{n−N−1}·fₙ·(qⁿ − q^{N−n}) with |g_m| ≤ C|q|^{m²}, |fₙ| ≤ C·Kⁿ|q|^{n²/2}, so |d_N| ≤ C′·|q|^{(|N|−1)²−|N|}, while |q^{N(N+3)/2}| = |q|^{(N²−3|N|)/2}; the quotient tends to 0. (iii): insert the definition of sₙ in (ii) and put n = j + k, k ≥ 0; the double sum converges absolutely, and Σⱼ zʲq^{j(j+3)/2+jk}(−q;q)ⱼ = B(zq^{k+1}). (iv): the substitution j = −m, k = m − 1; the evaluation of σ₋ is Rogers' false theta identity (checked to q³⁰⁰); the last statement holds because (−q;q)ⱼ → (−q;q)∞. ∎

**Of which type.**

- A level-one Appell–Lerch sum is A(z) = Σₙ (−1)ⁿq^{n(n+1)/2}zⁿ/(1 − xqⁿ) = Σ_{k≥0} xᵏ·ϑ(zqᵏ): a sum over the shifts of a theta function ϑ (which solves a first-order homogeneous equation), with A(z) − x·A(zq) = ϑ(z). By (iii) and Proposition C1, D is the same construction with B, a solution of the second-order equation L·y = 0, in place of ϑ, and with the theta-type weights (−1)ᵏzᵏq^{k(k+3)/2} in place of xᵏ. In the Laurent expansion, the geometric weights 1/(1 − xqⁿ) = Σⱼ xʲq^{nj} are replaced by the partial sums sₙ of Σ (−1)ʲ(−q;q)ⱼ. So D is of Appell–Lerch type in structure, one level up, and it is not a classical Appell–Lerch sum for general z.
- It is not a partial theta function: the coefficients are non-zero in both directions. For n → −∞ they carry the tails of a false theta series (iv); for n → +∞ they carry the partial sums of a divergent series of Kontsevich–Zagier type.
- This is the pattern of GGMW: there the Casoratians satisfy a first-order equation with a non-trivial series on the right (their (19), as quoted in `docs/related-ggmw.md`), and are called Appell–Lerch type sums.

**Corollary (two identities).** Since F(−1) = ψ₁\*, F(−q) = ψ₀\*/q, G(−1) = −f₁\*/2, G(−q) = −f₀\*/(2q), and F(1) = ψ₁, F(q) = ψ₀/q, G(1) = f₁/2, G(q) = f₀/(2q):

  f₀\*·ψ₁\* − f₁\*·ψ₀\* = 2q·D(−1) = −2q·Σ_{n∈ℤ} q^{n(n+3)/2}·sₙ,  f₁·ψ₀ − f₀·ψ₁ = 2q·D(1) = −2q·Σ_{n∈ℤ} (−1)ⁿ·q^{n(n+3)/2}·sₙ.

The second is an identity between Ramanujan's four fifth-order functions f₀, f₁, ψ₀, ψ₁; I have not looked for it in the literature.

---

## 2. The inverse of the 3 × 3 matrix

Let h₁, h₂ be two solutions of L·h = 0 and C(z) = h₁(z)h₂(zq) − h₁(zq)h₂(z) their Casoratian; C(z) = −zq²·C(zq) (Proposition C1 without the inhomogeneous term), so C is a constant multiple of Σ (−1)ⁿq^{n(n+3)/2}zⁿ times an elliptic function. For

  𝕁(z) = [[1, F(z), F(zq)], [0, h₁(z), h₁(zq)], [0, h₂(z), h₂(zq)]]

the first row of 𝕁(z)⁻¹ is (1, D_{h₂}(z)/C(z), −D_{h₁}(z)/C(z)) (block inversion; proved). This is the form of GGMW's Theorem 4: the first row of the inverse consists of the Casoratians of the trivial with the non-trivial series.

---

## 3. The dual side: the products H(·;p²)·𝒲 of Theorem C

Fix a and write, for m ∈ ℤ,

  Iₘ = I((a+m)/b),  wₘ = 𝒲(T·p^{−2m}),  Kₘ = −Iₘ₋₁·wₘ + T·p^{1−2m}·Iₘ·wₘ₋₁,

so that K₀ = −I((a−1)/b)·𝒲(T) + pT·I(a/b)·𝒲(p²T) is the state-integral part of Theorem Z of v18, and Kₘ is the same expression at a + m (T becomes T·p^{−2m}).

**Proposition C3.**

- (i) 𝒲 solves a second-order inhomogeneous equation: 𝒲(T) + Tp²·𝒲(p²T) − Tp³·𝒲(p⁴T) = 1.
- (ii) Iₘ solves the homogeneous three-term recursion Iₘ₋₁ = T·p^{−2m}·(Iₘ + Iₘ₊₁/p) ((B1) of v12).
- (iii) Kₘ + Kₘ₊₁ = −Iₘ.
- (iv) Kₘ = Σ_{j≥1} (−1)ʲ·Iₘ₋ⱼ.
- (v) Jₘ + Jₘ₊₁ = +Iₘ for Jₘ = J((a+m)/b), and Theorem Z reads F(z) = i·P·e^{−πia}·(J₀ + K₀).

Proof. (i): 𝒲(T; p) = 𝔉(−T, p; p²) with 𝔉(z, y; q) = Σ zⁿq^{n(n+1)/2}(y;q)ₙ, and (y;q)ₙ = (y;q)ₙ₋₁(1 − yq^{n−1}) gives 𝔉(z) = 1 + zq·𝔉(zq) − yzq·𝔉(zq²). (ii): v12, (B1), with aₘ = I(u − m/b). (iii): (i) at T·p^{−2m−2} reads wₘ₊₁ + T·p^{−2m}·wₘ − T·p^{1−2m}·wₘ₋₁ = 1. Hence Kₘ + Kₘ₊₁ = wₘ·(−Iₘ₋₁ + T·p^{−2m}·Iₘ + T·p^{−1−2m}·Iₘ₊₁) − Iₘ, and the bracket vanishes by (ii). (iv): iterate (iii) downwards, Kₘ = −Iₘ₋₁ − Kₘ₋₁; the series converges absolutely and Kₘ₋ⱼ → 0 by Theorem A of v12 (this is (B2) of v12 read backwards). (v): the kernels 1/(1 + e^{w}) and e^{w}/(1 + e^{w}) add up to 1; Theorem Z. ∎

**What this says.**

- **The products are cross terms of a Casoratian.** With I(u) = c·e^{−πiu²/2}·[H_e(ζ;q)·H(−ζ̃p;p²) + H_o(ζ;q)·H(ζ̃p;p²)] (Theorem C), the coefficient of each q-side sum H_e(ζ;q), H_o(ζ;q) in K₀ is a dual-side expression of the form

  −(dual homogeneous solution at a − 1)·𝒲(T) + pT·(dual homogeneous solution at a)·𝒲(p²T),

  where the dual homogeneous solutions are the sums H(∓ζ̃p;p²) with their Gaussian factors, and a ↦ a + 1 is the dual shift T ↦ T/p². These are the products H(·;p²)·𝒲 of Theorem C, and they are exactly the two cross terms of a Casoratian of the trivial series 𝒲 with a homogeneous solution. Theorem Z is therefore of the form of GGMW's Theorem 7:

  (integral of the trivial row) = (trivial series) + Σ (non-trivial q-series) × (Casoratian in the dual variable),

  read from the other side: F(z)·e^{πia}/(iP) = J₀ + K₀.
- **𝒲 is not the Appell–Lerch object; K is.** 𝒲 is the dual trivial series: by (i) it solves the dual second-order equation with right-hand side 1, as F does on the q-side. The object that solves a first-order equation with a homogeneous solution on the right, and is a sum over shifts, is K ((iii), (iv)). This answers the question of v18 in the form in which it has an answer: the Appell–Lerch type sum of the dual side is the Casoratian. Whether 𝒲(T) is in addition expressible by classical Appell–Lerch sums as a function of T is not decided here; at T = p^{±1} it is (next point).
- **At z = 1 the pair (K, J) is Zwegers' pair.** For a = 0, K₀ = I(0)·(φ₀(−p) − 1) + I(1/b)·φ₁(−p)/p (v18 §3.1; verified here), with I(0), I(1/b) pure powers (Proposition M of v24): the Casoratian is a combination of fifth-order mock theta functions, which are Appell–Lerch sums by the mock theta conjectures (Hickerson; quoted). And J is the Mordell integral (Theorem 17 of v24). The relations Jₘ + Jₘ₊₁ = Iₘ, Kₘ + Kₘ₊₁ = −Iₘ are then the two halves of the classical statement that the Appell–Lerch sum and the Mordell integral have the same failure of periodicity, with a theta function (here: a pure power) on the right.
- **For general z** the same two first-order equations hold with the non-modular Iₘ on the right. So along the family the kernel integral is a deformation of a Mordell integral (v24), and the dual Casoratian is the matching deformation of an Appell–Lerch sum.

---

## 4. Correction to v28 §3: the bilateral sum for general z

v28 §3 said that B(z) "is a theta function in z times a rank-one Nahm sum". That is right at z = −1 (identities (I), (II)) and wrong in general. By the formula in the proof of (I), (II) in v2, B(z) = ((−q;q)∞/(q;q)∞)·{A(z)·z⁻¹·r(1/z) + Θ₀(z)·(a Hecke-type double sum)}, with r(v) = Σ q^{n²+n}vⁿ/(q;q)ₙ, the Appell–Lerch sum A(z) = Σₙ zⁿq^{n(n−1)/2}/(1 + qⁿ) and Θ₀(z) = Σₙ zⁿq^{n(n−1)/2}. Correct is: B solves L·y = 0, so it is a combination of two Nahm-type solutions with coefficients that are elliptic functions of z.

**Verified:** if B were (multiplier)·r(1/z)/z with a multiplier that changes sign under z → zq, then B(z)·y(zq) + B(zq)·y(z) would vanish for y(z) = r(1/z)/z. Relative to its terms this quantity is 1.4, 0.43, 1.1 at three generic z, 1.4 at z = 1, and 9·10⁻⁴¹ at z = −1. So B is a single product at z = −1 and not at z = 1 (where B(1) = ψ₁ + ½f₁ is a sum of two theta products) or elsewhere.

Nothing else in v28 depends on this. The table in v28 §3 should read, in its last row: "a combination of Nahm-type solutions with elliptic coefficients".

---

## 5. Numerical checks (verified; `y1_casoratian.py`, `run_log.txt`)

| check | result |
|---|---|
| D by its two definitions; D(z) + zq²·D(zq) = −B(zq); the Laurent series (ii); the sum over shifts (iii): five pairs (z, q), relative | ≤ 2·10⁻⁴⁰ |
| the two identities of the Corollary at q = 0.3 + 0.2i, −0.4 + 0.3i, 0.6 | ≤ 2·10⁻⁴⁰ |
| σ₋ = Σ q^{n(3n+1)/2}(1 − q^{2n+1}) | exact to q³⁰⁰ |
| Proposition C3 (i), (iv), (v) and K₀ at z = 1: three x, four a | ≤ 3·10⁻²⁸ |
| Proposition C3 (ii), (iii), relative (cancellation at the smallest \|p\|) | 2·10⁻²¹, 5·10⁻²⁶ |
| §4: B(z)y(zq) + B(zq)y(z), relative | 1.4, 0.43, 1.1, 1.4 (z = 1), 9·10⁻⁴¹ (z = −1) |

---

## 6. Open

1. Whether D(z) (or K) is expressible by classical Appell–Lerch sums for z other than ±1, and whether 𝒲(T) is, as a function of T.
2. The matrix 𝕁 and its cocycle (question 4 of `docs/related-ggmw.md`): with §2 and §3 both sides of the matrix are now explicit.
3. The identity f₁ψ₀ − f₀ψ₁ = −2q·Σ (−1)ⁿq^{n(n+3)/2}sₙ: is it known, and what is its modular meaning (the left side is a mixed mock modular form of weight 1)?

---

## 7. Scripts (`psi_star/v29/`)

- `y1_casoratian.py [a b c] [N]`: the checks of §5 (3 minutes).
- `common.py`, `a1_formula.py`, `c2_z1.py` (from v18), `run_log.txt`, `README.md`.
