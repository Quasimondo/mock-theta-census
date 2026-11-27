# Mock theta search: addendum v17 (factorial constants; P·J is the Borel sum of Φ)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). The item left open in `psi-star-addendum-v16.md` §5: remainder bounds with factorial constants, and from them Borel summability. Scripts and output are in `psi_star/v17/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v17/check_v17.py`): the constants c₃ = 0.5816, L = 3.7195, ν = 0.5444, m_c = 1.7077, λ = 0.7904, A₀ = 3.22, B₀ = 15.91, B₁ = 238.6 and 1/(ν·c₁) = 45.9 recomputed; the infimum of |1 − e^{η}| on |η| ≥ 1, |Im η| ≤ π is 0.632121 = 1 − 1/e on a grid, above the proved 0.5816. Off the real axis (|t| = 0.2, arg t = −60°, own quadrature of P·J): |P·J − Σ_{n<N}(−1)ⁿτₙ| = 0.34, 0.16, 0.063, 0.022, 0.0070, 0.0019, 4.6·10⁻⁴, 9.3·10⁻⁵, 1.2·10⁻⁵ for N = 0…8, faster than 2^{−N}, then levelling at 1.0·10⁻⁵, the size |t|^{−1/2}·|e^{−S/t}| of the non-perturbative term (Corollary R only claims N ≤ 0.04/|t|). ãₖ = aₖ: `psi_star/v10/phi_psi1.py` computes Φ exactly as the termwise expansion of ½f₁\* = Σ(−1)ⁿτₙ and agrees with v13. The proofs were read; the estimates were not checked line by line, and the Nevanlinna–Sokal theorem is quoted, not checked against the original.

Labels: **proved** (a complete argument is written here; nobody else has checked it), **verified** (numerics), **quoted** (a theorem from the literature), **inference**.

---

## 0. Summary

1. **Theorem 4 (proved, §1–§4).** There are constants A and t₀ > 0 such that for all t with −210° ≤ arg t ≤ 10°, |t| ≤ t₀, and all K ≥ 0

   | P·J − Σ_{k<K} aₖ·tᵏ | ≤ A·K!·(240·|t|)ᴷ,

   for the ψ₁\* and the ψ₀\* integral. This includes the real x-axis (|q| = 1).
2. **Theorem 5 (proved, given the Nevanlinna–Sokal theorem; §5).** For every direction θ₀ with −120° ≤ θ₀ ≤ −80° the series Φ is Borel summable, and P·J(t) = 𝒮_{θ₀}[Φ](t) for Re(e^{iθ₀}/t) > 1/R. The Borel transform is analytic in a neighbourhood of the closed sector −120° ≤ arg ξ ≤ −80° and exponentially bounded there.
3. **Consequence.** The four exact formulas of v12 (kernel integral + bilinear terms) are now formulas "Borel sum of the perturbative series + bilinear terms", as theorems, for |q| < 1 with −90° < arg t < 10° and small |t|. This is the three-sector formula of v10, with one caveat about the name of the direction (§5.2).
4. **How.** Not through the window. Move the contour down across the first N lower kernel poles. The residues are τₙ = q^{n²+n}/(2(−q;q)ₙ), and the functional equation of Φ_b turns the moved integral back into the original integrand times an explicit product D_N, which is at most about 2^{−N}. With N = ⌊0.04/|t|⌋ the moved integral is exponentially small, and the finite residue sum is a sum of functions analytic in t near 0, for which Cauchy's estimate gives factorial constants directly.
5. **What the proof uses from v14–v16:** only the uniform bound |H| ≤ C·e^{−κ|s|} on the contour (v15 Lemma 3′, v16 §4) and the properties of Faddeev's integral (v16 Lemmas Φ, F, D). The Euler–Maclaurin and window expansions are not needed for this route.
6. **Numerically (verified):** the shift identity holds to 10⁻¹³ at real x and on both sides; the remainder after N residues falls like 2^{−N}; the ratios (|remainder|/(K!·|t|ᴷ))^{1/K} stay below 0.29 ≈ 1/|S|; P·J equals the Borel sums in the directions −120°, −100°, −80° (10⁻³¹ off the real axis, 4·10⁻¹⁴ on it).
7. **Not done:** Borel directions outside [−120°, −80°] (numerically the range is −176.5° … −17.1°); explicit values of A, t₀, R.

---

## 1. The shift identity

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

## 2. The size of D_N

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

## 3. The residues as functions of t

**Lemma T.** Let c = 0.3, m_c = 1 + e^{−c}·cos c = 1.7077, λ = e^{c}/m_c = 0.7904. Then τₙ is analytic for |t| < π/n and

  |τₙ(t)| ≤ (e^{c}/2)·λⁿ for |t| ≤ c/n,  hence |τ_{n,k}| ≤ (e^{c}/2)·λⁿ·(n/c)ᵏ for the Taylor coefficients τ_{n,k} = [tᵏ] τₙ.

The same holds for the ψ₀\* version.

Proof. For 1 ≤ j ≤ n and |t| ≤ c/n the number u = j·t has |u| ≤ c < π/2, so |1 + e^{−u}| ≥ Re(1 + e^{−u}) = 1 + e^{−Re u}·cos(Im u) ≥ m_c. And |q^{n²+n}| ≤ e^{|t|(n²+n)} ≤ e^{c(n+1)} (for ψ₀\*: e^{cn}). Cauchy's estimate on the circle |t| = c/n gives the coefficients. ∎

Consequently ãₖ := Σ_{n≥0} (−1)ⁿ·τ_{n,k} converges absolutely, and with Σ_{n≥0} nᵏ·λⁿ ≤ k!/(1 − λ)^{k+1} (because nᵏ ≤ (n+1)…(n+k)):

  |ãₖ| ≤ A₀·k!·B₀ᵏ,  A₀ = e^{c}/(2(1 − λ)) = 3.22,  B₀ = 1/(c·(1 − λ)) = 15.91.

---

## 4. Theorem 4

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

## 5. Theorem 5: Borel summability

### 5.1 Statement and proof

**Theorem (Nevanlinna; Sokal 1980; quoted).** Let f be analytic in the disc C_R = {z : Re(1/z) > 1/R} and suppose f(z) = Σ_{k<K} cₖzᵏ + R_K(z) with |R_K(z)| ≤ A·σᴷ·K!·|z|ᴷ for all K ≥ 0 and all z ∈ C_R. Then B(ζ) = Σ cₖζᵏ/k! converges for |ζ| < 1/σ, continues analytically to the strip {ζ : dist(ζ, ℝ₊) < 1/σ}, is bounded there by a constant times e^{|ζ|/R} (on every narrower strip), and f(z) = (1/z)·∫₀^∞ e^{−ζ/z}·B(ζ) dζ for z ∈ C_R.

Reference: A. D. Sokal, *An improvement of Watson's theorem on Borel summability*, J. Math. Phys. 21 (1980) 261–263, https://doi.org/10.1063/1.524408. I took the statement from quotations of it in later papers (for example arXiv:2511.04509 §4.2, arXiv:1011.6490 §4), not from the original. Watson's original version (a sector of opening greater than 180°) would also do here.

**Theorem 5.** Let B(ξ) = Σ aₖξᵏ/k! be the Borel transform of Φ (ψ₁\* or ψ₀\* version), and R = t₀. Then:

- (i) B is analytic in the set {ξ : dist(ξ, Σ₀) < 1/240}, Σ₀ = {−120° ≤ arg ξ ≤ −80°} ∪ {0}, and |B(ξ)| ≤ C·e^{|ξ|/R} on {dist(ξ, Σ₀) ≤ 1/241}.
- (ii) For every θ₀ ∈ [−120°, −80°] and every t with Re(e^{iθ₀}/t) > 1/R:

  P·J(t) = 𝒮_{θ₀}[Φ](t) := (1/t)·∫₀^{∞·e^{iθ₀}} e^{−ξ/t}·B(ξ) dξ.

Proof. Fix θ₀ and put f(z) = P·J(e^{iθ₀}·z). For z ∈ C_R: |z| < R and |arg z| < 90°, so t = e^{iθ₀}z lies in the open sector −210° < arg t < 10°, |t| < t₀, where P·J is analytic (Corollary of v16) and Theorem 4 holds. So the hypotheses hold with cₖ = aₖ·e^{ikθ₀} and σ = 240. The Borel transform of f is B(e^{iθ₀}ζ); substituting ξ = e^{iθ₀}ζ gives (ii), and the union of the strips over θ₀ gives (i). ∎

The union of the discs in (ii) contains every t with −210° < arg t < 10° and |t| < R·cos(dist(arg t, [−120°, −80°])).

### 5.2 What this says about the earlier formulas

- **|q| < 1 (proved).** For Im x > 0, −90° < arg t < 10° and t in a disc of (ii), the Corollary of v12 becomes

  ψ₁\*(q) = 𝒮_{θ₀}[Φ₁](t) + (κ/2)·[𝓑₂·p·W₃(p) − 𝓑₄·W(p)],  f₁\*(q) = 2·𝒮_{θ₀}[Φ₁](t) − κ·[𝓑₂·F_b(p) − 𝓑₄·F_a(p)],

  ψ₀\*(q) = 𝒮_{θ₀}[Φ₀](t) − i·q^{1/4}·(κ/2)·[𝓑₁·p·W₃(p) + 𝓑₃·W(p)],  f₀\*(q) = 2·𝒮_{θ₀}[Φ₀](t) + i·q^{1/4}·κ·[𝓑₁·F_b(p) + 𝓑₃·F_a(p)],

  with any θ₀ ∈ [−120°, −80°] within 90° of arg t, and Φ₁, Φ₀ the perturbative series of ψ₁\*, ψ₀\* at q → 1. This is the three-sector formula of v10 with 𝒯 = 𝒮_{θ₀}[Φ₁].
- **The caveat.** v10 defined 𝒯 as the Borel sum *along the direction of t* (for −90° < arg t < −17.1°). For −90° < arg t ≤ −80° that is literally what is proved. For −80° < arg t < −17.1° the two definitions agree if B has no singularity in −80° < arg ξ ≤ arg t and is exponentially bounded there. Numerically that is the case (the nearest singularities are S and its conjugate, at arg ∓17.1°; v13). It is not proved: Theorem 5 only reaches the directions −120° … −80°.
- **|q| = 1 and |q| > 1 (proved for P·J).** P·J is the Borel sum of Φ also for real x and for Im x < 0 (−210° < arg t ≤ −90°). What is missing on that side is the other half: Theorem A of v12 for |q| > 1.

---

## 6. Numerical checks (verified; `g1_shift.py`, `g2_checks.py`)

| check | result |
|---|---|
| Proposition 2: both sides for N = 1, 2, 5; ψ₁\* at (\|t\|, θ) = (0.3, −90°), (0.1, −90°), (0.2, −60°), (0.2, −150°); ψ₀\* at (0.2, −90°), (0.2, −120°) | differences 2·10⁻¹⁴ … 5·10⁻¹³ (double precision) |
| Corollary R at real x: \|P·J − Σ_{n<N} (−1)ⁿτₙ\| | \|t\| = 0.05: N = 1: 1.7·10⁻¹; 4: 2.1·10⁻²; 8: 1.4·10⁻³; 16: 9.5·10⁻⁶; 24: 1.4·10⁻⁷; 32: 8.3·10⁻⁹; 40: 4.0·10⁻⁹; 48: 1.3·10⁻⁸. Minimum near N = 2π/(3\|t\|) = 41.9, of size e^{−Cl₂(π/3)/\|t\|} = 1.5·10⁻⁹ up to a factor 3. Same picture at \|t\| = 0.1 (minimum 7.6·10⁻⁵ at N = 20) and 0.2 (8.1·10⁻³ at N = 8) |
| Corollary G: sup over the contour of log\|D_N\| + ν·N − L·c₁·\|s\|, θ = −210°, −150°, −90°, −30°, 10°; \|t\| = 0.02, 0.01, 0.004, 0.001 | negative in all 20 cases (−0.25 … −7.1) |
| Lemma G: sup of \|h′\| on the strips around the contours, α = −70° … 40° | 2.00 … 2.53 (L = 3.72) |
| infimum of \|1 − e^{η}\| on \|η\| ≥ 1, \|Im η\| ≤ π | 0.632121 = 1 − 1/e (proved bound: 0.5816) |
| Lemma T: max over \|t\| = c/n of \|τₙ\|/((e^{c}/2)·λⁿ), n = 1 … 80, both versions | ≤ 0.73, decreasing in n |
| ãₖ = aₖ: Σ_{n≤400} (−1)ⁿ[tᵏ]τₙ in exact rational arithmetic against `phi_exact.json`, k < 13, both versions | 1.3·10⁻⁶⁹ (the size of the omitted tail) |
| \|aₖ\|/(A₀·k!·B₀ᵏ), k < 13 | ≤ 0.10, decreasing fast |
| Gevrey ratios (\|P·J − Σ_{k<K} aₖtᵏ\|/(K!·\|t\|ᴷ))^{1/K}, K up to 16 (double precision; θ = −90°, −30°, −150°) and up to 24 (30 digits; θ = −200°, 5°), \|t\| = 0.2 and 0.1 | at most 0.290 in all 16 runs; e.g. θ = −90°, \|t\| = 0.2: 0.111, 0.124, 0.121, 0.103, 0.177, 0.209, 0.222, 0.224, 0.215, 0.224 for K = 1, 2, 3, 4, 6, 8, 10, 12, 14, 16. For comparison (\|aₖ\|/k!)^{1/k} = 0.271, 0.278, 0.280 at k = 50, 100, 150 and 1/\|S\| = 0.2905 |
| Theorem 5 (ii): P·J against 𝒮_{θ₀}[Φ], θ₀ = −120°, −100°, −80° (Padé [60/60]), \|t\| = 0.3 | real x: 4·10⁻¹⁴ in all three directions, ψ₁\* and ψ₀\*. θ = −150°: 2·10⁻³¹ (three directions). θ = −30°: 3·10⁻³¹ (two admissible directions). θ = −200°: 1·10⁻²⁴ (θ₀ = −120°). θ = 5°: 1·10⁻¹¹ (θ₀ = −80°; 85° from the direction, the Laplace integral is truncated) |

**An observation at real x (inference).** There |τₙ| = 1/(2·Π_{j≤n} 2|cos(j|t|/2)|), which decreases as long as j < 2π/(3|t|) and then grows. The least term is about exp(−(1/|t|)·∫₀^{2π/3} ln(2cos(u/2)) du) = e^{−Cl₂(π/3)/|t|} = |e^{−S/t}|, with S the Borel singularity of v13. So the residue series is itself an asymptotic series in N, and its optimal truncation error is the size of the first non-perturbative term.

---

## 7. What is open now

1. **More Borel directions.** Numerically 𝒮_{θ₀}[Φ] = P·J for −176.5° < θ₀ < −17.1° (v13). Proved: −120° ≤ θ₀ ≤ −80°. The proof of Theorem 4 needs, for another contour angle δ, only the uniform bound |H| ≤ C·e^{−κ|s|} with κ > L·c₁ and a lower bound for cos α; the scan of v14 (straight contours descend for −257.6° < arg t < 65.1°) suggests that varying δ reaches about −167° < θ₀ < −25°. The last 8° to 9° on each side need more than straight contours, or knowledge of the singularities of B.
2. **The Stokes jump** across arg t = −17.1° and the nature of the singularity at S (v10 §3).
3. **Asymptotics of 𝓑ⱼ** (v7 §3 is a derivation). With it and Theorem 5 the laws of v5 for the one-sided sums become theorems.
4. **Theorem A of v12 for |q| > 1.**
5. Identity (VII); W, W₃ and 𝒥 in the bounded direction.
6. **An independent reader.** v12 and v14–v17 are proofs written and checked by one author. The chain is: v12 (kernel identity), v15 Proposition 1 and Lemma 3′, v16 Lemmas Φ, F, D and §4 (uniform bound), v17 (this file). v14, the rest of v15 and Lemma S of v16 are not needed for Theorems 4 and 5 if one takes aₖ = Σ (−1)ⁿ[tᵏ]τₙ as the definition of Φ.

---

## 8. Scripts (`psi_star/v17/`)

- `g1_shift.py`: the shift identity, the remainder after N residues at real x, the constants (numpy; 20 seconds).
- `g2_checks.py [a b c d e]`: a: bound for D_N; b: bound for τₙ; c: aₖ from the residues in exact arithmetic; d: Gevrey ratios; e: Borel sums in three directions (a–c: seconds; d, e: about 5 minutes, mpmath for the points off the range of the fast integrator).
- `fad.py` (from v16), `borel.py`, `phi_exact.json` (from v13).
- `run_all.sh`, `run_log.txt`, `README.md`.
