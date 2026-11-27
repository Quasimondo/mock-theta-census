# Mock theta search: addendum v19 (asymptotics of the state integrals; the laws of v5 for ψ₁\*, ψ₀\*, f₁\*, f₀\* are theorems)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). The main gap named after v18: the growth and decay of 𝓑₁ … 𝓑₄ were a saddle-point derivation (v7 §3). Scripts and output are in `psi_star/v19/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v19/check_thm68.py`, own code without `common.py`): Theorem 8 tested directly against the q-series. ψ₁\* and ψ₀\* summed from their definitions at arg x = 1°, 5°, 10° and |x| = 1/50 … 1/400 (|ψ₁\*| up to 10²⁷), divided by the growing term A·t^{−1/2}·e^{5πi/12}·e^{−S̄/t}·W(p) (for ψ₀\*: −A·e^{iπ/12}·…): (ratio − 1)/t converges to 0.2507 + 0.0320i, i.e. the first correction s̄₁ = 1/4 + i√3/54 of v5, in all three directions, and for ψ₀\* to about −0.084 + 0.032i. Theorem 6 for the growing integral I(−1/(2b)) by direct quadrature with an own Φ_b at arg x = 10°, |x| = 1/6, 1/12, 1/24: (I/asymptotic − 1)/t stays at 0.07 … 0.10 + 0.04 … 0.09i, consistent with O(t) and with the 0.055·|t| of §7 after the factor c·p^{1/8}. The Γ₂ inequality of §5 holds with constant −0.082 (claim −0.015) on a grid of 41 values of τ ∈ [0, 0.2]. The contour estimates themselves were read, not rechecked line by line.

Labels: **proved** (a complete argument is written here or in the addenda cited; nobody else has checked it), **verified** (numerics), **inference**.

---

## 0. Summary

1. **Theorem 6 (proved; §3–§5).** As x → 0 with 0 ≤ arg x ≤ 11°, uniformly in the direction and including the real axis and tangential approach (horocycles):

   𝓑₄ = (4/3)^{1/4}·e^{−iπ/12}·e^{−i·Cl₂/t}·(1 + O(t)),  𝓑₃ = (4/3)^{1/4}·e^{+iπ/12}·e^{−i·Cl₂/t}·(1 + O(t))  (growing),

   𝓑₂ = (4/3)^{1/4}·e^{+iπ/12}·e^{+i·Cl₂/t}·(1 + O(t)),  𝓑₁ = (4/3)^{1/4}·e^{−iπ/12}·e^{+i·Cl₂/t}·(1 + O(t))  (decaying),

   with Cl₂ = Cl₂(π/3) and t = −2πix. These are the constants of v7 §3. For the decaying pair the statement holds for 0 ≤ arg x ≤ 60°.
2. **Theorem 8 (proved, by the chain v12, v15–v17, v19; §6).** For 0 < arg x ≤ 10° and |x| small, with S = π²/3 − i·Cl₂ and p the dual nome:

   ψ₁\*(q) = 𝒮[Φ₁](t) + (π/√3)^{1/2}·t^{−1/2}·[ e^{5πi/12}·e^{−S̄/t}·W(p)·(1 + O(t)) + e^{−5πi/12}·e^{−S/t}·p·W₃(p)·(1 + O(t)) ],

   and the same for f₁\* with (2𝒮, −2F_a, −2F_b), and for ψ₀\*, f₀\*. This is the forward law of v5 with its constants, now with a Borel sum, exact dual functions, and error terms. It holds however x approaches 0 inside the sector.
3. **Theorem 7 (proved; §4).** For every 0 < a < 1 the decaying state integral I(a/b) of the z-family has the saddle-point asymptotics at the root W in the lower half-plane (v18 §2.3), for 0 ≤ arg x ≤ ε(a).
4. **How.** Explicit contours on which the size of the integrand is given by closed formulas. Decaying case: the horizontal line through the saddle; the derivative of the height along it is a difference of two arguments of points on a ray, hence monotone. Growing case (z = −1): a ray from iπ, the imaginary axis from iπ to 0, and the positive real axis; on these the height is a quadratic, a Clausen function, and π·r. The bounds for Φ_b are those of v16.
5. **Not done:** the growing integral I((a−1)/b) for a ≠ ½ (the volume rule for general z stays a verified derivation); radial limits to roots of unity; arg x > 11° for the growing pair.

---

## 1. Setting

x → 0 with 0 ≤ arg x ≤ 60°, |x| ≤ ¼; t = −2πix, so arg t = −90° + ε with ε = arg x, τ = tan ε. Φ_b is Faddeev's integral (v16 §1); for real x > 0 the state integrals are defined with it. As in v18, with w = 2π·b·y and u = a′/b:

  I(a′/b) = (1/(2πb))·∫ e^{Q(w)/t}/Φ_b dw,  Q(w) = −w²/2 − 2πi·a′·w,  V(w) = Q(w) + Li₂(−e^{w}).

**Height.** |e^{V/t}| = exp(cos ε·H_τ(w)/|t|), H_τ = −Im V + τ·Re V. (Re(V/t) = Re(i·e^{−iε}·V)/|t|.)

V′(w) = −w − 2πi·a′ − log(1 + e^{w}), V″(w) = −(1 + 2e^{w})/(1 + e^{w}). Li₂(−e^{w}) is principal in the strip |Im w| ≤ π; for Re w < 0 it is analytic and 2πi-periodic.

---

## 2. Bounds for Φ_b (from v16)

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

## 3. Statement for the integrals

**Theorem 6.** Uniformly for |x| → 0:

- (decaying, 0 ≤ arg x ≤ 60°) I(1/(2b)) = (t/2π)^{1/2}·3^{−1/4}·e^{+iπ/12}·b^{−1}·exp((−5π²/12 + i·Cl₂)/t)·(1 + O(t)), and I(b + 1/(2b)) = e^{−2πi/3}·(the same)·(1 + O(t));
- (growing, 0 ≤ arg x ≤ 11°) I(−1/(2b)) = (t/2π)^{1/2}·3^{−1/4}·e^{−iπ/12}·b^{−1}·exp((−5π²/12 − i·Cl₂)/t)·(1 + O(t)), and I(b − 1/(2b)) = e^{+2πi/3}·(the same)·(1 + O(t)),

with principal roots. For Im x > 0 this is item 1 of the summary, by Theorem C of v12: 𝓑₂,₄ = I(±1/(2b))/(c·p^{1/8}), 𝓑₁,₃ = ±i·q^{1/4}·I(b ± 1/(2b))/(c·p^{1/8}), c·p^{1/8} = 2^{−1/2}·e^{−iπ/4}·e^{t/24}·e^{−5π²/(12t)}, and t^{1/2}/b = (2π)^{1/2}·e^{−iπ/4}.

**Theorem 7.** Let 0 < a < 1, W the root of W² + W = e^{−2πia} with Im W < 0, w₁ = Log W. There is ε(a) > 0 such that for 0 ≤ arg x ≤ ε(a), uniformly,

  I(a/b) = (1/(2πb))·(2πt·(1 + W)/(1 + 2W))^{1/2}·e^{V(w₁)/t}·(1 + O(t)),  −Im V(w₁) = −Vol(z) < 0,

with the principal root. For a = ½ one can take ε = 60°.

---

## 4. The decaying integrals (proof of Theorem 7 and of the first half of Theorem 6)

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

## 5. The growing integrals at z = −1 (proof of the second half of Theorem 6)

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

## 6. Consequences

### 6.1 Theorem 8: the forward laws (proved)

By v12, ψ₁\* = P·J(1/(2b)) + P·I(1/(2b))·p·W₃(p) − P·I(−1/(2b))·W(p) (since (κ/2)·𝓑 = P·I), with P = 2^{−1/2}·e^{π²/(12t) + t/6}/b. With Theorem 6 and 1/x = −2πi/t:

  P·I(1/(2b)) = +A·e^{−5πi/12}·t^{−1/2}·e^{−S/t}·(1 + O(t)),  P·I(−1/(2b)) = −A·e^{+5πi/12}·t^{−1/2}·e^{−S̄/t}·(1 + O(t)),

with A = (π/√3)^{1/2}, S = π²/3 − i·Cl₂ (checked numerically, §7). For 0 < arg x ≤ 10° the direction of t lies in (−90°, −80°], where P·J = 𝒮[Φ] by Theorem 5 of v17. Hence, for |t| < R:

- ψ₁\* = 𝒮[Φ₁] + A·t^{−1/2}·[ e^{5πi/12}·e^{−S̄/t}·W(p)·(1 + O(t)) + e^{−5πi/12}·e^{−S/t}·p·W₃(p)·(1 + O(t)) ],
- f₁\* = 2𝒮[Φ₁] − 2A·t^{−1/2}·[ e^{5πi/12}·e^{−S̄/t}·F_a(p)·(1 + O(t)) + e^{−5πi/12}·e^{−S/t}·F_b(p)·(1 + O(t)) ],
- ψ₀\* = 𝒮[Φ₀] − A·t^{−1/2}·[ e^{iπ/12}·e^{−S̄/t}·W(p)·(1 + O(t)) + e^{−iπ/12}·e^{−S/t}·p·W₃(p)·(1 + O(t)) ],
- f₀\* = 2𝒮[Φ₀] + 2A·t^{−1/2}·[ e^{iπ/12}·e^{−S̄/t}·F_a(p)·(1 + O(t)) + e^{−iπ/12}·e^{−S/t}·F_b(p)·(1 + O(t)) ],

The first exponential grows, the second decays; each O(t) is a factor of an exact coefficient and does not depend on p. The constants S̄ = π²/3 + i·Cl₂ and A·e^{5πi/12} are those of v5. The exponent S of the decaying term is the Borel singularity of Φ found in v13.

(For ψ₀\*: q^{1/2}·P·I(b ± 1/(2b)) = e^{∓2πi/3}·P·I(±1/(2b))·(1 + O(t)), and e^{∓2πi/3}·e^{∓5πi/12} = −e^{∓iπ/12}.)

### 6.2 What remains exact and what is asymptotic

The formulas of v12 are exact. Theorem 8 replaces the three coefficient functions by their leading behaviour. On a horocycle (p fixed) it gives the two-term law of v5 with relative error O(t) in each sector; radially (p → 0) it gives ψ₁\* = 𝒮 + A·e^{5πi/12}·t^{−1/2}·e^{−S̄/t}·(1 + O(t)).

### 6.3 The z-family

- Decaying sector: Theorem 7, all 0 < a < 1.
- Growing sector for a ≠ ½: not proved. The vertical line through the saddle is again the local descent path (d/dη Im V = −r₀ − ln|1 + e^{r₀+iη}| is monotone in η ∈ (0, π)), and the real axis works on the right (Im V = 2π|a′|·r). The left part needs a case analysis in a, because the left valley starts at Im w > 2π|a′|, which is above the strip for |a′| > ½, and the saddle is not on the imaginary axis.
- A statement at roots of unity (the volume rule as measured) needs in addition the radial limit of the formula of v18 and a bound for J(a/b) at real x for a ≠ ½.

---

## 7. Numerical checks (verified; `e1_checks.py`, `run_log.txt`)

| check | result |
|---|---|
| V(w₀) at a′ = −½ | −5π²/12 − i·Cl₂(π/3) to 12 digits |
| growing contour, τ = 0, 0.1, 0.2: sup of (H − H(w₀))/min(d², d) on the three pieces (d = distance to w₀) | Γ₁: −0.96, −0.89, −0.81; Γ₂: −0.31, −0.24, −0.16; Γ₃: −0.48, −0.33, −0.17 |
| closed formulas on the pieces: Im V = πr on Γ₃; Im V = −Cl₂(π − η), Re V = ¾η² − πη − π²/12 on Γ₂; Im V ≥ 0, Re V ≤ −π²/3 on Γ₁ | 10⁻¹⁵; minimum of Im V on Γ₁ is 10⁻⁵ (at ρ → 0), maximum of Re V is −π²/3 |
| ln(2cos(η/2))/(π − 3η/2) on [0, 2π/3) | minimum 0.220636 = ln 2/π at η = 0; limit 1/√3 at 2π/3 |
| decaying line, a = ½: the same sup for τ = 0, 0.2, 1, 1.7 | −0.38, −0.50, −0.99, −1.42 |
| decaying line, a = ¼, ⅓, 5/6, 0.03, 0.97 at τ = 0 (and τ = 0.05 or 0.01) | negative in all cases (−0.13 … −0.014 at τ = 0); critical point equation holds to 10⁻¹⁶ |
| 𝓑₄, 𝓑₃, 𝓑₂, 𝓑₁ computed exactly (Theorem C of v12) divided by the asymptotic value, minus 1; arg x = 0.5°, 5°, 11°; \|x\| = 1/100, 1/200, 1/400 | growing pair: 3.5·10⁻³, 1.7·10⁻³, 8.4·10⁻⁴; decaying pair: 3.2·10⁻³, 1.6·10⁻³, 8.2·10⁻⁴; the same in the three directions. That is 0.055·\|t\| and 0.051·\|t\|, while \|𝓑₄\| runs from 10⁷ to 10²⁸ and \|𝓑₂\| from 10⁻⁷ to 10⁻²⁸ |
| the coefficients of Theorem 8, P·I(±1/(2b)) and q^{1/2}·P·I(b ± 1/(2b)), exact against asymptotic; arg x = 1°, 10°; \|x\| = 1/100, 1/400 | ratio − 1 = 0.016 and 0.004 (ψ₁\*), 0.021 and 0.005 (ψ₀\* relative to ψ₁\*): signs and phases as stated, error O(t) |

---

## 8. Open

1. The growing integral I((a−1)/b) for a ≠ ½ (§6.3); then the volume rule for general z.
2. Radial limits of the exact formulas at roots of unity, and from them ψ₁\*(e^{2πi/X}) ≈ 3^{−1/4}·(X/2)^{1/2}·e^{Cl₂·X/2π} as a theorem.
3. Higher terms of the expansions of 𝓑ⱼ (the method gives them; only the leading term is written), and arg x > 11° for the growing pair (rotate Γ₂).
4. The items of v17 §7 and v18 §6 not touched here: more Borel directions, the Stokes jump, the trivial sector for a ≠ ½, |q| > 1, the Mordell integral at z = 1, the 3×3 matrix.
5. An independent reader. The chain for Theorem 8 is: v12 §2–§4; v15 Proposition 1 and Lemma 3′; v16 §1, Lemmas S, F, D, §4; v17; v19.

---

## 9. Scripts (`psi_star/v19/`)

- `e1_checks.py [a b c d e]`: a: heights on the growing contour; b: heights on the decaying lines; c: the concavity inequality; d: 𝓑ⱼ against their asymptotic values; e: the coefficients of Theorem 8 (110 digits; 40 seconds in all).
- `common.py` (from v12), `run_log.txt`, `README.md`.
