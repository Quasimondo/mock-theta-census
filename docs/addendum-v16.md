# Mock theta search: addendum v16 (the wedge around the real x-axis: P·J ~ Φ on a full sector, including |q| = 1)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). Items 1 and 2 of `psi-star-addendum-v15.md` §4: the wedge around the real x-axis, and the gluing of the two sides. Scripts and output are in `psi_star/v16/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v16/check_faddeev.py`, own implementation of Faddeev's integral from the formula in §1): (a) Lemma Φ (i): the integral equals the logarithm of the product formula to 5·10⁻²⁶ … 2·10⁻²⁵ at |x| = 0.2, arg x = ±20° (|q| < 1 and |q| > 1), two values of w each. (b) At real x = 0.0477 the inversion relation holds to 10⁻²⁶ … 10⁻²⁵. (c) An independent Padé–Laplace Borel sum of Φ at t = −0.3i agrees with the value of P·J at real x quoted in §6 to 4.6·10⁻¹⁴, the double-precision limit stated there. Lemmas S, F, D and Theorem 3 were read; their estimates were not checked line by line.

Labels: **proved** (a complete argument is written here; nobody else has checked it), **verified** (numerics), **inference**.

---

## 0. Summary

1. **Theorem 3 (proved, §1–§4).** For −150° ≤ arg t ≤ −30° (for x: |arg x| ≤ 60°, which contains the positive real x-axis, |q| = 1) and every K,

   | P·J − Σ_{k<K} aₖ·tᵏ | ≤ C_K·|t|ᴷ, uniformly in the direction,

   for the ψ₁\* and the ψ₀\* integral. Here Φ_b is defined by Faddeev's integral, which is one analytic function of b across the real axis.
2. **Corollary.** Together with Theorem 2 of v15: P·J is analytic on −210° < arg t < 10° and has the asymptotic expansion Φ uniformly on the closed sector, an opening of 220°. The gluing question of v15 is settled: the two product formulas are the two faces of the integral (Lemma Φ, with proof).
3. **How.** In the strip |Im w| < π the integral gives the expansion of Φ_b with a remainder that does not see whether x is real (Lemma S), and a bound for the leading term that holds up to the edges of the strip (Lemma F). Outside the strip, the functional equation in the direction i/b brings the argument back, at the cost of factors in the dual variable that are 1 up to e^{−c/|t|} on our contour (Lemma D). The leading exponent was already under control at real x (Proposition 1 of v15).
4. **Numerically (verified):** at real x, where no q-series exists, the kernel integral computed from Faddeev's integral equals the Borel sum of Φ to 4·10⁻¹⁴ (double precision), for ψ₁\* and ψ₀\*.
5. **What is left for Borel summability:** only the factorial growth of the constants C_K (item 3 of v14 §5). The sector is already wider than the 180° that Watson's theorem needs. §5 says why I expect the bookkeeping to work and what form it should take; it is not done.

---

## 1. Φ_b by Faddeev's integral

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

## 2. Two uniform lemmas in the strip

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

## 3. The contour

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

---

## 4. Theorem 3

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

---

## 5. What is left: factorial constants

Watson's theorem (or Nevanlinna's version) gives P·J = Borel sum of Φ once the constants satisfy C_K ≤ A·Bᴷ·K! on a sector wider than 180°. The sector is there. The constants are not yet tracked. Where they come from:

- **Lemma S** has factorial constants by construction: |a_k| is of size 2/π^{2k}, so |r_K(u)| ≤ C·(|u|/c₀)^{2K−1}, and the σ-integral of |σ|^{2K−2}·e^{(π−ε₁)|σ|}/|sinh πσ| is of size (2K−2)!/ε₁^{2K−1}. That is a remainder of size (2K)!·(C·|t|)^{2K−1}, the right kind.
- **The window step** is the delicate one. The crude version in v14 (split at |s| = |t|^{−1/2}, polynomials of degree 2m against e^{−κ|s|}) does not give factorial constants as written. The clean mechanism should be: for fixed s, the truncated integrand is analytic in t on a disc of radius about π/|s| (the condition |w| < π), with a bound C·e^{−κ|s|} there; Cauchy's estimate for the Taylor remainder in t gives (|t|·|s|/π)ᴷ, and ∫ |s|ᴷ·e^{−κ|s|} ds = 2·K!/κ^{K+1}. This needs the estimates of §2–§4 for complex t in a disc around directions, combined with the Gevrey remainder of Lemma S and an optimal truncation. I have not written it.

Until that is done, the statement "P·J is the Borel sum of Φ" remains verified (v10, v13, and §6 below at real x) and not proved.

Also still open: Theorem A of v12 for |q| > 1; the asymptotics of 𝓑ⱼ; the Stokes jump; identity (VII).

---

## 6. Numerical checks (verified; `z1_checks.py`)

| check | result |
|---|---|
| Lemma Φ (i): integral against the product formulas, θ = −70°, −80° (\|q\| < 1), −100°, −120° (\|q\| > 1), three w each | 10⁻²⁴ or better |
| Lemma Φ (ii) at real x: the two functional equations and the inversion relation of the integral, six cases | 2·10⁻²⁵ or better |
| fast double-precision Φ_b (strip integral, functional equation, inversion) against the products, 15 points inside and outside the strip | 2·10⁻¹² or better |
| Lemma S: defect/\|t\|^{2K−1}, K = 2, 3, \|t\| = 0.4 … 0.05, three w (one with Re w > 0), θ = −90°, −85°, −100° | constant in \|t\| to two digits, and the same in all three directions (e.g. 2.9·10⁻⁴ and 2.4·10⁻⁵ at w = −0.5 + 1.2i) |
| Lemma F: sup over Re w ∈ [−12, 0], \|Im w\| ≤ π, for θ = −90°, −75°, −120°, −60°, \|t\| = 0.4 … 0.05 | 0.347 to 0.351, attained at w = ±iπ |
| Lemma D: minimum of cos γ − \|sin α·sin 2β\| over \|β\| ≤ 30°; minimum of π/\|sin α\| | 0.383; 4.89 |
| uniform bound: sup over the contour (\|w\| ≤ 30, both halves) of \|H\|/e^{ψ/\|t\|}, θ = −90°, −89°, −91°, −85°, −95°, −60°, −120°, −150°, −30°, \|t\| = 0.4, 0.2, 0.1 | 0.70 to 0.75, not growing as \|t\| → 0 |
| double-precision integrator for P·J against the product-formula integrator of v13 at θ = −80°, \|t\| = 0.3 | 4·10⁻¹⁴ |
| **P·J at real x = 0.047746 (\|q\| = 1, \|t\| = 0.3) against the Borel sum of Φ** | ψ₁\*: 0.330580057167 − 0.033046552525i, difference 4·10⁻¹⁴; ψ₀\*: 0.332329469738 + 0.000288013873i, difference 4·10⁻¹⁴ |
| Theorem 3 at θ = −90°: \|(P·J − Σ_{k<K} aₖtᵏ)/(a_K·t^K) − 1\| at \|t\| = 0.2, 0.1, 0.05 | ψ₁\*: K = 1: 0.056, 0.028, 0.014; K = 3: 0.052, 0.020, 0.0085; K = 6: 0.47, 0.26, 0.13. ψ₀\*: K = 2: 0.19, 0.099, 0.050; K = 3: 0.25, 0.13, 0.064; K = 6: 0.33, 0.17, 0.088 |

The parts with real x use double precision throughout (about 10⁻¹²), because Φ_b at real b is only available through the integral.

---

## 7. Scripts (`psi_star/v16/`)

- `fad.py`: Faddeev's integral in the variable w (mpmath, in the strip; numpy, any w via functional equation and inversion); the product formulas.
- `z1_checks.py [a b c d e]`: the checks of §6 (a, b: mpmath, one minute; c, d, e: numpy and scipy, about one minute).
- `borel.py`, `phi_exact.json`: from v13.
- `run_log.txt`.
