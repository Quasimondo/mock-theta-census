# Mock theta search: addendum v26 (corrections after the stress test: analyticity of P·J, the e^{w}-versions, Theorem 7 with an explicit sector, smaller points)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). An answer to the findings F2, F3, F4 and the smaller points of `docs/psi-star-stress-test-v1.md` (Claude Opus 5.5). F5 (the ranges X ≥ 4, X ≥ 5) is `addendum-v25.md`. Scripts and output are in `psi_star/v26/`.

Labels: **proved** (a complete argument is written here; nobody else has checked it), **verified** (numerics), **measured** (numbers taken from the stress test, not recomputed here).

---

## 0. Summary

The stress test found no violated inequality. Its findings are places where a statement was incomplete. What changes:

1. **Analyticity of P·J (F3; proved, §1).** The proof of the Corollary of v16 said "the contour angle is fixed". That is true in the variable y, where the contour is the line L = e^{−10°·i}·ℝ for every t; in the variable s it turns with t, which is how the stress test read it. Written in y the argument is complete, and it needs neither a second contour angle nor independence of the contour.
2. **Lemma L and the e^{w}-versions (F2; proved, §2).** Lemma L is a statement for |t| ≤ t₁, and for the integrals with the extra factor e^{w} (those of ψ₀\*) the descent contours are admissible at real x only for |t| < π/3 (decaying) and |t| < π (growing). v19 did not say so. No theorem is affected: Theorem 6 is asymptotic, and the exact formulas at roots of unity use the line L.
3. **Theorem 7 with ε(a) = 60° for every 0 < a < 1 (F4; proved, §3).** On the horizontal line through the saddle, Re V has its maximum at the saddle (Lemma V, an elementary statement about λ²·|1 + λe^{iη}|²). The stress test had found this numerically.
4. **Smaller points (§4):** the numerical remark after Lemma F is corrected; the notation p^{n²/2} in Theorem C is fixed; three remarks on ranges and constants.
5. **Constants measured by the stress test (§5)** are recorded, since the proofs leave them implicit.
6. **Numerically (verified, §6):** Lemma V for 999 values of a; Theorem 7 at arg x = 30°, 45°, 60°; the growth rate 1 − (π/3)/|t| of the e^{w}-version; analyticity of P·J across the real x-axis.

---

## 1. Analyticity of P·J (the Corollary of v16 §4)

**Corollary (the statement of v16, with the range of |t| written out).** P·J is analytic in t on −210° < θ < 10°, 0 < |t| < t₀, where t₀ is the constant of the uniform bound.

Proof. With y = b·s the contour s = r·e^{−iγ}, γ = β + 10°, is y = |b|·r·e^{−10°·i}: the line L = e^{−10°·i}·ℝ of v12, for every t. Since w = 2πby and t = −2πib², e^{−w²/(2t)} = e^{−πiy²}, and

  P·J = e^{t/8}·(C₀/b)·∫_L e^{−πiy²}/(Φ_b(y)·2cosh(πy/b)) dy,  C₀ = 2^{−1/2}·e^{π²/(12t) + t/24},  b = (it/2π)^{1/2}.

- *The integrand is analytic in b for fixed y ∈ L.* Φ_b(y) is analytic in (b, y) for Re b > 0 (Lemma Φ). Its zeros y = −i·[(m + ½)·b + (n + ½)/b] have arguments in [−90° − |β|, −90° + |β|], the poles of the kernel y = i·b·(n + ½) have arguments ±90° + β. For |β| ≤ 60° neither set meets L, whose points have arguments −10° and 170°.
- *Domination.* The bound |H| ≤ C·e^{−κ|s|} (Lemma 3′ of v15 and v16 §4), uniform on the closed sector for |t| ≤ t₀, reads |integrand| ≤ C·e^{−κ|y|/|b|} in the variable y, and |b| is bounded above and below on compact subsets of the punctured sector.
- An integral over a fixed contour of functions analytic in a parameter, dominated locally uniformly, is analytic (Fubini and Morera). ∎

No contour other than L is used, so the questions raised in F3 (descent for δ ≠ 10°, independence of the contour for Im x ≤ 0) do not arise. The stress test's S3.7 (κ > 0 for δ between 5° and 15°) shows that the other route would work as well.

---

## 2. Lemma L and the integrals with the factor e^{w}

**Lemma L (v19 §2), completed.** Add to the hypotheses: every piece of Γ is a segment, a ray or a line, and |t| ≤ t₁ := c·cos ε₁/(2C), where ε ≤ ε₁ is the range of arg x, c is the constant of (L1) and C the constant in the factor e^{C|w|} of (L2). If (L2) holds without that factor there is no condition on |t| beyond the one of the bounds used.

Why. By (L1), (L2), on Γ: |E| ≤ C·e^{C|w\*|}·exp(cos ε·H_τ(w\*)/|t|)·exp(C·|w − w\*| − cos ε·c·min(|w − w\*|², |w − w\*|)/|t|). For |t| ≤ t₁ and |w − w\*| ≥ 1 the last exponent is at most −(c·cos ε/(2|t|))·|w − w\*|, which gives convergence and the tail estimate of the proof. Without the condition the integral over Γ need not converge.

**The e^{w}-versions at real x.** For I(b ± 1/(2b)) the integrand has the extra factor e^{w}, and (L2) holds with e^{|w|} (C = 1). At real x (t = −i|t|) the sizes on the right are:

- *Decaying integral (v19 §4), line Im w = η₁.* Between L and the line, and on the line, the size is e^{r·(1 + (2η + 2πa)/|t|)} with η ≤ η₁. So the integrand tends to 0 on the right, and the line is an admissible contour, exactly when |t| < −(2η₁ + 2πa). For a = ½ (η₁ = −2π/3) this is **|t| < π/3**, that is X ≥ 7 at x = 1/X.
- *Growing integral (v19 §5), contour Γ.* On Γ₃ and in the wedge between L and the real axis the size is at most e^{r·(1 − π/|t|)}. The condition is **|t| < π**, that is X ≥ 3.
- On the left the factor e^{w} only helps. For Im x > 0 the Gaussian factor gives decay for every t (v12 §2), and there is no condition.

v19 §4 proved the decay at infinity for the integrals without e^{w} and for the others only stated (L2). The condition is a missing hypothesis in "Moving the contour", as F2 says.

**Consequences.** None for the theorems: Theorem 6 is a statement for |x| → 0; Theorems 9, 11, 18 are stated with the integrals over L, which exist for every real x in their range. For numerical work: at small X the integrals I(b ± 1/(2b)) must be taken over L (the scripts of v20, v21 do this for X ≤ 8).

---

## 3. Theorem 7 with an explicit sector

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

## 4. Smaller corrections

1. **Lemma F of v16, the numerical remark.** It said that the supremum of |Λ(w) + Li₂(−e^{w})/t|·e^{−Re w/2} is 0.347, "the same for real x and for the neighbouring directions, and independent of |t|". Correct is: 0.347 = ½·ln 2 for real x and for |arg x| ≤ 10°; on the whole domain |arg x| ≤ 60°, |x| ≤ ¼ the supremum measured by the stress test is 0.511 (at |x| = ¼, arg x = ±60°, w = −0.441 ± iπ), and it depends on |t| and arg x. The lemma itself (existence of C₁) is not affected; C₁ = 0.52 works numerically.
2. **Notation in Theorem C (v12).** All fractional powers of q and p are defined through x: q^{α} = e^{2πixα}, p^{α} = e^{−πiα/x}. In D(u) the factor p^{n²/2} is e^{−πi·n²/(2x)}; with a principal power of the number p the odd-n terms can change sign.
3. **Lemma B (iii) of v19.** The factor 2^{N} is what the proof gives; numerically the factor 2 suffices (stress test, S2.8).
4. **(L1) on the contour Γ of v19 §5.** Proved for τ ≤ 0.2 (arg x ≤ 11.3°). Numerically it holds for τ < 0.3085 (arg x < 17.15°) (stress test, S6.3); the proof loses the difference in the chord bound for ln(2cos(η/2)).
5. **Lemma R2 of v20.** The bracket [1/(2K) − e^{−1/K + 2δ}/(2 − δ)] with K = 3/2 is positive for δ < 0.104, not for all δ; the proof uses it for δ → 0 only. (The argument gives e^{−1/K + δ}, for which the range is δ < 0.172.)
6. **Constants without slack.** Lemma A (equality at ρ = 0, θ = −30°), Lemma D (0.383022 at β = −30°) and Proposition 1 (iv) (κ = 0.185115 at θ = 10°) are sharp at an edge of their domains. The statements are correct as written.

---

## 5. Constants measured by the stress test

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

---

## 6. Numerical checks (verified; `u1_fixes.py`, `run_log.txt`)

| check | result |
|---|---|
| Lemma V: number of solutions of g = 1, and sup of Re V − Re V(w₁) on the line, for a = 0.001 … 0.999 (999 values) | one solution in every case; sup = 0 |
| local maximum of g at a = 0.9, 0.99 | 0.067, 0.063 (< 1) |
| Theorem 7′: (I(a/b)/formula − 1)/\|t\| at arg x = 30°, 45°, 60°, \|t\| = 0.4, 0.2, 0.1, for a = 0.1, 0.3, 0.5, 0.8, 0.95 (series of Theorem C, 90 digits) | 0.019, 0.040, 0.087, 0.051, 0.021 at most; flat in \|t\| |
| §2: slope of log\|integrand\| on the right of the line Im w = −2π/3 at X = 4, 6, 7, 8, 10, 12 against 1 − X/6 | agree to 4 decimals (+0.3333, 0, −0.1667, −0.3333, −0.6667, −1) |
| §2: I(b + 1/(2b)) over the line against L, X = 7, 8, 10, 12 | 4·10⁻¹⁰, 9·10⁻¹⁰, 5·10⁻⁹, 4·10⁻⁷ |
| §1: P·J over the fixed line L on circles \|t − t₀\| = 0.3\|t₀\| around x₀ = 0.1, 0.04 (arg x from −17.5° to 17.5°): negative Fourier modes c₋₁, c₋₂, c₋₃ | ≤ 3·10⁻¹⁵ (positive modes: 2·10⁻², 9·10⁻⁴, 6·10⁻⁵); mean value against P·J(t₀): 10⁻¹⁴ |

---

## 7. Scripts (`psi_star/v26/`)

- `u1_fixes.py [a b c d]`: the checks of §6 (2 minutes).
- `h1_roots.py` (from v20), `fad.py` (from v16), `common.py`, `a1_formula.py` (from v18), `run_log.txt`, `README.md`.
