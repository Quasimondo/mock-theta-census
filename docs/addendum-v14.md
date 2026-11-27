# Mock theta search: addendum v14 (first estimate: P·J ~ Φ is proved on a sector; straight contours descend on 320° of directions)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). Start of the estimate named in `psi-star-addendum-v13.md` §4: that the kernel integral P·J has the perturbative series Φ as its asymptotic expansion. Scripts and output are in `psi_star/v14/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v14/check_thm1.py`): U(0) = 0, U′(0) = −ln 2, U″ ≤ −1 on [−8, 8], the reflection formula (10⁻³¹), θ₀ = 12.4421° and max U = π²/60 at ln((√5 − 1)/2) confirmed. Theorem 1 at arg t = −45° with own quadrature: |ratio − 1| = 0.058, 0.028, 0.014 (K = 1) and 0.050, 0.014, 0.007 (K = 3) at |t| = 0.2, 0.1, 0.05, shrinking like |t|. The proof was read; its estimates were not checked line by line.

Labels: **proved** (a complete argument is written here; nobody else has checked it), **verified** (numerics), **exploration** (floating-point scan, no proof), **inference**.

---

## 0. Summary

1. **Theorem 1 (proved, §1–§2).** On every closed subsector of −90° < arg t < −12.44° (for x: 0° < arg x < 77.56°), and for every K,

   P·J(1/(2b)) = Σ_{k<K} aₖ·tᵏ + O(|t|ᴷ) as t → 0,

   with the coefficients aₖ of Φ (v3, exact in v13). The same holds for the ψ₀\* integral and its series. The bound −12.44° is −arctan(ln 2/π).
2. **With v12** this is the first proved asymptotic statement about the one-sided sums: on that sector, ψ₁\*(q) − (κ/2)·[𝓑₂·p·W₃(p) − 𝓑₄·W(p)] has the asymptotic expansion Φ(t).
3. **Why it works there.** On the contour where τ = i·t·s is real, the leading exponent of the integrand is (1/t) times a real, concave function U(τ) with U(0) = 0 and U′(0) = −ln 2. Concavity gives the decay on the whole contour in one line. The dual nome drops out once the integrand is written in its two forms on the two halves of the contour.
4. **This is not yet Borel summability.** Watson's theorem needs more than 180°, and the sector of Theorem 1 is 77.6° wide. Its two ends have different causes: at −12.44° the real contour stops descending; at −90° (real x) the product formulas lose uniformity.
5. **Exploration (§4).** At leading order, straight contours descend for −257° < arg t < 65°, a range of 322° that contains the real x-axis and the |q| > 1 side. On the grid the descent holds on exactly the set of contours allowed by the linear term at the origin. So the first limit in item 4 is a matter of proving one inequality for complex contours; the second needs a different representation of Φ_b near real b.

---

## 1. Statement

Notation of v12, v13: q = e^{−t}, t = −2πi·x = |t|·e^{iθ}, p = e^{−2π²/t}, b² = x, β = arg b = (θ + 90°)/2. Let θ₀ = arctan(ln 2/π) = 12.44°. Let aₖ and aₖ⁰ be the coefficients of Φ for ψ₁\* and ψ₀\* (v13 §1).

**Theorem 1.** Fix ε > 0 and K ≥ 0. There are C and t₀ > 0 such that for all t with −90° + ε ≤ θ ≤ −θ₀ − ε and |t| ≤ t₀:

- | P·J(1/(2b)) − Σ_{k<K} aₖ·tᵏ | ≤ C·|t|ᴷ,
- | P·q^{1/2}·J(b + 1/(2b)) − Σ_{k<K} aₖ⁰·tᵏ | ≤ C·|t|ᴷ.

**Corollary** (with the corollary of v12). On the same sector,

- ψ₁\*(q) − (κ/2)·[𝓑₂·p·W₃(p) − 𝓑₄·W(p)] ~ Σ aₖtᵏ, and ½·f₁\*(q) + (κ/2)·[𝓑₂·F_b(p) − 𝓑₄·F_a(p)] ~ Σ aₖtᵏ,
- the same for ψ₀\*, f₀\* with the formulas of v12.

This says nothing yet about the size of the bilinear terms themselves; that needs the asymptotics of 𝓑ⱼ (v7 §3, a derivation).

---

## 2. Proof

### 2.1 The contour and the integrand

In the sector, 0 < β < 45°, so δ = β is an admissible angle for the line L = e^{−iδ}·ℝ of v12 (the condition there is 0 < δ < 90° − β). J does not depend on the admissible δ: turning L keeps both ends in directions of Gaussian decay and crosses no pole.

On this L write y = b·s and s = −i·τ/t. Then τ = i·t·s runs over the real line, e^{2πby} = e^{τ}, and (v13 §1)

  P·J(1/(2b)) = e^{t/8} · ∫ H ds,  H = e^{−τ²/(2t)} · R(τ) · Π(τ) / (2·cos(πτ/t)),

  R(τ) = (−q;q)∞/(−q^{1/2}e^{τ}; q)∞,  Π(τ) = (−p·w; p²)∞/(p;p²)∞,  w = e^{2πs} = e^{−2πiτ/t}.

Two facts about the geometry, both from θ < 0: Re s = τ·|sin θ|/|t| and Im s = −τ·cos θ/|t|. For the ψ₀\* integral the integrand has the extra factor e^{τ} and the prefactor is e^{t/8 − t/2}.

### 2.2 Lemma 1 (Euler–Maclaurin with remainder)

Let u ∈ ℂ with 1 + u·e^{−ty} in the right half-plane or in the disk |z − 1| ≤ 1 − d for all y ≥ 0 (some d > 0), and let F(y) = log(1 + u·e^{−ty}) with the principal logarithm. Then for a ∈ (0, 1] and K ≥ 1

  Σ_{n≥0} F(n + a) = −Li₂(−u)/t + Σ_{k=1}^{K} B_k(a)·(−t)^{k−1}·Li_{2−k}(−u)/k! + ρ_K,  |ρ_K| ≤ C_K·|t|^{K−1}·|u|·(1 + |u|)^{K−1}/(cos θ·d^K).

Proof. Euler–Maclaurin with offset a: the sum is ∫₀^∞ F − Σ_{k≤K} B_k(a)·F^{(k−1)}(0)/k! plus a remainder bounded by (sup|B̃_K|/K!)·∫₀^∞ |F^{(K)}|. Here ∫₀^∞ F = −Li₂(−u)/t (substitute u′ = u·e^{−ty}; the path stays in the region where log(1 + u′) is the principal branch), and F^{(k)}(y) = −(−t)ᵏ·Li_{1−k}(−u·e^{−ty}). For k ≥ 1, Li_{1−k}(−v) is v times a polynomial of degree k − 1 in v, divided by (1 + v)ᵏ. So |F^{(K)}(y)| ≤ C_K·|t|ᴷ·|u|e^{−y·Re t}·(1 + |u|)^{K−1}/d^K, and ∫₀^∞ e^{−y·Re t} dy = 1/(|t|·cos θ). ∎

**The hypothesis holds for u = e^{τ} with τ ≤ τ₁ := (π/4)·cot|θ|, with d = d₀ := 1 − e^{−(π/4)·cot|θ|}.** The point u·e^{−ty} has modulus e^{τ − φ·cot|θ|} and argument φ = y·|Im t| ≥ 0. For φ ≤ 90° its real part is ≥ 0. For φ > 90° its modulus is below e^{τ₁ − (π/2)cot|θ|} = 1 − d₀. In the sector, cot|θ| ≥ tan ε, so d₀ is bounded below.

**Consequences.** With a = ½, u = e^{τ} for the denominator of R and a = 1, u = 1 for the numerator, for τ ≤ τ₁ and every K:

  (2.1) log R(τ) − τ²/(2t) = U(τ)/t + Σ_{k=1}^{K} (−t)^{k−1}·cₖ(τ)/k! + O(|t|ᴷ),

  **U(τ) = −τ²/2 + Li₂(−e^{τ}) + π²/12**,  cₖ(τ) = B_k(1)·Li_{2−k}(−1) − B_k(½)·Li_{2−k}(−e^{τ}),

uniformly in τ ≤ τ₁ and in the sector. (Apply the lemma with K + 1; c₁ = −½·ln 2.)

### 2.3 Lemma 2 (the other form of the integrand)

For all τ,

  e^{−τ²/(2t)}·R(τ)·Π(τ) = 2^{−1/2}·e^{−π²/(12t) + t/12} · e^{−τ²/t} · (−q^{1/2}e^{−τ}; q)∞/(−p/w; p²)∞.

Proof. R·Π = (−q;q)∞/((p;p²)∞·Φ_b(bs)). Insert 1/Φ_b(y) = e^{−πiy²}·Φ_b(−y)/Φ_b(0)² (v12, 2.2), e^{−πi·x·s²} = e^{−τ²/(2t)}, (−q;q)∞/(p;p²)∞ = 2^{−1/2}e^{π²/(12t) + t/24} (v6) and Φ_b(0)^{−2} = e^{t/24 − π²/(6t)}. ∎

By Lemma 1 with u = e^{−τ} (τ ≥ 0, so again d ≥ d₀), the logarithm of the right side is [−τ² − Li₂(−e^{−τ}) − π²/12]/t + O(1), and this bracket is U(τ) again, by the reflection formula Li₂(−e^{τ}) + Li₂(−e^{−τ}) = −π²/6 − τ²/2.

### 2.4 Lemma 3 (uniform bound on the whole contour)

There is C such that for all real τ and all t in the sector with |t| ≤ t₀,

  |H| ≤ C·exp(−c₁·|s|),  c₁ = π·|sin θ| − ln 2·cos θ > 0.

Proof.

- **The dual nome.** For τ ≤ τ₁: |p·w| = exp((2π/|t|)·(τ·|sin θ| − π·cos θ)) ≤ e^{−(3π²/2)·cos θ/|t|}, so Π = 1 + O(e^{−c/|t|}). For τ ≥ 0, in the form of Lemma 2: |p/w| ≤ |p|, so 1/(−p/w; p²)∞ = 1 + O(|p|). On each half the dual nome enters only through these factors.
- **The q-products.** By (2.1) with K = 1 for τ ≤ τ₁, and by Lemma 2 for τ ≥ 0: |e^{−τ²/(2t)}·R·Π| ≤ C·exp(cos θ·U(τ)/|t|). U is real, so Re(U/t) = cos θ·U/|t|.
- **The window.** 1/|2·cos ζ| ≤ C·e^{−|Im ζ|} for ζ = πτ/t on these rays, and |Im ζ| = π·|τ|·|sin θ|/|t|. (For |Im ζ| ≥ 1 this is |2cos ζ| ≥ e^{|Im ζ|}(1 − e^{−2}). For |Im ζ| ≤ 1 the ray gives |Re ζ| ≤ (π/ln 2)·|Im ζ|, a compact set on which cos has no zero.)
- **Concavity.** U(0) = 0, U′(τ) = −τ − log(1 + e^{τ}), so U′(0) = −ln 2, and U″ = −1 − e^{τ}/(1 + e^{τ}) < 0. Hence U(τ) ≤ −ln 2·τ for all τ. The total exponent is at most (1/|t|)·[−ln 2·cos θ·τ − π·|sin θ|·|τ|] ≤ −c₁·|τ|/|t| = −c₁·|s|. ∎

c₁ > 0 is exactly tan|θ| > ln 2/π. In the sector c₁ is bounded below.

### 2.5 Lemma 4 (expansion inside the window)

For |s| ≤ |t|^{−1/2} and every M ≥ 1,

  e^{t/8}·H = h₀(s)·[ Σ_{m<M} tᵐ·P_m(s) + E ],  |E| ≤ C_M·|t|ᴹ·(1 + |s|)^{2M},  h₀(s) = 2^{−1/2 − is}/(2·cosh πs),

with polynomials P_m of degree ≤ 2m that do not depend on t, P₀ = 1.

Proof. Here |τ| = |t·s| ≤ |t|^{1/2} ≤ τ₁. U and the cₖ are analytic for |τ| < π. Insert their Taylor polynomials at τ = 0 into (2.1) with K = M + 1 and put τ = i·t·s:

- U(τ)/t gives −i·s·ln 2 (the factor 2^{−is}) and terms t^{j−1}s^{j}, j ≥ 2; the Taylor remainder is O(|τ|^{M+2}/|t|) = O(|t|ᴹ|s|ᴹ) because |t|·|s|² ≤ 1;
- (−t)^{k−1}·cₖ(τ) gives terms t^{k−1+j}s^{j}; remainder O(|t|ᴹ|s|^{M−k+1});
- Π = 1 + O(e^{−c/|t|}), and 2·cos(πτ/t) = 2·cosh πs exactly.

So log(e^{−τ²/(2t)}·R·Π) = −½ln 2 − i·s·ln 2 + Σ_{m=1}^{M−1} tᵐ·Q_m(s) + O(|t|ᴹ(1 + |s|)^{M+1}) with deg Q_m ≤ m + 1. Each tᵐQ_m(s) is bounded for |s| ≤ |t|^{−1/2}, so the exponential can be expanded with the stated remainder. ∎

### 2.6 Conclusion

Split ∫ H ds at |s| = |t|^{−1/2}.

- Outside: by Lemma 3 the contribution is O(e^{−c₁|t|^{−1/2}}), smaller than every power of t.
- Inside: by Lemma 4 it is Σ_{m<M} tᵐ·∫ h₀·P_m ds + O(|t|ᴹ), where the integrals are first over |s| ≤ |t|^{−1/2} and may be extended to the whole ray at the same exponentially small cost (|h₀| ≤ C·e^{−c₁|s|} on the ray).

So P·J = Σ_{m<M} μ_m·tᵐ + O(|t|ᴹ) with μ_m = ∫ h₀·P_m ds over the ray.

**μ_m = a_m.** The P_m come from the formal series of (2.1) with τ = i·t·s. By the addition theorem B_k(½ − σ) = Σ_j C(k,j)·B_{k−j}(½)·(−σ)ʲ and Li_{n}(−e^{τ}) = Σ_j τʲ·Li_{n−j}(−1)/j!, this formal series is the one of v13 §1 with σ = i·s. The moments ∫ sʲ·h₀ ds over the ray equal those over the real line (no pole of 1/cosh between them, and decay on all rays in between because their slope stays below π/ln 2). These are the moments of v13. ∎

**ψ₀\*.** The extra factor e^{τ} = e^{i·t·s} is bounded by e^{c₁|s|/2} once |t| ≤ c₁/2, and inside the window it is a power series in t·s. The rest is unchanged.

---

## 3. Numerical checks (verified; `u1_checks.py`)

| check | result |
|---|---|
| U(0), U′(0) + ln 2, reflection formula, U ≤ −ln 2·τ | 0, 0, defect 10⁻²⁹, holds on the grid |
| maximum of U: at τ = ln((√5 − 1)/2), value π²/60 | defect 5·10⁻³² |
| Lemma 2 at three points | 3·10⁻³⁰ |
| Lemma 1: defect/\|t\|ᴷ for K = 1, 3, 5, \|t\| = 0.2 … 0.025, θ = −75°, −45°, −20°, τ = −0.6 and 0.25 | constant in \|t\| to three digits (0.0148, 8.1·10⁻⁵, 3.7·10⁻⁶ at τ = −0.6) |
| Lemma 3: sup over τ ∈ [−6, 6] of \|H\| divided by the bound, \|t\| = 0.4 … 0.05, θ = −80° … −14° | between 0.71 and 1.35, not growing as \|t\| → 0 |
| Theorem 1: (P·J − Σ_{k<K} aₖtᵏ)/(a_K·t^K) at \|t\| = 0.2, 0.1, 0.05 | see below |

The last check: modulus of the ratio minus 1, at three values of |t|, each half the one before.

- ψ₁\*, θ = −60°, |t| = 0.2, 0.1, 0.05: K = 1: 0.058, 0.028, 0.014; K = 3: 0.051, 0.017, 0.0076; K = 6: 0.61, 0.29, 0.14; K = 10: 0.63, 0.30, 0.14.
- ψ₁\*, θ = −30°, |t| = 0.2, 0.1, 0.05: K = 1: 0.059, 0.029, 0.014; K = 3: 0.046, 0.011, 0.0065; K = 6: 0.86, 0.34, 0.15; K = 10: 0.98, 0.34, 0.15.
- ψ₁\*, θ = −80°, |t| = 0.4, 0.2, 0.1: K = 1: 0.11, 0.056, 0.028; K = 3: 0.15, 0.052, 0.019; K = 6: 0.81, 0.51, 0.27; K = 10: 0.81, 0.52, 0.27.
- ψ₀\*, θ = −45°, |t| = 0.2, 0.1, 0.05: K = 2: 0.24, 0.11, 0.052; K = 3: 0.31, 0.14, 0.068; K = 6: 0.45, 0.20, 0.093; K = 10: 0.83, 0.21, 0.057.

The deviations shrink roughly in proportion to |t|, as they must if the next term of the series is the leading error.

---

## 4. Exploration: where straight contours descend (`u2_scan.py`)

On a ray s = r·e^{−iγ}, γ = δ + β, the leading-order size of the integrand is exp(ψ(ρ)/|t|) with ρ = |t|·|r| and

  ψ(ρ) = Re( e^{−iθ}·[−w²/2 + L(w) + π²/12] ) − π·ρ·cos γ,  w = i·t·s,

where L(w) = Li₂(−e^{w}) for Re w ≤ 0 and its continuation −π²/6 − w²/2 − Li₂(−e^{−w}) for Re w ≥ 0. On the real contour of §2 this is cos θ·U(τ) − π·|sin θ|·|τ|. The scan computes sup ψ(ρ)/ρ over both halves of the ray, for θ from −265° to 85° in steps of 5° and all admissible δ in steps of 0.5°.

- **Descent exists for −255° ≤ θ ≤ 60°** (grid values), none at −260° and at 65°.
- On the grid, descent holds for exactly those δ with 0 < δ < 90° − |β| and |δ + β| < 77.56° (0 exceptions among 6409 grid points, not counting points within 0.26° of that boundary). The second condition is the linear term at the origin: tan|γ| < π/ln 2. In other words, the behaviour away from the origin never costs anything. Solving the two conditions gives −257.6° < θ < 65.1°.
- At θ = −90° (real x, |q| = 1) the best value is sup ψ/ρ = −3.14, as good as anywhere.
- For the real contour (δ = β) the scan changes sign between θ = −15° (−0.145) and −10° (+0.136), in agreement with −12.44°.

This is floating point and leading order only. It says where a proof by straight contours can work, not that it does.

---

## 5. What is next

For Watson's theorem one needs Theorem 1 with factorial constants on a sector wider than 180°, for example −190° < θ < −5°. Three pieces are missing.

1. **The inequality sup ψ/ρ < 0 for complex contours.** On the real contour it was concavity. For general δ the potential is complex. Either an analytic argument (the scan suggests the statement is simple: it holds whenever the linear term allows it), or a rigorous interval computation on a compact range plus the two asymptotic ends.
2. **Uniformity across the real x-axis.** Lemma 1 uses the q-products and degenerates as θ → −90° (cos θ → 0 and d₀ → 0). Near real b one has to work with Φ_b itself: its expansion in the strip from Faddeev's integral representation, and the functional equation in the direction i/b outside the strip, which brings in finitely many factors of the dual kind. I have not started this.
3. **Factorial bounds.** The constants C_K of Lemma 1 are of size (K−1)!/(2π·d₀)ᴷ by their construction; carrying this through Lemma 4 is bookkeeping, not done.

Also open, as before: the upper end of the true range (72.9° against 65.1° for straight contours) is not needed for Watson.

---

## 6. Scripts (`psi_star/v14/`)

- `u1_checks.py [a b c d e]`: the checks of §3 (parts a, b: seconds; c: about 5 minutes; d: about 20 minutes, most of it the direction −80°; e: 3 minutes).
- `u2_scan.py`: the scan of §4 (numpy, scipy; 20 seconds).
- `borel.py`, `phi_exact.json`: from v13 (kernel integral, exact coefficients).
- `run_log.txt`.
