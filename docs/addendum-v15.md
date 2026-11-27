# Mock theta search: addendum v15 (the descent inequality is proved; P·J ~ Φ on both sides of the real axis)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). Item 1 of `psi-star-addendum-v14.md` §5: prove that the leading exponent decreases along complex straight contours. Scripts and output are in `psi_star/v15/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v15/check_prop1.py`, `check_thm2.py`): Proposition 1 at 800 random admissible points (own G, own quadrature of ℓ): the exact formula (i) holds to 1.7·10⁻¹⁷, the bound (iii) has no violations, and for δ = 10°, −60° ≤ β ≤ 50° the minimum of κ is 0.18512 at β = 50°. Theorem 2 with own quadrature (K = 1, |t| = 0.2, 0.1, 0.05): |ratio − 1| = 0.059, 0.029, 0.014 at θ = 0° and 0.053, 0.027, 0.014 at θ = −150° (|q| > 1), agreeing with §3. The proof was read; its estimates were not checked line by line.

Labels: **proved** (a complete argument is written here; nobody else has checked it), **verified** (numerics), **inference**.

---

## 0. Summary

1. **Proposition 1 (proved, §1).** Along any admissible straight contour, the leading exponent ψ(ρ) of the integrand is given by an exact formula in which the only non-elementary quantity is an average of log(1 + z) over points with |z| ≤ 1. Two one-line bounds on that average give

   ψ(ρ) ≤ −κ·ρ − c·ρ²·sin 2δ,  κ = (π/2)·cos γ − max(ln 2, −λ_α)·|sin γ|.

   With the single contour angle δ = 10°, κ ≥ 0.185 for all directions −210° ≤ arg t ≤ 10°. This range contains the real x-axis (arg t = −90°) and the sector −190° < arg t < −5° that Watson's theorem needs.
2. **The argument is not concavity.** The exponent is concave on only about three quarters of the admissible contours (and on none for arg t ≤ −200°), so the proof of v14 does not generalise as it stands. The averaging argument replaces it.
3. **Theorem 2 (proved, §2).** P·J has Φ as its asymptotic expansion, uniformly on closed subsectors of −90° < arg t ≤ 10° (|q| < 1) and of −210° ≤ arg t < −90° (|q| > 1), for ψ₁\* and ψ₀\*. This extends Theorem 1 of v14 past its limit −12.44°, through arg t = 0 (q → 1 radially), and gives the first statement on the |q| > 1 side.
4. **What is left for Borel summability** (§4): the wedge around the real x-axis, where the constants of the Euler–Maclaurin lemma blow up although the leading exponent is harmless there; and factorial constants. Item 1 is done; items 2 and 3 of v14 §5 remain.

---

## 1. The descent inequality

### 1.1 Setting

t = |t|·e^{iθ} with −270° < θ < 90°, β = (θ + 90°)/2 = arg b. A contour angle δ is admissible if 0 < δ < 90° − |β|. Put α = β − δ and γ = β + δ; both lie in (−90°, 90°). The contour is s = r·e^{−iγ} (r real), and

  w := i·t·s = |t|·r·e^{iα},  ρ := |w|.

So w runs over the line e^{iα}·ℝ; the right half of the contour (r > 0) has Re w > 0. Let

  G(w) = −w²/2 + Li₂(−e^{w}) + π²/12 for Re w ≤ 0,  G(w) = −w² − Li₂(−e^{−w}) − π²/12 for Re w ≥ 0

(the two expressions agree for |Im w| < π; G is the function U of v14 continued to complex argument). The leading exponent of the integrand of P·J is ψ(ρ)/|t| with

  ψ(ρ) = Re( e^{−iθ}·G(w) ) − π·ρ·cos γ.

The first term comes from the Gaussian and the quantum dilogarithm, the second from the kernel, |1/(2cosh πs)| ≈ e^{−π|Re s|}.

### 1.2 Proposition 1

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

## 2. Theorem 2

### 2.1 Statement

For Im x < 0 (|q| > 1) let Φ_b := Φ_{1/b}, with the right side given by the product formula of v12 (legitimate, because Im(1/b²) > 0). Explicitly, with t₁ = t for Re t > 0 and t₁ = −t for Re t < 0, q₁ = e^{−t₁}, p₁ = e^{−2π²/t₁}, w′ = e^{2πs}:

  1/Φ_b(b·s) = [ (−p₁·w′; p₁²)∞ / (−q₁^{1/2}·e^{w}; q₁)∞ ]^{±1},  + for Im x > 0, − for Im x < 0.

This is the continuation used in v13. All identities of v12 for Φ_b that were stated for Im x > 0 hold for Im x < 0 by applying them to 1/b; in particular Φ_b(y)·Φ_b(−y) = Φ_b(0)²·e^{πiy²} with Φ_b(0)² = e^{−t/24 + π²/(6t)}. Let C₀ = 2^{−1/2}·e^{π²/(12t) + t/24} and

  P·J := e^{t/8} · ∫ H ds,  H = C₀·e^{−w²/(2t)} / (Φ_b(b·s)·2·cosh πs),

over the contour of §1. For Im x > 0 this is P·J(1/(2b)) of v12 (v13 §1, with C₀ = (−q;q)∞/(p;p²)∞ by v6). For Im x < 0 it is the definition. The ψ₀\* version has the extra factor e^{w} and the prefactor e^{t/8 − t/2}.

**Theorem 2.** Fix ε > 0 and K ≥ 0. There are C and t₀ > 0 such that for all t with −210° ≤ θ ≤ 10°, |θ + 90°| ≥ ε and |t| ≤ t₀:

  | P·J − Σ_{k<K} aₖ·tᵏ | ≤ C·|t|ᴷ,

and the same for the ψ₀\* version with its coefficients.

In terms of x: ε ≤ arg x ≤ 100° and −120° ≤ arg x ≤ −ε.

### 2.2 Proof

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

**Lemma 4′ (expansion in the window)** and the **conclusion** are those of v14 §2.5–2.6 with the real variable τ replaced by the complex w, |w| ≤ |t|^{1/2} ≤ η/4: the functions G and Li_{2−k}(−e^{w}) are analytic for |w| < π, the Taylor expansion at 0 with w = i·t·s gives the same polynomials P_m(s), and the moments of h₀(s) = 2^{−1/2−is}/(2cosh πs) over the ray equal those over the real line because |γ| ≤ 60° < 77.56°. The formal series is the one of v13 on both sides of the real axis because (2.1′) is. For the ψ₀\* version, |e^{w}| ≤ e^{|t|·|s|}. ∎

---

## 3. Numerical checks (verified; `w1_descent.py`, `w2_checks.py`)

| check | result |
|---|---|
| Proposition 1 (i): exact formula against ψ computed from G, 3000 random points (β ∈ (−85°, 85°), all admissible δ, ρ ≤ 30, both halves) | max deviation 4·10⁻⁹ (floating point) |
| Proposition 1 (ii): bounds on X, Y at the same points | hold: X − ln 2 ≤ −6·10⁻⁴, X − λ_α ≥ 0.028, \|Y\| ≤ 0.70 |
| Proposition 1 (iii) at the same points | 0 violations |
| κ for δ = 10°, −60° ≤ β ≤ 50° in steps of 0.05° | minimum 0.1851, at β = 50° |
| actual sup ψ/ρ for δ = 10° at nine directions | between −3.07 and −0.97; always below −κ |
| inversion relation for \|q\| > 1, six points | 2·10⁻³⁰ |
| (2.1′) on the contour, θ = −200°, −150°, −110°, −70°, 0°; ρ = 0.5 and 3; K = 1, 3, 5; \|t\| = 0.2 … 0.025 | defect/\|t\|ᴷ constant in \|t\| to two or three digits in all 30 cases |
| Lemma 3′: sup over the contour of \|H\|/e^{ψ/\|t\|}, θ = −200°, −150°, −110°, −70°, −20°, 0°, 8°; \|t\| = 0.4, 0.2, 0.1 | between 0.70 and 0.88, not growing as \|t\| → 0 |
| Theorem 2: (P·J − Σ_{k<K} aₖtᵏ)/(a_K·t^K) at \|t\| = 0.2, 0.1, 0.05 | see below |

Theorem 2, modulus of the ratio minus 1 at |t| = 0.2, 0.1, 0.05 (contour angle 10°):

- ψ₁\*, θ = 0° (t real, beyond the range of v14): K = 1: 0.059, 0.029, 0.014; K = 3: 0.038, 0.0057, 0.0059; K = 6: 1.18, 0.36, 0.15; K = 10: 2.3, 0.35, 0.15.
- ψ₁\*, θ = −150° (|q| > 1): K = 1: 0.052, 0.027, 0.014; K = 3: 0.049, 0.021, 0.0095; K = 6: 0.36, 0.22, 0.12; K = 10: 0.36, 0.22, 0.12.
- ψ₁\*, θ = −200° (|q| > 1): K = 1: 0.052, 0.027, 0.014; K = 3: 0.049, 0.021, 0.0095; K = 6: 0.35, 0.21, 0.12; K = 10: 0.36, 0.21, 0.12.
- ψ₀\*, θ = −130° (|q| > 1): K = 2: 0.17, 0.092, 0.048; K = 3: 0.21, 0.12, 0.061; K = 6: 0.28, 0.16, 0.083; K = 10: 0.15, 0.051, 0.016.

The deviations shrink roughly in proportion to |t|. At θ = 0° the kernel integral is real up to an exponentially small imaginary part: 1.3·10⁻⁷ at |t| = 0.2 and 1.2·10⁻¹⁴ at |t| = 0.1, against e^{−π²/(3t)} = 7·10⁻⁸ and 5·10⁻¹⁵ (at |t| = 0.05 it is below the working precision). That is what a lateral Borel sum on the real axis should look like.

---

## 4. What is left

Target, as before: Theorem 2 with factorial constants on a full sector wider than 180°, for example −190° < arg t < −5°, which lies inside the range −210° … 10° of Theorem 2.

1. **The wedge |arg t + 90°| < ε around the real x-axis.** Proposition 1 holds there with κ ≥ 0.185, so the leading exponent is not the problem. What fails is Lemma 1′ and the dual factor: as arg t → −90°, cos θ → 0, the paths of Lemma A become vertical and hit the lattice (η → 0), and |p₁| → 1. Near real b one has to work with Φ_b itself: Faddeev's integral representation in the strip, the functional equation in the direction i/b outside it. This is item 2 of v14 §5, not started.
2. **Gluing.** Theorem 2 treats the two sides separately. Watson's theorem needs one analytic function on the sector, so one also needs that Φ_b, defined by the two product formulas, is analytic across real b. That is a standard property of Faddeev's dilogarithm which I have not looked up; the evidence here is v13 §3.
3. **Factorial constants** (item 3 of v14 §5): bookkeeping, not done.

---

## 5. Scripts (`psi_star/v15/`)

- `w0_concavity.py`: on how many admissible contours the exponent is concave (numpy; seconds).
- `w1_descent.py`: Proposition 1 (numpy, scipy; one minute).
- `w2_checks.py [a b c d]`: the checks for Theorem 2 (mpmath; a, b: seconds; c: about 5 minutes; d: about 15 minutes).
- `borel.py`, `phi_exact.json`: from v13.
- `run_log.txt`.
