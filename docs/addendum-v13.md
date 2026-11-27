# Mock theta search: addendum v13 (the perturbative series from the kernel integral; the kernel integral is the Borel sum on both sides of the real axis)

Session of 26–27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). First step towards the open item of `psi-star-addendum-v12.md` §7: that P·J is the Borel sum of the perturbative series Φ. Scripts and output are in `psi_star/v13/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v13/check_borel.py`): the first 41 coefficients in `phi_exact.json` are identical to an independent exact computation (`psi_star/v10/phi_psi1.py`, log-cosh method). The central claim of §3 was tested with own code (Padé [50/50] of the Borel transform at 300 digits, Laplace integral along the direction of t; own quadrature of J with Φ_b continued to |q| > 1 by inverted products): at |t| = 0.3, arg t = −80° (|q| = 0.949) P·J equals the Borel sum to 6·10⁻³⁷, and at arg t = −135° (Im x < 0, |q| = 1.236, where ψ₁\* does not converge) to 7·10⁻³⁹; the result does not depend on the contour angle (c = 0.5 and 0.65 agree exactly). Both values agree with the table of §3.

Nothing here is a proof. Labels: **derivation** (formal computation), **verified** (numerics), **inference**.

---

## 0. Summary

1. **Closed formula for Φ** (§1, derivation; checked exactly). Around y = 0 the kernel integral is a 1/cosh window, and expanding under it gives every coefficient of Φ as an explicit rational number built from Bernoulli polynomials. The formula reproduces all 151 coefficients of v3 for ψ₁\* and for ψ₀\* to the 120 digits stored there, and `phi_exact.json` now holds 301 exact coefficients of each.
2. **That agreement is weaker than I announced.** It turns out to be an identity: the window expansion is the termwise expansion of ½f₁\*, which is how v3 defined Φ. So it confirms the formula, not the asymptotic claim "P·J ~ Φ". The analytic part (the dual nome is negligible, tail estimates) is not done.
3. **The analytic claim was tested directly instead** (§3, verified). The kernel integral continues across the real x-axis into |q| > 1, where ψ₁\* does not exist, and equals the Borel sum of Φ there. At |t| = 0.3, over directions arg t from −220° to −60°: agreement to 10⁻²⁸…10⁻³¹. Over the wider range −250° to +30° the agreement is as good as the resummation allows (10⁻¹³ to 10⁻²⁰).
4. **Borel plane** (§2, verified for the first pair, otherwise indicative). The nearest singularities are S = π²/3 − i·Cl₂(π/3) and its conjugate. A second group lies near S − 2π² and its conjugate. Between them there are two sectors without singularities, each about 159° wide.
5. **The statement to prove** (§4, inference): P·J(x) = 𝒮_θ[Φ](t) on almost all of the cut x-plane, with any θ in the lower free sector. Watson's theorem needs factorial bounds on a sector of opening just over 180°, and the tested range leaves room for that.

---

## 1. The 1/cosh window and a closed formula for Φ (`s1_window.py`)

Notation of v12: q = e^{2πix} = e^{−t}, t = −2πi·x, p = e^{−πi/x} = e^{−2π²/t}, b² = x. Substituting y = b·s in P·J(1/(2b)), and using e^{πs}/(1 + e^{2πs}) = 1/(2·cosh πs), e^{2πxs} = q^{−is}:

  P·J(1/(2b)) = q^{−1/8} · ∫ e^{t·s²/2} / (2·cosh πs) · (−q;q)∞/(−q^{1/2−is}; q)∞ · (−p·e^{2πs}; p²)∞/(p;p²)∞ ds   (exact),

over the line e^{−iδ}·ℝ/b. Now, formally:

- drop the last factor (p is exponentially small as t → 0 along a ray);
- expand the ratio of q-products by Euler–Maclaurin: log(−q^a; q)∞ ~ π²/(12t) − Σ_{k≥1} B_k(a)·g^{(k−1)}(0)·t^{k−1}/k!, with g(u) = log(1 + e^{−u});
- integrate the resulting polynomials in s against 2^{−is}/(2·cosh πs).

With σ = i·s the result is

  **Φ(t) = 𝔼[ exp( t/8 − t·σ²/2 − Σ_{k≥2} [B_k(1) − B_k(½ − σ)]·g^{(k−1)}(0)·t^{k−1}/k! ) ],  𝔼[σʲ] = j!·[uʲ] e^{u/2}/(2 + e^u),**

where the exponential is expanded in powers of t and 𝔼 is applied to each coefficient, a polynomial in σ. For ψ₀\* (u = b + 1/(2b), prefactor P·q^{1/2}) the exponent gets the extra term t·(σ − ½).

- First coefficients: ψ₁\*: 1/3, 1/9, 5/162, 5/486, 29/17496, −1937/262440. ψ₀\*: 1/3, 0, 1/81, 1/81, 139/8748, 203/8748.
- Check: 151 coefficients each against `q2_series_v3.json`: worst relative deviation 4·10⁻¹²⁰ and 5·10⁻¹²⁰ (the precision of that file). The coefficient of t in the ψ₀\* series, stored there as 9·10⁻¹⁰¹⁸, is exactly 0.
- `phi_exact.json`: coefficients a₀…a₃₀₀ of both series as exact fractions (40 seconds; a₃₀₀ ≈ −2·10⁴⁵¹).

**Why the agreement is automatic.** 𝔼[σʲ] = ½·Σ_{n≥0} (−1)ⁿ·(n + ½)ʲ·2^{−n}: the moment functional evaluates at σ = n + ½, that is at the lower kernel poles y = −i·b·(n + ½). There the integrand is exactly q^{n²+n}/(2·(−q;q)ₙ), a term of ½f₁\*. So the formula is the termwise expansion of ½f₁\*(e^{−t}) at t = 0, and v3 computed Φ in exactly that way. The check confirms the Bernoulli formula and the bookkeeping of v12; it does not test whether P·J has Φ as its asymptotic series. For that one needs: the factor with p contributes only exponentially small terms, and the Euler–Maclaurin expansion holds uniformly up to |s| ~ |t|^{−1/2}, beyond which the Gaussian cuts off. Not done.

---

## 2. Singularities of the Borel transform (`s2_poles.py`)

Padé approximants [40/40] and [60/60] of the Borel transform Σ aₖξᵏ/k!, from the exact coefficients.

- **Nearest singularities** (verified): the poles accumulate at S = π²/3 − i·Cl₂(π/3) = 3.28987 − 1.01494i and at its conjugate, |S| = 3.4429, arg S = ∓17.145°. Distance of the nearest pole to S: 0.0027 ([40/40]), 0.0012 ([60/60]), the same for ψ₀\*. About 70% of all poles lie within 20° of the positive real axis, strung outward from these two points, which is how a Padé approximant shows a branch point.
- **Second group** (indicative): poles in the left half-plane starting near −16.5 ± 1.0i. S − 2π² = −16.449 − 1.015i. In three of the four approximants the nearest pole is within 0.09 of S − 2π² or its conjugate; in the fourth (ψ₀\*, [60/60]) it is at −16.17 on the real axis. A tower S + 2π²·k would fit p = e^{−2π²/t}: every power of the dual nome shifts the action by 2π² (inference).
- No poles closer than |ξ| = 97 in the directions 20° < |arg ξ| < 120°.

If the singularities are S + 2π²k and their conjugates (inference), the free sectors are 17.15° < |arg ξ| < 176.47°.

---

## 3. The kernel integral against the Borel sum, all directions (`s3_sectors.py`, `s3b_lateral.py`)

**Continuation.** With P = 2^{−1/2}·e^{π²/(12t) + t/6}/b (the elementary form of the prefactor, v6) and with the q-products continued by (a;q)∞ = 1/(a/q; 1/q)∞ for |q| > 1, the function P·J(x) is defined by the same integral on both sides of the real x-axis. The contour angle is δ = c·(π/2 − |arg b|) with c = 0.5 and 0.65; the two values agree to 10⁻³⁰ or better in every row.

**Comparison.** t = 0.3·e^{iφ}, x = i·t/(2π). 𝒮_θ is the Borel sum in direction θ (Padé–Laplace), with θ = φ where φ lies in the free sector and a lateral θ otherwise.

| φ = arg t | arg x | \|q\| | \|P·J − 𝒮_θ[Φ]\| | limited by |
|---|---|---|---|---|
| −250° | −160° | 1.108 | 3·10⁻¹³ | resummation: 5·10⁻⁹, 2·10⁻¹⁰, 5·10⁻¹¹, 2·10⁻¹², 3·10⁻¹³ for Padé orders 40, 60, 75, 90, 120 |
| −220° | −130° | 1.258 | 1.3·10⁻²⁸ | resummation |
| −190° | −100° | 1.344 | 1.5·10⁻³¹ | working precision (30 digits) |
| −160° | −70° | 1.326 | 2.2·10⁻³¹ | working precision |
| −135° | −45° | 1.236 | 1.9·10⁻³¹ | working precision |
| −110° | −20° | 1.108 | 2.5·10⁻³¹ | working precision |
| −100° | −10° | 1.054 | 2.0·10⁻³¹ | working precision |
| −80° | 10° | 0.949 | 1.0·10⁻³¹ | working precision |
| −60° | 30° | 0.861 | 4·10⁻²⁹ | resummation |
| −30° | 60° | 0.771 | 1.8·10⁻²¹ | 20 digits of P·J were carried into the lateral test |
| 0° | 90° | 0.741 | 3.8·10⁻²¹ | the same; θ = −45°, −60°, −75° |
| 30° | 120° | 0.771 | 4·10⁻²⁰ | θ = −45°, Padé order 75 |
| 60° | 150° | 0.861 | 3·10⁻⁵ | θ is 85° away from φ; not informative |

For ψ₀\* (P·q^{1/2}·J(b + 1/(2b)), 25 digits): φ = −190°, −135°, −80°, −60° give 2·10⁻²⁶, 8·10⁻²⁷, 3·10⁻²⁶, 1.6·10⁻²⁴.

Notes.

- The values at φ = −100° and −110° and at φ = −80° lie on the two sides of the real x-axis and change smoothly (0.3251 − 0.0316i, 0.3201 − 0.0294i against 0.3364 − 0.0335i); the Borel sum follows.
- At φ = −30°, 0°, 30° the Padé approximants of order 90 and 120 fail for θ = −25° (errors up to 5·10⁻³) and work for θ = −45° and below, as long as θ stays within about 80° of φ. I take that to be the nearness of θ = −25° to the cut starting at S; it is a property of the resummation, not of the integral.
- φ = 0° and 30° are beyond the direction of S (−17.1°). There P·J agrees with the Borel sum taken from below. This fits v10 §3: the Borel sum along the direction of t itself jumps at −17.1°, and P·J does not.

---

## 4. What follows, and what is next

**The statement** (inference from §2, §3): for −266° < arg t < 73°, which is almost all of the x-plane cut along the negative real axis,

  P·J(x) = 𝒮_θ[Φ](t),  t = −2πi·x,  for any θ in the lower free sector (−176.5°, −17.1°) with |arg t − θ| < 90°.

So the trivial sector 𝒯 of v10 is one analytic function on the cut plane. For Im x > 0 it is ψ₁\* minus the two bilinear terms (v12). For Im x < 0 the series ψ₁\*(q) and f₁\*(q) diverge, and the integral still exists.

**For a proof** by Watson's theorem one needs: P·J analytic on a sector in t of opening more than 180°, with Φ as asymptotic series and factorial bounds on the remainders, uniformly there. A sector such as −190° < arg t < −5° would do. It contains the positive real x-axis in its interior, which is natural for the integral and impossible for the q-series.

**Next steps, in order.**

1. The estimates named at the end of §1, on that sector: this would prove P·J ~ Φ.
2. Factorial bounds for the same expansion (Bernoulli polynomials grow factorially, so the formula of §1 is the right starting point).
3. Theorem A of v12 for |q| > 1: which q-series the residues give there. The sums that converge for |q| > 1 are ψ₁\*(1/q) and f₁\*(1/q), so I expect the mirror image of v12 (inference).
4. More coefficients are cheap now, so the second group of singularities can be located properly (Padé of order 150, or a direct fit of the large-order behaviour).

**Not checked.** Whether the continuation of Φ_b by inverted products is the standard one I assumed without looking it up; the agreement with the Borel sum at seven points with |q| > 1 is the evidence. Only |t| = 0.3 was used.

---

## 5. Scripts (`psi_star/v13/`)

- `s1_window.py [N]`: exact coefficients from the formula of §1; comparison with v3; writes `phi_exact.json`. Needs python-flint.
- `borel.py`: Padé–Laplace Borel sums in a given direction; P·J(x) on both sides of the real axis.
- `s2_poles.py`: poles of the Padé approximants.
- `s3_sectors.py`, `s3b_lateral.py`: the table of §3.
- `phi_exact.json` (301 exact coefficients for ψ₁\* and ψ₀\*), `q2_series_v3.json` (copy from v3), `run_log.txt`.
