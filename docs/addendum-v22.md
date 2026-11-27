# Mock theta search: addendum v22 (the growing state integral for every z; F(z) at even roots of unity)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). The open item of v19 §6.3: the growing integral I((a−1)/b) for a ≠ ½. Scripts and output are in `psi_star/v22/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v22/check_thm13.py`, own code): for a = 0.01 … 0.99 (15 values) the saddle W (Im W > 0) satisfies V′(w₀) = 0 and Im V(w₀) = −Vol(z) (Bloch–Wigner) to 10⁻²⁵, and condition (L1) at real x holds on all five pieces C₁ … C₅ with the dilogarithm continued as stated in §1: the worst sup of (H − H(w₀))/min(d², d) is −0.0036 (a = 0.99, on C₄ and C₅) and −0.0095 (a = 0.01, on C₁ and C₂), the values of §6. The claim of §5 is reproduced: along the horizontal ray from the growing saddle, −Im V rises by 1.20 and 7.36 at s = 2, 8 for a = ½ and by 3.07 and 16.4 for a = 0.2. For Theorem 14, on the path x = 1/8 + iε the dual series 𝒲(T), 𝒲(p²T) tend to 1 and p·T to z^X (errors 0.02–0.03 at ε = 10⁻⁴, falling about tenfold per decade of ε). Theorem 14 itself (the integrals at real b) and part (b) of Theorem 13 were not recomputed.

Labels: **proved** (a complete argument is written here or in the addenda cited; nobody else has checked it), **verified** (numerics), **inference**.

**Later remarks:** `addendum-v27.md` §2 (the sector ε(a) of Theorem 13 (b), measured), §5 (margins and constants on the contour C, measured); `addendum-v25.md` (Theorem 14 at X = 2).

---

## 0. Summary

1. **Theorem 13 (proved; §1–§3).** For every 0 < a < 1, as x → 0 on the real axis (and for 0 ≤ arg x ≤ ε(a) with some ε(a) > 0), uniformly:

   I((a−1)/b) = (1/(2πb))·(2πt/(−V″(w₀)))^{1/2}·e^{V(w₀)/t}·(1 + O(t)),

   with W the root of W² + W = 1/z in the upper half-plane, w₀ = Log W, V(w) = −w²/2 − 2πi(a−1)w + Li₂(−e^{w}), V″(w₀) = −(1 + 2W)/(1 + W), and |e^{V(w₀)/t}| = e^{Vol(z)/|t|} at real x. This is the saddle selected in v18 §2.2. With Theorem 7 of v19 both state integrals of the z-family now have proved asymptotics.
2. **The contour** is the same for all a except for one vertical line: the ray Im w = 2π from the far left to 2πi; the imaginary axis down to iπ; the top edge of the strip to the point above the saddle; the vertical line through the saddle; the real axis. On these pieces Im V is, in this order, 2πa·|Re w|, a Clausen function, linear, monotone on each side of the saddle, and linear. No computer-assisted step is needed.
3. **Theorem 14 (proved; §4).** For even X ≥ 4, x = 1/X, q₀ = e^{2πi/X} and every z = e^{2πia} with 0 < a < 1, the finite sum F(z; q₀) = Σ_{n<X/2} zⁿ q₀^{n(n+1)/2} (−q₀;q₀)ₙ equals

   F(z; q₀) = i·P·e^{−πia}·[ J(a/b) − I((a−1)/b) + z^X·I(a/b) ].

   The lemmas of v20 carry over. Checked to 10⁻⁹ … 10⁻¹³ for six values of z and X = 8 … 40.
4. **What is still missing for the volume rule as a theorem:** a bound for the kernel integral J(a/b) at real x for a ≠ ½ (numerically P·J(a/b) stays close to 1/(1 − 2z)). Given such a bound, Theorems 13, 14 and Theorem 7 of v19 give F(z; e^{2πi/X}) = (the formula of v18 §2.2)·(1 + O(1/X)) for even X.
5. **A correction to a suggestion.** A sketch from another session proposed a horizontal ray leaving the saddle to the left, tested with Re V. At real x the size of the integrand is governed by −Im V, and along that ray −Im V rises (by +7.4 at distance 8 for a = ½). §5 records this.

---

## 1. Setting and statement

Notation of v19 §1: x → 0 with 0 ≤ arg x ≤ 60°, |x| ≤ ¼, t = −2πix, τ = tan(arg x); w = 2π·b·y; for u = a′/b

  I(a′/b) = (1/(2πb))·∫ e^{Q(w)/t}/Φ_b dw,  Q(w) = −w²/2 − 2πi·a′·w,  V = Q + Li₂(−e^{w}),

and |e^{V/t}| = exp(cos(arg x)·H_τ/|t|) with H_τ = −Im V + τ·Re V. Here a′ = a − 1 ∈ (−1, 0). Put η\* = −2πa′ = 2π(1 − a) ∈ (0, 2π).

Li₂(−e^{w}) is principal in the closed strip |Im w| ≤ π (boundary values from inside), and for Re w ≤ 0 it is the 2πi-periodic function given by the power series. The two agree where both apply.

**The saddle.** Let W be the root of W² + W = e^{−2πia} with Im W > 0, w₀ = Log W = r₀ + iη₀, 0 < η₀ < π. Then arg W and arg(1 + W) lie in (0, π), their sum is ≡ −2πa (mod 2π) and lies in (0, 2π), so η₀ + arg(1 + W) = η\*; and |W(1 + W)| = 1. Hence w₀ + 2πia′ + Log(1 + W) = 0: V′(w₀) = 0. By v18 §2.1, Im V(w₀) = D(−W) = −Vol(z) < 0.

**Theorem 13.** Let 0 < a < 1.

- (a) At real x, as x → 0⁺: I((a−1)/b) = (1/(2πb))·(2πt/(−V″(w₀)))^{1/2}·e^{V(w₀)/t}·(1 + O(t)), with the root for which the quotient by −i has positive real part. Equivalently (v18 §2.3, where the branch was checked numerically): −i·P·e^{−πia}·I((a−1)/b) = (1 + W)·W^{1/2}·(π/(t(1 + 2W)))^{1/2}·e^{(V(w₀) + π²/12)/t}·(1 + O(t)) with principal roots and W^{1/2} = e^{w₀/2}.
- (b) There is ε(a) > 0 such that (a) holds uniformly for 0 ≤ arg x ≤ ε(a).

The constants are locally uniform in a. They degenerate as a → 0 or 1, where Vol(z) → 0.

---

## 2. The contour at real x and its heights (proof of (a))

**Contour C** = C₁ ∪ … ∪ C₅, from left to right:

- C₁: w = r + 2πi, r from −∞ to 0;
- C₂: w = iη, η from 2π down to π;
- C₃: w = r + iπ, r from 0 to r₀;
- C₄: w = r₀ + iη, η from π down to 0;
- C₅: w = r, r from r₀ to +∞.

**Heights (Im V; here τ = 0).**

- C₁ (left half-plane formula): e^{w} = e^{r} is real, so the dilogarithm is real, and Im Q = −2πr − 2πa′r = −2πa·r. Hence Im V = 2πa·|r| ≥ 0.
- C₂ (left half-plane formula on the axis): Q(iη) is real, and Im Li₂(−e^{iη}) = Cl₂(η − π). Hence Im V = Cl₂(η − π) ∈ [0, Cl₂(π/3)] for π ≤ η ≤ 2π, for every a.
- C₃ (boundary values from inside the strip): −e^{w} = e^{r}; for r ≤ 0 the dilogarithm is real and Im V = (η\* − π)·r; for r ≥ 0 it is Li₂(e^{r} − i0) with imaginary part −π·r, and Im V = −2πa·r. On C₃, r does not change sign, so Im V is linear between its value 0 at iπ and its value at the corner r₀ + iπ.
- C₄: d/dη Im V(r₀ + iη) = Re V′ = −r₀ − ln|1 + e^{r₀ + iη}| =: f(η). Since |1 + e^{r₀}e^{iη}|² = 1 + 2e^{r₀}cos η + e^{2r₀} decreases strictly on (0, π), f increases strictly, and f(η₀) = 0 with f′(η₀) > 0. So Im V decreases on (0, η₀], increases on [η₀, π), and has a nondegenerate minimum at w₀.
- C₅: the dilogarithm is real, Im V = η\*·r = 2π(1 − a)·r.

**(L1).** By C₄, the values at both ends of C₄ exceed Im V(w₀): the corner r₀ + iπ and the point r₀ on the real axis. So on C₃, Im V ≥ min(0, value at the corner) > Im V(w₀); on C₅, Im V ≥ η\*·r₀ > Im V(w₀) and grows linearly; on C₁, C₂, Im V ≥ 0 > −Vol, with linear growth on C₁. Together with the nondegenerate minimum: H(w) − H(w₀) ≤ −c·min(|w − w₀|², |w − w₀|) on C, for H = −Im V.

**(L2).** On C₃, C₄, C₅ (closed strip): Lemma B (i) of v19, and (ii) near w₀. On C₁, C₂ (Re w ≤ 0, π ≤ Im w ≤ 2π): one step of the functional equation, 1/Φ_b(y) = (1 + p·w′)/Φ_b(y − i/b), with |p| = 1 and |w′| = e^{Re w/x} ≤ 1 at real x, so the factor is at most 2; the shifted point has −π ≤ Im ≤ 0, Re ≤ 0 and the same e^{w}, where Lemma F of v16 applies. Hence |integrand| ≤ C·|e^{V/t}| on all of C. (C₂ runs through poles of Φ_b; these are zeros of the integrand.)

**Moving the contour** from the line L of v12 to C: the zeros of Φ_b lie on the negative imaginary axis below −iπ (real x), outside the two regions swept; in the left region η ≥ 2π > η\* eventually and the integrand is at most 2^{N}·C·e^{−|r|(η − η\*)/|t|}, in the right region it is at most C·e^{−r(η\* − 2η)/|t|} with η ≤ 0 (inversion relation), so it tends to 0 at infinity.

**Result.** Lemma L of v19 on C₄ with d = −i. ∎

---

## 3. Small arg x (proof of (b))

For arg x = φ > 0 the bound on C₁, C₂ fails on the imaginary axis (|p·w′| = e^{(Im w − π)·sin φ/|x|} there). Tilt these two pieces: fix φ₁ > 0, δ₁ = 90° − φ₁, P = iπ − ρ_P·e^{−iδ₁} with ρ_P·sin δ₁ = π (so Im P = 2π, Re P = −π·tan φ₁), and let

- C₁′: w = r + 2πi, r from −∞ to Re P;  C₂′: the segment from P to iπ.

- **Bound.** On C₂′, Lemma B (iii) of v19 applies for φ ≤ φ₁. On C₁′, |p·w′| = exp((−|r|·cos φ + π·sin φ)/|x|) ≤ 1 because |r| ≥ π·tan φ₁ ≥ π·tan φ.
- **Heights at τ = 0.** On C₁′ still Im V = 2πa·|r| ≥ 0. On C₂′ write w = r + iη, |r| = (η − π)·tan φ₁, φ′ = η − π ∈ [0, π]. Then Im V = |r|·(η − η\*) + Im Li₂(e^{−|r|}·e^{iφ′}), and since ∂/∂(ln R) Im Li₂(R·e^{iφ′}) = −arg(1 − R·e^{iφ′}) ∈ [0, (π − φ′)/2] for 0 < R ≤ 1: Im Li₂(e^{−|r|}e^{iφ′}) ≥ Cl₂(φ′) − |r|·(π − φ′)/2. Hence Im V ≥ −|r|·(η\* − η)₊ − |r|·(2π − η)/2 ≥ −(3π²/2)·tan φ₁. Choose φ₁ with (3π²/2)·tan φ₁ < Vol(z)/2. Then (L1) holds at τ = 0 on C′ = C₁′ ∪ C₂′ ∪ C₃ ∪ C₄ ∪ C₅.
- **τ > 0.** Re V − Re V(w₀) ≤ M·min(|w − w₀|², |w − w₀|) on C′ for some M: Re V is smooth along the pieces, V′(w₀) = 0, and Re V → −∞ on C₁′ and C₅. So (L1) holds for 0 ≤ τ ≤ c/(2M).

With ε(a) = min(φ₁, arctan(c/(2M))) Lemma L gives (b). ∎

The bound for φ₁ is crude. Numerically (§6) the tilt φ₁ = 5° works at τ = 0 for all a from 0.01 to 0.99, and arg x up to 5° works for a between about 0.1 and 0.75; near a = 0 and 1 the admissible arg x shrinks with the volume.

---

## 4. F(z) at even roots of unity

**Theorem 14.** Let X ≥ 4 be even, x = 1/X, b = X^{−1/2}, q₀ = e^{2πi/X}, 0 < a < 1, z = e^{2πia}. With the integrals taken with Faddeev's Φ_b:

  Σ_{n<X/2} zⁿ·q₀^{n(n+1)/2}·(−q₀;q₀)ₙ = i·P·e^{−πia}·[ J(a/b) − I((a−1)/b) + z^X·I(a/b) ].

Proof. For x = 1/X + iε the strip form of Theorem Z (v18 §1) holds: F(z) = i·P·e^{−πia}·[J(a/b) − I((a−1)/b)·𝒲(T) + p·T·I(a/b)·𝒲(p²T)], T = p·e^{2πia/x}. Let ε → 0⁺.

- *Left side.* Lemma R1 of v20 with the extra factor zⁿ of modulus 1: the proof bounds absolute values only. So F(z; q) tends to the finite sum.
- *Dual series.* p → 1 in a Stolz angle (v20 §1), and |T| = |p|^{1−2a}. The terms of 𝒲(T) are (−T)ᵏp^{k²+k}(p;p²)ₖ with |Tᵏp^{k²+k}| = |p|^{k² + (2−2a)k} ≤ |p|^{k²}, and those of 𝒲(p²T) have |p|^{k² + (4−2a)k}. The proof of Lemma R2 of v20 uses only |p|^{j²} and the bounds for (p;p²)ⱼ, so 𝒲(T) → 1 and 𝒲(p²T) → 1. Further p·T → e^{2πiaX} = z^X.
- *Integrals.* Lemma R3 of v20: the proof covers every u; the kernel e^{2πas}/(1 + e^{2πs}) is bounded on L for 0 ≤ a ≤ 1. ∎

For a = ½ this is Theorem 9 of v20.

**What follows, and what does not.** By Theorem 13 (a), −i·P·e^{−πia}·I((a−1)/b) is the saddle-point expression of v18 §2.2 up to 1 + O(1/X), of size X^{1/2}·e^{Vol(z)·X/2π}; by Theorem 7 of v19 the term with I(a/b) is of size X^{1/2}·e^{−Vol(z)·X/2π}. So

  F(z; e^{2πi/X}) = (1 + W)·W^{1/2}·(π/(t(1 + 2W)))^{1/2}·e^{(V(w₀) + π²/12)/t}·(1 + O(1/X)) + i·P·e^{−πia}·J(a/b) + (exponentially small).

The volume rule for even X is therefore equivalent to P·J(a/b) = o(X^{1/2}·e^{Vol(z)·X/2π}) at real x. That bound is not proved for a ≠ ½: it needs the descent inequality of v15 for the kernel e^{2πas}/(1 + e^{2πs}). Numerically |P·J(a/b)| is between 0.33 and 0.52 in all cases of §6, close to |1/(1 − 2z)|.

---

## 5. Why the horizontal ray does not work (verified)

A sketch from another session proposed, for the left part, the ray w = w₀ − s and tested Re(V(w₀ − s) − V(w₀)) ≤ −c·min(s², s). That inequality holds (reproduced), but it is not the relevant one:

- Near the real x-axis the size of e^{V/t} is e^{−Im V/|t|}; Re V enters only with the small factor τ.
- Along the horizontal line through the growing saddle, d/dr Im V = arg(1 + W) − arg(1 + e^{r−r₀}W), which is positive for r < r₀ and negative for r > r₀: Im V has a **maximum** at w₀ on that line. So the integrand grows away from the saddle in both horizontal directions. Measured along the proposed ray: −Im V rises above its saddle value by +1.2 at s = 2 and +7.4 at s = 8 for a = ½, and by +16.4 at s = 8 for a = 0.2, linearly without bound. At real x the integral over that ray does not converge.
- The saddle finder in the sketch (Newton from 2πi/3) converges for a = 0.1 to Im w = 1.86π with W in the lower half-plane: the other root on another sheet.

The horizontal line is the descent path of the **decaying** saddle (Theorem 7 of v19); for the growing one it is the vertical line. The computer-assisted approach suggested there is sound in itself; with the contour of §2 it is not needed at real x.

---

## 6. Numerical checks (verified; `m1_contour.py`, `m2_contour.py`, `n1_zroots.py`, `run_log.txt`)

| check | result |
|---|---|
| closed forms on C₁, C₂, C₃, C₅ against direct evaluation (a = 0.1) | 10⁻¹⁵ |
| monotonicity on C₄, a from 0.02 to 0.98 | holds |
| (L1) on C at τ = 0: sup of (H − H(w₀))/min(d², d) per piece, a = 0.01, 0.02, 0.05, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.98, 0.99 | negative on every piece for every a; worst −0.0036 (a = 0.99), −0.0095 (a = 0.01), −0.24 (a = 0.5) |
| the same on C′ with φ₁ = 5°, τ = 0 | negative for every a; worst −0.0036 (a = 0.99), −0.0146 (a = 0.01) |
| C′ with φ₁ = 5° at arg x = 3° and 5° | 3°: holds for 0.05 ≤ a ≤ 0.75, fails at a ≤ 0.02 and a ≥ 0.9; 5°: holds for 0.1 ≤ a ≤ 0.75 |
| first attempt (one straight ray from iπ at 30°, 60°, 75°, 85°; `m1_contour.py`) | works for a ≥ 0.3 at every angle; fails at a = 0.02 for every angle (the dilogarithm oscillates above Im w = 2π) |
| Theorem 14: finite sum against the three integrals, z = i, ω, −ω, e^{2πi·0.3}, e^{2πi·0.7}, −1; X = 8, 16, 24, 40 (double precision; the growing integral over C′) | relative difference 2·10⁻⁹ … 6·10⁻¹³; \|P·J(a/b)\| between 0.33 and 0.52 |
| Theorem 13 against the exact integral (v18 §2.3, `psi_star/v18/b2_dominant.py`) | ratio − 1 = O(t) at arg x = 1° for z = −1, i, ω, −ω |

---

## 7. Open

1. A bound for J(a/b) at real x for a ≠ ½ (the descent inequality with the kernel e^{2πas}/(1 + e^{2πs})). With it, the volume rule for F(z) at even roots of unity is a theorem, and the trivial sector for general z can be treated as in v17.
2. Odd X for z ≠ −1: by the census (README item 10) F(z) behaves differently there (it grows when z^X is close to 1); not looked at.
3. An explicit ε(a), and a wider sector near a = 0 and 1 (rotate C₄).
4. The remaining items of v19 §8, v20 §5, v21 §5.

---

## 8. Scripts (`psi_star/v22/`)

- `m1_contour.py`: the first attempt (straight ray from iπ) and the closed forms on the top edge, the vertical line and the real axis (4 minutes).
- `m2_contour.py`: the contour of Theorem 13, on the axis and tilted, at arg x = 0°, 3°, 5° (6 minutes).
- `n1_zroots.py`: Theorem 14 (4 minutes).
- `fad.py` (from v16), `run_all.sh`, `run_log.txt`, `README.md`.
