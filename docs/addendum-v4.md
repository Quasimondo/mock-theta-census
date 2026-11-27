# Mock theta search: addendum v4 (the W family and two more modularity laws)

Session of 26 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). This follows up §7 of `addendum-v3.md`, item 2: "find the partner of W and its identities". Scripts are in `psi_star/v4/`; `run_log.txt` there is their complete output.

**Independent check** (Claude Opus, same day, own code in `psi_star/v4/check_*.py`, `conv_trend.py`): (V), (VI), G(−q⁵) = (W − 1 + q²W₃)/q³ and G(−1/q) = 1 − W + qW₃ hold exactly to q³²⁹. The Stokes constant of §3 is confirmed: ψ(q)P₊(q)√t·e^{−π²/(30t)} → √(2π/√5)/2 to 10⁻¹⁸, with first correction (20+3√5)/120 to 15 digits. The converse law of §4 (W and W₃, four classes) and the mirror law of §5 (ψ₁\*, ψ₀\*, one class, ρ = 0; the run for the other classes did not finish) agree with the stated closed forms to 10⁻⁴ … 5·10⁻⁵ with Y, Z ≤ 2,000 … 3,900; for ρ = 0 the gap falls steadily to 1.5·10⁻⁶ at Y ≈ 5,800 with no plateau. So the closed forms are confirmed only to the precision those runs allow, which is below Fable's own 2·10⁻¹¹ and 10⁻⁹.

Status labels as in v3: **proved** = proof given here and machine-checked; **verified** = exact or high-precision computation, no proof; **inference** = my reading of the evidence.

Notation. ψ₁\*, ψ₀\*, F(z), W as in v3. φ = (1+√5)/2. ψ(q) = Σ_{n≥0} q^{n(n+1)/2} = (q²;q²)∞/(q;q²)∞ (Ramanujan's ψ, no star).

---

## 0. Summary

1. **W has one companion, W₃, and the pair has identities of type (I)/(II)** (§1, proved, checked to q³¹⁹). With G(z) = Σ zʲq^{j²}(q;q²)ⱼ, W = G(−q) and W₃ = G(−q³). Completing each to a bilateral sum gives ψ(q) times a Rogers–Ramanujan-type sum with a half-integer shift.
2. **Even-root rule for W** (§2, verified to 5·10⁻³⁷ on all even orders ≤ 40). W diverges at roots of even order, and its geometric-series value there equals minus the "other half" of the bilateral sum. This is the exact analogue of the odd-root rule for ψ\*.
3. **Resurgence of W at ξ → 1** (§3, verified, 29 digits). Borel singularity at −π²/30, Stokes constant √(2π/√5)/2π, first correction (20+3√5)/120.
4. **Converse law** (§4, verified to 2·10⁻¹¹). Near x = 0 the pair (W, W₃) goes back to (ψ₀\*, ψ₁\*) at the dual root, through a closed-form 2×2 matrix with golden-ratio entries and the Rogers–Ramanujan phases e^{2πi·Y/120}, e^{−2πi·11Y/120}. No exponential growth in this direction.
5. **Mirror law at q = −1** (§5, verified to 10⁻⁹). Near x = ½ the pair (ψ₁\*, ψ₀\*) goes to itself at the image point under x ↦ (x−1)/(2x−1), again through a golden-ratio matrix. This is the exact form of what the v1/v2 fits at q = −1 were approximating.
6. **The third sector** (§6, inference). Both sides are three-dimensional: (1, ψ₁\*, ψ₀\*) and (1, W, W₃), each closed under its q-difference equation. The two bounded laws fill a 2×2 block; the constant 1 is the third direction. In the growing direction of v3 the companion W₃ does **not** appear, and the decaying sector is still not measured. I could not complete the 3×3 matrix.
7. **Precision caveat.** Items 4 and 5 are limits obtained by extrapolation, and they are below the brief's 15-digit rule: 2·10⁻¹¹ and 10⁻⁹. The closed forms were found by fitting, then compared.

---

## 1. The family G(z), the companion W₃, identities (V) and (VI)

Define, for |q| < 1,

  G(z) = Σ_{j≥0} zʲ q^{j²} (q;q²)ⱼ,  W = G(−q),  W₃ = G(−q³).

**q-difference equation (proved, checked to q³¹⁹).** From (q;q²)ⱼ = (q;q²)ⱼ₋₁(1 − q^{2j−1}):

  G(z) = 1 + zq·G(zq²) − zq²·G(zq⁴).

So every G(−q^{2k+1}), k ∈ ℤ, is a combination of 1, W, W₃ with coefficients in ℤ[q, q⁻¹]. Examples: G(−q⁵) = (W − 1 + q²W₃)/q³ and G(−1/q) = 1 − W + qW₃. The W family is therefore the three-dimensional module ⟨1, W, W₃⟩, and **W₃ is the only companion**.

The same holds on the ψ\* side: F(z) = 1 + zq·F(zq) + zq²·F(zq²), so ψ₁\* = 1 − ψ₀\* − q²F(−q²), and every qʲF(−qʲ) lies in ⟨1, ψ₁\*, ψ₀\*⟩ (`f_relation.py`).

**Theorem (proved, checked to q³¹⁹; `w_identities.py`).** Let

  P₊(q) = Σ_{k≥0} q^{2k²+k}/(q²;q²)ₖ,  P₋(q) = Σ_{k≥0} q^{2k²−k}/(q²;q²)ₖ.

Then

  **(V)**  W(q) + Σ_{n≥1} q^{2n²−n}/(q;q²)ₙ = ψ(q)·P₊(q),

  **(VI)**  q·W₃(q) + Σ_{n≥1} q^{2n²−3n+1}/(q;q²)ₙ = ψ(q)·P₋(q).

P± are the two-variable Rogers–Ramanujan function Σ zᵏQ^{k²}/(Q;Q)ₖ at Q = q², z = q^{±1}. They are not functions of q² alone, and they are not modular (see §3: their asymptotic series do not terminate).

*Proof.* Extend the sum for W to negative j using (q;q²)₋ₙ = (−1)ⁿq^{n²}/(q;q²)ₙ. The term j = −n of Σ(−1)ʲq^{j²+j}(q;q²)ⱼ becomes q^{2n²−n}/(q;q²)ₙ, so the left side of (V) is the bilateral sum Σ_{j∈ℤ}(−1)ʲq^{j²+j}(q;q²)ⱼ. Write (q;q²)ⱼ = (q;q²)∞/(q^{2j+1};q²)∞ and expand, with Q = q²,

  1/(yQ;Q)∞ = (1/(Q;Q)∞)·Σ_{k≥0} (−1)ᵏQ^{k(k+1)/2} / ((Q;Q)ₖ(1 − yQ^{k+1})),  y = q^{2j−1}.

This leaves the inner sums Aₖ = Σ_{j∈ℤ}(−1)ʲq^{j²+j}/(1 − q^{2j+2k+1}). Shifting j by k and using x^{−k}/(1−x) = 1/(1−x) + x⁻¹ + … + x^{−k} with x = q^{2j+1}, every extra term is a multiple of Σ_{j∈ℤ}(−1)ʲq^{j²+j} = 0 (the terms j and −1−j cancel). Hence Aₖ = (−1)ᵏq^{k²}·A₀. Finally A₀ = Σ_{j∈ℤ}(−1)ʲq^{j(j+1)}/(1 − q^{2j+1}) = (q²;q²)∞²/(q;q²)∞² = ψ(q)², by the partial-fraction expansion of the reciprocal theta function. Collecting factors gives (V). For (VI) the same steps with q^{j²+3j} = q^{j²+j}·x/q give Bₖ = (−1)ᵏq^{k²−2k−1}A₀. ∎

Compare (I), (II) of addendum-v2: there the prefactor is θ₄ and the sums are sign-twisted Rogers–Ramanujan sums; here the prefactor is ψ and the sums carry a half-integer shift. The one-sided halves Σ_{n≥1} q^{2n²−n}/(q;q²)ₙ and Σ_{n≥1} q^{2n²−3n+1}/(q;q²)ₙ play the role of f₁\*, f₀\*. I did not find them among the classical mock theta functions (short look only).

---

## 2. Even-root rule for W (verified)

At a root ξ of odd order W terminates. At even order 2s it does not: (ξ;ξ²)ₛ = 2, so the terms grow like 2^{j/s} and block s+j is −2 times block j. The geometric-series value is W_geo(ξ) = ⅓·Σ_{j<s}(−1)ʲξ^{j²+j}(ξ;ξ²)ⱼ. This is the value used in the law of v3 §5.

**Statement.** For every root ξ of even order 2s,

  W_geo(ξ) = −F_a(ξ),  F_a(ξ) = Σ_{n≥1} ξ^{2n²−n}/(ξ;ξ²)ₙ = (2/3)·Σ_{n=1}^{s} ξ^{2n²−n}/(ξ;ξ²)ₙ,

and W(q) → W_geo(ξ) as q → ξ radially.

- The finite identity ⅓·(first s terms of W) = −(2/3)·(first s terms of F_a): **verified** for all primitive roots of even order ≤ 40, largest error 5·10⁻³⁷ (`w_checks.py`, part a). F_a converges at ξ as a geometric series with ratio −½.
- The radial limit follows from (V) once ψ(q)·P₊(q) → 0 at even-order roots. **Verified** numerically (for example at ξ = −1: |ψP₊| = 1.5·10⁻⁵ at t = 1/8 and 2·10⁻²² at t = 1/32, while (V) itself holds to 10⁻⁴⁶ … 10⁻⁶⁴). Not proved here. A proof needs the decay of ψ at even-order roots, which is standard, and a bound on P₊ there, which I did not write down.

---

## 3. Resurgence of W at ξ → 1 (verified; `w_checks.py`, parts b and c)

W(e^{−t}) = Σ aₖtᵏ with k!·aₖ = 1, −1, 11, −253, 10859, −745981, 75045131, −10399927693, … (integers; checked for k ≤ 30).

Large order: aₖ = (−1)ᵏ Σⱼ cⱼ Γ(k+½−j) A^{−(k+½−j)} with

- **A = π²/30** (ratio test with Richardson: 0.3289868144 against 0.3289868134; the fit with A fixed leaves a residual of 4·10⁻³⁷);
- **π·c₀² = 1/(2√5)**, to 4·10⁻²⁹. Equivalently c₀ = √(2π/√5)/2π = 0.26678878226364336;
- first coefficient of the attached series **s₁ = (20+3√5)/120** = 0.222568366104161 (PSLQ on 18 digits).

The attached series is that of ψ·P₊: directly, ψ(q)P₊(q)·√t·e^{−π²/(30t)} → π·c₀ at q = e^{−t} (relative error 5·10⁻²⁰), and its first correction equals s₁ to 15 digits. So the singularity of W at −π²/30 is the bilateral side of (V), in the same way that the singularities of ψ\* in v2/v3 are the bilateral sides of (I), (II). Since s₁ ≠ 0 and the series does not stop, ψ·P₊ is not a modular form.

---

## 4. The converse law (verified to 2·10⁻¹¹)

Let Y → ∞ through rationals with fixed fractional part, τ = −2πi/Y, q′ = e^{−πiY}, c = √(2π/√5) = 1.676283357. Then

  W(e^{2πi/Y})·√τ → c·[ −e^{iπY/60}·ψ₀\*(q′) + i·e^{−11iπY/60}·ψ₁\*(q′) ],

  W₃(e^{2πi/Y})·√τ → c·[ −φ·e^{iπY/60}·ψ₀\*(q′) − (i/φ)·e^{−11iπY/60}·ψ₁\*(q′) ].

Values at roots are the exact ones: v3 §5.1 for ψ\*, §2 above for W, W₃.

How it was found (`wflint.py`, `conv_limits.py`, `conv_fit.py`):

- Limits were computed for 8 fractional parts ρ ∈ {0, ½, ⅓, ⅔, ¼, ¾, ⅕, ⅖} and Y = ρ + M₀ + 120k, M₀ = 0…9, k up to 181 (Y up to about 21,700), by Richardson extrapolation in Y^{−1/2}. Estimated extrapolation error 4·10⁻¹¹.
- Step 1: in every class the M₀-dependence is exactly two frequencies, e^{iπY/60} and e^{−11iπY/60} (residual 5·10⁻¹²).
- Step 2: the two amplitudes, as functions of ρ, were fitted on (1, ψ₁\*, ψ₀\*) at several candidate dual points. Only e^{2πi(−Y/2)} fits (relative residual 5·10⁻¹⁴ to 10⁻¹³; the five other candidates leave ≥ 0.84). The coefficient of 1 and the two "wrong" coefficients come out as zero within the error.
- The closed forms above then deviate from the 16 measured amplitudes by at most 1.9·10⁻¹¹.

Remarks.

- The phases are Rogers–Ramanujan: e^{iπY/60} = e^{π²/(30τ)} and e^{−11iπY/60} = e^{−11π²/(30τ)}. The first is the singularity −π²/30 of §3, and at Y even the law reduces to the Stokes constant of §3 (c/2π = c₀, with ψ₀\*(1) = ⅓). That is an independent 29-digit check of the constant c in one class. The second predicts a Borel singularity of W at +11π²/30 (inference; not tested).
- The matrix is c·[[i, −1], [−i/φ, −φ]] on (ψ₁\*, ψ₀\*) after the phases, with determinant −i·c²·(φ + 1/φ) = −2πi.
- Direction of the map: x = 1/Y ↦ −Y/2 = −1/(2x). The forward law of v3 is x = 1/X ↦ −X/2. Both are the Fricke involution of Γ₀(2), in opposite directions (inference). They are not symmetric: forward, ψ\* grows like e^{0.1615·X} and only W appears; converse, nothing grows and both ψ₀\*, ψ₁\* appear.

---

## 5. The mirror law at q = −1 (verified to 10⁻⁹)

Let x = ½ + 1/Z, so q = −e^{2πi/Z}, with Z → ∞ through rationals of fixed fractional part. Let t = −2πi/Z, q̂ = −e^{−iπZ/2} = e^{2πi(½ − Z/4)}, κ = √(π/√5) = 1.185310. Then

  ψ₁\*(q)·√t → κ·[ −(1/φ)·e^{iπZ/120}·ψ₀\*(q̂) + iφ·e^{−11iπZ/120}·ψ₁\*(q̂) ],

  ψ₀\*(q)·√t → κ·[ φ·e^{iπZ/120}·ψ₀\*(q̂) − (i/φ)·e^{−11iπZ/120}·ψ₁\*(q̂) ].

How it was found (`jflint.py`, `mirror_limits.py`, `mirror_fit.py`):

- Limits for ρ ∈ {0, ½, ⅓, ⅔}, Z = ρ + M₀ + 240k, M₀ = 0…15, k up to 40 (Z up to about 9,600). Estimated extrapolation error 2·10⁻⁷ in the worst class.
- Fit on (1, U₁, U₂) with two frequencies, for 16 candidates: ψ\* or W at eight dual points. Only ψ\* at e^{2πi(½ − Z/4)} fits, with relative residual 4·10⁻¹⁰ and 1.2·10⁻⁹; every other candidate leaves ≥ 0.75. In particular **W does not appear at q = −1**.
- The four nonzero fitted coefficients are −0.7325626885, 1.917874017·i, 1.917874017, −0.7325626882·i, against κ/φ = 0.7325626885 and κφ = 1.9178740173. The other eight are zero within 2·10⁻⁷.

Remarks.

- x ↦ (x−1)/(2x−1) is in Γ₀(2). It sends ½ to ∞ (the same cusp class of Γ₀(2)) and ½ + 1/Z to ½ − Z/4. So near q = −1 the pair transforms **into itself**: ψ\* is a vector-valued quantum modular object for this element, with weight ½ and no exponential growth.
- This agrees with addendum-v2 §3: the Borel singularity −π²/60 found there at q = −1 is the first phase here, e^{iπZ/120} = e^{π²/(60t)}. The second phase predicts another singularity at +11π²/60 (inference; not tested).
- The matrix is κ·[[iφ, −1/φ], [−i/φ, φ]] on (ψ₁\*, ψ₀\*), determinant iκ²(φ² − φ⁻²) = iπ.

---

## 6. What this says about the third sector

**What is established.**

- Each side is a rank-3 object: ⟨1, ψ₁\*, ψ₀\*⟩ and ⟨1, W, W₃⟩ (§1, proved). Each has three sectors at its cusp: for ψ\* at q → 1 the actions 0, S, S̄ of v3 §3.2; for W at ξ → 1 the actions 0, −π²/30 and (predicted) +11π²/30.
- Bounded 2×2 blocks, in closed form: (W, W₃) → (ψ₀\*, ψ₁\*) under Fricke (§4), and (ψ₁\*, ψ₀\*) → itself at the cusp ½ (§5). In both, the coefficient of the constant 1 at order τ^{−1/2} is zero: the constant is the trivial sector and enters at order τ⁰, as the perturbative series itself.
- Growing direction (v3): ψ₁\* and ψ₀\* both go to W alone, with ratio e^{2πi/3}.

**What is not.** (inference from here on)

- The hope behind this follow-up was that the companion would supply the missing row of the forward law. It does not do so at leading order: the v3 fit on W alone is good to 6·10⁻²¹, so W₃ has no room in the growing sector. If W₃ enters the forward law, it is in the decaying sector e^{−S/t} or in the O(1) part.
- Those parts cannot be read off on the real line: they sit under a divergent series times e^{0.16·X}. The first corrections s̄₁ of ψ₁\* and ψ₀\* differ, so no fixed combination of the two cancels the growth.
- The natural picture is a 3×3 matrix between (1, ψ₁\*, ψ₀\*) at x and (1, W, W₃) at −1/(2x), with one row per sector. Known so far: the full bounded block in one direction, and one row (the growing one) in the other.

**How to get the rest.** Move off the real line. For x in the upper half-plane both ψ\*(e^{2πix}) and W(e^{2πi(−1/2x)}) converge, and identities (I), (II), (V), (VI) split each into a bilateral part and a one-sided part. The test is whether (bilateral part of ψ\*) and (bilateral part of the W family) are exchanged by the 2×2 matrix of §4 up to a function that extends analytically across the real axis. That needs no extrapolation and would give 30+ digits, including the sector that is invisible on the real line. This is the method of Garoufalidis–Zagier for the knot case.

---

## 7. Negative and untested

- W, W₃ at ξ → −1: not tested. By §2 and the forward law one expects growth with a Clausen-type action there, since that is where ψ\* near x = 0 lands.
- F₀\*, F₁\*, T₀\*, F(i), F(−ω): not touched in this round.
- Proofs of the laws in §4, §5 and of the radial limit in §2: not attempted.
- Literature for G(z), P±, (V), (VI): not searched beyond the standard mock theta lists.

---

## 8. Scripts (`psi_star/v4/`)

Run each with `python3 <name>` from the unpacked folder. Needs mpmath and python-flint.

| script | what it checks | time |
|---|---|---|
| `w_identities.py` | (V), (VI) and the G relations, exact to q³¹⁹ | seconds |
| `f_relation.py` | ψ₁\* = 1 − ψ₀\* − q²F(−q²), exact to q³¹⁹ | seconds |
| `w_checks.py` | §2 even-root rule; §3 series, singularity, constants | ~1 min |
| `conv_limits.py 0,1/2,1/3,2/3,1/4,3/4,1/5,2/5 long` | limits for §4 → `conv_limits_hi.json` (included) | ~25 min |
| `conv_fit.py conv_limits_hi.json` | §4 fits and closed-form comparison | seconds |
| `mirror_limits.py 0,1/2,1/3,2/3` | limits for §5 → `mirror_limits.json` (included) | ~15 min |
| `mirror_fit.py` | §5 fits and closed-form comparison | seconds |
| `qs.py`, `jval.py`, `wval.py`, `jflint.py`, `wflint.py` | helpers: exact series; values of ψ\*, W, W₃ at roots of unity | |
