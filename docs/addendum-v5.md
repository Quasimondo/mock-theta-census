# Mock theta search: addendum v5 (the laws inside the unit disk; the third sector)

Session of 26 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). This carries out the step proposed in `addendum-v4.md` §6: move off the real line. Scripts are in `psi_star/v5/`; `run_log.txt` there is their complete output.

**Independent check** (Claude Opus, same day, `psi_star/v5/check_*.py`): the integer coefficients of M in §5 are exactly those of 2(q̃W₃F_a − W·F_b) (exact series); ψ(q) = q^{−1/8}√(π/2τ)θ₄(q′) holds to 10⁻⁴¹ at complex Y; the perturbative series Φ_W, recomputed in exact rational arithmetic, matches v4 §3. The converse law of §2 for W, with both first-correction constants, was tested inside the disk on two horocycles (b = 1, 0.6) at |Y| ≈ 240 … 1,930: the relative residual falls by a factor 4.0 each time |Y| doubles (1.8·10⁻⁵ → 2.9·10⁻⁷), i.e. like τ², as it must if the leading terms and both first corrections are right. The forward law and the third sector were not rechecked (too costly at this scale).

Everything in this addendum is **verified** (high-precision numerics, no proof) unless marked **inference**. Notation as in v3/v4: ψ₁\*, ψ₀\*, f₁\*, f₀\*, R₁, R₀, θ₄ from addendum-v2; W, W₃, P₊, P₋ and ψ(q) from v4; φ = (1+√5)/2. In addition, the one-sided halves of (V), (VI):

  F_a(q) = Σ_{n≥1} q^{2n²−n}/(q;q²)ₙ,  F_b(q) = Σ_{n≥1} q^{2n²−3n+1}/(q;q²)ₙ,

so that W + F_a = ψ·P₊ and q·W₃ + F_b = ψ·P₋.

---

## 0. Summary

1. **The third sector is found** (§5). In the forward law the decaying sector e^{−S/t} carries the companion: its coefficient function is **q̃·W₃(q̃)**, determined up to adding the growing-sector function W. The measured quantity is a Wronskian-type series with integer coefficients, matched on 32 points to 5·10⁻¹⁸.
2. **The forward law is now complete to three sectors** (§4, §5): trivial sector with coefficient exactly 1, growing sector with W, decaying sector with q̃W₃. The same holds for f\* with (2, −2F_a, −2F_b) and for the bilateral part θ₄R with (0, 2ψP₊, 2ψP₋).
3. **All three laws hold inside the unit disk, as statements about holomorphic functions** (§2–§4). Precision went from 10⁻⁹…10⁻¹¹ (v4) to about 10⁻²⁹. Inside the disk the laws distinguish ψ\* from ½f\*, which the real line cannot: W pairs with ψ\*, F_a with f\*, and the bilateral parts with each other.
4. **First corrections** (§2–§5): in every sector the next order is the same function of the dual nome times a new constant. All constants were identified, in ℚ(√5) for the bounded laws.
5. **A by-product** (§6): the sign-twisted Rogers–Ramanujan sums R₀, R₁ and the half-shifted ones P₊, P₋ are dual to each other in both directions, with explicit matrices.

No proofs yet. The method that made this possible is in §1.

---

## 1. Method

Take the limit along a **horocycle**: X = X₀ + (period)·k − i·b with b > 0 fixed and k → ∞. Then q = e^{2πi/X} is inside the disk and tends to 1, while the dual nome q̃ (a function of e^{−πiX}) stays fixed with |q̃| = e^{−πb} < 1. All series converge on both sides. Three things follow.

- **No half-integer powers.** The trivial sector is subtracted exactly, using the perturbative series at optimal truncation (error e^{−c·K}, far below everything else). What remains has corrections in integer powers of t, so polynomial extrapolation in t reaches 10⁻²⁵ with ten values.
- **The dual side is a free complex variable.** Fits are on functions of q̃ inside the disk, sampled at 8–32 points on two circles, instead of on values at roots of unity.
- **The growing sector can be cancelled exactly** (§5), because two different functions on the q side carry the same divergent series.

Evaluation: `horo.py`. The recurrences run in midpoint arithmetic at up to 1,900 bits, since ball radii blow up over 10⁵ complex multiplications. The values are therefore not rigorous enclosures; stability was checked by rerunning at higher precision, and identities (I), (II), (V), (VI) hold at the evaluation points to 10⁻⁴¹ or better.

---

## 2. Converse law inside the disk (`e1_limits.py`, `e1_fit.py`)

Y = Y₀ + 120k − ib, q = e^{2πi/Y}, τ = −2πi/Y, q′ = e^{−πiY}. Write u₀ = q′^{−1/60}, u₁ = q′^{11/60}, c = √(2π/√5). Φ_W, Φ_{W₃} are the perturbative series of W, W₃ at q = e^{−τ} (exact rational coefficients, `pert.py`). 12 points (Y₀, b), k = 10…46, extrapolation error ≤ 1.5·10⁻²⁵.

**One-sided sums** (relative residual 2·10⁻³⁰ and 4·10⁻³⁰; coefficients of f₀\*, f₁\* in the fit are zero to 10⁻²⁸):

  W(q) = Φ_W(τ) + τ^{−1/2}·c·[ −u₀ψ₀\*(q′)·(1 + (20+3√5)/120·τ + …) + i·u₁ψ₁\*(q′)·(1 + (20−3√5)/120·τ + …) ]

  W₃(q) = Φ_{W₃}(τ) + τ^{−1/2}·c·[ −φ·u₀ψ₀\*(q′)·(1 + (152−9√5)/120·τ + …) − (i/φ)·u₁ψ₁\*(q′)·(1 + (152+9√5)/120·τ + …) ]

**Bilateral parts** (relative residual 1.4·10⁻³⁰ and 2·10⁻³⁰):

  P₊(q) = 5^{−1/4}·[ u₀R₀(q′)·(1 + (5+3√5)/120·τ + …) + i·u₁R₁(q′)·(1 + (5−3√5)/120·τ + …) ]

  P₋(q) = 5^{−1/4}·[ φ·u₀R₀(q′)·(1 + (17−9√5)/120·τ + …) − (i/φ)·u₁R₁(q′)·(1 + (17+9√5)/120·τ + …) ]

**The other halves** follow from these with ψ(q) = q^{−1/8}·√(π/2τ)·θ₄(q′), which is exact, and (I), (II):

  F_a(q) = −Φ_W(τ) + τ^{−1/2}·(c/2)·[ u₀f₀\*(q′) − i·u₁f₁\*(q′) ]·(1 + O(τ)).

Remarks.

- The leading constants agree with c, cφ, c/φ and 5^{−1/4}, φ·5^{−1/4}, 5^{−1/4}/φ to all 19–20 digits printed. The first-correction constants agree with the closed forms to 10⁻²² (`constants.py`); they were guessed from the two sectors being Galois conjugates.
- The fits of the first correction on the same functions leave residuals 10⁻²⁵ to 10⁻²⁶. So the law has the Garoufalidis–Zagier shape Σ_σ (power series in τ)·(function of q′), at least to this order.
- This is v4 §4 with three improvements: 29 digits instead of 11, the dual function is ψ\* and not ½f\*, and the trivial sector is exactly Φ_W times 1.

---

## 3. Mirror law inside the disk (`e3_mirror.py`)

Z = Z₀ + 240k − ib, q = −e^{2πi/Z} = −e^{−t}, t = −2πi/Z, q̂ = −e^{−πiZ/2}, û₀ = e^{iπZ/120}, û₁ = e^{−11iπZ/120}, κ = √(π/√5). Φ₁, Φ₀ are the perturbative series of ψ₁\*, ψ₀\* at q = −1. 8 points, k = 13…31, extrapolation error ≤ 5·10⁻²⁵.

  ψ₁\*(q) = Φ₁(t) + t^{−1/2}·κ·[ −(1/φ)·û₀ψ₀\*(q̂)·(1 + (23+6√5)/60·t + …) + iφ·û₁ψ₁\*(q̂)·(1 + (23−6√5)/60·t + …) ]

  ψ₀\*(q) = Φ₀(t) + t^{−1/2}·κ·[ φ·û₀ψ₀\*(q̂)·(1 + (−13+6√5)/60·t + …) − (i/φ)·û₁ψ₁\*(q̂)·(1 + (−13−6√5)/60·t + …) ]

Relative residuals 4·10⁻³¹ and 1.3·10⁻³⁰ at leading order, 2·10⁻²⁶ for the first corrections; κ/φ and κφ reproduced to 20 digits; coefficients of f₀\*, f₁\* zero to 10⁻²⁸. Again the dual function is ψ\* itself, not ½f\*.

---

## 4. Forward law inside the disk: the growing sector (`e2_forward.py`, `e2_fit.py`)

X = X₀ + 2m − ib, q = e^{2πi/X}, t = −2πi/X, q̃ = e^{−πiX}. S̄ = π²/3 + i·Cl₂(π/3), C̄ = √(π/√3)·e^{5πi/12} as in v3. 12 points, K = 800…3,680, extrapolation error ≤ 2.4·10⁻²⁵. Fits on (W, W₃, F_a, F_b, 1) at q̃:

  ψ₁\*(q) ≃ C̄·t^{−1/2}·e^{−S̄/t}·**W(q̃)**·(1 + (1/4 + i/(18√3))·t + …)

  ψ₀\*(q) ≃ e^{2πi/3}·C̄·t^{−1/2}·e^{−S̄/t}·**W(q̃)**·(1 + (−1/12 + i/(18√3))·t + …)

  f₁\*(q) ≃ C̄·t^{−1/2}·e^{−S̄/t}·**(−2F_a(q̃))**·(same series as ψ₁\*)

  f₀\*(q) ≃ e^{2πi/3}·C̄·t^{−1/2}·e^{−S̄/t}·**(−2F_a(q̃))**·(same series as ψ₀\*)

Relative residual 5·10⁻²⁹; the coefficients come out as 1.000… and −2.000…, all others zero to 10⁻²⁶. First corrections: residual 10⁻²⁴.

Consequences.

- By (II), θ₄(q)R₁(q) = 2ψ₁\* − f₁\* ≃ C̄t^{−1/2}e^{−S̄/t}·2ψ(q̃)P₊(q̃). The growing sector is the triple (W, −2F_a, 2ψP₊) for (ψ\*, f\*, θ₄R), tied together by identity (V).
- W₃ and F_b are absent from the growing sector, for ψ₁\* and ψ₀\* alike. This confirms v4 §6 at 29 digits.
- On the real axis W_geo = −F_a at even-order roots (v4 §2), which is why the real-line law could not tell ψ\* ↔ W from f\* ↔ −2F_a.

---

## 5. The third sector (`e4_decay.py`, `e4_run.py`, `e4_fit.py`, `wronskian.py`)

**Idea.** ψ₁\* and f₁\* carry the same divergent series in the growing sector, with coefficient functions W(q̃) and −2F_a(q̃). Since q̃ is fixed along the horocycle, the combination

  D₁ = 2F_a(q̃)·ψ₁\*(q) + W(q̃)·f₁\*(q) − 2ψ(q̃)P₊(q̃)·Φ(t)

has no growing sector at all. Φ is the perturbative series of ψ₁\* at q → 1 (v3; that of f₁\* is 2Φ), truncated at 140 terms. If the trivial sector is exactly Φ times 1, then D₁ is the decaying sector alone.

**Result.**

- D₁ is tiny: at K = 1,000, |ψ₁\*| = 1.8·10⁷¹ and |D₁| = 2.4·10⁻⁶⁹. So the growing series cancel to all orders and the trivial sector has coefficient exactly 1 (to 140 orders of magnitude).
- D₁·√t·e^{S/t}/C converges, with S = π²/3 − i·Cl₂(π/3) and C = √(π/√3)·e^{−5πi/12}, the complex conjugates of S̄, C̄. Call the limit M(q̃). The same construction for ψ₀\*, f₀\* with the factor e^{−2πi/3} gives the same M.
- Sampling M on the circle |q̃| = e^{−2π/5} at 20 points and taking Fourier coefficients gives a power series with **integer coefficients**:
  M(q̃) = −2 − 2q̃ + 2q̃² − 2q̃³ − 2q̃⁴ + 0·q̃⁵ − 8q̃⁶ − 2q̃⁷ + 0·q̃⁸ − 4q̃⁹ − 4q̃¹⁰ − 4q̃¹¹ + 0·q̃¹² − 6q̃¹³ …
- This is exactly **M = 2·(q̃W₃F_a − W·F_b)(q̃)**. On 32 points (two circles) the largest relative deviation is 5·10⁻¹⁸, against an extrapolation error of 2·10⁻¹⁶.
- First correction: M·(1/4 − i/(18√3)) for ψ₁\* and M·(−1/12 − i/(18√3)) for ψ₀\*, the conjugates of §4 (12 digits, spread 6·10⁻¹⁷).

**Reading (inference).** If the decaying sector of ψ₁\* has coefficient J_ψ and that of f₁\* has J_f, then M = 2F_a·J_ψ + W·J_f. The solution of 2F_aJ_ψ + WJ_f = 2q̃W₃F_a − 2WF_b that mirrors the growing sector is

  **J_ψ = q̃·W₃(q̃),  J_f = −2F_b(q̃),**  and for the bilateral part 2J_ψ − J_f = 2ψ(q̃)P₋(q̃).

It is determined only up to adding g(q̃)·(W, −2F_a), the growing-sector pair times a function g; the natural case is g constant. That freedom is expected: the real axis is the Stokes line between S̄ and S, so "the decaying sector" is defined only relative to a resummation of the growing one.

**The forward law, three sectors.** For X → ∞ with Im X = −b:

  ψ₁\*(q) = Φ(t)·1 + C̄t^{−1/2}e^{−S̄/t}·A(t)·W(q̃) + C·t^{−1/2}e^{−S/t}·Ā(t)·[q̃W₃(q̃) + λW(q̃)],

with A(t) = 1 + (1/4 + i/(18√3))t + …, Ā its conjugate series, λ undetermined. For ψ₀\* the exponential terms get e^{±2πi/3} and A = 1 + (−1/12 + i/(18√3))t + …. Replacing (1, W, q̃W₃) by (2, −2F_a, −2F_b) gives f\*, and by (0, 2ψP₊, 2ψP₋) gives θ₄R₁.

So the answer to the question behind this follow-up is yes: **W's companion is the coefficient of the third sector.** In terms of the two rank-3 modules of v4 §1, the three sectors (0, S̄, S) of ⟨1, ψ₁\*, ψ₀\*⟩ at q → 1 pair with the basis (1, W, q̃W₃) of the W module.

---

## 6. By-product: R₀, R₁ and P₊, P₋ are dual

Dividing the θ₄R statements by the exact modular factor θ₄(q) = 2√(π/t)·q̃^{1/8}·ψ(q̃) removes ψ and leaves a closed system between the four Rogers–Ramanujan-type sums (inference from §2, §4, §5; the P± → R direction is measured directly):

- q → 1 along the horocycle, growing and decaying sector: R₁(q) ≃ (C̄/√π)·e^{−(π²/12 + i·Cl₂(π/3))/t}·P₊(q̃) + (C/√π)·e^{−(π²/12 − i·Cl₂(π/3))/t}·[P₋(q̃) + λP₊(q̃)], and the same for −R₀ with the factors e^{±2πi/3}.
- In the other direction (§2): (P₊, P₋)(q) → 5^{−1/4}·[[1, i], [φ, −i/φ]]·(u₀R₀, u₁R₁)(q′).

So "Σ zⁿq^{n²}/(q)ₙ at z = −1, −q" and "the same function at Q = q², z = q^{±1}" are exchanged by the Fricke involution. Both are solutions of the homogeneous q-difference equations, which is why they have no trivial sector.

---

## 7. Open

- **Proofs.** With (I), (II), (V), (VI) the laws reduce to statements about the bilateral parts (the system in §6) plus statements about the one-sided halves. The bilateral system looks provable by standard means (the R and P sums have Nahm-type integral representations).
- **λ**, and whether g in §5 is constant. This needs a Borel resummation of the growing series next to its Stokes line, or the exact method of the previous bullet.
- **Higher orders.** Only the first correction was checked in each sector. In §5 the full series Ā(t) should be the conjugate of A(t).
- **W, W₃ near ξ = −1.** Not tested. By duality one expects growth and a law that sends the W family to itself (inference).
- **Mirror and converse laws:** a third (constant) sector is the trivial one; there is no hidden sector left in those two directions at the measured precision.
- F₀\*, F₁\*, T₀\*, F(i), F(−ω): not touched.

---

## 8. Scripts (`psi_star/v5/`)

| script | section | what it does |
|---|---|---|
| `horo.py` | §1 | evaluators inside the disk; self-test with (I), (II), (V), (VI) |
| `pert.py` | §2, §3 | exact perturbative series of W, W₃ at q → 1 and of ψ₁\*, ψ₀\* at q → −1 |
| `e1_limits.py`, `e1_fit.py` | §2 | converse law: limits, fits, first corrections |
| `e3_mirror.py` | §3 | mirror law: limits, fits, first corrections |
| `e2_forward.py`, `e2_fit.py` | §4 | forward law, growing sector |
| `e4_decay.py`, `e4_run.py`, `e4_fit.py` | §5 | decaying sector: cancellation, Taylor coefficients, comparison with 2(q̃W₃F_a − WF_b) |
| `wronskian.py` | §5 | exact q-series of q·W₃F_a − W·F_b |
| `constants.py` | §2, §3 | closed forms of the first-correction constants |

Limit files (`*.json`) are included; each limit script takes one to two minutes. Needs mpmath and python-flint.
