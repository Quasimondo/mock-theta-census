# Mock theta search: addendum v8 (the law for F₀\*, F₁\*)

Session of 26 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). This does item 4 of `addendum-v7.md` §4: the two-sector law for F₀\*, F₁\* predicted by the state-integral structure, compared with direct numerics. Scripts and output are in `psi_star/v8/`.

**Independent check** (Claude Opus, same day, `psi_star/v8/check_VII.py`, `check_F0_law.py`): identity (VII) holds exactly to r³²⁹, the two forms of 𝒦 agree, and the quoted integer coefficients of 𝓔, 𝒥, 𝒦 are those of the stated series. The one-sided law F₀\* ≃ A(t)·𝒥(r)·(1 + (1/8 + i/(9√3))t) was computed directly on two horocycles (|r| = 0.21, 0.39) at |𝒳| ≈ 200 … 1,600: the deviation is the same at both and falls by a factor 4 per doubling (1.6·10⁻⁴ → 2.4·10⁻⁶), i.e. like t². The other three functions (G₀, qF₁\*, −G₁) were not rechecked.

Labels: **verified** (numerics or exact series check, no proof), **inference**.

Notation (v3 §4). F₀\* = Σ(−1)ⁿq^{2n²}/(q;q²)ₙ, F₁\* = Σ(−1)ⁿq^{2n²+2n}/(q;q²)ₙ₊₁, and the other halves of (III), (IV):

  G₀ = Σ_{n≥1} q^{n²}(q;q²)ₙ,  G₁ = Σ_{n≥0} q^{(n+1)²}(q;q²)ₙ,

so F₀\* + G₀ = θ₄(q²)(−q;q)∞·X(q²) and q·F₁\* − G₁ = θ₄(q²)(−q;q)∞·Y(q²).

---

## 0. Summary

1. **The prediction holds** (§1). The growing sector of the bilateral parts of F₀\*, F₁\*, computed from the four state integrals of v6/v7 with no free parameter, agrees with direct numerics: leading function, constant and first correction, on 20 points to 10⁻¹⁷.
2. **The one-sided law is found too** (§2). F₀\* and q·F₁\* go to 𝒥(r) = Σ r^{2m²+m}/(−r²;r²)ₘ, and G₀, −G₁ go to 𝒦(r) = 2Σ r^{(n+1)²}(−r²;r²)ₙ, at the dual nome r. Both have integer coefficients; matched to 10⁻¹⁷.
3. **A new bilateral identity (VII)** (§2, exact to r²⁹⁹): 𝒥 + 𝒦 = Σ_{n∈ℤ} r^{n²}(−1;r²)ₙ = ψ(r)·(r²;r⁴)∞·Σ r^{m(m+1)/2}/(r²;r²)ₘ.
4. **Pattern** (§3). The dual of the sign-twisted F₀\*, F₁\* is the *unsigned* fifth-order f/ψ family at a half-integer shift, just as the dual of the sign-twisted ψ\* is the W family. Sign on one side, half shift on the other.

---

## 1. The bilateral parts: prediction and check (`fstar.py`)

Setting: 𝒳 = 𝒳₀ + 4m − ib, q = e^{2πi/𝒳}, t = −2πi/𝒳, dual nome **r = e^{−πi𝒳/2}**, fixed along the sequence.

**Prediction.** X and Y enter (III), (IV) at argument q², so v6 §2 applies with t replaced by 2t and dual nome r. With the exact modular factors θ₄(q²) = 2√(π/2t)·r^{1/8}·ψ(r) and (−q;q)∞ = 2^{−1/2}e^{π²/(12t)+t/24}(r²;r⁴)∞, the growing sector is

  F₀\* + G₀ ≃ A(t)·𝓔(r),  q·F₁\* − G₁ ≃ e^{2πi/3}·A(t)·𝓔(r),

  A(t) = √(π/(2√3))·e^{−iπ/12}·t^{−1/2}·exp(−(π²/24 + (i/2)·Cl₂(π/3))/t),  𝓔(r) = ψ(r)·(r²;r⁴)∞·E₊′(r),

with E₊′(r) = Σ r^{m(m+1)/2}/(r²;r²)ₘ from v6. The growth is e^{Cl₂(π/3)·K/4π} = e^{0.0808·K}, half the rate of ψ\*. First corrections: (1/8 + i/(9√3))·t for the first line, by the same substitution.

**Check.** Each of F₀\*, G₀, q·F₁\*, −G₁ was computed directly on the horocycle (K = 800…2,000), divided by A(t) or e^{2πi/3}A(t), and extrapolated at 20 points on the circle |r| = e^{−2π/5}. Extrapolation error 2·10⁻¹⁶.

- Taylor coefficients of the sum (F₀\* + G₀)/A in r: 1, 2, 0, 1, 2, −1, 2, 1, 0, 1, 1, 3, −1, 1, … These are the coefficients of 𝓔(r). The same for (q·F₁\* − G₁)/(e^{2πi/3}A). So the constant, the phase e^{2πi/3} and the function are all as predicted.
- First correction of F₀\* and G₀: 0.125 + 0.06415003·i = 1/8 + i/(9√3), as predicted (12 digits, same at all 20 points). For q·F₁\* and −G₁ it is −13/24 + i/(9√3).

This is a test of v6/v7 on functions that were not used to find them.

---

## 2. The one-sided law (`fstar.py`, `fstar_check.py`)

The Taylor coefficients of the four limits separately are integers:

- F₀\*/A and q·F₁\*/(e^{2πi/3}A): 1, 0, 0, 1, 0, −1, 0, 1, 0, −1, 1, 1, −1, −1, …
- G₀/A and −G₁/(e^{2πi/3}A): 0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 2, 0, 2, …

They are

  **𝒥(r) = Σ_{m≥0} r^{2m²+m}/(−r²;r²)ₘ,  𝒦(r) = 2·Σ_{n≥0} r^{(n+1)²}(−r²;r²)ₙ = Σ_{n≥1} r^{n²}(−1;r²)ₙ.**

Largest relative deviation of the 80 measured limits from 𝒥, 𝒦: 1.2·10⁻¹⁷.

**Law (growing sector).** F₀\*(q) ≃ A(t)·𝒥(r)·(1 + (1/8 + i/(9√3))t + …),  q·F₁\*(q) ≃ e^{2πi/3}A(t)·𝒥(r)·(1 + (−13/24 + i/(9√3))t + …), and the same with 𝒦 for G₀ and −G₁.

**Identity (VII)** (exact to r²⁹⁹). 𝒥 is the n ≤ 0 half and 𝒦 the n ≥ 1 half of one bilateral sum, since (−1;r²)₋ₘ = r^{m(m+1)}/(−r²;r²)ₘ:

  Σ_{n∈ℤ} r^{n²}(−1;r²)ₙ = 𝒥(r) + 𝒦(r) = ψ(r)·(r²;r⁴)∞·Σ_{m≥0} r^{m(m+1)/2}/(r²;r²)ₘ.

I did not write a proof; the method of (V), (VI) should apply.

Not measured: the decaying sector (prediction for the bilateral part: the same with E₋′ in place of E₊′ and the conjugate exponential) and the trivial sector (G₀ has an ordinary perturbative series at q → 1; F₀\* should carry minus it).

---

## 3. The pattern (inference)

Write s = r². Then 𝒥 = Σ s^{m²+m/2}/(−s;s)ₘ and 𝒦 = 2Σ s^{(n+1)²/2}(−s;s)ₙ. Ramanujan's unsigned fifth-order functions are f₀ = Σ s^{n²}/(−s)ₙ, f₁ = Σ s^{n²+n}/(−s)ₙ, ψ₀ = Σ s^{(n+1)(n+2)/2}(−s)ₙ, ψ₁ = Σ s^{n(n+1)/2}(−s)ₙ. So 𝒥 sits halfway between f₀ and f₁, and 𝒦/2 halfway between ψ₁ and ψ₀.

| near q → 1 (base q) | dual function | dual nome |
|---|---|---|
| ψ₁\*, ψ₀\* (sign-twisted ψ) | W, then q̃W₃: Σ(−1)ʲp^{j²+j}(p;p²)ⱼ | p = e^{−πi𝒳} |
| f₁\*, f₀\* (sign-twisted f) | −2F_a, −2F_b: Σ p^{2n²−n}/(p;p²)ₙ | p |
| G₀, G₁ (halves of F\*) | 𝒦: unsigned ψ at half shift | r = e^{−πi𝒳/2} |
| F₀\*, F₁\* (sign-twisted F) | 𝒥: unsigned f at half shift | r |

In the variable z of G(z) = Σ zʲq^{j²}(q;q²)ⱼ (v4 §1): the W family is G at z = −q, −q³, and G₀, G₁ are G at z = 1, q² (G₀ = G(1) − 1, G₁ = q·G(q²)). So W and the F\* halves are the same module at z and −z, and the two rows of the table for ψ\* and for F\* are dual to each other in the sense of v4 §4: W near 1 goes to ψ\*, and G₀ near 1 goes to an unsigned ψ. The whole census of sign-twisted fifth-order functions is closed under these dualities, up to half shifts.

---

## 4. Open

- Proof of (VII); decaying and trivial sectors for F\*.
- The unsigned 𝒥, 𝒦 at half shift: their own behaviour at roots of unity and their dual should close the loop back to F\* (not tested).
- The one-sided halves in general (v7 §4, item 2) remain numerical.

## 5. Scripts

`fstar.py` (limits and Taylor coefficients), `fstar_check.py` (identity (VII) exactly; comparison with 𝒥, 𝒦; first corrections). Helpers from v5.
