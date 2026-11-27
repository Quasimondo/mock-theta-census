# Mock theta search: addendum v27 (after the second stress test: Theorem 15 with X·Vol as the parameter, the sector of Theorem 13 (b), remarks on Lemma J and on the closing lemma)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). An answer to the findings F1, F2, F3, F5 of `docs/psi-star-stress-test-v2.md` (Claude Opus 5.5), which tested addenda v22–v26 and found no violation. Scripts and output are in `psi_star/v27/`.

Labels: **proved** (a complete argument is written here or in the addenda cited; nobody else has checked it), **verified** (numerics here), **measured** (numbers taken from the stress test, not recomputed here).

---

## 0. Summary

1. **Theorem 15 with κ = X·Vol(a)/2π as the parameter (F1; proved, §1).** The error of Theorem 15 has three parts, and only one of them is of order 1/X:

   |F/Π − 1| ≤ C(a)/X + T_J + T_dec·(1 + C′(a)/X),  T_J = e^{−κ}·|(1 + 2W)/(1 + W)|^{1/2}/(2·X^{1/2}·sin πa),  T_dec = |W|·e^{−2κ}.

   T_J comes from the kernel integral (Lemma J) and T_dec is the ratio of the decaying to the growing state integral. Near z = 1 the volume is small, and the error is of order 1/X only once κ exceeds about ½·ln X + ln(1/sin πa). Before that it is T_dec: at a = 0.001, X = 800 one has κ = 0.385 and X·T_dec = 599, against the stress test's X·|F/Π − 1| = 650. That C(a) stays bounded as a → 0, 1 is measured (|K| ≤ 1.6), not proved.
2. **Theorem 13 (b): the sector ε(a) (F2; measured, §2).** The angle measured by the stress test is 20° near a = 0.3 and 17° at a = ½, and it is proportional to the volume near both ends: ε(a)/Vol(a) is about 26° per unit of volume for a ≤ 0.3 and falls to about 7° at a = 0.999. A proof of ε(a) ≥ c·Vol(a) is open; the proof of v22 gives only a tilt of at most arctan(Vol/(3π²)) (2° at a = ½) and an unspecified second condition.
3. **Lemma J (F3; proved, §3).** The turn of the contour is made for each fixed x, with a constant that grows as x → 0; that constant does not enter the result, so the bound of the lemma is uniform. The lemma holds for 0 < x ≤ 1, not only for x ≤ ¼.
4. **The closing lemma of v12 (F5; §4).** The remark "times a factor between 1 and 1.5" referred to the modulus, not to the exponent, and is reworded.
5. **By-product (proved, §5):** with Lemma R3′ of v25, Theorem 17 of v24 (kernel integral = Mordell integral) holds on the unit circle, for real 0 < x ≤ 1.
6. **Numerically (verified, §6):** the inequality of item 1 at 27 pairs (a, X) from a = 0.001 to 0.999; the two exact moduli; Lemma J for x up to 1.

---

## 1. Theorem 15 with X·Vol as the parameter

Notation of v22, v23: 0 < a < 1, z = e^{2πia}, W the root of W² + W = 1/z with Im W > 0 (the other root is −1 − W), w₀ = Log W, Vol = Vol(z) = −Im V(w₀), X even, t = −2πi/X, and

  Π = (1 + W)·W^{1/2}·(π/(t(1 + 2W)))^{1/2}·e^{(V(w₀) + π²/12)/t},  κ = X·Vol(z)/2π.

**Theorem 15′.** F(z; e^{2πi/X}) = Π·(1 + ρ₁) + E_J + E_dec, where

- (i) |Π| = (X/2)^{1/2}·|(1 + W)/(1 + 2W)|^{1/2}·e^{κ} exactly;
- (ii) Π·(1 + ρ₁) = −i·P·e^{−πia}·I((a−1)/b), and |ρ₁| ≤ C(a)/X (Theorem 13 (a));
- (iii) E_J = i·P·e^{−πia}·J(a/b), and |E_J| ≤ 1/(2√2·sin πa) for every X (Lemma J);
- (iv) E_dec = i·P·e^{−πia}·z^X·I(a/b), and |E_dec| = (X/2)^{1/2}·(|W|/|1 + 2W|)^{1/2}·e^{−κ}·(1 + ρ₂), |ρ₂| ≤ C′(a)/X (Theorem 7).

Hence

  |F/Π − 1| ≤ C(a)/X + T_J + T_dec·(1 + C′(a)/X),  T_J = e^{−κ}·|(1 + 2W)/(1 + W)|^{1/2}/(2·X^{1/2}·sin πa),  T_dec = |W|·e^{−2κ}.

Proof. The decomposition is Theorem 14 of v22. (i): |W|·|1 + W| = |1/z| = 1, so |1 + W|·|W|^{1/2} = |1 + W|^{1/2}; (π/|t|)^{1/2} = (X/2)^{1/2}; 1/t = iX/2π, so |e^{V(w₀)/t}| = e^{−Im V(w₀)·X/2π} = e^{κ} and |e^{π²/(12t)}| = 1. (iv): Theorem 7 with the lower root W₁ = −1 − W: |P| = 2^{−1/2}/b, b² = 1/X, |t| = 2π/X give |P·I(a/b)| = (X/2)^{1/2}·|(1 + W₁)/(1 + 2W₁)|^{1/2}·e^{−κ}·(1 + ρ₂), and 1 + W₁ = −W, 1 + 2W₁ = −(1 + 2W). The quotients: |E_J|/|Π| ≤ T_J, and |E_dec|/|Π| = (|W|/|1 + W|)^{1/2}·e^{−2κ}·(1 + ρ₂) = |W|·e^{−2κ}·(1 + ρ₂) because |W|/|1 + W| = |W|². ∎

**What this says about uniformity (F1).**

- T_J and T_dec are explicit and uniform: |W| lies between 0.618 and 1.618, |(1 + 2W)/(1 + W)|^{1/2} between 1.17 and 1.91. They are below 1/X when κ ≥ ½·ln X + ln(1/sin πa) + 1, say. Near z = 1, Vol ≈ 2π·ln φ·min(a, 1 − a) (v18), so this needs X of the order (1/a)·ln(1/a).
- For smaller κ the growing term has not yet separated from the others, and the relative error is T_dec (and T_J), not 1/X. This is the stress test's observation: at a = 10⁻³, X = 800 it finds X·|F/Π − 1| = 650, and X·T_dec = 599 there.
- The constants C(a), C′(a) of Theorems 13 (a) and 7 are not controlled as a → 0, 1 by the proofs. **Measured** (stress test, T15.1): the limit K(a) = lim X·(F/Π − 1) exists for every a, 1.15 ≤ |K(a)| ≤ 1.6 on (0, 1), and K → −1.152i at both ends. So the expected uniform statement is |F/Π − 1| ≤ C/X + T_J + C·T_dec with absolute C. It is open; it needs Theorems 13 (a) and 7 with constants uniform in a, where the saddles approach the edge Im w = ±π of the strip.
- Theorem 15 as stated in v23 ("(1 + O(1/X)) for each a, constants depending on a") is correct; Theorem 15′ is the same statement with the error written out.

---

## 2. The sector of Theorem 13 (b)

Theorem 13 (b) of v22 says: there is ε(a) > 0 such that the asymptotic formula of (a) holds uniformly for 0 ≤ arg x ≤ ε(a). The proof gives ε(a) = min(φ₁, arctan(c/(2M))) with a tilt φ₁ < arctan(Vol(z)/(3π²)) and constants c, M that are not computed.

**Measured sector** (stress test, T13b.4: the largest arg x for which (L1) holds on the tilted contour C′, optimised over the tilt; T13b.5: the formula of (a) holds up to that angle with a flat O(t) constant ≤ 0.10):

| a | 0.001 | 0.005 | 0.01 | 0.02 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 0.8 | 0.9 | 0.95 | 0.99 | 0.999 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Vol(a) | 0.0030 | 0.0151 | 0.0302 | 0.0604 | 0.151 | 0.299 | 0.575 | 0.804 | 1.015 | 0.804 | 0.575 | 0.299 | 0.151 | 0.0302 | 0.0030 |
| ε(a), degrees | 0.08 | 0.39 | 0.78 | 1.55 | 3.9 | 7.8 | 16.2 | 20.1 | 17.2 | 9.6 | 6.2 | 3.0 | 1.4 | 0.25 | 0.02 |
| ε/Vol, degrees | 26 | 26 | 26 | 26 | 26 | 26 | 28 | 25 | 17 | 12 | 11 | 10 | 9 | 8 | 7 |
| tilt allowed by the proof, degrees | 0.006 | 0.03 | 0.06 | 0.12 | 0.29 | 0.58 | 1.1 | 1.6 | 2.0 | 1.6 | 1.1 | 0.58 | 0.29 | 0.06 | 0.006 |

Remarks.

- ε(a) is proportional to the volume near both ends, with different constants: about 26° per unit of volume for a ≤ 0.3, and about 7° to 10° for a ≥ 0.9 (the last two entries of ε have one or two digits). On the measured range ε(a) ≥ 6°·Vol(a). The value c ≈ 20° suggested in the stress test fits a ≤ 0.3 only.
- For a ≥ 0.9 the limit is set by the pieces C₄, C₅, not by the tilt (stress test); turning C₄ is the place to gain.
- Open: a proof of ε(a) ≥ c·Vol(a). It needs the margin of (L1) on C to be bounded below by a multiple of Vol; the stress test measures sup (H − H(w₀))/min(d², d) ≈ −0.315·Vol as a → 0 and −0.118·Vol as a → 1.

---

## 3. Lemma J: the contour turn, and the range of x

Two additions to Lemma J of v23 §1.

**The turn is pointwise in x.** The proof bounds 1/|Φ_b| in the left wedge and |Φ_b| there (for the right wedge) by a constant C, from Lemmas F and D of v16. Their proofs apply on every ray of the wedge (the margin of Lemma D on a ray at angle δ′ < δ is at least that on L). C depends on x and on δ: for δ = 10° the stress test finds sup 1/|Φ_b| = e^{1.77} at x = 0.005, and it grows as x → 0, because 1/|Φ_b| ≈ exp(−Im Li₂(−e^{w})/|t|) and Im Li₂ < 0 in part of the wedge. This is harmless: the turn is made for one x at a time, and C does not enter the result, since on the real line only |Φ_b| = 1 is used. So the bound of Lemma J is uniform in x, although the estimate that justifies the turn is not.

**Lemma J holds for 0 < x ≤ 1** (v23 stated x ≤ ¼). The restriction came from the range of Lemmas F and D; by the proof of Lemma R3′ (v25 §1) they hold at every real x₀ ≤ 1, with constants depending on x₀, which is all the turn needs. **Measured** (stress test, J.4): the largest value of |P·J(a/b)|·2√2·sin πa over 0 < x ≤ 1 and 10⁻³ ≤ a ≤ 1 − 10⁻³ is 0.9655 (a = 0.284, x = 1/16), so the bound is sharp to 3.5 %; it tends to 0 linearly at both ends of (0, 1).

So in Theorem 15′ the bound for E_J holds for every even X ≥ 2.

---

## 4. The closing lemma of v12 §4.1: the remark on the factor

The sentence in the proof, "(The second value belongs to the point on C; `t5_bounds.py` confirms that the supremum on the line is this bound times a factor between 1 and 1.5.)", is to be read as follows, and is replaced by:

"(The second value belongs to the point on C. Numerically, for the product of the two q-factors and the Gaussian, without the linear factor, the supremum of the modulus on the line is e^{−c·σ_N²} times a factor that stays bounded in N: between 1.0 and 1.5 in the two cases of `t5_bounds.py`. In terms of the exponent E with the linear factor included: sup E/(−c·σ_N²) → 1 as N → ∞, from below (0.94 at N = 2, 0.98 to 1.00 at N = 8, 0.9998 at N = 32), and sup E + c·σ_N² ≤ 0.141·(σ_N + 1) in the samples of the stress test. So c is the exact leading constant, and the difference is the lower-order term e^{O(σ_N)} that the proof allows.)"

The proof is not affected.

---

## 5. Further facts from the second stress test

**Theorem 17 on the unit circle (proved, given the results quoted in v24 §4).** For real 0 < x ≤ 1, with J over L and Faddeev's Φ_b, and t = −2πix:

  J(1/b) = e^{−πi/4}·(√15·t/π)·q^{−1/60}·p^{1/24}·L(2/5, 10t),  J(b + 1/b) = e^{−πi/4}·(√15·t/π)·q^{−19/60}·p^{1/24}·L(1/5, 10t),

where L(r, α) for α on the imaginary axis is the integral over the ray arg x = −(arg α)/2.

Proof. For Im x > 0 this is Theorem 17. The left sides are continuous on {Im x ≥ 0} near every real x₀ ∈ (0, 1] by Lemma R3′ of v25 (u = 1/b and u = b + 1/b are continuous and bounded there). The right sides are analytic in t for |arg t| < 180°: on the ray x = ρ·e^{−i(arg α)/2} one has αx² = |α|ρ² and αx = |α|ρ·e^{i(arg α)/2}, so the Gaussian is e^{−3|α|ρ²/2} and cosh(3αx/2) has no zero; for Re α > 0 the turned integral equals the original one, because the poles of 1/cosh(3αx/2) have arguments ±90° − arg α, outside the wedge between the two rays, and Re(αx²) > 0 in it. ∎

This is item 2 of v24 §6 as far as the kernel integrals are concerned. **Measured** (stress test, M.5, M.6): the identity holds at x = 1/X, X = 1 … 100, to 1.2·10⁻¹¹, and J(x + iε) − J(x) is proportional to ε.

**Measured, for the record** (none changes a statement):

| fact | where | value |
|---|---|---|
| margin of (L1) on the contour C, sup (H − H(w₀))/min(d², d) | v22 §2 | negative for a from 10⁻⁵ to 1 − 10⁻⁵; ≈ −0.315·Vol (a → 0), −0.118·Vol (a → 1) |
| constant in 1/\|Φ_b\| ≤ C·\|e^{Li₂/t}\| on C | v22 §2 (L2) | exactly 2, attained at the corner 2πi |
| O(t) constant of Theorem 13 (a) | v22 | ≤ 0.13 for 0.25 ≤ a ≤ 0.98 once κ ≳ 2 |
| C₁ of Lemma F on x ∈ [¼, 1] | v25 §1 | 0.39 |
| \|E\| in Lemma D for \|t\| ≤ 2π·1.01 | v25 §1 (claims 1.2) | ≤ 0.17 |
| (R3), (R4) of v24 | v24 §4 | hold to q¹⁰⁰⁰ (still not proved here) |
| λ₋ and g(λ₋) in Lemma V | v26 §3 (claims < 0.16) | λ₋ ∈ [0.500, 0.693], g(λ₋) ≤ 0.083 |
| Theorems 9, 11, 14 at X = 1, 2, 3 in 25-digit arithmetic | v25 | 10⁻¹⁵ or better at X = 1, 2; 10⁻¹² at X = 3 (limit of the stress test's quadrature, F4) |

---

## 6. Numerical checks (verified; `w1_round2.py`, `run_log.txt`)

| check | result |
|---|---|
| Theorem 15′ (i): \|Π\| against its closed form, 27 pairs (a, X), a = 0.001 … 0.999, X = 10 … 6400 | ≤ 6·10⁻³⁰ |
| Theorem 15′ (iv): \|P·I(a/b)\| against (X/2)^{1/2}·(\|W\|/\|1 + 2W\|)^{1/2}·e^{−κ}, a = 0.1, 0.3, 0.5, 0.8, X = 16, 32, 64 | relative difference ≤ 0.011, 0.006, 0.003: O(1/X) |
| the inequality: \|F/Π − 1\|/(1.7/X + T_J + 1.1·T_dec) at the 27 pairs | ≤ 0.952 |
| small κ: a = 0.001, X = 800 (κ = 0.385): X·\|F/Π − 1\|, X·T_dec, X·T_J | 648.6, 599.4, 5827 |
| large κ: X·\|F/Π − 1\| at a = 0.01 (X = 3200), 0.3 (800), 0.5 (200), 0.99 (3200) | 1.152, 1.451, 1.591, 1.152 |
| Lemma J for x = 0.3, ⅓, ½, ¾, 1 and a = 0.05 … 0.95: \|P·J(a/b)\|·2√2·sin πa | ≤ 0.866 |
| ε/Vol from the stress test's table | 6.6° … 28.2° |

The constant 1.7 in the third line stands for C(a); it is above every value of X·|F/Π − 1| found at large κ (the largest here is 1.62, at a = ½, X = 50).

---

## 7. Scripts (`psi_star/v27/`)

- `w1_round2.py [a b c]`: the checks of §6 (1 minute).
- `n1_zroots.py` (from v22), `fad.py` (from v16), `run_log.txt`, `README.md`.
