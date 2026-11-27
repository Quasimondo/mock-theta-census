# Mock theta search: addendum v23 (the kernel integral at real x is bounded; the tetrahedron-volume rule at even roots of unity is a theorem)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). The one bound left open in `addendum-v22.md` §4. Scripts and output are in `psi_star/v23/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v23/check_lemmaJ.py`, own code including an own Faddeev integral): Lemma U holds to 10⁻¹⁵ (|Re log Φ_b| at 50 real points, X = 8, 40, 200). |P·J(a/b)| computed over the real line agrees with every value of §3 tested, to all six digits (a = ¼ and ½ at X = 8, 40, 200; a = 0.7 at X = 8; a = 0.05 and 0.95 at X = 40), lies below 1/(2√2·sin πa), and approaches |1/(1 − 2z)|. The case a = 0.05 at X = 8 was not reached, because its integration range is too long for the simple Faddeev quadrature used. The two analytic steps, |Φ_b| = 1 by the conjugation symmetry and ∫ e^{2πas}/(1 + e^{2πs}) ds = 1/(2 sin πa), are standard and were checked by hand. Theorem 15 against the exact sums at X = 400 … 3200: X·(F/prediction − 1) converges, and Richardson extrapolation gives −0.21657 − 1.34534i (a = ¼), −0.18821 − 1.50019i (a = ⅓), 0.09094 − 1.15995i (a = 5/6) and −0.05388 − 1.15907i (a = 1/20): the prediction holds with an O(1/X) correction.

Labels: **proved** (a complete argument is written here or in the addenda cited; nobody else has checked it), **verified** (numerics).

**Later remarks:** `addendum-v27.md` §1 (Theorem 15 with its error written out in κ = X·Vol/2π; it is not uniform near z = 1), §3 (Lemma J: the contour turn is pointwise in x; the lemma holds for x ≤ 1).

---

## 0. Summary

1. **Lemma J (proved; §1).** For real x > 0 and 0 < a < 1: |P·J(a/b)| ≤ 1/(2√2·sin πa). The reason is elementary: for real b Faddeev's dilogarithm has modulus 1 on the real line, the contour of J can be turned onto the real line, and there the integrand has modulus e^{2πas}/(1 + e^{2πs}). No descent inequality is needed.
2. **Theorem 15 (proved; §2).** For every z = e^{2πia} with 0 < a < 1, as X → ∞ through even integers:

   F(z; e^{2πi/X}) = (1 + W)·W^{1/2}·(π/(t(1 + 2W)))^{1/2}·e^{(V(w₀) + π²/12)/t}·(1 + O(1/X)),  t = −2πi/X,

   with W the root of W² + W = 1/z in the upper half-plane, w₀ = Log W, V as in v22. In particular |F(z; e^{2πi/X})| = e^{Vol(z)·X/2π + O(log X)}: the growth rate is the volume of the ideal tetrahedron with shape W, divided by 2π. This is the rule found numerically in `census3/zfamily/` and v18, now a theorem for even X.
3. **A simplification for z = −1.** At a = ½ Lemma J gives |P·J(1/(2b))| ≤ 2^{−3/2} at real x. Theorem 10 of v20 and Theorem 12 of v21 used the uniform bound of v15/v16 for this; they no longer need it.
4. **Numerically (verified):** |Φ_b| = 1 on the real line to 10⁻¹²; the integral over the real line equals the one over the tilted line to 10⁻⁹ and respects the bound in all 21 cases; X·(F/prediction − 1) converges for six values of z.
5. **Not covered:** odd X for z ≠ −1; z = 1; the expansion of P·J(a/b) beyond its boundedness.

---

## 1. The kernel integral at real x

**Lemma U (unitarity).** For real b > 0 and real y: |Φ_b(y)| = 1.

Proof. log Φ_b(y) = ∫_{ℝ+i0} f(ξ) dξ with f(ξ) = e^{−2iyξ}/(4ξ·sinh(bξ)·sinh(ξ/b)) (v16 §1). For real b, y the conjugate of f(ξ) is −f(−ξ̄). The map ξ → −ξ̄ sends the contour ℝ + i0 to itself with the orientation reversed, so the conjugate of the integral is −∫_{ℝ+i0} f(ξ′) dξ′. Hence log Φ_b(y) is purely imaginary. ∎

**Lemma J.** Let x > 0 (b = √x real), x ≤ ¼, 0 < a < 1, and J(a/b) the integral of v12 over L = e^{−iδ}·ℝ, 0 < δ < 90°, with Faddeev's Φ_b. Then

  J(a/b) = ∫_ℝ e^{−πiy² + 2πay/b}/(Φ_b(y)·(1 + e^{2πy/b})) dy,  |J(a/b)| ≤ b/(2·sin πa),  |P·J(a/b)| ≤ 1/(2√2·sin πa).

Proof.

- *Turning the contour.* Between L and ℝ lie the two wedges −δ ≤ arg y ≤ 0 and π − δ ≤ arg y ≤ π. The kernel poles y = i·b·(n + ½) are on the imaginary axis and the zeros of Φ_b on the negative imaginary axis, so the integrand is analytic in the wedges. In both wedges sin(2·arg y) ≤ 0, so |e^{−πiy²}| ≤ 1. In the left wedge 1/|Φ_b| ≤ C (Lemmas F and D of v16, as in Lemma R3 of v20); in the right wedge 1/Φ_b(y) = e^{−πiy²}·Φ_b(−y)/Φ_b(0)² with −y in the left wedge, where |Φ_b| ≤ C as well, so 1/|Φ_b(y)| ≤ C·|e^{−πiy²}| ≤ C. With s = y/b the remaining factor e^{2πas}/(1 + e^{2πs}) is at most C·e^{−2π·min(a, 1−a)·|Re s|} on the arcs |y| = R of the wedges, and |Re s| ≥ (R/b)·cos δ there. So the arcs contribute nothing as R → ∞.
- *On the real line.* |e^{−πiy²}| = 1 and |Φ_b(y)| = 1 (Lemma U), so the integrand has modulus e^{2πas}/(1 + e^{2πs}), and ∫_ℝ e^{2πas}/(1 + e^{2πs})·b ds = b/(2·sin πa).
- *Prefactor.* P = 2^{−1/2}·e^{π²/(12t) + t/6}/b with t = −2πix purely imaginary, so |P| = 2^{−1/2}/b. ∎

Remarks.

- The bound holds for every real x in the range, not only asymptotically. For a = ½ it is 0.3536, and the true value tends to 1/3.
- Off the real axis the argument does not apply as it stands (|Φ_b| ≠ 1 on the real line for complex b); there the bound for a = ½ is the one of v15/v16.

---

## 2. The volume rule at even roots of unity

**Theorem 15.** Let 0 < a < 1, z = e^{2πia}, W the root of W² + W = 1/z with Im W > 0, w₀ = Log W, V(w) = −w²/2 − 2πi(a − 1)w + Li₂(−e^{w}), Vol(z) = −Im V(w₀) = |D(−W)| > 0. For even X and t = −2πi/X:

  F(z; e^{2πi/X}) = (1 + W)·W^{1/2}·(π/(t(1 + 2W)))^{1/2}·e^{(V(w₀) + π²/12)/t}·(1 + O(1/X)) + E,  |E| ≤ 1/(2√2·sin πa) + O(X^{1/2}·e^{−Vol(z)·X/2π}),

with principal roots and W^{1/2} = e^{w₀/2}. The first term has modulus |1 + W|·|W|^{1/2}·(X/(2|1 + 2W|))^{1/2}·e^{Vol(z)·X/2π}, so E is exponentially smaller, and

  lim_{X→∞, X even} (1/X)·log|F(z; e^{2πi/X})| = Vol(z)/2π.

Proof. Theorem 14 of v22: F(z; q₀) = i·P·e^{−πia}·[J(a/b) − I((a−1)/b) + z^X·I(a/b)] at x = 1/X.

- −i·P·e^{−πia}·I((a−1)/b) is the first term, by Theorem 13 (a) of v22.
- |i·P·e^{−πia}·J(a/b)| ≤ 1/(2√2·sin πa) by Lemma J.
- |P·I(a/b)| = O(X^{1/2}·e^{−Vol(z)·X/2π}) by Theorem 7 of v19 at real x (its saddle value has Im V = +Vol(z)), and |z^X| = 1. ∎

So for every z on the unit circle other than 1, the finite sums F(z; e^{2πi/X}) at even X grow at the rate of an ideal tetrahedron's volume, with the shape W determined by z·W·(1 + W) = 1; the constant and the phase are those verified in v18 §2.2. The constants in the O-terms depend on a and degenerate as z → 1, where Vol(z) → 0 and the bound for E blows up.

**For z = −1.** Lemma J at a = ½ replaces the reference to v16 (Theorem 3) in the proofs of Theorem 10 (v20) and Theorem 12 (v21). The chain for those two theorems is then: v12 (formula), v16 Lemmas Φ, F, D and v19 (state integrals at real x), v20/v21 (radial limits), and Lemma J; the descent inequality of v15 and the uniform bound for the kernel integrand are needed only for Borel summability (v17) and the law inside the disk (v19 Theorem 8).

---

## 3. Numerical checks (verified; `p1_unitary.py`, `run_log.txt`)

| check | result |
|---|---|
| Lemma U: \|Re log Φ_b(y)\| at 200 random real points, x from 0.005 to 0.25 | at most 8·10⁻¹³ (double precision) |
| Lemma J: integral over the real line against the integral over L (δ = 10°), a = 0.05, ¼, ⅓, ½, 0.7, 5/6, 0.95; X = 8, 40 | differences 10⁻¹¹ … 5·10⁻⁹ |
| Lemma J: \|P·J(a/b)\| against 1/(2√2·sin πa), the same a, X = 8, 40, 200 | below the bound in all 21 cases (for example a = ½: 0.3264, 0.3330, 0.3333 against 0.3536; a = ¼: 0.468, 0.471, 0.452 against 0.5; a = 0.05: 0.97, 1.82, 0.93 against 2.26) |
| Theorem 15: X·(F/prediction − 1) at X = 400, 800, 1600, 3200 for a = ¼, ⅓, 5/6, 0.3, 0.05, 1/√2 | converges in each case (for example a = ¼: −0.2201 − 1.3440i, …, −0.2170 − 1.3452i) |

The extrapolated agreement of F with the prediction to 10⁻²⁶ … 10⁻³⁶ is in v18 §2.2.

---

## 4. Open

1. Odd X for z ≠ −1. By the census (README item 10) F(z) then grows only when z^X is close to 1; the dual nome tends to −1 and the dual series 𝒲(T) need the treatment of v21 with the extra variable T.
2. The asymptotic expansion of P·J(a/b) (expected 1/(1 − 2z) + …) and its Borel summability for a ≠ ½.
3. z → 1: the constants degenerate; the mock theta case z = 1 itself is of a different kind (v18 §3).
4. The remaining items of v19 §8, v20 §5, v21 §5, v22 §7.

---

## 5. Scripts (`psi_star/v23/`)

- `p1_unitary.py [a b c]`: the three checks of §3 (3 minutes).
- `n1_zroots.py` (from v22), `fad.py` (from v16), `run_log.txt`, `README.md`.
