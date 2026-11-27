# Mock theta search: addendum v7 (the bilinear combinations are state integrals)

Session of 26 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). This does step 1 of `addendum-v6.md` §4: write the state integral explicitly and compare it with 𝓑₁…𝓑₄. Script and output are in `psi_star/v7/`.

**Independent check** (Claude Opus, same day, `psi_star/v7/check_state_integral.py`, own implementation of Φ_b, the contour and the q-series): the identity of §1 holds for 𝓑₁ and 𝓑₃ at x = 1/(3 − 0.4i), 0.3 + 0.2i, 1.2 + 0.5i, 0.11 + 0.03i, to 10⁻²⁶ … 10⁻³¹ (one case, 𝓑₃ at x = 0.11 + 0.03i, to 10⁻¹¹, limited by the quadrature mesh). Note for anyone reimplementing: q^{1/2} and q̃^{1/2} must be taken as e^{πix} and e^{−πi/x}, not as principal square roots, or the check fails whenever |Re x| or |Re 1/x| exceeds ½. The saddle computation of §3 (saddles z = ±2πi/3, V = −π²/12 ∓ i·Cl₂(π/3), V″ = ∓i√3/2) was redone by hand and is correct.

Labels: **proved**, **verified** (numerics), **derivation** (formal computation given here, confirmed numerically, analytic details not written out), **inference**.

---

## 0. Summary

1. **𝓑₁…𝓑₄ are state integrals** (§1). For every x in the upper half-plane,
   𝓑ⱼ = (explicit elementary factor)⁻¹ · ∫ Φ_b(y)·e^{−πi·y²/2 + 2π·uⱼ·y} dy,  b² = x,
   with Faddeev's quantum dilogarithm Φ_b and four explicit values uⱼ. Derivation by residues; verified at 4 values of x for all four j, to 10⁻²⁷ … 10⁻³¹.
2. **Consequence** (§2, inference from 1 and known properties of Φ_b). The integral does not see the dual nome and is analytic in b across the positive real axis. That is the one analytic claim that v6 left open. With it, the two-sector law for R₀, R₁, and so for the bilateral parts θ₄R₀, θ₄R₁ of ψ\*, follows from the exact inversion formulas of v6 §2.
3. **The saddle point of the integral reproduces the measured constants** (§3, derivation): actions ±i·Cl₂(π/3) and modulus (4/3)^{1/4}.
4. **Still open:** a written-out justification of the contour closing; the phase e^{∓iπ/12} from the saddle; the one-sided halves f\*, F_a, F_b.

---

## 1. The identity (`state_integral.py`)

Notation. b² = x with Im x > 0; q = e^{2πix}, q̃ = e^{−2πi/x}, p = q̃^{1/2} = e^{−πi/x}; c_b = i(b + b⁻¹)/2;

  Φ_b(y) = (−q^{1/2}e^{2πby}; q)∞ / (−q̃^{1/2}e^{2πy/b}; q̃)∞.

In v5/v6 notation x = 1/𝒳, and p is the dual nome there. Let

  I(u) = ∫_C Φ_b(y)·exp(−πi·y²/2 + 2π·u·y) dy,

with C running from ∞·e^{i(π−δ)} through 0 to ∞·e^{iδ} (any small δ > 0; the integrand decays like a Gaussian at both ends, and C passes below all poles of Φ_b).

**Claim.** With 𝓑₁…𝓑₄ as in v6 §2,

  I(uⱼ) = −i·b·((q;q)∞/(q̃;q̃)∞)·exp(−πi·c_b²/2 + 2π·uⱼ·c_b)·𝓑ⱼ,

  u₁ = b/4,  u₂ = −b/4,  u₃ = b/4 − 1/(2b),  u₄ = −b/4 − 1/(2b).

**Derivation.** The poles of Φ_b are simple, at y = c_b + i·b·k + i·n/b with k, n ≥ 0, with residue

  −(b/2π)·((q)∞/(q̃)∞)·(−1)ⁿq̃^{n(n+1)/2} / ((q;q)ₖ(q̃;q̃)ₙ).

At such a pole the exponential factor is exp(−πic_b²/2 + 2πuc_b) times q^{(k²+k)/4}·q̃^{−(n²+n)/4}·i^{k+n}·(−1)^{kn}·ζᵏ·ζ̃ⁿ, where ζ = e^{2πiub}, ζ̃ = e^{2πiu/b}. Closing C in the upper half-plane and summing residues gives

  I(u) = −i·b·((q)∞/(q̃)∞)·e^{−πic_b²/2 + 2πuc_b}·Σ_{k,n≥0} (−1)^{kn}·Zᵏq^{k²/4}/(q)ₖ · Z̃ⁿq̃^{n²/4}/(q̃)ₙ,  Z = iζq^{1/4}, Z̃ = −iζ̃q̃^{1/4}.

The sign (−1)^{kn} splits the sum by the parity of k:

  Σ = H_even(Z; q)·H(Z̃; q̃) + H_odd(Z; q)·H(−Z̃; q̃),

with H(z;q) = Σ zᵐq^{m²/4}/(q)ₘ as in v6 §1. For u = b/4: Z = i·q^{1/2} and Z̃ = q̃^{1/4}, so the sum is N₀·E₊′(p) + i·q^{−1/4}·Y·E₊(p) = 𝓑₁. The other three values of u give 𝓑₂, 𝓑₃, 𝓑₄ in the same way. ∎ (formal: the closing of the contour is not justified here)

**Check.** At x = 1/(3 − 0.4i), 0.3 + 0.2i, 1.2 + 0.5i, 0.11 + 0.03i, numerical quadrature of I(uⱼ) against the q-series side: |ratio − 1| between 10⁻³¹ and 3·10⁻²⁷ in all 16 cases (working precision 30 digits). The four u-values and the prefactor were derived first and matched on the first run; nothing was fitted.

The double series converges absolutely for |q|, |q̃| < 1, since both factors are Nahm sums with A = ½. This is the case A = ½, B = 1 of the family ∫Φ_b^B·e^{−Aπiy²}dy that Garoufalidis–Kashaev treat for integer A; the only new feature of the half-integer exponent is the parity split, which is what produces two terms and hence a 2×2 structure.

---

## 2. What follows (inference)

- **Independence of the dual nome.** I(u) is defined without reference to q or q̃ separately. Φ_b is analytic in b for b² outside (−∞, 0], in particular for real b > 0, where it is given by Faddeev's integral formula instead of the product. So each 𝓑ⱼ, times its elementary factor, extends analytically from the upper half-plane across the positive real x-axis.
- **Asymptotics.** As x → 0 inside the domain of analyticity, the integral has a stationary-phase expansion in powers of x (§3). The horocycle limit of v5/v6 is x → 0 tangentially to the positive real axis, which is interior to that domain. So the expansion holds there, with coefficients that cannot depend on p. This is what v6 §2 measured to 10⁻¹⁷.
- **The law.** The inversion formulas of v6 §2 are exact identities. Combined with the two points above they give, to all orders in t,
  R₀(q) = 2^{−1/2}·e^{−π²/(12t) − t/24}·[𝓑₃·P₊(p) + 𝓑₁·P₋(p)],
  and the same for R₁ with 𝓑₂, 𝓑₄, where now 𝓑ⱼ stands for its asymptotic series. Multiplying by the modular form θ₄ gives the two exponential sectors of the bilateral parts of ψ₀\*, ψ₁\*.

So for the bilateral parts the modularity law is reduced to two standard facts about state integrals: the residue evaluation (§1) and stationary phase near b = 0. I have not written either out rigorously.

---

## 3. Saddle point (derivation)

With z = 2πby and ħ = 2πx, Φ_b(y) ≈ exp(Li₂(−e^z)/(iħ)), and the integrand is exp(V(z)/(iħ) + O(1)) with

  V(z) = Li₂(−e^z) + z²/4.

- Saddles: V′ = 0 gives 1 + e^z = e^{z/2}, so e^{z/2} = e^{±iπ/3}, z = ±2πi/3.
- Values: V = −π²/12 ∓ i·Cl₂(π/3). In the variable t = −2πi·x of v5/v6 this is a factor e^{(π²/12 ± i·Cl₂)/t}.
- The elementary factor in §1 is e^{π²/(12t)} to leading order. Dividing, 𝓑ⱼ ~ e^{±i·Cl₂(π/3)/t}, as measured.
- Hessian: V″ = ∓i·√3/2 at the two saddles, so the Gaussian integral contributes modulus (2/√3)^{1/2} = (4/3)^{1/4}, as measured. I did not track the phase; the measured value is e^{∓iπ/12}.

Which saddle dominates for which uⱼ is decided by the linear term: u₁, u₂ pick the decaying exponential on the horocycle and u₃, u₄ the growing one (verified in v6, not derived here).

---

## 4. Open

1. Justify the contour closing for A = ½ (the Garoufalidis–Kashaev argument for integer A should carry over), and track the phase and the higher orders in §3.
2. The one-sided halves. ψ\* = ½f\* ± ½θ₄R, and W = −F_a + ψP₊: the laws for f\*, F_a, F_b (v5) are not covered by the homogeneous module. A natural guess is a state integral with one more factor or a different contour (inference).
3. The converse and mirror directions for the bilateral parts: by §2 they should follow from the same four integrals after exchanging the roles of q and q̃ and using Γ₀(2). Not checked.
4. F₀\*, F₁\*: their bilateral parts are X, Y, which already appear in 𝓑₁…𝓑₄. The inversion formulas of v6 give their two-sector law at once; it has not been compared with direct numerics.

---

## 5. Scripts (`psi_star/v7/`)

`state_integral.py`: q-series side of 𝓑₁…𝓑₄, the four integrals by tanh-sinh quadrature on the contour C (δ = 0.4), and the comparison. Needs mpmath only. About 3 minutes at 30 digits.
