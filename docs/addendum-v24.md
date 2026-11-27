# Mock theta search: addendum v24 (z = 1: the kernel integral is a Mordell integral, and the formula of the family is the classical transformation law of ψ₀, ψ₁, f₀, f₁)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). The item "Mordell integral at z = 1" of `addendum-v18.md` §3.2. Scripts and output are in `psi_star/v24/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v24/check_thm17.py`, own code: product formula for Φ_b, own quadratures, own q-series, own Mordell integral). At x = 0.3 + 0.2i, 0.3i and 0.11 + 0.08i (25 digits):

- **Proposition M.** All four state integrals I(0), I(1/b), I(b), I(b + 1/b), computed by direct quadrature instead of the series of Theorem C, agree with the closed forms to 10⁻²¹ … 10⁻²⁶.
- **Theorem 17.** M₁ and M₀, computed by quadrature of the kernel integrals, agree with (60t/π)^{1/2}·L(2/5, 10t) and (60t/π)^{1/2}·L(1/5, 10t) to 10⁻²⁶. M₁ = 0.7247535755 at x = 0.3i, as in §5.
- **The four laws of Theorem 16, with M_j replaced by the Mordell integrals,** hold to 10⁻²⁵ … 10⁻²⁷. This uses only Ramanujan's series and the integral L, so it is independent of all the quoted identities.

On the literature point of §0.5: Garvan (*New fifth and seventh order mock theta function identities*) states that Watson studied the third-order functions under the modular group but could not find similar transformation properties for the fifth and seventh order. Andersen (Res. Number Theory 2016) quotes the fifth-order transformation law as Proposition 4.10 of Zwegers' thesis. The correction is therefore right. arXiv:2604.20034 was not looked up.

Labels: **proved** (a complete argument is written here; nobody else has checked it), **proved, given …** (the argument is complete once the named results from the literature are granted; these are listed in §4, each with the check made here), **verified** (numerics), **quoted**.

Notation: Im x > 0, b = √x, t = −2πix (Re t > 0), q = e^{2πix} = e^{−t}, p = e^{−πi/x} = e^{−2π²/t}, r = p^{1/2}; powers q^α = e^{2πixα}, p^α = e^{−πiα/x}. Φ_b, L, I(u), J(u), P, c, D(u), F(z), F̌, 𝒲, 𝒢₀, 𝒢₁ as in v12/v18. G, H are the Rogers–Ramanujan functions; s₁ = sin(π/5), s₂ = sin(2π/5); c₋ = (2/√5)·s₁·(π/t)^{1/2}, c₊ = (2/√5)·s₂·(π/t)^{1/2}, that is c± = (π(5 ± √5)/(10t))^{1/2}.

**Later remarks:** `addendum-v27.md` §5 (Theorem 17 on the unit circle, for real 0 < x ≤ 1).

---

## 0. Summary

1. **Proposition M (proved, given classical identities; §1).** At z = 1 and z = q the state integrals are pure powers times entries of the modular S-matrix of the Rogers–Ramanujan functions:

   I(0) = e^{−πi/4}·(2/√5)·s₂·e^{−πi(b² + b⁻²)/30},  I(b + 1/b) = −e^{−πi/4}·(2/√5)·s₂·e^{−19πi(b² + b⁻²)/30},

   I(1/b) = e^{−πi/4}·(2/√5)·s₁·e^{−πi(b² + 19b⁻²)/30},  I(b) = e^{−πi/4}·(2/√5)·s₁·e^{−πi(19b² + b⁻²)/30}.

2. **Theorem 16 (proved, given Proposition M; §2).** With M₁ := 2i·q^{11/60}·P·J(1/b) and M₀ := 2i·q^{29/60}·P·J(b + 1/b), the formula of the z-family (Theorem Z of v18) reads at z = 1 and z = q

   q^{11/60}·f₁(q) = M₁ + 2c₊·p^{−1/120}·(F₀(p) − 1) − 2c₋·p^{71/120}·F₁(p),
   2q^{11/60}·ψ₁(q) = −M₁ + 2c₊·p^{−1/120}·φ₀(−p) + 2c₋·p^{−49/120}·φ₁(−p),
   q^{−1/60}·f₀(q) = M₀ + 2c₋·p^{−1/120}·(F₀(p) − 1) + 2c₊·p^{71/120}·F₁(p),
   2q^{−1/60}·ψ₀(q) = −M₀ + 2c₋·p^{−1/120}·φ₀(−p) − 2c₊·p^{−49/120}·φ₁(−p).

3. **Theorem 17 (proved, given in addition the known transformation law of χ₀, χ₁ and two series identities that are checked to q³⁰⁰ but not proved here; §3).** The kernel integrals are Mordell integrals:

   M₁ = (60t/π)^{1/2}·L(2/5, 10t),  M₀ = (60t/π)^{1/2}·L(1/5, 10t),

   L(r, α) = ∫₀^∞ e^{−3αx²/2}·[cosh((3r − 2)αx) + cosh((3r − 1)αx)]/cosh(3αx/2) dx.

   So the four formulas of Theorem 16 are the classical transformation laws of Ramanujan's fifth-order functions f₀, f₁, ψ₀, ψ₁, term by term: state integrals ↔ S-matrix entries, dual series ↔ φ₀, φ₁, F₀, F₁, kernel integral ↔ Mordell integral.
4. **Numerically (verified):** everything above to 10⁻³⁰ at six points x; the series identities used to order q³⁰⁰; Theorem 17 also on the unit circle (x = 1/X, X = 5 … 40) to 10⁻¹² or better.
5. **A correction to v18.** v18 and the README spoke of "Watson's law" for the fifth-order functions. Watson (1937) proved Ramanujan's relations between these functions but did not obtain their transformation laws (so the introduction of Zwegers' thesis); those are in Zwegers' thesis (2002, Propositions 4.10 and 4.13) and, with explicit Mordell integrals, in the survey of Gordon and McIntosh (2012).
6. **Not done:** a proof of Theorem 17 inside this framework (a direct evaluation of J(1/b)); it would be a new proof of the fifth-order laws.

---

## 1. The state integrals at z = 1 and z = q

**Proposition M.** For Im x > 0, with ε = e^{−πi/4}:

  I(0) = ε·(2/√5)·s₂·q^{−1/60}·p^{1/30},  I(1/b) = ε·(2/√5)·s₁·q^{−1/60}·p^{19/30},

  I(b) = ε·(2/√5)·s₁·q^{−19/60}·p^{1/30},  I(b + 1/b) = −ε·(2/√5)·s₂·q^{−19/60}·p^{19/30}.

(In terms of b these are the formulas of §0, since q = e^{2πib²}, p = e^{−πi/b²}. The first one is the evaluation ∫_L e^{−πiy²}/Φ_b(y) dy = e^{−πi/4}·(2/√5)·sin(2π/5)·e^{−πi(b² + b⁻²)/30}, symmetric under b → 1/b as it must be. I did not look for it in the literature on state integrals.)

Proof.

*Step 1: bilinear forms.* By Theorem C of v12, I(u) = c·e^{−πiu²/2}·D(u), c = 2^{−1/2}·ε·q^{−1/24}·p^{1/12}, D(u) = Σ_{k,n≥0} (−1)^{kn+n}·e^{πiu(bk + n/b)}·q^{k²/4}·p^{n²/2+n}/((q;q)_k·(p²;p²)_n). At the four points e^{πiu(bk + n/b)} is 1, (−1)^k·p^{−n}, q^{k/2}·(−1)ⁿ, q^{k/2}·(−1)^{k+n}·p^{−n}. Summing over n for k even and k odd separately, with (S1), (S2) at the nomes ±r, and then over k with (S3)–(S6):

  (q;q²)∞·(−p;p)∞·D(0) = G(q⁴)·H(−r) + q^{1/4}·H(−q)·H(r),
  (q;q²)∞·(−p;p)∞·D(1/b) = G(q⁴)·G(−r) − q^{1/4}·H(−q)·G(r),
  (q;q²)∞·(−p;p)∞·D(b) = G(−q)·H(r) + q^{3/4}·H(q⁴)·H(−r),
  (q;q²)∞·(−p;p)∞·D(b + 1/b) = G(−q)·G(r) − q^{3/4}·H(q⁴)·G(−r).

(The first two are in v18 §3.1.)

*Step 2: transformation formulas.* Put 𝔤(τ) = e^{−2πiτ/60}·G(e^{2πiτ}), 𝔥(τ) = e^{22πiτ/60}·H(e^{2πiτ}), v = (𝔤, 𝔥)ᵀ, ζ = e^{2πi/60}. Classical (Q1): v(τ + 1) = diag(ζ⁻¹, ζ¹¹)·v(τ) and v(−1/τ) = 𝒮·v(τ), 𝒮 = (2/√5)·[[s₂, s₁], [s₁, −s₂]]. Let σ = −1/(4x), so r = e^{2πiσ}.

- −1/σ = 4x gives

  (T1) G(q⁴) = q^{1/15}·(2/√5)·[s₂·r^{−1/60}·G(r) + s₁·r^{11/60}·H(r)],
  (T2) H(q⁴) = q^{−11/15}·(2/√5)·[s₁·r^{−1/60}·G(r) − s₂·r^{11/60}·H(r)].

- σ + ½ is obtained from τ′ = x + ½ by τ′ → τ′ − 1 → −1/(τ′ − 1) → −1/(τ′ − 1) − 2 → its negative reciprocal (indeed (x − ½)/(2x) = ½ − 1/(4x)); all four points are in the upper half-plane. Applying the two rules of (Q1) along this chain: v(σ + ½) = 𝒮·T⁻²·𝒮·T⁻¹·v(x + ½) with T = diag(ζ⁻¹, ζ¹¹), hence v(x + ½) = T·𝒮·T²·𝒮·v(σ + ½), and the matrix product is (with 4s₁s₂ = √5)

  T·𝒮·T²·𝒮 = (2/√5)·[[s₁, s₂·ζ⁻⁶], [s₂·ζ⁶, −s₁]].

  Since 𝔤(x + ½) = e^{−πi/60}·q^{−1/60}·G(−q), 𝔥(x + ½) = e^{11πi/60}·q^{11/60}·H(−q), and the same with (σ, r), the phases cancel (ζ^{∓6}·e^{±12πi/60} = 1):

  (T3) G(−q) = q^{1/60}·(2/√5)·[s₁·r^{−1/60}·G(−r) + s₂·r^{11/60}·H(−r)],
  (T4) H(−q) = q^{−11/60}·(2/√5)·[s₂·r^{−1/60}·G(−r) − s₁·r^{11/60}·H(−r)].

*Step 3: the products collapse.* Insert (T1)–(T4) in the four bilinear forms. In each, the two terms with equal letters cancel (H(r)H(−r) in the first and third, G(r)G(−r) in the second and fourth; the powers of q agree because ¼ − 11/60 = 1/15 and ¾ − 11/15 = 1/60), and the rest is a multiple of G(r)H(−r) + G(−r)H(r) = 2/(r²;r⁴)∞² = 2/(p;p²)∞² (R5):

  numerators = (4/√5)·(p;p²)∞⁻²·{ s₂·q^{1/15}·p^{−1/120}, s₁·q^{1/15}·p^{11/120}, s₁·q^{1/60}·p^{−1/120}, s₂·q^{1/60}·p^{11/120} }.

*Step 4: eta.* (−p;p)∞ = 1/(p;p²)∞, and (q;q²)∞·(p;p²)∞ = √2·q^{1/24}·p^{1/24} (from η(−1/τ) = (−iτ)^{1/2}·η(τ): (q;q²)∞ = q^{1/24}·η(x)/η(2x) and (p;p²)∞ = p^{1/24}·η(−1/(2x))/η(−1/x) = √2·p^{1/24}·η(2x)/η(x)). So D = numerator·(p;p²)∞/(q;q²)∞ = (2√2/√5)·{ s₂·q^{1/40}·p^{−1/20}, s₁·q^{1/40}·p^{1/20}, s₁·q^{−1/40}·p^{−1/20}, s₂·q^{−1/40}·p^{1/20} }. Multiplying by c·e^{−πiu²/2}, where e^{−πiu²/2} = 1, p^{1/2}, q^{−1/4}, −q^{−1/4}·p^{1/2} at the four points, gives the statement. ∎

So at volume zero the two state integrals of the family stop being transcendental functions of t: they are the numbers (2/√5)·sin(π/5), (2/√5)·sin(2π/5) times pure powers of q and p. For comparison, at z = −1 they are (4/3)^{1/4}·e^{±iπ/12}·e^{∓iCl₂/t}·(1 + O(t)) with a full power series in t (v19).

---

## 2. The formula of the family at z = 1 and z = q

**Theorem 16.** For Im x > 0 the four identities of §0.2 hold, with M₁ = 2i·q^{11/60}·P·J(1/b), M₀ = 2i·q^{29/60}·P·J(b + 1/b).

Proof. Theorem Z of v18 (for every a ∈ ℂ) at a = 0 and at a = x. In both cases T = p, and 𝒲(1/p) = φ₀(−p), p·𝒲(p) = −φ₁(−p), 𝒢₀(p) = F₀(p), 𝒢₁(p) = F₁(p) (v18 §3.1).

- a = 0 (z = 1): ψ₁ = F(1) = i·P·[J(0) + I(0)·(φ₀(−p) − 1) + I(1/b)·φ₁(−p)/p] and f₁ = F̌(1) = −2i·P·[J(0) − I(0)·F₀(p) + I(1/b)·F₁(p)].
- a = x (z = q, e^{πia} = q^{1/2}, a/b = b): ψ₀(q) = q·F(q) = i·P·q^{1/2}·[J(b) + I(b)·(φ₀(−p) − 1) + I(b + 1/b)·φ₁(−p)/p] and f₀(q) = F̌(1/q) = −2i·P·q^{1/2}·[J(b) − I(b)·F₀(p) + I(b + 1/b)·F₁(p)].
- The kernels 1/(1 + e^{w}) and e^{w}/(1 + e^{w}) add up to 1, so J(u) + J(u + 1/b) = I(u). Replace J(0) by I(0) − J(1/b) and J(b) by I(b) − J(b + 1/b).
- P = 2^{−1/2}·q^{−1/6}·p^{−1/24}/b and 1/b = ε·(2π/t)^{1/2}, so by Proposition M (i·ε² = 1): 2i·P·I(0) = 2c₊·q^{−11/60}·p^{−1/120}, 2i·P·I(1/b) = 2c₋·q^{−11/60}·p^{71/120}, 2i·P·q^{1/2}·I(b) = 2c₋·q^{1/60}·p^{−1/120}, 2i·P·q^{1/2}·I(b + 1/b) = −2c₊·q^{1/60}·p^{71/120}. ∎

Remarks.

- **The mixing matrix.** The coefficients form (2/√5)·[[s₁, s₂], [s₂, −s₁]], the pattern of the matrix in Zwegers' Proposition 4.10 for (f₀, f₁) against (F₀ − 1, F₁) at the nome q^{1/2} of −1/τ, which is our p (quoted; I compared the pattern, not the normalisation). It appears here as the values of two state integrals.
- **Which integral.** Adding the laws: q^{11/60}·(f₁ + 2ψ₁) and q^{−1/60}·(f₀ + 2ψ₀) contain no integral. They are modular (the known theta-product formulas for f_j + 2ψ_j). In the formula of the family this is J(u) + J(u + 1/b) = I(u): the two choices a and a + 1 of the strip differ by a state integral, which at z = 1 is a pure power. For z ≠ 1 the corresponding difference I((a − 1)/b) is the exponentially large state integral of v22.

---

## 3. The kernel integral is a Mordell integral

**Theorem 17.** For Re t > 0: M₁ = (60t/π)^{1/2}·L(2/5, 10t) and M₀ = (60t/π)^{1/2}·L(1/5, 10t). Equivalently

  J(1/b) = e^{−πi/4}·(√15·t/π)·q^{−1/60}·p^{1/24}·L(2/5, 10t),  J(b + 1/b) = e^{−πi/4}·(√15·t/π)·q^{−19/60}·p^{1/24}·L(1/5, 10t).

Proof. Subtract the laws of Theorem 16 in pairs and use (R1), (R2) in the form χ₀(p) − 2 = 2(F₀(p) − 1) − φ₀(−p), p^{71/120}·χ₁(p) = 2p^{71/120}·F₁(p) + p^{−49/120}·φ₁(−p):

  q^{11/60}·(ψ₁ − f₁) = −(3/2)·M₁ − c₊·p^{−1/120}·(χ₀(p) − 2) + c₋·p^{71/120}·χ₁(p),
  q^{−1/60}·(ψ₀ − f₀) = −(3/2)·M₀ − c₋·p^{−1/120}·(χ₀(p) − 2) − c₊·p^{71/120}·χ₁(p).

By (R3), (R4) the left sides are q^{71/60}·χ₁(q²) and q^{−1/60}·(χ₀(q²) − 2). The transformation law (Q2) of χ₀, χ₁, taken at α = 2t (so its nome is q², its dual nome q₁⁴ is p, and its coefficients (π(5 ± √5)/(5α))^{1/2} are c±), says

  q^{−1/60}·(χ₀(q²) − 2) = −c₋·p^{−1/120}·(χ₀(p) − 2) − c₊·p^{71/120}·χ₁(p) − (135t/π)^{1/2}·L(1/5, 10t),
  q^{71/60}·χ₁(q²) = −c₊·p^{−1/120}·(χ₀(p) − 2) + c₋·p^{71/120}·χ₁(p) − (135t/π)^{1/2}·L(2/5, 10t).

Comparing: (3/2)·M_j = (135t/π)^{1/2}·L, and (2/3)·135^{1/2} = 60^{1/2}. The formulas for J follow with 1/(2i·P) = ½·ε·(t/π)^{1/2}·q^{1/6}·p^{1/24}. ∎

Remarks.

- **Read in the other direction:** given Theorem 16 and (R1)–(R4), Theorem 17 is equivalent to the law of χ₀, χ₁. A proof of Theorem 17 from the integral itself would prove the fifth-order laws from the formula of the family. I do not have one. In variables in which both kernels have their poles at multiples of i, the Gaussian of the Mordell integral is e^{−3πiv²/(10x)} and that of J is e^{−πixs²}; so the two integrals are related by something like a Fourier transform, in which 1/Φ_b has to supply the factor 5.
- **The trivial sector along the family.** P·J(a/b) is: at a = ½ the Borel sum of the perturbative series of ψ₁\* (v17); for 0 < a < 1 bounded on the unit circle by 1/(2√2·sin πa) (v23); at a = 1 a Mordell integral, and at a = 0 a Mordell integral plus a pure power. The kernel e^{2πas}/(1 + e^{2πs}) is the only thing that changes.
- **Another form of L.** L(r, α) = ∫₀^∞ e^{−3αx²/2}·2cosh(3(r − ½)αx)/(2cosh(αx) − 1) dx (from cosh A + cosh B and cosh 3y = cosh y·(2cosh 2y − 1)); the kernel has its poles at αx = ±πi/3 + 2πin. For M₁ this is (60t/π)^{1/2}·∫₀^∞ e^{−15tx²}·2cosh(3tx)/(2cosh(10tx) − 1) dx. I did not rewrite it in Zwegers' h(z; τ).
- **On the unit circle** both sides of Theorem 17 are defined (J over L with Faddeev's integral; L(r, 10t) on a turned ray). The identity holds there numerically (§5); its proof needs the continuity of J(1/b), J(b + 1/b) up to real x, which Lemma R3 of v20 gives for other values of u by an argument that applies here too (not written out).

---

## 4. What is taken from the literature

| | statement | status here |
|---|---|---|
| (Q1) | v(τ + 1), v(−1/τ) for the Rogers–Ramanujan pair 𝔤, 𝔥 | classical, quoted from memory; its consequences (T1)–(T4) are checked to 10⁻³⁰ at six points |
| (S1)–(S6) | Σ qⁿ²/(q⁴;q⁴)ₙ = G(q)/(−q²;q²)∞, Σ q^{n²+2n}/(q⁴;q⁴)ₙ = H(q)/(−q²;q²)∞, Σ q^{m²}/(q)₂ₘ = G(q⁴)/(q;q²)∞, Σ q^{m²+m}/(q)₂ₘ₊₁ = H(−q)/(q;q²)∞, Σ q^{m²+m}/(q)₂ₘ = G(−q)/(q;q²)∞, Σ q^{m²+2m}/(q)₂ₘ₊₁ = H(q⁴)/(q;q²)∞ | Rogers; in Slater's list (numbers not looked up); exact to q³⁰⁰ |
| (R5) | G(q)H(−q) + G(−q)H(q) = 2/(q²;q⁴)∞² | one of Ramanujan's forty identities; exact to q³⁰⁰ |
| (R1), (R2) | χ₀(q) = 2F₀(q) − φ₀(−q), q·χ₁(q) = 2q·F₁(q) + φ₁(−q) | Ramanujan, proved by Watson (1937); exact to q³⁰⁰ |
| (R3), (R4) | ψ₁(q) − f₁(q) = q·χ₁(q²), ψ₀(q) − f₀(q) + 2 = χ₀(q²) | exact to q³⁰⁰; they should follow from Ramanujan's relations in Watson (1937), which I did not look up, so within this addendum they are **verified**, not proved |
| (Q2) | the law of χ₀, χ₁ with L(1/5, 5α), L(2/5, 5α) | Gordon–McIntosh (2012), p. 118, quoted as printed in arXiv:2604.20034, eqs. (5), (6a), (6b); tested here on its own (q-series against the Mordell integrals) to 10⁻³⁰ at six points |

So: Proposition M and Theorem 16 are proved given (Q1), (S1)–(S6), (R5); Theorem 17 is proved given in addition (R1)–(R4) and (Q2). The laws of Theorem 16 with M_j replaced by the Mordell integrals are also verified directly, without any of these (§5, lines b and c together).

References: G. N. Watson, The mock theta functions (2), Proc. London Math. Soc. 42 (1937); S. Zwegers, Mock theta functions, thesis, Utrecht 2002, §4.4; B. Gordon, R. J. McIntosh, A survey of classical mock theta functions, in Partitions, q-series and modular forms, Springer 2012; arXiv:2604.20034, "On uniqueness of mock theta functions", §3.1.

---

## 5. Numerical checks (verified; `run_log.txt`)

Six points x = 0.3 + 0.2i, 0.11 + 0.03i, 0.2 + 0.5i, 0.3i, 0.05 + 0.4i, 0.6 + 0.25i; 30 digits.

| check | largest deviation |
|---|---|
| a. the four state integrals (series of Theorem C) against Proposition M, relative | 1.9·10⁻³⁰ |
| a. I(0) by quadrature against its closed form, relative | 2.1·10⁻³¹ |
| b. the four laws of Theorem 16 against the q-series | 2.2·10⁻³⁰ |
| c. M₁, M₀ against (60t/π)^{1/2}·L(2/5, 10t), L(1/5, 10t) | 2.2·10⁻³⁰ |
| c. J(1/b), J(b + 1/b) against the closed forms of Theorem 17, relative | 3.4·10⁻³⁰ |
| d. the quoted law (Q2), both lines | 6.0·10⁻³¹ |
| e. the four bilinear forms; (T1)–(T4); the matrix T𝒮T²𝒮; the eta identity | 5.9·10⁻³¹; 1.2·10⁻³⁰; 9.9·10⁻³¹; 2.0·10⁻³¹ |
| series identities (R1)–(R5), (S1)–(S6) | exact to q³⁰⁰ |
| Theorem 17 on the unit circle, x = 1/X, X = 5, 8, 20, 40 (double precision) | 1.9·10⁻¹⁴ … 4.5·10⁻¹² |

Sample values: at t = 0.6π (x = 0.3i), M₁ = 0.72475357548025069091, M₀ = 1.163932047808188582.

---

## 6. Open

1. A proof of Theorem 17 from the integral (see the first remark of §3).
2. The laws on the unit circle as theorems (continuity of J at u = 1/b, b + 1/b; then radial limits of ψ₁, f₁ at roots of unity from this formula).
3. General z: is P·J(a/b) expressible by Mordell-type integrals, and are the dual series 𝒲(T; p), 𝒢(T; p) Appell–Lerch sums in T? At T = p^{±1} they are.
4. z = q^m for other integers m (T = p again; three-term recursion in m).
5. Whether the evaluation of I(0) is in the state-integral literature.
6. The remaining items of v19 §8, v20 §5, v21 §5, v22 §7, v23 §4.

---

## 7. Scripts (`psi_star/v24/`)

- `s1_series.py [N]`: (R1)–(R5), (S1)–(S6) to order N = 300 (7 seconds).
- `s2_law.py [a b c d e] [digits]`: the checks a–e of §5 (2.5 minutes).
- `s3_real.py`: Theorem 17 at real x (1 minute).
- `common.py`, `a1_formula.py`, `c2_z1.py` (from v18), `n1_zroots.py` (from v22), `fad.py` (from v16), `run_all.sh`, `run_log.txt`, `README.md`.
