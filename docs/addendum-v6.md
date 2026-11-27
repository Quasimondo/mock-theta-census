# Mock theta search: addendum v6 (towards a proof: the state-integral form of the law)

Session of 26 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). This explores the proof route named at the end of `addendum-v5.md`: the duality between R₀, R₁ and P₊, P₋. Scripts are in `psi_star/v6/`; `run_log.txt` there is their complete output.

**Independent check** (Claude Opus, same day, `psi_star/v6/check_casoratian.py`, `check_B1.py`): the reformulations N₀ = (−q;q)∞R₀, N₁ = (−q;q)∞R₁, E± = (q;q²)∞P± and both Casoratian identities X·N₀ + N₁·Y = 1, E₊′E₋ + E₋′E₊ = 2 hold exactly to q³²⁹; the proof of the lemma was read and is correct. For 𝓑₁ with the stated c and first correction, the ratio to the predicted law was computed at three dual nomes (|p| = 0.017, 0.043, 0.111) and |𝒳| ≈ 200 … 1,600: the deviation is the same at all three to three digits and falls by a factor 4 per doubling of |𝒳| (3.6·10⁻⁵ → 5.8·10⁻⁷), i.e. like t², as it must if the law is p-independent and both constants are right. 𝓑₂, 𝓑₃, 𝓑₄ were not rechecked.

Labels: **proved** (proof given, machine-checked), **verified** (numerics, no proof), **inference**.

---

## 0. Summary

1. **The law is now an exact identity plus one analytic claim** (§2). There are four explicit bilinear combinations 𝓑₁…𝓑₄ of q-series in q and in the dual nome. Solving for R₀ gives, exactly,
   R₀(q) = 2^{−1/2}·e^{−π²/(12t) − t/24}·[𝓑₃·P₊(p) + 𝓑₁·P₋(p)].
   The whole content of the modularity law is that 𝓑₁…𝓑₄ **do not depend on the dual nome p** in their asymptotics. That is verified to 10⁻¹⁷ at two orders, with closed-form constants.
2. **The "up to W" freedom of v5 is gone.** The coefficient of the decaying sector is a single function 𝓑₁ of t; nothing depends on a choice of resummation any more.
3. **Both sides are one function** (§1, proved). R₀, R₁, P₊, P₋ are all special values of H(z;q) = Σ zᵐq^{m²/4}/(q;q)ₘ. Its second solutions are the sums X, Y of v3's identities (III), (IV), so the bilateral parts of ψ\* and of F\* together form one fundamental matrix. Two exact Casoratian identities come out: **X·N₀ + N₁·Y = 1** and **E₊′E₋ + E₋′E₊ = 2**.
4. **A proof route exists in the literature** (§3, inference). The 𝓑ⱼ have the form of factorised one-dimensional state integrals of Faddeev's quantum dilogarithm, and the equation of H is a q-hypergeometric equation of the kind Garoufalidis–Wheeler prove to be modular. I have not checked their hypotheses against this case.

---

## 1. One function behind both sides (proved, exact to q³⁰⁹; `exact_checks.py`)

Let 𝓡(z;q) = Σ zⁿq^{n²}/(q;q)ₙ, so R₀ = 𝓡(−1;q), R₁ = 𝓡(−q;q), P± (q) = 𝓡(q^{±1};q²).

**Lemma.** With Θ(z) = Σ_{n∈ℤ} zⁿq^{n²},

  (q;q)∞·𝓡(z;q) = Θ(z)·Σⱼ z^{−j}q^{j²+j}/(q;q)₂ⱼ − q·Θ(zq)·Σⱼ z^{−j}q^{j²+2j}/(q;q)₂ⱼ₊₁.

*Proof.* Write 1/(q)ₙ = (q^{n+1};q)∞/(q)∞ = (1/(q)∞)·Σₘ(−1)ᵐq^{m(m+1)/2+nm}/(q)ₘ. The factor (q^{n+1})∞ vanishes for n < 0, so the sum over n can run over ℤ and equals Θ(zqᵐ). Then Θ(zq^{2j}) = z^{−j}q^{−j²}Θ(z) and Θ(zq^{2j+1}) = z^{−j}q^{−j²−j}Θ(zq). ∎

At z = −1, −q one of the two theta factors vanishes and the other is θ₄. At z = q^{±1} in base q² both are multiples of ψ(q). This gives

  (−q;q)∞·R₀ = **N₀** := Σ(−1)ⁿq^{n²+n}/(q;q)₂ₙ,  (−q;q)∞·R₁ = **N₁** := Σ(−1)ⁿq^{n²+n}/(q;q)₂ₙ₊₁,

  (q;q²)∞·P₊ = **E₊** := Σ(−1)ᵐq^{m(m+1)/2}/(q²;q²)ₘ,  (q;q²)∞·P₋ = **E₋** := Σ(−1)ᵐq^{m(m+3)/2}/(q²;q²)ₘ.

All four are values of **H(z;q) = Σₘ zᵐq^{m²/4}/(q;q)ₘ**, the Nahm sum with A = ½:

- N₀ and i·q^{−1/4}·Y are the even and odd parts of H(i·q^{1/2}; q); X and i·q^{1/4}·N₁ those of H(i; q). Here X = Σ(−1)ⁿq^{n²}/(q)₂ₙ and Y = Σ(−1)ⁿq^{(n+1)²}/(q)₂ₙ₊₁ are the sums of identities (III), (IV) in v3.
- E₊ = H(−Q^{1/4}; Q), E₋ = H(−Q^{3/4}; Q) with Q = q². Their partners are E₊′ = H(Q^{1/4}; Q), E₋′ = H(Q^{3/4}; Q), the same sums without the sign.

**Equation and Casoratian.** H(z) − H(zq) = z·q^{1/4}·H(zq^{1/2}). If ε(zq^{1/2}) = −ε(z), then ε(z)H(−z) solves the same equation, and the Casoratian of the two solutions is invariant under z ↦ zq^{1/2} and analytic on ℂ\*, hence constant:

  H(z)·H(−zq^{1/2}) + H(zq^{1/2})·H(−z) = 2 for all z.

At z = Q^{1/4} and z = i this reads

  **E₊′·E₋ + E₋′·E₊ = 2,  X·N₀ + N₁·Y = 1.**

H is also self-dual under q ↦ 1/q: H(z; 1/q) = H(−zq^{1/2}; q). So the family makes sense on both sides of the unit circle, which is what a holomorphic quantum modular form needs.

---

## 2. The four bilinear combinations (verified; `hfam.py`, `e5_state.py`)

Setting as in v5: X-plane variable 𝒳 = 𝒳₀ + 2m − ib, q = e^{2πi/𝒳}, t = −2πi/𝒳, dual nome p = e^{−πi𝒳} (called q̃ in v5). Define

  𝓑₁ = N₀(q)·E₊′(p) + i·q^{−1/4}·Y(q)·E₊(p),   𝓑₃ = N₀(q)·E₋′(p) − i·q^{−1/4}·Y(q)·E₋(p),

  𝓑₂ = X(q)·E₊(p) + i·q^{1/4}·N₁(q)·E₊′(p),   𝓑₄ = X(q)·E₋(p) − i·q^{1/4}·N₁(q)·E₋′(p).

By the Casoratian these invert **exactly**:

  N₀ = ½(𝓑₁E₋ + 𝓑₃E₊),  i·q^{−1/4}·Y = ½(𝓑₁E₋′ − 𝓑₃E₊′),  X = ½(𝓑₂E₋′ + 𝓑₄E₊′),  i·q^{1/4}·N₁ = ½(𝓑₂E₋ − 𝓑₄E₊).

With (−q;q)∞ = 2^{−1/2}·e^{π²/(12t) + t/24}·(p;p²)∞ (exact), the first one is the formula for R₀ in the summary; the last one gives R₁ in the same way.

**Result.** Along the horocycle, with Cl₂ = Cl₂(π/3) and c = (4/3)^{1/4}·e^{−iπ/12}:

  𝓑₁ = c·e^{+i·Cl₂/t}·(1 + (−1/24 − i/(18√3))·t + …),   𝓑₂ = c̄·e^{+i·Cl₂/t}·(1 + (+1/24 − i/(18√3))·t + …),

  𝓑₃ = c̄·e^{−i·Cl₂/t}·(1 + (−1/24 + i/(18√3))·t + …),   𝓑₄ = c·e^{−i·Cl₂/t}·(1 + (+1/24 + i/(18√3))·t + …),

**with no dependence on p**. Evidence:

- 𝓑₁ and 𝓑₂ are the decaying ones. Their two terms are each of size 10²⁸ at K = 400 and the sum is 8.9·10⁻²⁹: the growing sector cancels over 56 orders of magnitude. Only the sign choice shown does this; the other three unit factors do not (`hfam.py`).
- Limits at six values of p on two circles: the constants agree with each other to 10⁻¹⁹ (𝓑₁, 𝓑₂) and 10⁻¹⁷ (𝓑₃, 𝓑₄), and with c = 1.03795484930204250 − 0.27811916365044995·i to 17 digits. Extrapolation error 2·10⁻¹⁶.
- The coefficient of t is also the same at all six points (15 digits), and its ratio to the leading constant matches ±1/24 ± i/(18√3) as shown to 11 digits in all four cases.

**Consistency with v5.** Inserting these into the exact formula gives R₀ ≃ 3^{−1/4}e^{iπ/12}·e^{−(π²/12 + i·Cl₂)/t}·P₊(p) + 3^{−1/4}e^{−iπ/12}·e^{−(π²/12 − i·Cl₂)/t}·P₋(p), which is v5 §6 with both constants reproduced. In v5 the P₋ term was known only up to adding a multiple of P₊; here its coefficient is the single function 𝓑₁(t), and the coefficient of P₊ is 𝓑₃(t).

**What remains to prove.** Exactly this: each 𝓑ⱼ, as a function of x = 1/𝒳 in the upper half-plane, extends analytically across the positive real axis. Then its asymptotics at x → 0 cannot depend on p, and the laws of v5 for R₀, R₁, hence for the bilateral parts θ₄R, follow from the exact inversion formulas.

---

## 3. Why this should be provable (inference)

- **Shape of the 𝓑ⱼ.** A one-dimensional state integral ∫Φ_b(y)·e^{−πi·A·y² + 2πuy}dy, evaluated by residues at the poles of Faddeev's quantum dilogarithm, gives a double sum over (m, n) that factorises into a q-series times a q̃-series. For A = ½ the cross term is (−1)^{mn}, so the sum splits by the parity of m: (even part of H)·H̃(z̃) + (odd part of H)·H̃(−z̃). That is precisely the form of 𝓑₁…𝓑₄, including the pairing of N₀ with E′ and of Y with E. The saddle of such an integral has action ±i·Cl₂(π/3), as observed. I derived this form on paper and then found the cancellation numerically; I did not evaluate the integral itself.
- **Literature.** Garoufalidis–Kashaev, "From state integrals to q-series" (Math. Res. Lett. 2017, arXiv:1304.2705) prove the factorisation for the family ∫Φ_b(x)^B·e^{−Aπix²}dx with integer A, B. Garoufalidis–Wheeler, "Modular q-holonomic modules" (arXiv:2203.17029) state that the generalised q-hypergeometric equation is modular, meaning its fundamental matrices define a cocycle with improved analyticity. Wheeler's thesis (Bonn 2023) has sections on dualities between modules of Nahm sums and on deformations of the Rogers–Ramanujan identities. I read abstracts and tables of contents only. The case here is A = ½ with the deformation variable at z = i, i·q^{1/2}; whether it is literally covered needs checking.
- **What a proof would still need for ψ\* itself.** The bilateral parts are θ₄R (for ψ\*) and, by (III), (IV), the X, Y sums (for F\*). The one-sided halves f\* and F_a, F_b are not solutions of the homogeneous equation; their laws (v5 §2, §4) need a separate argument, probably the inhomogeneous version of the same module.

---

## 4. Next steps

1. Check the hypotheses of Garoufalidis–Wheeler against the equation H(z) − H(zq) = zq^{1/4}H(zq^{1/2}), or write the A = ½ state integral explicitly and compare it numerically with 𝓑₁.
2. Test analyticity directly: evaluate 𝓑₁ just below the real axis, using H(z;1/q) = H(−zq^{1/2};q), and compare with the values just above.
3. The same construction on the converse side (P± near 1 against R₀, R₁ at the dual nome) and at q = −1.
4. Repeat for F₀\*, F₁\*: by §1 their bilateral parts are the second row of the same matrix, so their law should come for free. Not tested.

---

## 5. Scripts (`psi_star/v6/`)

| script | section | what it does |
|---|---|---|
| `exact_checks.py` | §1 | the four reformulations and the two Casoratian identities, exact to q³⁰⁹ |
| `hfam.py` | §1, §2 | evaluators for N₀, N₁, X, Y and E±, E±′; which combinations cancel the growing sector |
| `e5_state.py` | §2 | limits and first corrections of 𝓑₁…𝓑₄ at six points |
