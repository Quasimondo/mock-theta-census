# Mock theta search: addendum to handoff v1

Session of 21 Nov 2026 (Claude Fable 5.2 with Mario Klingemann), working only from `mock-theta-handoff-v1.md`. The v4 archive was not available, so everything here was recomputed independently. Scripts are in `mock-theta-addendum-scripts.zip`.

Status labels: **proved** = short proof given here, also machine-checked; **verified** = exact computation, no proof; **inference** = my reading of the evidence.

---

## 0. Summary

1. **H(x) = (1−x²)·g₂(x,q) + 1 − x is proved.** It telescopes (§1).
2. **Two new exact identities for ψ₀\*, ψ₁\*** (§2, proved, checked to q³⁰⁰). They tie the pair to the sign-twisted f₀, f₁ and to the Rogers–Ramanujan series at z = −1:
   - 2ψ₀\*(q) = f₀\*(q) − θ₄(q)·R₀(q)
   - 2ψ₁\*(q) = f₁\*(q) + θ₄(q)·R₁(q)
3. **The "missing components" are R₀, R₁** (§3, verified to 18 digits). The exact asymptotic series of ψ₀\*, ψ₁\* at q = −1 has integer coefficients, and its large-order behaviour is governed by R₀, R₁, including exact constants in Q(√5).
4. **Consequence (inference):** (ψ₁\*, ψ₀\*) is not a classical weight-½ vector-valued quantum modular form with constant multipliers. The pair law with weight ≈ 0.508 in v1 is an effective fit to something of a different kind: a resurgent object of Nahm-sum type (Garoufalidis–Zagier setting). Hypotheses (a)–(d) of v1 §7 should be replaced by this.

---

## 1. Proof of the H(x) identity

Let hₙ = q^{n(n−1)/2}(−q;q)ₙ / ((xq;q)ₙ(q/x;q)ₙ) and tₙ = (−q;q)ₙ q^{n(n+1)/2} / ((x;q)ₙ₊₁(q/x;q)ₙ₊₁), so H = Σhₙ and g₂ = Σtₙ.

Put Uₙ = (1 − xqⁿ)·hₙ. Then

  Uₙ − Uₙ₊₁ = hₙ − (1 − x²)·tₙ.

Check: tₙ = qⁿhₙ / ((1−x)(1−qⁿ⁺¹/x)) and hₙ₊₁/hₙ = xqⁿ(1+qⁿ⁺¹) / ((1−xqⁿ⁺¹)(x−qⁿ⁺¹)); both sides reduce to hₙ·(x − qⁿ⁺¹ − xqⁿ − x²qⁿ)/(x − qⁿ⁺¹).

Summing over n ≥ 0 with U₀ = 1 − x and Uₙ → 0 gives H(x) − (1−x²)g₂(x) = 1 − x. ∎

The step was also checked symbolically for n = 0..5 (`s11_checks.py`). The identity is elementary; it is probably known or implicit in Gordon–McIntosh, but I did not find it stated.

---

## 2. New identities for ψ₀\*, ψ₁\*

Notation:
- f₀\*(q) = Σ (−1)ⁿ q^{n²}/(−q;q)ₙ, f₁\*(q) = Σ (−1)ⁿ q^{n²+n}/(−q;q)ₙ (sign-twisted fifth-order f₀, f₁)
- R₀(q) = Σ (−1)ⁿ q^{n²}/(q;q)ₙ, R₁(q) = Σ (−1)ⁿ q^{n²+n}/(q;q)ₙ (Rogers–Ramanujan series Σ zⁿq^{n²}/(q)ₙ at z = −1 and z = −q)
- θ₄(q) = Σ (−1)ⁿ q^{n²} = (q;q)∞/(−q;q)∞

**Identities (proved; verified to q³⁰⁰ in `s1_identities.py`):**

  (I) 2ψ₀\*(q) = f₀\*(q) − θ₄(q)·R₀(q)
  (II) 2ψ₁\*(q) = f₁\*(q) + θ₄(q)·R₁(q)

**Proof sketch.** Extend the sum to all n ∈ ℤ using (−q)₋ₘ = q^{m(m−1)/2}/(2(−q)ₘ₋₁). The negative half of Σ(−1)ⁿq^{n(n+1)/2}(−q)ₙ is −½f₁\*, so the bilateral sum is B = ψ₁\* − ½f₁\*. Write (−q)ₙ = (−q)∞/(−qⁿ⁺¹)∞ and use the partial fraction 1/(−xq)∞ = (1/(q)∞) Σⱼ (−1)ʲ q^{j(j+1)/2} / ((q)ⱼ(1 + xq^{j+1})), valid for all x. With A(w) = Σₙ wⁿ q^{n(n−1)/2}/(1+qⁿ) and Θ(w) = Σₙ wⁿ q^{n(n−1)/2}, the general bilateral sum is

  Σₙ zⁿ q^{n(n+1)/2}(−q)ₙ = ((−q)∞/(q)∞) · { A(z)·z⁻¹ Σⱼ q^{j(j+1)} z⁻ʲ/(q)ⱼ + Θ(z)·(a Hecke-type double sum) }.

At z = −1 two things happen: Θ(−1) = 0, and A(−1) is an Appell sum at a point where its theta function vanishes, so it is modular: A(−1) = −½θ₄². This gives B = ½θ₄R₁, which is (II). The case z = −q gives (I). ∎

**Context.** For the unsigned functions the bilateral series are known to be modular (e.g. Σ_{n∈ℤ} q^{n²}/(−q)ₙ = f₀ + 2ψ₀; see "Bilateral series and Ramanujan's radial limits", arXiv:2202.12141, and Mc Laughlin, arXiv:1906.11997). (I) and (II) are the sign-twisted analogues: the bilateral series is no longer modular, but equals θ₄ times a Rogers–Ramanujan series at z = −1. At z = +1 the same formula has the modular Rogers–Ramanujan H(q) in that slot.

**Double-sum forms (verified to q³⁰⁰, `s11_checks.py`):**

  (−q)∞·ψ₁\* = Σ_{k even, n≥0} (−1)ⁿ q^{k(k+1)/2 + n² + n + kn}/(q)ₖ
  (−q)∞·ψ₀\* = Σ_{k odd, n≥0} (−1)ⁿ q^{k(k+1)/2 + n² + kn}/(q)ₖ

Equivalently (−q)∞·ψ₁\* = Σ_{j≥m≥0} (−1)^{m+j} q^{m²+j²+j}/(q)₂ₘ: a Nahm-type sum weighted by tails of a false theta function. This explains why no eta-quotient multiplier makes ψ\* sparse (v1 §7, fact 4): the Bailey pair with βₙ = (−1)ⁿ has no short αₙ (`s5` check, not included).

---

## 3. Exact asymptotics at q = −1 and what governs them

**Perturbative series (verified, exact).** At q = −e^{−t} each term of ψ\* is O(t^{⌈n/2⌉}), so the expansion is a well-defined formal power series, computed in exact rational arithmetic to 120 terms (`s2_expand.py`):

  ψ(−e^{−t}) = Σₖ cₖ tᵏ/k!

- ψ₁\*: cₖ = 1, −1, 27, −1237, 106587, −14667181, 2954135907, −819380500597, …
- ψ₀\*: cₖ = −1, 4, −74, 3370, −286886, 39242914, −7873061654, 2177801727490, …

All 121 computed cₖ are integers. (Worth an OEIS search on the absolute values.)

**Large-order behaviour (verified numerically, constants identified by PSLQ, `s6`, `s7`).** With A = π²/60 and aₖ = cₖ/k!:

  aₖ ≈ K·(−1)ᵏ Γ(k+½) A^{−k−½} · [1 + s₁A/(k−½) + s₂A²/((k−½)(k−3⁄2)) + …]

| | K²·π | s₁ | s₂ |
|---|---|---|---|
| ψ₀\* (K < 0) | (5+3√5)/40 | (13−6√5)/60 | (781−300√5)/7200 |
| ψ₁\* (K > 0) | (3√5−5)/40 | −(23+6√5)/60 | (1141+420√5)/7200 |

So the leading Borel singularity sits at −π²/60 (the "1/60" of the fifth-order world), with exponent ½, and the series attached to it, 1 + s₁t + s₂t² + …, is not a pure exponential. Its coefficients keep growing (0.0153, −0.0223, 0.0368, −0.0799, 0.2098 for ψ₀\*).

**The attached series is R₀, R₁ (verified, `s10_R.py`).** Direct evaluation at 600–700 digits, t = 1/60 … 1/83:

  Rᵢ(−e^{−t}) ~ ρᵢ · e^{π²/(60t)} · (1 − s₁t + s₂t² − …)

with ρ₀² = (5+3√5)/10, ρ₁² = (3√5−5)/10, i.e. **ρᵢ = 2√π·|Kᵢ| exactly**, and the same s₁ (18 digits) and s₂ (15 digits) as in the table. The signs of K match the signs of ∓½θ₄Rᵢ in (I), (II), since θ₄(−e^{−t}) ≈ √(π/t).

In words: the divergence of the perturbative series of ψᵢ\* at q = −1 is exactly the term ∓½θ₄Rᵢ of identities (I), (II). The two sides of the identity are the resurgent "perturbative" and "non-perturbative" sectors of one object.

---

## 4. Interpretation (inference)

- For a classical quantum modular form (a vector of q-series with constant multipliers and an analytic cocycle, e.g. false theta functions or the Kontsevich–Zagier function), the series attached to each Borel singularity is a pure exponential e^{αt}. Here it is a genuinely divergent series with coefficients in Q(√5). So **no finite-dimensional law of the kind fitted in v1 §7 can be exact**, whatever the weight, and adding log terms or a third q-series component with constant multipliers will not fix it. The fit to 10⁻⁹ is real but effective: it approximates a matrix cocycle whose entries are non-trivial smooth functions.
- The right framework is the one for Nahm sums and knot invariants: radial asymptotics at each root of unity are finite sums of e^{V/t}·Φ(t) with divergent Φ, and modularity holds for a matrix of such objects (Garoufalidis–Zagier, "Asymptotics of Nahm sums at roots of unity", arXiv:1812.07690; their "Knots, perturbative series and quantum modularity"; Wheeler's thesis on modular q-difference equations). ψ₀\*, ψ₁\* play the role of the Kashaev-like row (terminating at even roots, integer perturbative coefficients); R₀, R₁ and f₀\*, f₁\* are the q-series companions.
- Structural reason: F(z) = Σ zⁿq^{n(n+1)/2}(−q)ₙ satisfies F(z) = 1 + zq·F(zq) + zq²·F(zq²), an inhomogeneous second-order q-difference equation. ψ₁\* = F(−1), ψ₀\* = q·F(−q), and the v1 relation 1 − ψ₁\* − ψ₀\* = F(−1/q) is this equation at z = −1/q. A rank-2 equation plus inhomogeneity suggests a 3×3 matrix-valued object with first row (1, ψ₀\*, ψ₁\*), as for the Kashaev invariant and its descendants.
- ψ\* being bounded at every root is then a cancellation between two exponentially large Nahm-type pieces (f\* and θ₄R), not false-theta behaviour.
- c8438 and c9542 (the halves ½ψ + ½ψ\*) are therefore mock only in Ramanujan's sense, as v1 suspected: classical mock theta function plus a bounded non-mock piece.

Not checked: the same analysis at other roots (q → 1, i, …); whether ψ\*(−e^{−t}) equals the Borel sum of its series; the exact 3×3 cocycle.

---

## 5. Suggested next steps

1. **Repeat §3 at q = i and at q → 1** (odd roots need numerical radial limits). Predict the exponents from the twisted Nahm equation for R(−1) and compare.
2. **Test the 3×3 matrix law**: columns built from (1, ψ₀\*, ψ₁\*), the homogeneous solutions θ₄R₀, θ₄R₁ (and f₀\*, f₁\*), under x ↦ x/(2x+1) and x ↦ x+1, with a cocycle required only to be smooth, not of the form ε(cx+d)^{−k}.
3. **Literature**: look for R(−1) = Σ(−1)ⁿq^{n²}/(q)ₙ in the Garoufalidis–Zagier/Wheeler line and in Andrews' Lost Notebook papers on the "entire Rogers–Ramanujan function"; look for (I), (II) in the bilateral-series papers above. My search was brief and found neither ψ\* nor (I), (II).
4. **Who to ask**: this now points less to the mock-theta specialists and more to Garoufalidis, Zagier, Campbell Wheeler (resurgence and Nahm sums), or Lovejoy/Osburn (quantum q-series identities).
5. **Classifier**: R₀, R₁ have rational κ = 1/60 at q = −1 but no exact theta law. The census could flag this class explicitly ("rational exponent, divergent series"), since it is where ψ\*-like objects come from.

---

## 6. Scripts (`mock-theta-addendum-scripts.zip`)

- `s1_identities.py N`: exact series; checks (I), (II) to q^N.
- `s2_expand.py K`: exact expansion at q = −e^{−t} to order K → `exp_m1.json`.
- `s6_largeorder.py`, `s7_id.py`: large-order fit and PSLQ identification (run `s2_expand.py 120` first).
- `s10_R.py`: asymptotics of R₀, R₁ at q = −e^{−t} and comparison with the table.
- `s11_checks.py`: telescoping step for H(x); double-sum forms.

Requirements: `pip install python-flint mpmath sympy`. Paths inside the scripts point to `/home/claude/work`; adjust before running.
