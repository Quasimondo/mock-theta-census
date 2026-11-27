# Mock theta search: addendum v28 (what kind of object is f\*? Not mock modular; a one-sided Nahm sum)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). The first open question of the README: can f₁\*, f₀\* (and with them ψ₁\*, ψ₀\*) be written as Nahm sums, generalised Nahm sums or indefinite theta series? Scripts and output are in `psi_star/v28/`.

Labels: **proved** (a complete argument is written here or in the addenda cited; nobody else has checked it), **proved, given …** (complete once the named standard facts are granted), **verified** (numerics), **interpretation** (a way of reading the results, not a theorem), **from memory** (literature not rechecked).

---

## 0. Summary

1. **Not mock modular (Theorem 19; proved, given standard growth bounds; §1).** Along every ray arg x = ε near the unit circle, ψ₁\*, ψ₀\* grow like exp((Cl₂·cos ε − (π²/3)·sin ε)/|t|) and f₁\*, f₀\* like exp((Cl₂·cos ε − (7π²/3)·sin ε)/|t|) (from Theorem 8). A modular form, a mock modular form, a false theta function, or any finite sum of products of such, grows at most like exp(D·sin ε/|x|): its exponent vanishes with ε. Since Cl₂ > 0, none of the four functions is of that kind, for any weight and level. In particular none is an indefinite theta series of Zwegers' type divided by a modular form. The argument uses only Cl₂ ≠ 0.
2. **A one-sided Nahm sum (Proposition N; proved; §2).** With M = (−q;q)∞ and Q₁(n, k) = n² + n + nk + k(k+1)/2:

   M·f₁\* = Σ_{n,k≥0} (−1)ⁿ·q^{Q₁}/(q;q)ₖ,  M·θ₄R₁ = Σ_{n,k≥0} (−1)^{n+k}·q^{Q₁}/(q;q)ₖ,  M·ψ₁\* = the terms with k even.

   So f₁\* and θ₄R₁ are the two characters of one double sum, over a positive definite form of determinant 1, in which only the variable k carries a Pochhammer symbol. Identity (II), ψ₁\* = ½(f₁\* + θ₄R₁), is the splitting by the parity of k. If n runs over ℤ instead of n ≥ 0, the sum is θ₄ times a rank-one Nahm sum. So f\* is a Nahm-type sum restricted to a half-lattice, as a partial theta function is a theta function restricted to n ≥ 0.
3. **A generalised Nahm sum with an indefinite form (proved; §2).** ψ₁\* = (q²;q²)∞·Σ_{n,j≥0} (−1)ⁿ·q^{n(n+1)/2 + 2nj + 2j}/((q;q)ₙ(q²;q²)ⱼ). Its Nahm equations reduce to X² + X + 1 = 0, the tetrahedron equation of ψ₁\*. This is the wide sense of "Nahm sum" used for knot invariants, not the positive definite one.
4. **The picture (interpretation; §3).** In the variable z, F(z) and −F̌(1/z)/(2z) both solve the inhomogeneous equation y(z) − zq·y(zq) − zq²·y(zq²) = 1, and their difference, the bilateral sum, solves the homogeneous one and is a theta function times a rank-one Nahm sum. At z = 1 this is: fifth-order mock theta functions against Rogers–Ramanujan functions. At z = −1: ψ₁\*, f₁\* against the Nahm sum R₁ with volume Cl₂. So f\* is to the non-modular Nahm sum R what a mock theta function is to a modular form.
5. **Answer to the question.** f\* is not an indefinite theta series and not a Nahm sum of the family of R; it is the one-sided companion of that Nahm sum. Objects of this kind are known for knots (the series of the trivial connection, from memory); the instance, and its continuous link to Ramanujan's mock theta functions through z, I have not seen. Not excluded by a proof: a positive definite Nahm sum of higher rank equal to f\*.

---

## 1. Growth at the cusp: f\* and ψ\* are not mock modular

Notation: x in the upper half-plane, q = e^{2πix}, t = −2πix, so |t| = 2π|x| and Im(−1/x) = sin(arg x)/|x|.

**Definition.** A function f holomorphic in the upper half-plane has *modular growth at 0* if there is D ≥ 0 such that for every ε ∈ (0°, 90°) there are C, B with

  |f(x)| ≤ C·|x|^{−B}·exp(D·sin(arg x)/|x|)  for ε ≤ arg x ≤ 180° − ε, |x| ≤ 1.

**Lemma (standard).** The following have modular growth at 0, and the class is closed under sums and products:

- (a) weakly holomorphic modular forms of any real weight, for any subgroup of finite index in SL₂(ℤ), with a multiplier system of modulus 1;
- (b) mock modular forms: holomorphic parts f = F − f⁻ of harmonic weak Maass forms F whose non-holomorphic part f⁻ = Σ_{n>0} c(n)·Γ(1 − k, 4πny/h)·e^{−2πinx/h} has coefficients c(n) = O(n^{A}) (shadow holomorphic at the cusps);
- (c) partial and false theta functions Σ_{n≥0} χ(n)·n^{ν}·q^{(n+α)²/m}, χ bounded.

Sketch. (a) g(x′) = (x′)^{−k}·f(−1/x′) is weakly holomorphic for the conjugate group, with an expansion Σ_{n≥−n₀} aₙ·e^{2πi(n+κ)x′/h}; so |g(x′)| ≤ C·e^{2πn₀·Im x′/h} for Im x′ ≥ 1, and |f(x)| = |x|^{−k}·|g(−1/x)|. The part of the region with Im(−1/x) ≤ 1 is compact in the upper half-plane. (b) The same for F, by the growth condition in the definition of a harmonic weak Maass form. For f⁻: the n-th term has modulus Γ(1 − k, u)·e^{u/2} with u = 4πny/h, which is at most C·(1 + u^{−B₀})·e^{−u/2} for some B₀; with c(n) = O(n^{A}) the sum is O(y^{−B}) for 0 < y ≤ 1, and y ≥ |x|·sin ε in the cone. (c) |f| ≤ Σ n^{ν}·e^{−2πy(n+α)²/m} ≤ C·y^{−(ν+1)/2}. D is the order of the pole at the cusp (times 2π/h) and does not depend on ε. ∎

**Theorem 19.** Let λ₁(ε) = Cl₂·cos ε − (π²/3)·sin ε and λ₂(ε) = Cl₂·cos ε − (7π²/3)·sin ε, with Cl₂ = Cl₂(π/3) = 1.0149…

- (a) For 0 < ε ≤ 10° and x → 0 along arg x = ε: |ψ₁\*(q)| = A·|t|^{−1/2}·e^{λ₁(ε)/|t|}·(1 + o(1)), and the same for ψ₀\*.
- (b) For 0 < ε ≤ 2.5°: |f₁\*(q)| = 2A·|t|^{−1/2}·e^{λ₂(ε)/|t|}·(1 + o(1)), and the same for f₀\*.
- (c) None of ψ₁\*, ψ₀\*, f₁\*, f₀\* has modular growth at 0. So none of them is a weakly holomorphic modular form, a mock modular form in the sense of the Lemma, a false theta function, or a finite sum of products of these (a mixed mock modular form), of any weight and level.

Proof. (a), (b): Theorem 8 of v19, valid for 0 < arg x ≤ 10°. On the ray, 1/t = e^{i(90° − ε)}/|t|, so with S̄ = π²/3 + i·Cl₂ and S = π²/3 − i·Cl₂: |e^{−S̄/t}| = e^{λ₁(ε)/|t|} and |e^{−S/t}| = e^{−(Cl₂·cos ε + (π²/3)·sin ε)/|t|}. Further |p| = e^{−2π²·sin ε/|t|} → 0, so W(p) → 1, p·W₃(p) → 0, F_a(p) = p·(1 + O(p)), F_b(p) → 1, and the Borel sums are bounded. In the law for ψ₁\* the term with e^{−S̄/t}·W(p) dominates when λ₁(ε) > 0, which holds for ε < 17.1°. In the law for f₁\* that term carries F_a(p), of modulus e^{−2π²·sin ε/|t|}·(1 + o(1)), which changes the exponent to λ₂(ε); it dominates the bounded terms when λ₂(ε) > 0, that is ε < 2.52°. (c): modular growth would give limsup |t|·log|f(x)| ≤ 2πD·sin ε on every ray. By (a), (b) the limit is λⱼ(ε), which tends to Cl₂ > 0 as ε → 0, while 2πD·sin ε → 0. ∎

Remarks.

- **Indefinite theta series.** Zwegers' indefinite theta functions (signature (r − 1, 1), rational data) are holomorphic parts of real-analytic modular forms whose corrections are of the type in (b) of the Lemma (quoted: Zwegers' thesis, chapter 2), and the mock theta functions are such series divided by modular forms. Theorem 19 (c) excludes any expression of f\* or ψ\* as a finite sum of products of such functions with modular forms having the stated growth. It does not address series over cones with irrational data.
- **What the argument uses.** Only that the exponent has a component proportional to cos ε, the imaginary part ±Cl₂ of S. Irrationality of Cl₂/π² is not needed.
- **The same holds for θ₄R₁ = 2ψ₁\* − f₁\*** and for R₀, R₁, as it must for Nahm sums with a non-torsion solution. As far as growth is concerned, f\* and ψ\* sit with those, not with modular or mock modular objects.
- **For comparison:** Ramanujan's ψ₁ grows like |p|^{−1/120} = e^{(π²/60)·sin ε/|t|} and f₁ is bounded (Theorem 16 of v24): modular growth, as the Lemma says.

---

## 2. Double sums

M = (−q;q)∞, θ₄ = (q;q)∞/(−q;q)∞, R₀ = Σ (−1)ⁿq^{n²}/(q;q)ₙ, R₁ = Σ (−1)ⁿq^{n²+n}/(q;q)ₙ, and

  Q₁(n, k) = n² + n + nk + k(k+1)/2,  Q₀(n, k) = n² + nk + k(k+1)/2.

Both have the quadratic part ½·(n, k)·[[2, 1], [1, 1]]·(n, k)ᵀ, positive definite with determinant 1.

**Proposition N.**

- (N1) M·f₁\* = Σ_{n,k≥0} (−1)ⁿ·q^{Q₁(n,k)}/(q;q)ₖ, and M·f₀\* = the same with Q₀.
- (N2) M·θ₄·R₁ = Σ_{n,k≥0} (−1)^{n+k}·q^{Q₁(n,k)}/(q;q)ₖ, and M·θ₄·R₀ = the same with Q₀.
- (N3) M·ψ₁\* = the part of (N1) with k even; M·ψ₀\* = the part of (N1) for Q₀ with k odd.
- (N4) With n over all of ℤ: Σ_{n∈ℤ, k≥0} (−1)ⁿ·q^{Q₁}/(q;q)ₖ = −θ₄·Σ_{m≥0} (−1)ᵐ·q^{m²+m}/(q;q)₂ₘ₊₁, and for Q₀: +θ₄·Σ_{m≥0} (−1)ᵐ·q^{m²+m}/(q;q)₂ₘ.
- (N5) ψ₁\* = (q²;q²)∞·Σ_{n,j≥0} (−1)ⁿ·q^{n(n+1)/2 + 2nj + 2j}/((q;q)ₙ·(q²;q²)ⱼ).

Proof. (N1): 1/(−q;q)ₙ = (−q^{n+1};q)∞/M and Euler's (−x;q)∞ = Σ_k q^{k(k−1)/2}xᵏ/(q;q)ₖ at x = q^{n+1}. (N2): M·θ₄·R₁ = Σₙ (−1)ⁿq^{n²+n}·(q^{n+1};q)∞ and Euler's (x;q)∞ = Σ_k (−1)ᵏq^{k(k−1)/2}xᵏ/(q;q)ₖ. (N3): identities (I), (II) of v2, 2ψ₀\* = f₀\* − θ₄R₀ and 2ψ₁\* = f₁\* + θ₄R₁, added to (N1), (N2). (N4): for fixed k, Σ_{n∈ℤ} (−1)ⁿq^{n² + n(k+1)} is 0 for k even (the terms n and −n − k − 1 cancel) and (−1)^{(k+1)/2}·q^{−(k+1)²/4}·θ₄ for k odd; k = 2m + 1 gives the exponent m² + m. For Q₀ the roles of even and odd k are exchanged. (N5): (−q;q)ₙ = (q²;q²)ₙ/(q;q)ₙ and (q²;q²)ₙ = (q²;q²)∞·Σⱼ q^{(2n+2)j}/(q²;q²)ⱼ. ∎

What they say.

- **f\* and θ₄R are the two characters of one lattice sum.** The sums (N1), (N2) differ by the character (−1)ᵏ. With it, the sum over k is an infinite product and the double sum collapses to the rank-one Nahm sum R; without it, it does not collapse, and one gets f\*. Identities (I), (II) are the decomposition by the parity of k.
- **f\* is one-sided.** By (N4) the same sum with n over ℤ is θ₄ times a rank-one Nahm sum (the sums with denominators (q;q)₂ₘ, (q;q)₂ₘ₊₁ are the even and odd parts of H(z;q) = Σ zᵏq^{k²/4}/(q;q)ₖ of v6, the member A = ½). The variable n has no Pochhammer symbol and runs over half of its lattice: f\* stands to that Nahm sum as a partial theta function stands to a theta function.
- **(N5) as a generalised Nahm sum.** The quadratic form is ½n² + 2nj, with matrix [[1, 2], [2, 0]], indefinite but non-negative on the cone n, j ≥ 0; the symmetrisers are (1, 2). The Nahm equations (ratios of consecutive terms equal to 1) are 1 − X = −XY and 1 − Y = X², with X = qⁿ, Y = q^{2j}; eliminating Y gives X² + X + 1 = 0. So the sum has the single pair of solutions X = e^{±2πi/3} of the regular ideal tetrahedron, in agreement with the saddle W² + W + 1 = 0 of v18 and with the Garoufalidis–Zagier recipe (`psi_star/knot/gz_recipe.py`). Indefinite forms with numerator symbols traded in this way are the type of the Kashaev-invariant sums, not of the positive definite Nahm sums of the modularity conjecture.

---

## 3. The picture in the variable z (interpretation)

Let L·y(z) = y(z) − zq·y(zq) − zq²·y(zq²).

- **L·F = 1** (from (−q;q)ₙ = (−q;q)ₙ₋₁·(1 + qⁿ)), and **L·[F̌(1/z)/z] = −2** (from F̌(v) + F̌(vq) − vq²·F̌(vq²) = 2). So F(z) and −F̌(1/z)/(2z) solve the same inhomogeneous equation, one regular at z = 0 and the other at z = ∞. (proved; elementary)
- Their difference B(z) = F(z) + F̌(1/z)/(2z) = Σ_{n∈ℤ} zⁿq^{n(n+1)/2}(−q;q)ₙ solves L·y = 0. It is a theta function in z times a rank-one Nahm sum with the Jacobi variable z (v2, proof of (I), (II); v18, bilateral identity). The homogeneous equation is the one of the Rogers–Ramanujan sums Σ vⁿq^{n²+n}/(q;q)ₙ up to the multiplier θ(−v)/θ(v).
- The power-series solution of L·y = 1 at z = 0 is unique (cₙ = (qⁿ + q^{2n})·cₙ₋₁, c₀ = 1), so F is not a rational function of z plus a solution of the homogeneous equation: the extension of the Nahm module by the constant is not trivial in that sense.

| z | homogeneous solutions (bilateral sum) | inhomogeneous solutions (one-sided sums) |
|---|---|---|
| 1 | Rogers–Ramanujan functions times theta functions: modular (ψ₁ + ½f₁ is a theta product) | ψ₁, f₁: Ramanujan's fifth-order mock theta functions |
| −1 | θ₄·R₁: a Nahm sum with the non-torsion solution e^{±iπ/3}, volume Cl₂ | ψ₁\*, f₁\* |
| general, \|z\| = 1 | theta times Nahm sum, volume Vol(z) | F(z), F̌(1/z) |

So the same construction that produces a mock theta function from a modular form at z = 1 produces f\* and ψ\* from a non-modular Nahm sum at z = −1. In the exact formulas this is visible as well: at z = 1 the kernel integral is a Mordell integral (v24), at z = −1 the Borel sum of a divergent series (v17).

**Knots (from memory).** For a hyperbolic knot the q-series attached to the non-trivial flat connections form a fundamental solution of a linear q-difference equation, and Garoufalidis–Gu–Mariño–Wheeler (arXiv:2111.04763) extend that matrix by a row for the trivial connection, whose series has a perturbative expansion with constant term 1 and Borel singularities at the other connections. ψ₁\* with its three sectors (trivial with Φ₁, growing, decaying) has the same shape, with R₀, R₁ in the role of the non-trivial block. I did not compare the constructions in detail.

---

## 4. What is decided and what is not

- **Indefinite theta series / mock modular: no** (Theorem 19, given the standard growth bounds of the Lemma).
- **Nahm sum of the family of R: no.** f\* solves the inhomogeneous equation. It is that Nahm sum restricted to a half-lattice (Proposition N).
- **Generalised Nahm sum: yes, in the wide sense** with an indefinite form (N5), as for knot invariants.
- **Positive definite Nahm sum of some higher rank: not decided.** No argument here excludes it. Against it speak: the trivial sector (f₁\* → 2/3 at q → 1 with a divergent series and no exponential), which a Nahm sum has only if one of its Nahm solutions has volume and Chern–Simons value zero; and that ψ₁\* is a finite sum at every even root of unity.
- **New or known:** the kind of object is known from knots; f\* is a new instance of it, the first I know of that is tied to Ramanujan's mock theta functions by a parameter. This rests on my memory of the literature and one search, not on a survey.

---

## 5. Numerical checks (verified; `x1_kind.py`, `run_log.txt`)

| check | result |
|---|---|
| (N1)–(N5) and the versions with Q₀ | exact to q³⁰⁰ |
| Theorem 19 (a): \|t\|·log(\|ψ₁\*\|/(A\|t\|^{−1/2})) at ε = 1°, 2°, 5°, 10°, \|t\| = 0.05, against λ₁(ε) | 0.95747, 0.89961, 0.72449, 0.42843 against 0.95737, 0.89951, 0.72435, 0.42824 |
| Theorem 19 (b): the same for f₁\* with 2A at ε = 1°, 2°, against λ₂(ε) | 0.61299, 0.21069 against 0.61287, 0.21062 |
| f₁\* at ε = 5°, 10° (λ₂ < 0) | bounded: \|f₁\*\| = 0.668, 0.669 (limit 2/3) |
| Ramanujan's ψ₁ on the same rays: \|t\|·log\|ψ₁\| at \|t\| = 0.2, 0.1, 0.05 | decreases: about 0.26, 0.17, 0.10 for every ε (limit (π²/60)·sin ε ≤ 0.03); \|f₁\| bounded |
| L·F = 1, L·[F̌(1/z)/z] = −2, L·B = 0 at four pairs (z, q) | 10⁻³⁰ |

---

## 6. Scripts (`psi_star/v28/`)

- `x1_kind.py [a b c] [N]`: the checks of §5 (5 seconds).
- `run_log.txt`, `README.md`.
