# ψ₁\* and the z-family compared with Garoufalidis–Gu–Mariño–Wheeler

27 Nov 2026. Claude Opus 5.5 at the request of Mario Klingemann, from a reading of the full paper.

**The paper:** S. Garoufalidis, J. Gu, M. Mariño, C. Wheeler, *Resurgence of Chern–Simons theory at the trivial flat connection*, arXiv:2111.04763v2 (18 Apr 2026); Comm. Math. Phys. 406 (2025) 20. Below, "GGMW (n)" means their equation (n). Everything is about their first example, the figure-eight knot 4₁, Sections 2 and 3. Their 5₂ section was only skimmed.

**Summary.** ψ₁\* fits the structure GGMW describe for 4₁ almost point by point:
- a 3 × 3 matrix with a trivial row that solves an inhomogeneous q-difference equation and a 2 × 2 non-trivial block that solves the homogeneous one;
- a state integral with a tanh-type kernel for the trivial row;
- factorisation into q- and q̃-series by residues;
- a trivial perturbative series whose Borel singularities are the actions of the non-trivial connections;
- a Jacobi variable that deforms the hyperbolic structure.

Two differences matter for what is new here:
1. For 4₁, the asymptotics and Borel summability of the trivial-row state integral are **conjectures** (GGMW (74), (75); Conjectures 5, 6, 14). For ψ₁\* the corresponding statements are **proved** (Theorems II, IV, V of `docs/psi-star-theorems-v1.md`).
2. GGMW's objects come from knots. Ours is a rank-one example from a sign-twisted fifth-order mock theta function, and its Jacobi variable z ends, at z = 1, in Ramanujan's mock theta functions and a Mordell integral.

---

## 1. Dictionary

| GGMW, 4₁ | ψ₁\* and the z-family | remark |
|---|---|---|
| 2 × 2 block G^{(0)}_m, G^{(1)}_m (GGMW (6), (8)), solutions of f_m + (q^{m+1} − 2)f_{m+1} + f_{m+2} = 0 (GGMW (5)) | the bilateral sum B(z) = Σ_{n∈ℤ} zⁿq^{n(n+1)/2}(−q;q)ₙ, a theta function times the rank-one Nahm sum H, solution of y(z) − zq·y(zq) − zq²·y(zq²) = 0 (v28 §3); at z = −1 this is θ₄R₁ | second-order equations; our shift z → zq plays the role of the descendant m → m + 1 |
| trivial row G^{(2)}_m (GGMW (9), (11)), solution of the same equation with right-hand side 1; 3 × 3 matrix J_m (GGMW (12), (13)) | F(z), with L·F = 1 (v28 §3); ψ₁\* = F(−1); and F̌(1/z), with L·[F̌(1/z)/z] = −2 | same structure. Theirs is the ε²-coefficient of a deformation (GGMW (10)), with harmonic sums E₁, E₂; ours is a plain q-series |
| state integral Z(τ) = −(i/2b³)∫_C tanh(πy/b)·Φ_b(−y + i/2b)/Φ_b(y − i/2b) dy (GGMW (73)), descendants (76) | P·J(u), J(u) = ∫_L e^{−πiy² + 2πuy}/(Φ_b(y)(1 + e^{2πy/b})) dy (v11, v12) | 1/(1 + e^{2πy/b}) = ½(1 − tanh(πy/b)), so J = ½I − ½·(tanh-kernel integral). Same kernel; their integrand is a ratio of two dilogarithms, ours one dilogarithm with a Gaussian. GGMW trace the kernel to Kashaev [Kas97] (a factor 1/sinh) |
| poles: the trivial-row q-series comes from triple poles (double pole of the dilogarithm ratio times the tanh pole) (proof of GGMW Theorem 7) | ψ₁\* comes from simple poles of the kernel (v12, Theorems A–C) | hence no harmonic sums in ψ₁\* |
| factorisation Z_{λ,μ} = q^{λ/2}(−1)^μ[G^{(2)}_λ(q) + τ⁻¹G^{(1)}_λ(q)L^{(0)}_μ(q̃) − τ⁻²G^{(0)}_λ(q)L^{(1)}_μ(q̃)] + … (GGMW Theorem 7, (78)) | P·J(1/(2b)) = ψ₁\*(q) − P·I(1/(2b))·p·W₃(p) + P·I(−1/(2b))·W(p), with I bilinear in H(·;q) and H(·;p²) (v12 Theorem C) | in both, the trivial-row integral is the trivial q-series plus (non-trivial q-series) × (q̃-series) |
| q̃-series of the trivial row: L^{(j)}_μ = Casoratians G^{(j)}_{μ+1}G^{(2)}_μ − G^{(j)}_μG^{(2)}_{μ+1} (GGMW (16)), Appell–Lerch type sums (GGMW (17)); first row of J_m⁻¹ (Theorem 4) | the dual side: H(·;p²) times the dual series W, W₃ (= 𝒲(T;p) at T = 1, p²), with 𝒲 in the same family as F (v18 §1) | in GGMW the q̃-factors are trivial × non-trivial products combined into Casoratians; ours are such products not yet combined. See §3, question 1 |
| perturbative series φ^{(σ₀)} with constant term 1 (GGMW (29)); singularities at the actions of σ₁, σ₂ (GGMW (32), (33)) | Φ₁ with a₀ = 1/3 (v13); singularities at S, S̄ = π²/3 ∓ i·Cl₂(π/3) (v13) | in both, the trivial series "knows" the non-trivial actions |
| Z(τ) ∼ Φ^{(σ₀)}(τ) and Z = s(Φ^{(σ₀)}) − (i/2)τ^{−3/2}s(Φ^{(σ₂)}): **conjectured** (GGMW (74), (75), Remark 9: Z = s_med(Φ^{(σ₀)}) for τ > 0) | P·J ∼ Φ₁ with factorial bounds on a 220° sector, and P·J = Borel sum of Φ₁ for −120° ≤ θ₀ ≤ −80°: **proved** (v14–v17, checked against Sokal 1980) | the main difference in status |
| Jacobi variable x, "the monodromy of the meridian" (Section 3.2); C_m(x,q) = Σ (−1)ᵏq^{k(k+1)/2+km}/((x⁻¹;q)_{k+1}(x;q)_{k+1}) (GGMW (94)), related to the Gukov–Manolescu series F_{4₁} (GGMW (95)); inhomogeneous equation with A(x, qᵐ, q) (GGMW (97)) | z = e^{2πia}, F(z) = 𝔉(z, −q; q); shape equation zW(1 + W) = 1 for a single tetrahedron; Vol(z) (v18, v22–v23) | both deform the hyperbolic structure. Their x enters through denominators (meromorphic in x), our z through zⁿ (entire in z) |
| trivial series at x ≈ 1: φ^{(σ₀)}(x) = −1/(x⁻¹ − 3 + x) − … (GGMW (92)), the inverse Alexander polynomial | P·J(a/b) ≈ 1/(1 − 2z) numerically (v22, v23); termwise Σ(2z)ⁿ at q = 1 | the analogue of the inverse Alexander polynomial is 1/(1 − 2z); at z = −1 it is a₀ = 1/3 |
| resurgence in x: Conjecture 14, (116); Remark 16 (C₀(x,q) involves the Borel sums of all three series) | Theorem Z (exact, every z); the volume rule at even roots of unity (Theorem 15, proved); Borel summability of the trivial sector for a ≠ ½ open | |
| Kashaev invariant ⟨4₁⟩_N from the same integral over a contour A_N around N poles (GGMW (69)); colored Jones from the (u, τ) integral (GGMW (139)–(142)) | finite sums at roots of unity equal the three integrals at real b (Theorems 9, 11, 14 of v20–v22) | the same mechanism: residue sums of the state integral give the finite state sums |

---

## 2. What is different

- **Origin.** GGMW's series come from a knot: colored Jones, Kashaev invariant, Habiro ring. ψ₁\* comes from Ramanujan's ψ₁ with a sign twist. It is not the Kashaev invariant of 4₁: it has weight ½ instead of 3/2, and a different perturbative series (v18 §2.5, v20 corollary). It has the volume and torsion of one tetrahedron, not two.
- **Rank.** Our state integrals have one Faddeev dilogarithm and a Gaussian (rank one, one tetrahedron). Theirs have a ratio of two dilogarithms. This is why our trivial q-series comes from simple poles and has no harmonic sums. It also makes the analytic estimates (v14–v17, v19) tractable.
- **The Jacobi variable at its end point.** In GGMW, x near 1 is the regime studied (the complete hyperbolic structure). Our z runs over the whole unit circle. At z = 1 the volume vanishes, the state integrals become entries of the Rogers–Ramanujan S-matrix, the dual series become Ramanujan's φ₀, φ₁, F₀, F₁, and the kernel integral becomes a Mordell integral (v24, Theorems 16, 17). GGMW mention mock theta functions only in the title of a reference [Zwe01], and Appell–Lerch sums in Theorem 4 and §4.3.
- **Status.** GGMW's asymptotic and resurgence statements for the trivial row are conjectures supported by numerics (Conjectures 5, 6, 14 and (74), (75)). For ψ₁\* the analogous statements are proved, given standard results: v12, v14–v17, v19–v21, v25, collected in `docs/psi-star-theorems-v1.md`. No referee has checked these proofs yet.
- **A coincidence to keep in mind.** For 5₂, Rogers–Ramanujan functions appear as entries of a 2 × 2 block (GGMW Conjecture 19, (179), (180)). At z = 1 the Rogers–Ramanujan S-matrix appears in our state integrals (v24). Whether this is more than a coincidence was not looked at.

---

## 3. Questions this suggests

1. **Casoratians and Appell–Lerch sums.** In GGMW the first row of J_m⁻¹ consists of the Casoratians L^{(j)}_m of the trivial and the non-trivial series, which are Appell–Lerch type sums (Theorem 4) and satisfy a first-order equation L_{m−1} − L_m = G_m (GGMW (19)). For our equation L·y = 1, a short computation gives the analogue: for a homogeneous solution h and D(z) = h(z)F(zq) − h(zq)F(z),

   D(z) + zq²·D(zq) = −h(zq).

   Questions:
   - What is D for h = B, the bilateral sum, and is it an Appell–Lerch or partial-theta-type sum in z?
   - Does it give the inverse of our 3 × 3 matrix?
   - Does that inverse explain the products H(·;p²)·W(p) in Theorem C (v12)? It would answer the open question whether 𝒲 is an Appell–Lerch sum (v18 §6): in GGMW the Appell–Lerch object is the Casoratian, not the trivial series itself.
2. **|q| > 1.** GGMW extend all entries to |q| > 1 by Jₘ(q⁻¹) = (constant matrix)·J_{−m−1}(q)·(constant matrix) (GGMW (7), (14), (98)). The analogue would give Theorem A of v12 for |q| > 1, which is open (`docs/psi-star-theorems-v1.md` §0.7). Is F̌(1/z) the "inverted" partner of F in the sense of GGMW's remark after (75) (Habiro series against inverted Habiro series, the "upside-down cake" of Garoufalidis–Zagier)?
3. **Transfer of the proofs.** Our proof that P·J is the Borel sum of Φ₁ (descent of the leading exponent, uniform bound, factorial bounds, Nevanlinna–Sokal) uses the special form of the integrand only through Φ_b and the kernel. Could it prove GGMW (74), (75) for 4₁, with the ratio Φ_b(−y + i/2b)/Φ_b(y − i/2b) in place of 1/Φ_b? The main new difficulty would be the descent inequality (Part 3) for a two-dilogarithm potential.
4. **The matrix and its cocycle.** GGMW define W_{S,λ,μ}(τ) = J_λ(q̃)⁻¹·diag(τ^{3/2}, τ^{1/2}, τ^{−1/2})·J_μ(q) (GGMW (83)) and show that it extends to a holomorphic function on ℂ \ (−∞, 0] (Theorem 10). The analogue for ψ₁\* would be the 3 × 3 matrix-valued quantum modular form expected in v18 §4. Writing it down, with the rows (F, two homogeneous solutions) and the dual side, is now a concrete task with a template.
