# Handoff: ψ\* modularity project (state after addendum v21)

Prepared 26 Nov 2026 for a fresh session; updated after addenda v12–v21 (26–27 Nov). Project owner: Mario Klingemann. Public repo with the earlier material: https://github.com/Quasimondo/mock-theta-census (`docs/addendum-v2.md`, `docs/twists.md`, `census3/`).

This folder holds everything produced in the sessions that wrote addenda v3–v21:

- `addenda/` — `addendum-v2.md` (copy from the repo, defines (I), (II)) and `psi-star-addendum-v3.md` … `v21.md`.
- `scripts/v3` … `scripts/v21` — the scripts of each addendum, each with a `README.md`, a `run_log.txt` (complete output), and cached `*.json` limit files.

The theorems and their proofs in logical order are collected in `docs/psi-star-theorems-v1.md` (built by `psi_star/theorems/assemble.py`); start there for the proofs. For the history, read in this order: this file, then v21 and v20 (odd and even roots of unity), v19 (asymptotics of the state integrals; the forward laws as theorems), v18 (the z-family: everything below is its special case z = −1), v17 (Borel summability; its §7 lists what is open), v16 and v15 (the uniform bound it rests on), v14 and v13 (background for these), v12 (the proof of the kernel identity), v11, v10, v7, v6 (the route to it), then v5 (the laws), then v8, v9 (the F\* family). v3, v4 are background.

---

## 1. The task to continue

**The task of v11 is done (v12):** the identity

  ψ₁\*(q) = P·J(1/(2b)) + (κ/2)·[𝓑₂·p·W₃(p) − 𝓑₄·W(p)],  J(u) = ∫_L e^{−πiy² + 2πuy}/(Φ_b(y)(1 + e^{2πy/b})) dy,

is proved, together with its analogues for ψ₀\*, f₁\*, f₀\* and a general-u version (Theorems A, B, C of v12). Mechanism: residues at the kernel poles on one side give the one-sided sum; on the moved contour the kernel is a geometric series, whose terms I(u ∓ m/b) (kernel-free integrals) obey a three-term recursion that sums to W, W₃ (or F_a, F_b); and I(u) is the integral of v7 after two Fourier transforms of Φ_b. v7's contour closing is justified in v12 §4.1.

**Next tasks (1 and 2 are done: v17 for the Borel directions −120° … −80°, v19–v21 for the forward laws; what is left of each is listed in the item; 3 and 4 are open):**

1. Prove that P·J(1/(2b)) is the Borel sum of the perturbative series Φ. State after v13: Φ has a closed formula from the 1/cosh window of the integral (exact coefficients to order 300 in `scripts/v13/phi_exact.json`); numerically P·J equals the lateral Borel sum 𝒮_θ[Φ], θ in (−176.5°, −17.1°), for arg t from −250° to +30° at |t| = 0.3, including |q| > 1 (10⁻²⁸…10⁻³¹ for −220°…−60°). Borel singularities: S = π²/3 − i·Cl₂ and conjugate; a second group near S − 2π². State after v16: P·J is analytic on −210° < arg t < 10° and has the asymptotic expansion Φ uniformly on the closed sector, including the real x-axis (|q| = 1). Ingredients: descent inequality with one contour angle δ = 10° (v15, Proposition 1); products away from the real axis (v15, Theorem 2); near the real axis Φ_b by Faddeev's integral (v16: Lemma Φ = integral equals both product formulas, with proof; Lemma S = expansion in the strip |Im w| < π uniform in x; Lemma F = leading term up to the strip edge; Lemma D = dual factors outside the strip are 1 + O(e^{−c/|t|}) on the contour). Numerically at real x the kernel integral equals the Borel sum to 4·10⁻¹⁴. State after v17: **done.** Theorem 4: |P·J − Σ_{k<K} aₖtᵏ| ≤ A·K!·(240·|t|)ᴷ on −210° ≤ arg t ≤ 10°, |t| ≤ t₀, all K (ψ₁\* and ψ₀\* versions). Mechanism: shift the contour down across N = ⌊0.04/|t|⌋ lower kernel poles; exact identity P·J = Σ_{n<N} (−1)ⁿτₙ + (−1)^N·e^{t/8}·∫ H·D_N ds with τₙ = q^{n²+n}/(2(−q;q)ₙ) and D_N = q^{N²}e^{−2Nw}/(−q^{1/2}e^{−w};q)_N (functional equation of Φ_b); |D_N| ≤ e^{−0.544·N + 0.149·|s|} against the uniform bound |H| ≤ C·e^{−0.185·|s|}; Cauchy estimates for τₙ on |t| ≤ 0.3/n. Theorem 5 (with the Nevanlinna–Sokal theorem): P·J = 𝒮_{θ₀}[Φ] for every Borel direction −120° ≤ θ₀ ≤ −80°; the Borel transform is analytic and exponentially bounded near that closed sector. So the formulas of v12 hold with the Borel sum in place of the kernel integral, for |q| < 1, −90° < arg t < 10°, small |t|. Left over here: (a) Borel directions outside [−120°, −80°] (numerically −176.5° … −17.1°; v10 defined 𝒯 along the direction of t, which is proved only for −90° < arg t ≤ −80°; route: repeat the uniform bound of v15/v16 with other contour angles δ, v17 §7); (b) the constants A, t₀, R are not explicit; (c) Theorem A of v12 for |q| > 1.
1b. **The z-family (v18).** F(z) = Σ zⁿq^{n(n+1)/2}(−q;q)ₙ, z = e^{2πia}. Proved (v12 at u = a/b): for 0 < Re a < 1, F(z) = i·P·e^{−πia}·[J(a/b) − I((a−1)/b)·𝒲(T) + p·T·I(a/b)·𝒲(p²T)], T = e^{πi(2a−1)/x}; 𝒲(T;p) = Σ(−T)ᵏp^{k²+k}(p;p²)ₖ; I(a/b) = c·e^{−πia²/(2x)}·D(a/b), bilinear in H(·;q) and H(·;p²). Volume rule (derivation, verified 10⁻²⁶…10⁻³⁶ at q = e^{2πi/X}): potential V = −w²/2 − 2πi(a+k)w + Li₂(−e^{w}), critical points W² + W = 1/z, Im V = D(−W) on every branch (proved); selected saddle = root W in the upper half-plane with w = Log W, carried by I((a−1)/b); the other root is the decaying sector I(a/b). At z = 1 (a = 0): dual functions are the fifth-order mock theta functions φ₀(−p), φ₁(−p), F₀(p), F₁(p); I(0), I(1/b) are Rogers–Ramanujan products; Vol ≈ 2π·ln φ·|a|. Open: steepest descent for I((a−1)/b) (would prove the volume rule); trivial sector for a ≠ ½ (v13–v17 with the kernel e^{2πas}/(1+e^{2πs})); P·J(0) versus Watson's Mordell integral; the 3×3 matrix (quantum Jacobi question). Pitfall: the finite sum at a root of unity needs 0.07·X extra digits.
2. **Done in v19** for the forward laws. Theorem 6: 𝓑₄, 𝓑₃ = (4/3)^{1/4}e^{∓iπ/12}e^{−i·Cl₂/t}(1 + O(t)), 𝓑₂, 𝓑₁ = (4/3)^{1/4}e^{±iπ/12}e^{+i·Cl₂/t}(1 + O(t)), uniformly for 0 ≤ arg x ≤ 11° (decaying pair: ≤ 60°), including tangential approach. Contours: decaying = horizontal line through the saddle −2πi/3 (derivative of the height is a difference of arguments, monotone; works for the z-family for every 0 < a < 1, Theorem 7); growing = ray from iπ at 30°, imaginary axis from iπ to 0, positive real axis (heights: quadratic, −Cl₂(π − η), π·r). Theorem 8: ψ₁\* = 𝒮[Φ₁] + (π/√3)^{1/2}t^{−1/2}[e^{5πi/12}e^{−S̄/t}W(p)(1+O(t)) + e^{−5πi/12}e^{−S/t}p·W₃(p)(1+O(t))], S = π²/3 − i·Cl₂, for 0 < arg x ≤ 10°; same for f₁\*, ψ₀\*, f₀\*. Left: the growing integral I((a−1)/b) for a ≠ ½ (volume rule for general z); converse and mirror laws of v5 (other functions, not started). **v20:** at q₀ = e^{2πi/X}, X ≥ 4 even, x = 1/X: ψ₁\*(q₀) = P·[J(1/(2b)) + I(1/(2b)) − I(−1/(2b))] exactly (Theorem 9; radial limit lemmas R1–R3; W = p·W₃ = 1 at p = 1), and ψ₁\*(q₀) = 3^{−1/4}(X/2)^{1/2}e^{Cl₂X/2π}e^{2πi/3 − iπX/6}(1 + O(1/X)), ψ₀\* the same times e^{2πi/3} (Theorem 10). Figure-eight ratio 2/X follows with the Kashaev asymptotics (quoted). **v21:** odd X ≥ 5: radial limit A/3 = P·[J − ⅓·I(1/(2b)) − ⅓·I(−1/(2b))] (Theorem 11; uses the odd-root rule of v3 §1 and the new Lemma R2′: W(p) → ⅓, p·W₃(p) → −⅓ as p → −1, via (V), (VI)), and A/3 = ⅓·3^{−1/4}(X/2)^{1/2}e^{Cl₂X/2π}e^{2πi/3 − iπX/6}(1 + O(1/X)) (Theorem 12). So the forward law on the unit circle holds for every X with W(1) = 1, W(−1) = ⅓. Open: roots e^{2πih/k} with h ≠ 1; X = 3; the second-order term −2πi·s̄₁/X as a theorem. Numerical note: at real x integrate the state integrals over the descent contours of v19, not over L (cancellation).
3. The same construction for W, W₃ in the bounded direction (guess: exchange the roles of q and p), for 𝒥 (v9); the Stokes jump across arg t = −17.1° (v10 §3); identity (VII) (v8).
4. Someone other than the author should check the proofs. The chain for "one-sided sum = Borel sum + bilinear terms" is: v12 §2–§4 (kernel identity); v15 Proposition 1 and Lemma 3′; v16 Lemmas Φ, F, D and §4 (uniform bound); v17.

## 2. Status table

| statement | status | where |
|---|---|---|
| (I), (II): 2ψ₀\* = f₀\* − θ₄R₀, 2ψ₁\* = f₁\* + θ₄R₁ | proved | v2 |
| odd-root rule; (III), (IV) for F₀\*, F₁\* | proved | v3 |
| (V), (VI) for W, W₃; q-difference modules ⟨1, ψ₁\*, ψ₀\*⟩, ⟨1, W, W₃⟩ | proved | v4 |
| reformulation through H(z;q); X·N₀ + N₁·Y = 1; E₊′E₋ + E₋′E₊ = 2 | proved | v6 |
| identity (VII) | exact to r²⁹⁹, no proof | v8 |
| 𝓑₁…𝓑₄ = state integrals of Φ_b | proved (closing lemma in v12 §4.1) | v7, v12 |
| two-sector law for R₀, R₁ and the bilateral parts | follows from the line above + stationary phase (not written out) | v6, v7 |
| forward law ψ\* ↔ (1, W, p·W₃); f\* ↔ (2, −2F_a, −2F_b) | numerics 10⁻¹⁸…10⁻²⁹; now proved, see the v19 line below | v5 |
| converse law W, W₃ → ψ₀\*, ψ₁\*; mirror law at q = −1 | numerics 10⁻²⁹…10⁻³⁰ | v5 |
| F\* ↔ 𝒥, 𝒦 and 𝒥 → F\* | numerics 10⁻¹⁷, 10⁻²⁸ | v8, v9 |
| exact three-sector formula for ψ₁\* with Borel sum | numerics 10⁻⁵⁵ at 5 points | v10 |
| exact formulas ψ₁\*, ψ₀\*, f₁\*, f₀\* = kernel integral + bilinear terms; Theorems A, B, C for general u | proved, numerics 10⁻²⁸…10⁻³⁰ | v12 |
| kernel integral = Borel sum of Φ in the directions −120° … −80° (Theorem 5), from factorial remainder bounds on −210° ≤ arg t ≤ 10° (Theorem 4); on both sides of the real x-axis and on it | proved (using the Nevanlinna–Sokal theorem); numerics 10⁻³¹ off the real axis, 4·10⁻¹⁴ on it | v17 |
| the same for the other Borel directions in (−176.5°, −17.1°), in particular 𝒯 of v10 along the direction of t for −80° < arg t < −17.1° | numerics only (10⁻²⁸…10⁻³¹ over 160° of directions) | v10, v13 |
| shift identity P·J = Σ_{n<N} (−1)ⁿτₙ + moved integral; N residues approximate P·J up to C·e^{−0.544·N} for N ≤ 0.04/\|t\| | proved; numerics 10⁻¹³ | v17 |
| exact formula for F(z), all z (Theorem Z), both routes; dual functions 𝒲(T;p), 𝒢(T;p); z = 1: ψ₁, f₁ with φ₀, φ₁, F₀, F₁ of the dual nome and Rogers–Ramanujan coefficients | proved (specialisation of v12); numerics 10⁻²⁷…10⁻³¹; series identities exact to order 300 | v18 |
| volume rule \|F(z; e^{2πi/X})\| ~ e^{Vol(z)·X/2π}, with constant and phase: saddle W in the upper half-plane, w = Log W, of the term I((a−1)/b) | derivation; Im V = D(−W) proved; numerics 10⁻²⁶…10⁻³⁶ for ten z | v18 |
| asymptotics of 𝓑₁…𝓑₄ with constants (4/3)^{1/4}e^{∓iπ/12}, error O(t), uniform for 0 ≤ arg x ≤ 11° (Theorem 6); decaying integral I(a/b) for all 0 < a < 1 (Theorem 7) | proved; numerics: relative error 0.05·\|t\| in three directions | v19 |
| forward law of v5 for ψ₁\*, f₁\*, ψ₀\*, f₀\* with Borel sum, exact dual functions and O(t) errors, 0 < arg x ≤ 10° (Theorem 8) | proved (chain v12, v15–v17, v19) | v19 |
| ψ₁\*, ψ₀\* at even roots of unity: exact three-integral formula (Theorem 9) and asymptotics with constant and phase (Theorem 10); figure-eight ratio 2/X | proved (the last given the Kashaev asymptotics, quoted); numerics 10⁻⁸…10⁻¹³ for X = 4…60 | v20 |
| ψ₁\*, ψ₀\* at odd roots e^{2πi/X}, X ≥ 5: radial limit = one third of a period = P·[J − ⅓I(+) − ⅓I(−)] (Theorem 11), asymptotics (Theorem 12); radial limits of W, p·W₃ at p = −1 (Lemma R2′) | proved; numerics 10⁻⁹…10⁻¹³ for X = 5…61 | v21 |
| closed formula for the coefficients of Φ (Bernoulli polynomials under a 1/cosh window) | derivation; exact check of 151 coefficients | v13 |
| P·J analytic on −210° < arg t < 10° with asymptotic expansion Φ uniformly on the closed sector (includes real x, \|q\| = 1), for ψ₁\* and ψ₀\* | proved | v14, v15, v16 |
| Faddeev's integral = both product formulas; functional equations and inversion for all Re b > 0 | proved (classical; proof in v16 §1) | v16 |
| descent of the leading exponent along straight contours, margin κ ≥ 0.185 for −210° ≤ arg t ≤ 10° including the real x-axis | proved | v15 |

The exact formulas for the one-sided sums are proved (v12), and so is the Borel-sum statement for the kernel integral (v17, directions −120° … −80°). The forward asymptotic laws of v5 are proved as well (v19). The converse and mirror laws, and the F\* family, are numerics only.

---

## 3. Notation and the key formulas

**Variables.** x in the upper half-plane; b² = x; q = e^{2πix} = e^{−t}, t = −2πi·x. Often x = 1/𝒳 with Im 𝒳 < 0. Dual nome p = e^{−πi/x} = e^{−πi𝒳} (called q̃ in v5); q̃ = p² is the S-dual nome of q. Cl₂ = Cl₂(π/3) = 1.0149416064096536. φ = (1+√5)/2.

**q-series.**

- ψ₁\* = Σ(−1)ⁿq^{n(n+1)/2}(−q;q)ₙ, ψ₀\* = Σ(−1)ⁿq^{(n+1)(n+2)/2}(−q;q)ₙ
- f₀\* = Σ(−1)ⁿq^{n²}/(−q;q)ₙ, f₁\* = Σ(−1)ⁿq^{n²+n}/(−q;q)ₙ
- R₀ = Σ(−1)ⁿq^{n²}/(q;q)ₙ, R₁ = Σ(−1)ⁿq^{n²+n}/(q;q)ₙ; θ₄ = (q;q)∞/(−q;q)∞
- W = Σ(−1)ʲq^{j²+j}(q;q²)ⱼ, W₃ = Σ(−1)ʲq^{j²+3j}(q;q²)ⱼ
- F_a = Σ_{n≥1}q^{2n²−n}/(q;q²)ₙ, F_b = Σ_{n≥1}q^{2n²−3n+1}/(q;q²)ₙ
- P₊ = Σq^{2k²+k}/(q²;q²)ₖ, P₋ = Σq^{2k²−k}/(q²;q²)ₖ; ψ(q) = Σ_{n≥0}q^{n(n+1)/2}
- (V) W + F_a = ψ·P₊, (VI) q·W₃ + F_b = ψ·P₋
- N₀ = Σ(−1)ⁿq^{n²+n}/(q)₂ₙ = (−q)∞R₀, N₁ = Σ(−1)ⁿq^{n²+n}/(q)₂ₙ₊₁ = (−q)∞R₁
- X = Σ(−1)ⁿq^{n²}/(q)₂ₙ, Y = Σ(−1)ⁿq^{(n+1)²}/(q)₂ₙ₊₁
- E₊ = Σ(−1)ᵐq^{m(m+1)/2}/(q²;q²)ₘ = (q;q²)∞P₊, E₋ = Σ(−1)ᵐq^{m(m+3)/2}/(q²;q²)ₘ = (q;q²)∞P₋; E₊′, E₋′ the same without the sign
- H(z;q) = Σzᵐq^{m²/4}/(q;q)ₘ; H(z) − H(zq) = z·q^{1/4}·H(zq^{1/2})

**Faddeev dilogarithm** (convention of Garoufalidis–Kashaev): Φ_b(y) = (−q^{1/2}e^{2πby}; q)∞ / (−p·e^{2πy/b}; q̃)∞; c_b = i(b + 1/b)/2. Poles at c_b + i·b·k + i·n/b, zeros at the negatives (k, n ≥ 0).

**Bilinear combinations (v6) = state integrals (v7).**

- 𝓑₁ = N₀(q)E₊′(p) + i·q^{−1/4}Y(q)E₊(p), 𝓑₃ = N₀(q)E₋′(p) − i·q^{−1/4}Y(q)E₋(p)
- 𝓑₂ = X(q)E₊(p) + i·q^{1/4}N₁(q)E₊′(p), 𝓑₄ = X(q)E₋(p) − i·q^{1/4}N₁(q)E₋′(p)
- ∫_C Φ_b(y)·e^{−πiy²/2 + 2πuⱼy}dy = −i·b·((q)∞/(q̃)∞)·e^{−πic_b²/2 + 2πuⱼc_b}·𝓑ⱼ, with u₁ = b/4, u₂ = −b/4, u₃ = b/4 − 1/(2b), u₄ = −b/4 − 1/(2b); C from ∞e^{i(π−δ)} through 0 to ∞e^{iδ}.
- Asymptotics on horocycles: 𝓑₁, 𝓑₂ ~ (4/3)^{1/4}e^{∓iπ/12}·e^{+i·Cl₂/t} (decaying), 𝓑₃, 𝓑₄ ~ (4/3)^{1/4}e^{±iπ/12}·e^{−i·Cl₂/t} (growing).

**Exact three-sector formula (v10), κ = θ₄(q)·(p;p²)∞ / (2i·q^{1/4}·(−q;q)∞·ψ(p)):**

  ψ₁\*(q) = 𝒯(x) + (κ/2)·[𝓑₂·p·W₃(p) − 𝓑₄·W(p)],  f₁\*(q) = 2𝒯(x) − κ·[𝓑₂·F_b(p) − 𝓑₄·F_a(p)].

**Trivial sector (v10, v11):** 𝒯 = Borel sum of the perturbative series of ψ₁\* along the direction of t (for −90° < arg t < −17.1°), and

  𝒯(x) = (−q;q)∞/(b·q^{1/8}·(p;p²)∞) · ∫_{e^{−iδ}ℝ} e^{−πiy² + πy/b} / (Φ_b(y)·(1 + e^{2πy/b})) dy.

**Proved in v17:** this integral (P·J of v12, defined for −210° < arg t < 10° with δ = 10° and Φ_b from Faddeev's integral) equals the Borel sum 𝒮_{θ₀}[Φ] in every direction −120° ≤ θ₀ ≤ −80°, for Re(e^{iθ₀}/t) > 1/R. For −90° < arg t ≤ −80° this is 𝒯 as defined above; for −80° < arg t < −17.1° the equality with the Borel sum along the direction of t is numerical. Residue form: P·J = Σ_{n<N} (−1)ⁿ·q^{n²+n}/(2(−q;q)ₙ) + O(e^{−0.544·N}) for N ≤ 0.04/|t|.

**Constants of the laws (v5):** forward: S̄ = π²/3 + i·Cl₂, C̄ = √(π/√3)e^{5πi/12}, first corrections 1/4 + i/(18√3) (ψ₁\*), −1/12 + i/(18√3) (ψ₀\*). Converse: c = √(2π/√5), exponents −1/60, 11/60 of the dual nome. Mirror: κ_mirror = √(π/√5). F\* family (v8, v9): dual nome r = e^{−πi𝒳/2}, exponents −1/120, 71/120.

---

## 4. Numerical methods that worked

- **Horocycle limits** (v5 §1): 𝒳 = 𝒳₀ + (period)·k − i·b, so q → 1 inside the disk while the dual nome stays fixed. Subtract the trivial sector by optimal truncation of the perturbative series, then extrapolate polynomially in t (`e1_limits.extrap`). Sampling the dual nome on a circle and taking Fourier coefficients gives the Taylor coefficients of the unknown dual function (`e4_run.py`, `fstar.py`).
- **Exact cancellation of the growing sector** (v5 §5, v6 §2): combine two q-side functions with dual-side functions so that the growing sector cancels identically; what remains is the decaying sector.
- **Prony fit** for unknown phases (v9 `jloop.py`).
- **State integrals**: mpmath `quad` on the tilted contour with Φ_b from `qp` products (`scripts/v7/state_integral.py`, `scripts/v11/kernel_integral.py`). Needs Im x not too small (the q-product in Φ_b needs about 70/(2π·Im x) factors).
- **Φ_b at real b and fast kernel integrals** (v16 `fad.py`, v17 `g1_shift.py`): Faddeev's integral in the variable w by the trapezoid rule, continued by functional equation and inversion; double precision, valid for |arg x| < 90°. `g1_shift.integral(t, x, N, case)` gives P·J (N = 0) or the integral moved across N poles.
- **Borel sum**: Padé [60/60] of the Borel transform at 150 digits, Laplace integral along the ray of t (`scripts/v10/exact3.py`). The 151 coefficients of the series at q → 1 are in `q2_series_v3.json` (from `scripts/v3/q2_series.py`).

Pitfalls met:

- python-flint balls blow up over 10⁴–10⁵ complex multiplications (wrapping effect). The long recurrences use `.mid()` after every step; values are then not rigorous, and precision was checked by rerunning with more bits (done explicitly only for `horo.near`).
- `arb(fmpq)` is evaluated at the precision in force when it is created: set `ctx.prec` first.
- Sums near the unit circle have intermediate terms of size e^{0.16·K} (ψ\*, R) or e^{0.08·K} (W, F\*, P): choose the precision accordingly (see `horo.py`, `hfam.py`).
- Stop criteria must use the size of the terms, not of q^{n²} alone.
- θ₄ and (−q;q)∞ near q = 1: use the modular formulas, not the products.

Environment used: Python 3 with mpmath and python-flint. Each scripts folder is self-contained; run from inside it.

---

## 5. Ground rules the addenda follow

From the owner's brief: label every claim **proved / verified / inference** (plus "derivation" for formal computations); exact identities checked to at least q³⁰⁰; constants to at least 15 digits, or state the achieved precision; name the script for each claim; negative results with a clear test are welcome. Addenda are written summary first, one section per topic, scripts listed at the end, with a zip containing the scripts and a complete `run_log.txt`.

## 6. Things not checked

- Literature: only abstracts, tables of contents or search snippets were read for Garoufalidis–Kashaev (arXiv:1304.2705), Garoufalidis–Wheeler (arXiv:2203.17029) and Wheeler's thesis. Of Garoufalidis–Gu–Mariño–Wheeler (arXiv:2111.04763v2) §2.6–2.7 (pp. 19–22) were read for v12: their trivial-connection integral for 4₁ has the same kernel (tanh(πy/b)) and their Theorem 7 the same shape; their proof (closing downward) does not apply here. Whether this case follows from results elsewhere in that paper is open.
- The other session's figure-eight comparison (|ψ₁\*(e^{2πi/X})| / J_{4₁}(X/2) → 2/X) is reproduced in v18 §2.5, and follows from the asymptotic formula of v18 §2.2. `census3/zfamily/` of the repository selected its saddle by |Im V| alone, which cannot distinguish the two roots; v18 §2.4 has the translation.
- No search for (V), (VI), (VII), W, 𝒥, 𝒦 in the literature beyond standard mock theta lists.
