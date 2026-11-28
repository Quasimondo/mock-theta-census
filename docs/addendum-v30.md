# Mock theta search: addendum v30 (the 3 × 3 matrix: state integrals = dual frame⁻¹ · Jacobi factor · q-frame; the q-frame under Γ₀(2); the relations (R3), (R4))

Session of 28 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). Two tasks: first the short one, a proof of the series identities (R3), (R4) of v24; then question 4 of `docs/related-ggmw.md` §3, "the matrix and its cocycle". Scripts and output are in `psi_star/v30/`.

Labels: **proved** (a complete argument is written here or in an earlier addendum; nobody else has checked it), **verified** (numerics), **derivation** (a computation whose analytic steps are not all justified), **inference**, **from memory** (a statement of the literature that I did not look up).

Notation as in v18 and v29: Im x > 0 unless stated, b = √x, q = e^{2πix}, p = e^{−πi/x}; fractional powers q^α = e^{2πixα}, p^α = e^{−πiα/x}. For a ∈ ℂ: ζ = e^{πia}, z = ζ², ζ̃ = e^{πia/x}, T = p·ζ̃². F(z) = Σ zⁿq^{n(n+1)/2}(−q;q)ₙ, 𝒲(T) = Σ (−T)ᵏp^{k²+k}(p;p²)ₖ, H(y; r) = Σ yᵐr^{m²/4}/(r;r)ₘ with even and odd parts H_e, H_o in y. I(u), J(u) are the state integrals of v12, P = 2^{−1/2}·q^{−1/6}·p^{−1/24}/b.

---

## 0. Summary

1. **(R3), (R4) (proved, given four relations of Ramanujan's last letter; §1).** ψ₁ − f₁ = q·χ₁(q²) and ψ₀ − f₀ + 2 = χ₀(q²) follow in two lines from f₀(q) + 2F₀(q²) − 2 = φ₀(−q²) + ψ₀(q), f₁(q) + 2q·F₁(q²) = ψ₁(q) − φ₁(−q²)/q and the two relations for χ₀, χ₁ that v24 already used. So Theorem 17 of v24 (the kernel integral is a Mordell integral) no longer rests on an identity that was only checked to q³⁰⁰. The four relations are stated from memory; their form is confirmed to q³⁰⁰.
2. **The matrix (proved; §2, Theorem 20).** With the six state integrals collected in

     𝐒(a, x) = [ 1, 𝒥(a), 𝒥(a+x) ; 0, 𝒾₀(a), 𝒾₀(a+x) ; 0, 𝒾₁(a), 𝒾₁(a+x) ],

   𝒥(a) = i·P·e^{−πia}·J(a/b), 𝒾₀(a) = i·P·e^{−πia}·I(a/b), 𝒾₁(a) = i·P·e^{−πia}·I((a−1)/b), one has

     **𝐒(a, x) = 𝐉̃(ζ̃; p)⁻¹ · Δ(a, x) · 𝐉(ζ; q).**

   𝐉 is a matrix of q-series (first row 1, F(z), F(zq); the other rows from H(±ζ; q)), 𝐉̃ a matrix of p-series (first row 1, p·T·𝒲(p²T), −𝒲(T); the other rows from H(·; p²)), and Δ = diag(1, δ, δ) with the Jacobi-type factor δ = e^{πi/4}·q^{−5/24}·p^{1/24}·e^{−πia − πia²/(2x)}/(2b). This is the analogue of (83) of Garoufalidis–Gu–Mariño–Wheeler. Its first row is Theorem Z, the other two are Theorem C of v12. det 𝐒 is elementary.
3. **Continuation in x (proved; §3).** Every entry of 𝐒 is holomorphic in (a, x) for x ∈ ℂ \ (−∞, 0]. This is the analogue of GGMW's Theorem 10: the product of the three factors, defined for Im x > 0, continues across the positive real axis.
4. **The q-frame is a matrix-valued quantum Jacobi form on Γ₀(2) (proved; §4, Theorem 21).** 𝐉(a; x+1) = C_T·𝐉(a; x) with a constant matrix, and for γ = (1 0; 2 1), x′ = x/(1+2x), a′ = a/(1+2x):

     𝐉(a′; x′) = Λ(a, x) · 𝐉(a; x) · Ω_γ(a, x),  Ω_γ = 𝐒(a, x)⁻¹·𝐒(a′, x′),

   with Λ diagonal and elementary, and Ω_γ holomorphic on ℂ \ (−∞, 0]. T and γ generate Γ₀(2). The first row of 𝐉 is (1, F(z), F(zq)), so this is the precise form of the statement expected in v18 §4. The dual frame belongs to a different function (the second argument of 𝔉 is q^{1/2} instead of −q), which is why the group is Γ₀(2) and not SL₂(ℤ).
5. **The natural domain is ℂ \ (−∞, −½] (verified, not proved; §5).** Ω_γ continues across the interval (−½, 0), where neither 𝐒(a, x) nor 𝐒(a′, x′) is defined: the jump across the axis goes to 0 linearly at three points, and the negative Fourier modes on a circle around −0.2 are 10⁻⁵ of the positive ones. Across (−∞, −½) the jump stays of order 1. This is the domain that Garoufalidis–Zagier's definition of a holomorphic quantum modular form asks for (cx + d ∉ (−∞, 0]).
6. **|q| > 1 (§6).** The matrix 𝐒 exists there, and the equations L·𝒥 = 1, L·𝒾 = 0 hold there (proved by continuation, verified). The 2 × 2 block is given by series in q⁻¹ and p⁻² (proved: Φ_b = Φ_{1/b}). The trivial row is different: for |q| > 1 the equation L·y = 1 has no solution that is a convergent Laurent series in z (proved), and neither residue sum of the kernel integral converges (derivation). So below the real axis the kernel integral is not a q-series plus bilinear terms; it is an integral of Mordell type throughout. This differs from the figure-eight knot, where the series converge on both sides (GGMW (7), (14), as summarised in `docs/related-ggmw.md`).

---

## 1. The relations (R3), (R4) of v24

v24 used four series identities: (R1) χ₀(q) = 2F₀(q) − φ₀(−q), (R2) q·χ₁(q) = 2q·F₁(q) + φ₁(−q), (R3) ψ₁(q) − f₁(q) = q·χ₁(q²), (R4) ψ₀(q) − f₀(q) + 2 = χ₀(q²). (R3) and (R4) were only verified to q³⁰⁰.

**Relations of the last letter (from memory; form verified to q³⁰⁰, `r1_letter.py`).** In his last letter Ramanujan gave, for each of the two fifth-order families, a chain of relations, which Watson proved (G. N. Watson, *The mock theta functions (2)*, 1937). In the form that the exact series confirm:

- (A0) f₀(q) + 2F₀(q²) − 2 = φ₀(−q²) + ψ₀(q)  (both equal θ₄(0, −q)·G(−q)),
- (A1) f₁(q) + 2q·F₁(q²) = ψ₁(q) − φ₁(−q²)/q  (both equal θ₄(0, −q)·H(−q)),

together with (R1), (R2). I did not compare these with the printed text. Andrews–Garvan list Watson's identities as their (2.12)–(2.21); a comparison with that list would turn "from memory" into "quoted".

**Proposition (proved, given (A0), (A1), (R1), (R2)).** (R3) and (R4) hold.

Proof. (A0) reads ψ₀(q) − f₀(q) + 2 = 2F₀(q²) − φ₀(−q²), and the right side is χ₀(q²) by (R1) at q². (A1) reads ψ₁(q) − f₁(q) = 2q·F₁(q²) + φ₁(−q²)/q, and the right side is q·χ₁(q²) by (R2) at q² divided by q. ∎

**Consequence.** In v24 §6 the line for (R3), (R4) can read "proved from the relations of Ramanujan's last letter (Watson 1937)". Theorem 17 of v24 is then proved given the law of χ₀, χ₁ (Gordon–McIntosh, quoted there) and classical identities only.

For orientation, `r1_letter.py` also confirms the other members of the two chains: 2φ₀(−q²) − f₀(q) = θ₄(0, q)·G(q), f₁(q) + 2φ₁(−q²)/q = θ₄(0, q)·H(q), and the Lost Notebook relations ψ₀(q) − F₀(q²) + 1 = q·ψ(q²)·H(q⁴), ψ₁(q) − q·F₁(q²) = ψ(q²)·G(q⁴).

---

## 2. The matrix

### 2.1 Two frames and the Jacobi factor

**Lemma W (v6 §1; proved).** For |r| < 1 and a fixed choice of r^{1/4}: H(y) − H(yr) = y·r^{1/4}·H(y·r^{1/2}), and

  H_e(y)·H_e(y·r^{1/2}) − H_o(y)·H_o(y·r^{1/2}) = 1,  equivalently  H(y)·H(−y·r^{1/2}) + H(−y)·H(y·r^{1/2}) = 2.

Proof. The first equation is the recursion of 1/(r;r)ₘ. Splitting it by parity gives H_e(y) − H_e(yr) = y·r^{1/4}·H_o(y·r^{1/2}) and the same with e and o exchanged. Hence the left side W(y) of the second equation satisfies W(y·r^{1/2}) = W(y). W is entire in y, so its Taylor coefficients wₙ satisfy wₙ·(r^{n/2} − 1) = 0, and W = W(0) = 1. ∎

Write H(y) for H(y; q) with q^{1/4} = e^{πix/2}, and H̃(y) for H(y; p²) with (p²)^{1/4} = p^{1/2}.

**The q-frame.** With μ = q^{−3/4}·ζ⁻¹:

  𝐉(a; x) = 𝐉(ζ; q) = [ 1, F(z), F(zq) ; 0, H(ζ), μ·H(ζq^{1/2}) ; 0, H(−ζ), −μ·H(−ζq^{1/2}) ].

Its entries are series in q^{1/4} whose coefficients are Laurent polynomials in ζ (the third column carries the factor q^{−3/4}). By Lemma W, det 𝐉 = −2μ.

**The dual frame.**

  𝐉̃(a; x) = 𝐉̃(ζ̃; p) = [ 1, p·T·𝒲(p²T), −𝒲(T) ; 0, H̃_e(ζ̃p²), ζ̃⁻¹p^{−1/2}·H̃_o(ζ̃p) ; 0, H̃_o(ζ̃p²), ζ̃⁻¹p^{−1/2}·H̃_e(ζ̃p) ].

Its entries are series in p^{1/2} whose coefficients are Laurent polynomials in ζ̃. By Lemma W (at y = ζ̃p, r = p²), det 𝐉̃ = ζ̃⁻¹p^{−1/2}, and the lower block is the inverse of

  M̃ = [ H̃_e(ζ̃p), −H̃_o(ζ̃p) ; −ζ̃p^{1/2}·H̃_o(ζ̃p²), ζ̃p^{1/2}·H̃_e(ζ̃p²) ].

**The Jacobi factor.** Δ(a, x) = diag(1, δ, δ),

  δ(a, x) = i·P·c·e^{−πia − πia²/(2x)} = (e^{πi/4}/(2b))·q^{−5/24}·p^{1/24}·e^{−πia − πia²/(2x)},

where c = 2^{−1/2}e^{−πi/4}q^{−1/24}p^{1/12} is the constant of Theorem C (v12). P is elementary because (p;p²)∞ = √2·p^{1/24}q^{1/24}·(−q;q)∞ (the transformation of η(2x)/η(x)); checked to 10⁻³¹.

**In which sense these are frames.** L·y = y(a) − zq·y(a+x) − zq²·y(a+2x) is the operator of v28 in the shift a → a + x. With γ(a) = e^{−πia − πia²/(2x)}:

- y₊(a) = γ(a)·H(ζ) and y₋(a) = ζ̃·γ(a)·H(−ζ) satisfy L·y = 0 (the equation of H at y = ±ζ; note γ(a+x) = μ·γ(a), γ(a+2x) = q⁻²ζ⁻²·γ(a), and ζ̃ changes sign under a → a + x).
- So diag(1, γ, ζ̃γ)·𝐉 = [ 1, F(a), F(a+x) ; 0, y₊(a), y₊(a+x) ; 0, y₋(a), y₋(a+x) ] is a fundamental matrix of the inhomogeneous equation, as in v29 §2. The Gaussian e^{−πia²/(2x)} is not a function of z and q. It is the automorphy factor of a Jacobi form, and it is put into Δ so that 𝐉 and 𝐉̃ contain q-series and p-series only.

### 2.2 Theorem 20

**Theorem 20 (proved).** For Im x > 0 and every a ∈ ℂ:

  𝐉̃(a; x)·𝐒(a, x) = Δ(a, x)·𝐉(a; x),  that is  𝐒 = 𝐉̃⁻¹·Δ·𝐉.

Moreover det 𝐒 = 𝒾₀(a)·𝒾₁(a+x) − 𝒾₀(a+x)·𝒾₁(a) = −(i/(2x))·q^{−7/6}·p^{7/12}·e^{−3πia + πia/x − πia²/x}.

Proof.

- *Rows 2 and 3.* Theorem C of v12 in the form of Theorem Z: I(a/b) = c·e^{−πia²/(2x)}·[H_e(ζ)·H̃(−ζ̃p) + H_o(ζ)·H̃(ζ̃p)] = c·e^{−πia²/(2x)}·[H(ζ)·H̃_e(ζ̃p) − H(−ζ)·H̃_o(ζ̃p)]. The substitutions are:
  - a → a − 1: ζ → −ζ, ζ̃ → ζ̃p, e^{−πia²/(2x)} → e^{−πia²/(2x)}·ζ̃p^{1/2};
  - a → a + x: ζ → ζq^{1/2}, ζ̃ → −ζ̃, γ → μ·γ.

  They give (𝒾₀(a), 𝒾₀(a+x)) = δ·[H̃_e(ζ̃p)·ρ₂ − H̃_o(ζ̃p)·ρ₃] and (𝒾₁(a), 𝒾₁(a+x)) = δ·ζ̃p^{1/2}·[−H̃_o(ζ̃p²)·ρ₂ + H̃_e(ζ̃p²)·ρ₃], where ρ₂, ρ₃ are the last two entries of rows 2 and 3 of 𝐉. So the lower 2 × 2 block of 𝐒 is δ·M̃ times that of 𝐉.
- *Row 1.* The strip form of Theorem Z at a and at a + x (T does not change): F = 𝒥 − 𝒲(T)·𝒾₁ + p·T·𝒲(p²T)·𝒾₀ in both columns. So (1, 𝒥(a), 𝒥(a+x)) = (1, F(z), F(zq)) − p·T·𝒲(p²T)·(row 2 of 𝐒) + 𝒲(T)·(row 3 of 𝐒).
- *Together.* 𝐒 = [ 1, v·M̃ ; 0, M̃ ]·Δ·𝐉 with v = (−p·T·𝒲(p²T), 𝒲(T)), and the inverse of the first factor is 𝐉̃.
- *Determinant.* det 𝐒 = δ²·det 𝐉/det 𝐉̃ = δ²·(−2q^{−3/4}ζ⁻¹)·ζ̃p^{1/2}. ∎

**Corollary (the converse of Theorem Z; proved).** The first row of 𝐉̃ = Δ·𝐉·𝐒⁻¹ is

  (1, p·T·𝒲(p²T), −𝒲(T)) = (1, F(z), F(zq))·𝐒(a, x)⁻¹.

Written out, with d = det 𝐒:

  𝒲(T) = { [F(z) − 𝒥(a)]·𝒾₀(a+x) − [F(zq) − 𝒥(a+x)]·𝒾₀(a) }/d,  p·T·𝒲(p²T) = { [F(z) − 𝒥(a)]·𝒾₁(a+x) − [F(zq) − 𝒥(a+x)]·𝒾₁(a) }/d.

So the dual trivial series is recovered exactly from F at z and zq and the six state integrals.

Remarks.

- **Comparison with GGMW.** Their (83) is W(τ) = J(q̃)⁻¹·diag(τ^{3/2}, τ^{1/2}, τ^{−1/2})·J(q), with the state integrals as entries of W (their Theorem 7 for the first row). Here the diagonal matrix has the weight factor 1/b (weight ½) and the Jacobi factor e^{−πia²/(2x)}, and the two frames are different functions (item 4 of the summary).
- **The inverse frames.** By v29 §2 the first row of 𝐉⁻¹ consists of 1 and the Casoratians of F with the two H-rows, divided by det 𝐉. The first row of 𝐉̃⁻¹ is (1, v·M̃): the dual Casoratians K of v29 §3.
- **Shifts.** 𝐒(a + x, x) = 𝐒(a, x)·𝐀(z), 𝐀 = [ 1, 0, −1/(zq²) ; 0, 0, 1/(zq²) ; 0, 1, −1/q ], because every row of 𝐒 satisfies L·y = 1 or 0 (the dual coefficients do not change under a → a + x).

---

## 3. Continuation in x

For 0 < δ₀ < π/2 let L_{δ₀} = e^{−iδ₀}·ℝ, and β = arg b.

**Lemma E (bounds on the line; proved, given Lemma Φ of v16).** Let K be a compact subset of {x ∈ ℂ \ (−∞, 0] : |β| < π/2 − δ₀}. There is C with

  |1/Φ_b(y)| ≤ C for y = s·e^{−iδ₀}, s ≤ 0,  |e^{πiy²}/Φ_b(y)| ≤ C for s ≥ 0,  |1/(1 + e^{2πy/b})| ≤ C,  for all x ∈ K.

Proof.

- *Zeros and poles.* The zeros of Φ_b are y = −i[(m+½)b + (n+½)/b], the poles are their negatives (m, n ≥ 0), the poles of the kernel are i·b·(n+½), n ∈ ℤ. Their arguments lie in ±π/2 + [−|β|, |β|]. The line has arguments −δ₀ and π − δ₀. For |β| < π/2 − δ₀ the line meets none of them.
- *The second bound follows from the first.* By the inversion relation, e^{πiy²}/Φ_b(y) = Φ_b(−y)/Φ_b(0)². So it is enough to bound Φ_b(y)^{±1} for s ≤ 0.
- *Im x > 0.* For s ≤ 0, Re(by) = |b|·s·cos(β − δ₀) ≤ 0 and Re(y/b) = (s/|b|)·cos(β + δ₀) ≤ 0. In the product formula Φ_b(y) = (−q^{1/2}e^{2πby}; q)∞/(−p·e^{2πy/b}; p²)∞ both exponentials have modulus ≤ 1, so |Φ_b(y)|^{±1} ≤ (−|q|^{1/2};|q|)∞·(−|p|;|p|²)∞/((|q|^{1/2};|q|)∞·(|p|;|p|²)∞). This is locally uniform on Im x > 0.
- *Im x < 0.* Φ_b = Φ_{1/b} (Lemma Φ (i)), and the hypotheses are symmetric under b → 1/b.
- *Near a point x₀ > 0.* Put w = 2πby and choose the integer N with |Im w − 2πN| ≤ π; then |N| ≤ |b|·|s| + ½. By the functional equation in the direction i/b, used |N| times, Φ_b(y) = Φ_b(y − iN/b)·Π (1 + p^{±(2j+1)}·e^{2πy/b})^{∓1}, 0 ≤ j < |N|.
  - The first factor is exp of Faddeev's integral ∫ e^{−i(w − 2πiN)σ}/(4σ·sinh(πσ)·sinh(πxσ)) dσ over the line Im σ = c₀, with a fixed 0 < c₀ < min(1, Re(1/x)). There |e^{−i(w−2πiN)σ}| ≤ e^{π|Re σ|} (because Re w ≤ 0), and |σ·sinh(πσ)·sinh(πxσ)| ≥ c′·e^{π(1 + Re x)|Re σ|} locally uniformly. So the integral is bounded.
  - In the product, |p^{±(2j+1)}·e^{2πy/b}| ≤ e^{π|Im(1/x)|}·exp(−|s|·[2π·cos(β + δ₀)/|b| − 2π|b|·|Im(1/x)|]). At x₀ the bracket is 2π·cos δ₀/b₀ > 0, so on a neighbourhood every term is ≤ 2e^{−κ|s|} with κ > 0. For |s| ≥ s₀ the product of at most |b||s| + 1 such factors is bounded above and below.
  - For |s| ≤ s₀ the function Φ_b(y) is continuous and without zeros or poles on a compact set.
- *The kernel.* 1 + e^{2πy/b} has no zero on the line, and tends to 1 or to ∞ at its two ends, locally uniformly. ∎

**Proposition E (proved).** For every u ∈ ℂ and every admissible δ₀ the integrals

  I(u) = ∫_{L_{δ₀}} e^{−πiy² + 2πuy}/Φ_b(y) dy,  J(u) = ∫_{L_{δ₀}} e^{−πiy² + 2πuy}/(Φ_b(y)·(1 + e^{2πy/b})) dy

converge absolutely, do not depend on δ₀, and are holomorphic in (u, x) on ℂ × (ℂ \ (−∞, 0]). For Im x > 0 they are the integrals of v12. Hence 𝐒(a, x) is holomorphic in (a, x) for x ∈ ℂ \ (−∞, 0], and det 𝐒 is given by the formula of Theorem 20 there; in particular 𝐒 is invertible.

Proof. |e^{−πiy²}| = e^{−π·sin(2δ₀)·s²} on the line. By Lemma E the integrands are bounded by C·e^{−π·sin(2δ₀)·s² + 2π|u||s|}, locally uniformly in (u, x). Φ_b(y) is holomorphic in b for Re b > 0 (Lemma Φ (ii)), so the integrands are holomorphic in (u, x) for fixed y, and the integrals are holomorphic (Fubini and Morera). Two admissible lines bound two sectors that contain no zero or pole, and the bounds hold uniformly for the lines in between, so the integrals agree (Cauchy). Every x in the cut plane has |β| < π/2, hence an admissible δ₀. ∎

Remarks.

- For real x > 0 this is the statement, used since v20, that the state integrals at real b are the limits of those in the upper half-plane; Lemma E gives it for every x₀ > 0 and every u at once.
- v12 §1 already takes any δ₀ with 0 < δ₀ < π/2 − β, so Theorem Z and Theorem C hold on all of Im x > 0. The 10° line of the later addenda is one admissible choice, for arg x < 160°.
- The negative real axis is a true boundary for 𝐒: as β → ±π/2 the zeros and poles close in on the line from both sides.

---

## 4. The q-frame under Γ₀(2)

**Theorem 21 (proved).** Let Im x > 0, a ∈ ℂ.

- (T) 𝐉(a; x + 1) = C_T·𝐉(a; x),  C_T = [ 1, 0, 0 ; 0, (1+i)/2, (1−i)/2 ; 0, (1−i)/2, (1+i)/2 ].
- (γ) Let x′ = x/(1 + 2x), a′ = a/(1 + 2x). Then

    𝐉(a′; x′) = Λ(a, x)·𝐉(a; x)·Ω_γ(a, x),  Ω_γ(a, x) := 𝐒(a, x)⁻¹·𝐒(a′, x′),  Λ = diag(1, λ, −λ),

    λ(a, x) = δ(a, x)/δ(a′, x′) = (1 + 2x)^{−1/2}·e^{πi/12}·e^{−5πi·xx′/6}·e^{−πi(a − a′)}·e^{−πi·aa′}.
- (continuation) Ω_γ(a, x) is holomorphic in (a, x) for x ∈ ℂ \ (−∞, 0].
- T and γ generate Γ₀(2)/{±1}. So for every g ∈ Γ₀(2) there is a matrix Ω_g with 𝐉(g·(a, x)) = Λ_g·𝐉(a; x)·Ω_g(a, x), Λ_g of the block form diag(1, 2 × 2) and elementary, Ω_T = 1, and Ω_{g₁g₂}(a, x) = Ω_{g₂}(a, x)·Ω_{g₁}(g₂·(a, x)). Each Ω_g is a product of matrices 𝐒^{±1} at points of the orbit.

Proof.

- (T) F does not change. q^{m²/4} is multiplied by i^{m²}, which is 1 for even and i for odd m, so H(ζ) → H_e(ζ) + i·H_o(ζ) = ((1+i)/2)·H(ζ) + ((1−i)/2)·H(−ζ). In the third column q^{−3/4} → i·q^{−3/4} and q^{1/2} → −q^{1/2}, and i·(H_e − i·H_o) = ((1+i)/2)·(H_e + H_o) − ((1−i)/2)·(H_e − H_o) gives the same matrix.
- (γ) 1/x′ = 1/x + 2 and a′/x′ = a/x. So p, ζ̃ and T are the same at (a′, x′) and at (a, x), and p^{1/2} changes sign. Then H̃(y) → H̃(−y), that is H̃_e → H̃_e, H̃_o → −H̃_o, and one reads off 𝐉̃(a′; x′) = E·𝐉̃(a; x) with E = diag(1, 1, −1). Theorem 20 at both points gives Δ(a′, x′)·𝐉(a′; x′)·𝐒(a′, x′)⁻¹ = E·Δ(a, x)·𝐉(a; x)·𝐒(a, x)⁻¹. For λ: b′/b = (1 + 2x)^{−1/2} (principal values; arg x′ = arg x − arg(1 + 2x) ∈ (0, π)), x − x′ = 2xx′, 1/x − 1/x′ = −2, a²/x − a′²/x′ = 2aa′.
- (continuation) x′ ∈ (−∞, 0] exactly for x ∈ (−½, 0], and 1 + 2x ≠ 0. So for x in the cut plane both (a, x) and (a′, x′) are in the domain of Proposition E, and 𝐒(a, x) is invertible.
- (group) Γ₀(2)/{±1} is generated by T and an element of order 2, for which one can take E₂ = γT⁻¹ = (1 −1; 2 −1). The cocycle rule is the associativity of the definition. ∎

Remarks.

- **What this says about F.** The first rows: (1, F(z′; q′), F(z′q′; q′)) = (1, F(z; q), F(zq; q))·Ω_γ(a, x), where z′ = e^{2πia′}, q′ = e^{2πix′}. The matrix Ω_γ is not a q-series; it consists of state integrals. At real x the frames themselves have no meaning (H(·; q) has no limit at roots of unity), but Ω_γ has: this is the quantum-modular statement.
- **A check: E₂² = −1.** −1 acts as a → −a. Applying (T) and (γ) twice gives 𝐉(−a; x) = Λ(a₁, x₁−1)·C_T⁻¹·Λ(a, x−1)·C_T⁻¹·𝐉(a; x)·Ω_γ(a, x−1)·Ω_γ(a₁, x₁−1) with x₁ = (x−1)/(2x−1), a₁ = a/(2x−1). In particular F(1/z) = (1, F(z), F(zq))·(second column of the product of the two Ω's). Verified to 2·10⁻²⁹.
- **Why Γ₀(2) (derivation).** In the family 𝔉(z, y; q) = Σ zⁿq^{n(n+1)/2}(y;q)ₙ of v18, F has y = −q, the dual series 𝒲 has y = (nome)^{1/2}, and x → x + 1 turns that into y = −(nome)^{1/2}. These three are permuted by SL₂(ℤ) like the three even theta characteristics, and Γ₀(2) is the stabiliser of the first. A formula of the type of Theorem Z for general y would give the full group at once; it was not looked for.
- **Weight and index.** λ contains (1 + 2x)^{−1/2} and e^{−πi·a²/(1+2x)}: weight ½ in the convention of Δ, and index −¼ in the Jacobi variable a.

---

## 5. Across the negative real axis (numerical)

Theorem 21 gives Ω_γ on ℂ \ (−∞, 0]. Garoufalidis–Zagier's definition of a holomorphic quantum modular form asks for the domain {x : cx + d ∉ (−∞, 0]} (from memory), which is ℂ \ (−∞, −½] for γ. On (−½, 0) neither 𝐒(a, x) nor 𝐒(a′, x′) exists (x′ < 0 there), so the question is whether the quotient continues.

**Observation (verified, not proved).** Ω_γ(a, x) continues holomorphically across (−½, 0), and not across (−∞, −½).

How it was tested (a = 0.3). Above the axis Ω_γ = 𝐉(a; x)⁻¹·Λ⁻¹·𝐉(a′; x′) by series. Below the axis it is 𝐒(a, x)⁻¹·𝐒(a′, x′) with the integrals on an admissible line (δ₀ = (90° − |β|)/2), Φ_b from the product formula for Im x < 0; the 2 × 2 block was also computed by series (§6).

| test | result |
|---|---|
| jump max\|Ω_γ(x₀ + iε) − Ω_γ(x₀ − iε)\|, first row, ε = 0.08, 0.04, 0.02, at x₀ = −0.15 | 0.42, 0.21, 0.105 |
| the same at x₀ = −0.25 | 0.73, 0.36, 0.175 |
| the same at x₀ = −0.35 | 2.25, 0.88, 0.385 |
| the same at x₀ = −0.7 | 2.2, 1.3, 1.4 (no decrease) |
| the block alone, by series on both sides, ε = 0.004, 0.002, 0.001, at x₀ = −0.1 … −0.45 | halves with ε at every point |
| the block alone at x₀ = −0.6 | 1.7, 1.66, 1.58 |
| Fourier modes of the full Ω_γ on the circle \|x + 0.2\| = 0.1 (16 points): negative modes / mode 0 | 3·10⁻⁶ … 2·10⁻⁵ (aliasing level) |
| the block alone on \|x + 0.2\| = 0.15 (64 points, series) | ≤ 3·10⁻¹⁴ |
| the same on circles around −0.7 and −0.75 | negative modes of order 1 |

A proof would need a representation of Ω_γ that makes sense for x in (−½, 0); I do not have one.

---

## 6. The lower half-plane: |q| > 1

For Im x < 0 one has |q| > 1 and |p| > 1.

**Proposition Q.**

- (i) (proved) L·𝒥 = 1, L·𝒾₀ = L·𝒾₁ = 0 and 𝐒(a + x, x) = 𝐒(a, x)·𝐀(z) hold for all x ∈ ℂ \ (−∞, 0]. So do 𝒥(a) − 𝒥(a+1) = 𝒾₀(a) and the recursion of the 𝒾 in the shift a → a + 1 (v29 §3).
- (ii) (proved) I(u; x) = I(u; 1/x) for fixed u. Hence for Im x < 0

    I(a/b) = c·e^{−πia²/(2x)}·[ H(ζ̃; p⁻²)·H_e(ζq^{−1/2}; q⁻¹) − H(−ζ̃; p⁻²)·H_o(ζq^{−1/2}; q⁻¹) ],

  with (p⁻²)^{1/4} = p^{−1/2} and (q⁻¹)^{1/4} = q^{−1/4}: the 2 × 2 block of 𝐒 is bilinear in series in p⁻² and q⁻¹. In the block, passing from |q| < 1 to |q| > 1 replaces H(±ζ; q) by H_{e,o}(ζq^{−1/2}; q⁻¹) and H_{e,o}(ζ̃p; p²) by H(±ζ̃; p⁻²).
- (iii) (proved) For |q| > 1 the equation y(z) − zq·y(zq) − zq²·y(zq²) = 1 has no solution that is a Laurent series in z converging on some annulus.
- (iv) (derivation) For Im x < 0 neither residue sum of the kernel integral converges: the residues at the upper poles y = i·b·(k+½) are the terms of Σ zᵏq^{k(k+1)/2}(−q;q)ₖ, those at the lower poles the terms of F̌(1/z), and both grow like a Gaussian in k when |q| > 1.

Proof. (i) Both sides are holomorphic in x on the cut plane (Proposition E) and agree for Im x > 0. (ii) Φ_b = Φ_{1/b} as functions of y, and a line admissible for b is admissible for 1/b; so the integral for x equals the one for 1/x, and Theorem C applies at 1/x with â = a/x. The constant c and the Gaussian are given by the same expressions. (iii) For y = Σ cₙzⁿ the equation is cₙ = qⁿ(1 + qⁿ)·cₙ₋₁ + δ_{n,0}. If cₙ ≠ 0 for some n ≥ 0, then |cₙ| grows like |q|^{n²} as n → +∞. Otherwise c₋₁ = −½ and c₋ₘ₋₁ = c₋ₘ·qᵐ/(1 + q⁻ᵐ), so |c₋ₘ| grows like |q|^{m²/2}. In both cases the series diverges for every z ≠ 0. ∎

**What follows (inference).** In GGMW the three rows of the matrix are q-hypergeometric series that converge for |q| < 1 and for |q| > 1, and the cocycle is a quotient of series on both sides. Here this is true for the 2 × 2 block only. The trivial series F(z) = Σ zⁿq^{n(n+1)/2}(−q;q)ₙ diverges for |q| > 1, the equation has no other series solution there, and the first row of 𝐒 is given by the integral alone: 𝒥(a) = i·P·e^{−πia}·J(a/b), with J(u; x) equal to the integral with the mirrored kernel 1/(1 + e^{2πyb̂}) at x̂ = 1/x. So the "|q| > 1 side" of the open list has this answer for the block, and for the trivial row the answer that the state integral is the only object there. Whether it has another closed form (for instance as a sum over a lattice in the dual variable) is open.

---

## 7. Numerical checks (verified; scripts in `psi_star/v30/`, output in `run_log.txt`)

| check | result |
|---|---|
| §1: (A0), (A1), (R1)–(R4) and four further relations (`r1_letter.py`) | exact to q³⁰⁰ |
| Theorem 20: max\|𝐉̃·𝐒 − Δ·𝐉\|, 𝐒 by quadrature, 3 points x (one with Re x < 0), a = ¼, ½, 0.3 + 0.1i, 0.8 (`m1_matrix.py`) | ≤ 7·10⁻²⁹ |
| the same at a = 1.7 (the quadrature is less accurate) | ≤ 2·10⁻²⁰ |
| det 𝐒 against the formula; det 𝐉 = −2μ; P elementary | ≤ 5·10⁻²⁹ (a = 1.7: 2·10⁻²¹); 4·10⁻³¹; 2·10⁻³¹ |
| the Corollary: 𝒲(T) and p·T·𝒲(p²T) from F(z), F(zq) and the six state integrals | ≤ 10⁻³⁰ (three points) |
| Lemma W and the equation of H, for r = q and r = p² | ≤ 4·10⁻³¹ |
| Proposition E: integrals in double precision, Φ_b from Faddeev's integral, against the series of Theorem 20, Im x > 0 (`m2_extend.py`) | ≤ 5·10⁻¹⁴ |
| Fourier modes of 𝐒(a, x) on circles around x₀ = 0.5, 0.3, 1.3: negative modes k = 1 … 4 | ≤ 10⁻¹¹, ≤ 10⁻⁸, ≤ 2·10⁻¹⁵ (positive modes of order 1) |
| Proposition Q (i): L·𝒥 − 1, L·𝒾₀, L·𝒾₁ at x = 0.5 − 0.15i, 0.3 − 0.1i, 0.8 − 0.3i (\|q\| up to 6.6) | ≤ 3·10⁻¹² |
| Proposition Q (ii): I(u; x) = I(u; 1/x); the block by integrals against the series in q⁻¹, p⁻² | ≤ 3·10⁻¹⁴; ≤ 4·10⁻¹¹ |
| Theorem 21 (T) (`m3_modular.py`) | ≤ 2·10⁻²⁹ |
| Theorem 21 (γ): max\|𝐉(a′;x′) − Λ·𝐉·Ω_γ\|, four points; λ against δ/δ′ | ≤ 2·10⁻²⁹; 4·10⁻³¹ |
| the relation E₂² = −1, two points; F(1/z) from F(z), F(zq) | ≤ 2·10⁻²⁹; 15 digits shown |
| §5 (`m4_block.py`, `m5_jump.py`, `m6_circle.py`) | see the table there |

A correction to a tool: `fad.py` (v23–v27) evaluates Faddeev's integral on the line Im σ = ½, which is only valid for Re(1/x) > ½. The copy in `psi_star/v30/` takes the line Im σ = min(1, Re(1/x))/2 with a step proportional to it. The earlier addenda used it for 0 < x ≤ 1 only, where nothing changes.

---

## 8. Open

1. A proof that Ω_γ continues across (−½, 0) (§5).
2. The formula of the type of Theorem Z for 𝔉(z, y; q) with general y, which would replace Γ₀(2) by SL₂(ℤ) acting on (a, the logarithm of y, x).
3. Whether the entries of 𝐒⁻¹, which are Casoratians of state integrals, are single state integrals at the dual modulus.
4. A closed form of the trivial row for Im x < 0 other than the integral (§6).
5. The comparison of (A0), (A1) with the printed relations (§1).
6. Unchanged from the list of 28 Nov: odd and general roots of unity for z ≠ −1; Borel summability of the trivial sector for z ≠ −1; constants uniform near z = 1; a direct proof of Theorem 17.

---

## 9. Scripts (`psi_star/v30/`)

- `r1_letter.py [N]`: §1, exact series (under a second).
- `frames.py`: the two frames, Δ, and 𝐒 by quadrature (mpmath). `m1_matrix.py`: Theorem 20, determinants, Lemma W (3 minutes).
- `npint.py`: the state integrals for Re x > 0 in double precision, Φ_b from Faddeev's integral (`fad.py`). `m2_extend.py`: Proposition E and Proposition Q (2 minutes).
- `m3_modular.py`: Theorem 21 and the relation E₂² = −1 (25 minutes).
- `m4_block.py`: the 2 × 2 block of Ω_γ on both half-planes by series; Fourier modes and jumps (2 minutes).
- `lowint.py`: the state integrals for Im x < 0 in double precision, Φ_b from the product formula. `m5_jump.py`, `m6_circle.py`: §5 (5 minutes together).
- `common.py`, `a1_formula.py` (from v18; `common.py` with a longer line for small δ₀), `fad.py` (from v27, corrected as said in §7), `run_log.txt`, `README.md`.
