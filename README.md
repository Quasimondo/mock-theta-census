# mock-theta-census

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24443194.svg)](https://doi.org/10.5281/zenodo.24443194)

A complete numerical census of short Ramanujan-style q-series, classified as **modular**, **mock theta**, **false/quantum type** or **neither** by their radial asymptotics at roots of unity. Along the way it found one pair that does not appear in the literature we checked:

  **ψ₁\*(q) = Σₙ≥0 (−1)ⁿ q^{n(n+1)/2} (−q;q)ₙ**, **ψ₀\*(q) = Σₙ≥0 (−1)ⁿ q^{(n+1)(n+2)/2} (−q;q)ₙ**

These are Ramanujan's fifth-order mock theta functions ψ₁, ψ₀ with a sign twist (−1)ⁿ. They are bounded at every root of unity tested and terminate at even-order roots. At q = −1 they have a divergent asymptotic series with integer coefficients, whose large-order behavior is controlled by Ramanujan's entire Rogers–Ramanujan function at a = 1. Their structure is now largely understood: each is half of a sign-twisted fifth-order function plus or minus θ₄ times Ramanujan's entire Rogers–Ramanujan function (identities (I), (II)), and near q = 1 they obey a quantum modularity law of Kashaev-invariant type whose dual is a new finite sum W (item 5). As far as we could find, neither the functions, the identities nor W appear in the literature. **We'd be glad to hear if they do.**

Two later findings extend this. The same sign twist gives bounded, non-trivial companions for several other mock theta functions (fifth-order F₀, F₁, eighth-order T₀); for the third-order ψ it gives a known function, −σ\*/2 of Andrews–Dyson–Hickerson. And a census of 6,327 double sums found one mock-theta-like double sum, **D(q) = Σ_{m,n≥0} q^{m²−mn+n²} / ((−q;q)ₘ (−q;q)ₙ)**, which grows exactly like Ramanujan's third-order f(q) but has golden-ratio asymptotics at q = 1. It is now explained by a proved identity: **D = θ₄·f + 4R**, with R = Σ_{d>k≥0} (−1)^{d+k+1} q^{d²}/(−q;q)ₖ a sum of theta tails.

**Interactive pages:** [quasimondo.github.io/mock-theta-census](https://quasimondo.github.io/mock-theta-census/)
- [Explorer](https://quasimondo.github.io/mock-theta-census/explorer-v4.html): mock theta functions on the unit disk, with a radial probe toward any root of unity, including ψ₁\* and ψ₀\* and the double sum D
- [Fingerprints](https://quasimondo.github.io/mock-theta-census/fingerprints-v3.html): the classifier on 51 known functions, both censuses, the sign twists, and the double sum D
- [Formula sheet](docs/formulas.md): every formula in LaTeX, Unicode, plain text and Mathematica

All results here are **numerical evidence and conjectures** unless marked as proved. Computations used Arb ball arithmetic (python-flint) at 50+ significant digits, and exact integer arithmetic where possible.

## Main results

1. **A validated classifier.** Along q = e^{−t}ζ, fit log|g| = c/t + α log t + β + γt (+ exponentially small terms) at 19 roots of unity ζ, and test whether κ = c·m²/4π² is rational and whether the fit is exact.
   - On a labeled catalog of 51 known functions, it classifies every one correctly except a sum of two modular forms with different q-offsets, which is a documented limitation.
   - It recovers known normalizations (q^{−1/24}f, q^{−1/60}f₀, q^{−1/168}F₀ of order 7) and dilogarithm values of non-modular Nahm sums to 11–12 digits.
2. **A census.** All 9,936 formulas Σ ±q^{(An²+Bn)/2} × (one or two q-product factors from a fixed list) are classified. They reduce to 8,869 distinct functions: **25 mock, 28 modular, 1,357 bounded at all 19 roots**, and the rest rejected.
   - 23 of the 25 mock hits are known functions, rescalings, or specializations of universal mock theta functions. This includes the elementary identity **H(x) = (1−x²)·g₂(x,q) + 1 − x** (proved by telescoping), which identifies two hits as g₂(∓ω, q).
   - The other two hits are ½ψ + ½ψ\*, i.e. mock only in Ramanujan's original sense.
3. **Identities for ψ\*** (verified exactly to q³⁹⁹; proof sketch in `docs/addendum-v2.md`):
   - 2ψ₀\* = f₀\* − θ₄·R₀
   - 2ψ₁\* = f₁\* + θ₄·R₁

   with f₀\* = Σ(−1)ⁿq^{n²}/(−q;q)ₙ, f₁\* = Σ(−1)ⁿq^{n²+n}/(−q;q)ₙ, θ₄ = Σ(−1)ⁿq^{n²}, R₀ = Σ(−1)ⁿq^{n²}/(q;q)ₙ, R₁ = Σ(−1)ⁿq^{n²+n}/(q;q)ₙ. These are the z = −1 analogues of the classical bilateral-series identities for f₀ + 2ψ₀ and f₁ + 2ψ₁; the method of relating fifth-order functions to modular bilateral series is that of Bajpai, Kimport, Liang, Ma and Ricci (*Bilateral series and Ramanujan's radial limits*, Proc. AMS 2015).
4. **Asymptotics at q = −1** (numerical, 15–18 digits).
   - ψ₁\*(−e^{−t}) = Σ cₖ tᵏ/k! with integer cₖ = 1, −1, 27, −1237, 106587, −14667181, …
   - Borel singularity at −π²/60, with Stokes data π·K² = (5 ± 3√5)/40 and subleading coefficients such as (13 − 6√5)/60.
   - These match the exponential asymptotics of R₀, R₁ at q → −1 exactly.
   - The same holds at q = i (singularity at −π²/240, constants in ℚ(√5, i)). At q → 1 the singularities are at π²/3 ∓ i·Cl₂(π/3), with closed-form Stokes constants (`docs/addendum-v3.md` §3).
   - Radial limits at every root of unity are finite. At even order m the series terminates; at odd order m the terms repeat with a factor −2 every m steps, and the limit is one third of the first m terms. **Proved** in `docs/addendum-v3.md` §1 (from (II): f₁\* converges at the root with ratio −½, θ₄R₁ → 0, and reversing the finite sum gives the stated one; the analytic estimates are a sketch).
5. **Transformation behavior.** Under x ↦ x/(2x+1) (Γ₀(2)), ψ₁\* and ψ₀\* transform *together*, to ~10⁻⁹ on held-out roots. Unrelated pairs reach only 10⁻⁴ to 10⁻¹, and the Kontsevich–Zagier control reaches 10⁻¹⁰. But the best-fit weight is ≈ 0.508, not ½, so this is not a classical weight-½ quantum modular form with constant multipliers. `docs/addendum-v3.md` §5 explains why: under x ↦ 1/x the pair does not go to itself but to a **new finite sum W(ξ) = Σ(−1)ʲ ξ^{j²+j}(ξ;ξ²)ⱼ at the dual root**. As X → ∞, ψ₁\*(e^{2πi/X}) = C̄·t^{−1/2}·e^{−S̄/t}·W(e^{−πiX})·(1 + s̄₁t + …) + O(1), with t = −2πi/X, S̄ = π²/3 + i·Cl₂(π/3), C̄ = √(π/√3)·e^{5πi/12}, s̄₁ = 1/4 + i√3/54 (verified on 56 residue classes to 6·10⁻²¹, independently rechecked). This is quantum modularity of the Kashaev-invariant / Garoufalidis–Zagier type, with 2·Cl₂(π/3), the volume of the figure-eight knot complement, as the growth rate. `docs/addendum-v4.md` adds the other side: W = G(−q) has exactly one companion W₃ = G(−q³) in the family G(z) = Σ zʲq^{j²}(q;q²)ⱼ, and the pair satisfies bilateral identities (V), (VI) of the same type as (I), (II), with Ramanujan's ψ(q) = Σq^{n(n+1)/2} as prefactor (proved). In the converse direction (W, W₃) goes back to (ψ₀\*, ψ₁\*) at the dual root through a 2×2 matrix with golden-ratio entries and Rogers–Ramanujan phases, and near q = −1 the pair (ψ₁\*, ψ₀\*) goes to itself under x ↦ (x−1)/(2x−1) ∈ Γ₀(2), with weight ½ (both verified numerically, 10⁻⁹ to 10⁻¹¹, independently rechecked to ~10⁻⁵). `docs/addendum-v5.md` tests all three laws inside the unit disk (along horocycles, where everything converges) to about 10⁻²⁹ and finds the third sector: in the forward law the decaying term e^{−S/t} has coefficient function q̃·W₃(q̃), so the three sectors (trivial, growing, decaying) of ψ\* pair with the basis (1, W, q̃W₃). It also finds that R₀, R₁ and P₊, P₋ are exchanged by the Fricke involution. All verified numerically. `docs/addendum-v6.md` turns this into a proof programme: R₀, R₁, P₊, P₋ (and the X, Y of (III), (IV)) are all values of one Nahm-type sum H(z;q) = Σ zᵐq^{m²/4}/(q;q)ₘ, with exact Casoratian identities X·N₀ + N₁·Y = 1 and E₊′E₋ + E₋′E₊ = 2 (proved). The law for R₀, R₁ is then an exact formula plus one analytic claim: four bilinear combinations 𝓑₁…𝓑₄ of q-series in q and the dual nome do not depend on the dual nome (verified to 10⁻¹⁷). `docs/addendum-v7.md` shows that they are exactly one-dimensional state integrals ∫Φ_b(y)·e^{−πiy²/2 + 2πu·y}dy of Faddeev's quantum dilogarithm, with b² = x and four explicit values of u (residue derivation; verified to 10⁻²⁷ … 10⁻³¹ and rechecked independently). Since such an integral does not see the dual nome and is analytic across the real axis, this reduces the law for the bilateral parts θ₄R₀, θ₄R₁ of ψ\* to two standard facts about state integrals: justifying the contour closing (the A = ½ analogue of Garoufalidis–Kashaev) and stationary phase near b = 0. Neither is written out rigorously yet. `docs/addendum-v8.md` tests this structure on F₀\*, F₁\*, which were not used to find it: the growing sector of their bilateral parts follows from the same four integrals with no free parameter (verified), and the one-sided halves go to 𝒥(r) = Σ r^{2m²+m}/(−r²;r²)ₘ and 𝒦(r) = 2Σ r^{(n+1)²}(−r²;r²)ₙ, unsigned fifth-order-type sums at a half shift, with a new bilateral identity (VII) 𝒥 + 𝒦 = ψ(r)(r²;r⁴)∞·Σ r^{m(m+1)/2}/(r²;r²)ₘ (exact to r³²⁹, not yet proved). `docs/addendum-v9.md` closes the loop: near r = 1, 𝒥 goes back to the pair (q^{−1/120}(F₀\* − 1), q^{71/120}F₁\*) at the dual nome, the exponents of Ramanujan's unsigned fifth-order pair (F₀ − 1, F₁) with the sign-twisted functions in their place (verified, rechecked independently). So both families, ψ\* ↔ W and F\* ↔ 𝒥, are closed two-way dualities. `docs/addendum-v10.md` turns the one-sided law into an exact formula at finite 𝒳: ψ₁\* = 𝒮[Φ](t) + (κ/2)·[𝓑₂·p·W₃(p) − 𝓑₄·W(p)], with 𝒮[Φ] the Borel sum of the perturbative series and 𝓑₂, 𝓑₄ the state integrals (verified to ~10⁻⁵⁰, rechecked independently). `docs/addendum-v11.md` gives that handle: the trivial-sector term is itself a state integral, (−q;q)∞/(b·q^{1/8}(p;p²)∞)·∫ e^{−πiy² + πy/b}/(Φ_b(y)(1 + e^{2πy/b}))dy, whose naive residue sum would be ψ₁\* itself (verified to ~10⁻²⁰ and rechecked independently to 10⁻³⁶ … 10⁻⁴¹). So every term of the three-sector formula is an explicit integral times a q-series in the dual nome. `docs/addendum-v12.md` **proves** this exact formula, and its analogues for ψ₀\*, f₁\*, f₀\*: moving the contour across the kernel poles gives the one-sided sum, the rest is a geometric series of kernel-free integrals obeying a three-term recursion whose generating sums are W, W₃ (or F_a, F_b), and those integrals are the state integrals of v7 after two Fourier transforms of Φ_b (which also completes the proof of v7). The key identities were rechecked independently to 10⁻³¹; the proof has been read but not refereed. What is still open: that the kernel integral is the Borel sum of Φ, and hence the asymptotic laws of v5 as theorems. `docs/addendum-v13.md` gives Φ in closed form (Bernoulli polynomials, 301 exact coefficients) and shows that the kernel integral continues across the real axis into |q| > 1, where ψ₁\* does not exist, and equals the Borel sum of Φ on almost all of the cut plane (verified to 10⁻²⁸ … 10⁻³¹, rechecked independently at two points to 10⁻³⁷). `docs/addendum-v14.md` and `docs/addendum-v15.md` **prove** the asymptotic statement P·J ~ Φ on closed subsectors of −210° ≤ arg t ≤ 10° on both sides of the real x-axis (the key step is an elementary descent inequality along complex straight contours), and so the first proved asymptotic law for the one-sided sums. `docs/addendum-v16.md` closes the wedge around the real x-axis using Faddeev's integral representation of Φ_b, which is one analytic function across real b: P·J is analytic and has the expansion Φ uniformly on the closed sector −210° ≤ arg t ≤ 10°, an opening of 220° that includes |q| = 1. `docs/addendum-v17.md` supplies the factorial constants (|P·J − Σ_{k<K} aₖtᵏ| ≤ A·K!·(240|t|)ᴷ on that sector, by moving the contour across N ≈ 0.04/|t| kernel poles) and concludes, via the Nevanlinna–Sokal theorem, that **P·J is the Borel sum of Φ** in the directions −120° … −80°. So the exact formulas of v12 become "Borel sum of the perturbative series + state-integral terms" as theorems. Still open: Borel directions outside [−120°, −80°], the asymptotics of 𝓑ⱼ (needed to turn this into the asymptotic laws), Theorem A for |q| > 1, the Stokes jump, identity (VII). The whole chain (v12, v14–v17) is written and checked by one author and rechecked numerically here; it has not been refereed. A handoff summarising the state is in `docs/psi-star-handoff-v17.md`.

   **Relation to the figure-eight knot** (`psi_star/knot/`): ψ₁\*(e^{2πi/X}) has exactly the exponential growth and the torsion 3^{−1/4} of the Kashaev invariant of 4₁ at N = X/2 (|ψ₁\*|/J_{X/2} → 2/X), but weight ½ instead of 3/2 and a different perturbative series, so it is not that invariant in disguise. Whether it is one of the other q-series attached to 4₁ is open.

6. **Sign twists: a pattern beyond ψ\*.** Inserting (−1)ⁿ into the other catalog mock theta functions gives more companions that are bounded at all 19 roots and are not a known function at ±q: the 5th-order F₀, F₁ and the 8th-order T₀, and the family F(z) at z = i, −ω. The twist of the 3rd-order ψ is also bounded, but it is a known function: −σ\*(q)/2, the Andrews–Dyson–Hickerson companion of σ(q) (ℚ(√6), Maass forms, quantum modular). Three more twists turn out to be false theta functions. None of the dense bounded twists becomes sparse after multiplying by an eta quotient, so none of them is a single Hecke-type double sum over such a product (`census3/hecke/`). F₀\*, F₁\* do follow the ψ\* pattern, with a different ingredient (proved, `docs/addendum-v3.md` §4): (III) F₀\* + Σₙ≥₁ q^{n²}(q;q²)ₙ = θ₄(q²)(−q;q)∞·X(q²) and (IV) qF₁\* − Σₙ≥₀ q^{(n+1)²}(q;q²)ₙ = θ₄(q²)(−q;q)∞·Y(q²), with X = Σ(−1)ⁿq^{n²}/(q;q)₂ₙ, Y = Σ(−1)ⁿq^{(n+1)²}/(q;q)₂ₙ₊₁. Details, the full table and first checks: [`docs/twists.md`](docs/twists.md).
7. **The double sum d1151, and an identity for it.** The rank-2 census left one candidate,

     D(q) = Σ_{m,n≥0} q^{m²−mn+n²} / ((−q;q)ₘ (−q;q)ₙ).

   - **Observed:** the growth fingerprint of Ramanujan's third-order f(q) (κ = 1/24 at all 12 even-order roots), and at q = 1 an exact expansion whose large-order behaviour carries golden-ratio constants (π²/60, π·K² = 8(5+√5)/25, s₁ = (35+6√5)/300; numerical).
   - **Proved:** D = θ₄(q)·f(q) + 4R(q) with R = Σ_{d>k≥0} (−1)^{d+k+1} q^{d²}/(−q;q)ₖ. Each diagonal m − n = d of D is ±f plus a finite sum of 1/(−q)ₖ (a three-line telescoping lemma), and the signs are those of θ₄. The same argument gives an identity for the whole family Σ q^{m²−mn+n²}/((zq)ₘ(q/z)ₙ) in terms of the rank generating function, a theta series and tail sums. So D is not a new mock theta function in the strict sense; whether the general-z identity is known is open.

   Details: [`docs/d1151.md`](docs/d1151.md) (how it was found), [`docs/d1151-identity.md`](docs/d1151-identity.md) (the identity and its proof). Census results: `census2/results/`.
8. **One-parameter families at roots of unity** (`census3/`). F(z), g₂(x,q) and g₃(x,q) at all z, x of order ≤ 12 (136 functions).
   - g₂, g₃ behave like mock theta functions everywhere, with κ = (1−4a)²/16 and (1−6a)²/24 at q = 1 for x = e^{2πia} below a = 1/4, 1/6.
   - F(z) is bounded at even-order roots. At odd order m it grows exactly when z^m is closer than a₀ = 0.107366 of a turn to 1, so it is bounded everywhere only for z = −1, ±i, −ω, −ω̄. a₀ is where a dilogarithm saddle-point value changes sign; its predicted κ matches the fits to 4–5 digits (heuristic).
   - Negative results: ψ₀\*, ψ₁\* have no q → 1/q partner (like σ ↔ σ\*) among 1,110 testable bounded census functions (`census3/qm/`).

## Open questions

- Are ψ₀\*, ψ₁\*, the family F(z) = Σ zⁿq^{n(n+1)/2}(−q;q)ₙ (which satisfies F(z) = 1 + zqF(zq) + zq²F(zq²)), or the identities in (3) known?
- The full matrix of the quantum modularity law: the first row (dual W) is found; the decaying third sector, the asymptotics of W itself and a proof by stationary phase are open (`docs/addendum-v3.md` §5.6, §7).
- Rigorous proof for the bilateral parts: the contour closing for the A = ½ state integral and the stationary-phase expansion with its phase (addendum-v7 §4). The one-sided halves f\*, F_a, F_b need a separate (inhomogeneous) argument, possibly a state integral with one more factor; the converse and mirror directions should follow from the same four integrals but are not checked; proof of (VII) and the decaying sector for F\* are open.
- A proof of the F(z) threshold rule (item 8), including the reduction from a root of order m to z^m.
- How do X, Y of identities (III), (IV) behave at roots of unity, and do F₀\*, F₁\* have a dual like W? Identities for the other bounded twists T₀\*, F(i), F(−ω) are still missing.
- Is the general-z identity for Σ q^{m²−mn+n²}/((zq)ₘ(q/z)ₙ) (`docs/d1151-identity.md`, section 3) known?

## Layout

| Folder | Contents |
|---|---|
| `census/` | Evaluator and labeled catalog (`catalog.py`), per-root probe (`probe.py`), generator and dedupe (`search.py`), staged search (`stage.py`, `filt.py`, `dclass.py`), q→1 dilogarithm pre-filter (`dilog.py`), identification tools (`ident.py`, `shadow.py`, `companions.py`, …), plus all result data (`cands.json`, `labels.json`, `stageC_summary.json`, `*.jsonl`) |
| `census2/` | Rank-2 (double-sum) census: generator (`gen2.py`), evaluators (`sum2_fast.py`, `sum2.py`), staged runner (`run2.py`), report (`report2.py`); results of the first run in `census2/results/`, the analysis of d1151 in `census2/d1151/`; sign-twist scripts are in `census/twist.py` and `census/twist_checks/` |
| `psi_star/` | Quantum-modularity tests with Kontsevich–Zagier control (`qmf*.py`); exact identities, expansion at q = −1 and large-order analysis (`s1`–`s11`) |
| `census3/` | One-parameter families F(z), g₂(x,q), g₃(x,q) at all roots of unity of order ≤ 12 (`zfam.py`, `zreport.py`, `zextra.py`, results and `README.md`); `hecke/`: sparsity test of the bounded twists against eta quotients and an f_{a,b,c} matcher; `qm/`: radial limits of the 1,357 bounded census functions and the q → 1/q duality search |
| `docs/` | `d1151.md` (the double sum d1151, how it was found), `d1151-identity.md` (its identity and proof), `d1151-brief.md` (the brief that led to it); `twists.md` (sign twists of the catalog functions); GitHub Pages site (`index.html`, `explorer-v4.html`, `fingerprints-v3.html`); `formulas.md` (every formula in LaTeX, Unicode, plain text and Mathematica, with equivalent forms and coefficients), `handoff-v1.md` (full account of method and results), `addendum-v2.md` (identities, proofs, resurgence analysis), `addendum-v3.md` (odd-root rule, identities (III), (IV), resurgence at q = i and q → 1, the law with W), `addendum-v4.md` (the W family, identities (V), (VI), converse and mirror laws), `addendum-v5.md` (the laws inside the disk, the third sector), `addendum-v6.md` (one function H behind both sides, Casoratians), `addendum-v7.md` (the bilinear combinations as state integrals), `addendum-v8.md` (the law for F₀\*, F₁\*, identity (VII)), `addendum-v9.md` (𝒥 back to F₀\*, F₁\*; status table), `addendum-v10.md` (the one-sided law as an exact formula with a Borel sum), `addendum-v11.md` (the trivial sector as a state integral with a kernel), `addendum-v12.md` (proof of the exact formula), `addendum-v13.md` (closed form of Φ, the kernel integral as Borel sum on both sides of the real axis), `addendum-v14.md`, `addendum-v15.md`, `addendum-v16.md` (P·J ~ Φ proved on a 220° sector), `addendum-v17.md` (factorial constants; P·J is the Borel sum of Φ), `psi-star-handoff-v17.md` (status of v3–v17 for a fresh session), `psi-star-brief.md`, `psi-star-brief-v2.md` (the briefs that led to them) |
| `psi_star/formulas_check.py` | Checks every equivalent form in `docs/formulas.md` to 50 digits |
| `psi_star/s12_odd_limits.py` | Radial limits of ψ₁\*, ψ₀\* at odd-order roots vs. the closed form (first period)/3 |
| `psi_star/v17/` | Scripts of `docs/addendum-v17.md`, plus the independent check `check_v17.py` |
| `psi_star/v16/` | Scripts of `docs/addendum-v16.md` (`fad.py`: Faddeev's integral), plus the independent check `check_faddeev.py` |
| `psi_star/v15/` | Scripts of `docs/addendum-v15.md`, plus independent checks (`check_prop1.py`, `check_thm2.py`) |
| `psi_star/v14/` | Scripts of `docs/addendum-v14.md`, plus the independent check `check_thm1.py` |
| `psi_star/v13/` | Scripts of `docs/addendum-v13.md` (`phi_exact.json`: 301 exact coefficients), plus the independent check `check_borel.py` |
| `psi_star/v12/` | Scripts of `docs/addendum-v12.md` (`run_all.sh`, `t1`–`t5`, `common.py`), plus the independent check `check_thmAC.py` |
| `psi_star/v11/` | Script of `docs/addendum-v11.md`, plus the independent check `check_kernel.py` (with `tfun.py`) |
| `psi_star/v10/` | Scripts of `docs/addendum-v10.md`, plus independent checks (`phi_psi1.py`, `check_exact.py`) |
| `psi_star/knot/` | Comparison of ψ₁\* with the Kashaev invariant of the figure-eight knot (`fig8.py`, `README.md`) |
| `psi_star/v9/` | Scripts of `docs/addendum-v9.md`, plus independent checks (`phiJ.py`, `check_Jloop.py`) |
| `psi_star/v8/` | Scripts of `docs/addendum-v8.md`, plus independent checks (`check_VII.py`, `check_F0_law.py`) |
| `psi_star/v7/` | Script of `docs/addendum-v7.md` (`state_integral.py`), plus the independent check `check_state_integral.py` |
| `psi_star/v6/` | Scripts of `docs/addendum-v6.md`, plus independent checks (`check_casoratian.py`, `check_B1.py`) |
| `psi_star/v5/` | Scripts of `docs/addendum-v5.md`, plus independent checks (`check_M.py`, `check_psi_mod.py`, `check_converse_disk.py`) |
| `psi_star/v4/` | Scripts of `docs/addendum-v4.md`, plus independent checks (`check_V_VI.py`, `check_laws.py`, `conv_trend.py`, `check_W_stokes.py`) |
| `psi_star/v3/` | Scripts of `docs/addendum-v3.md` (`run_log.txt` is their output), plus the independent checks `check_III_IV.py` and `check_law.py` |

## Reproducing

```bash
pip install -r requirements.txt
cd census
python3 catalog.py                       # identity checks of the catalog formulas
python3 probe.py f theta3 "Rogers-Ramanujan G" "Nahm A=4,B=0" | python3 fmt.py
cd ../psi_star
python3 s1_identities.py 400             # identities (I), (II) exactly to q^400
python3 s2_expand.py 120 && python3 s6_largeorder.py && python3 s7_id.py
python3 qmf_scale.py                     # transformation test with Kontsevich control
cd v3
python3 check_III_IV.py && python3 check_law.py   # identities (III), (IV) and the law with W
```

The full census takes a few hours on 2 cores. `stage.py` is resumable, and `census/labels.json` holds the final label for every candidate, so nothing needs re-running to inspect results. Scripts expect to run from their own folder.

## How this was made

This was an exploratory project by Mario Klingemann, an artist and creative technologist, not a professional mathematician. The computations and much of the code were developed in conversation with AI assistants (Anthropic's Claude Opus 5.5 and Claude Fable 5.2). Every claim in the documents is backed by a script in this repository. Corrections are very welcome.

## How to cite

Mario Klingemann, *mock-theta-census: a numerical census of Ramanujan-style q-series and the sign-twisted fifth-order functions ψ₀\*, ψ₁\**, Zenodo, 2026. https://doi.org/10.5281/zenodo.24443194 (all versions; each release also has its own version DOI on Zenodo)

GitHub's "Cite this repository" button (from `CITATION.cff`) gives the same reference in BibTeX and APA.

## License

Code: MIT (see `LICENSE`). Data and documents: CC BY 4.0.
