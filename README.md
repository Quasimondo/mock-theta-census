# mock-theta-census

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24443194.svg)](https://doi.org/10.5281/zenodo.24443194)

A complete numerical census of short Ramanujan-style q-series, classified as **modular**, **mock theta**, **false/quantum type** or **neither** by their radial asymptotics at roots of unity. Along the way it found one pair that does not appear in the literature we checked:

  **ψ₁\*(q) = Σₙ≥0 (−1)ⁿ q^{n(n+1)/2} (−q;q)ₙ**, **ψ₀\*(q) = Σₙ≥0 (−1)ⁿ q^{(n+1)(n+2)/2} (−q;q)ₙ**

These are Ramanujan's fifth-order mock theta functions ψ₁, ψ₀ with a sign twist (−1)ⁿ. They are bounded at every root of unity tested and terminate at even-order roots. At q = −1 they have a divergent asymptotic series with integer coefficients, whose large-order behavior is controlled by Ramanujan's entire Rogers–Ramanujan function at a = 1. As far as we could find, neither the functions nor the identities below appear in the literature. **We'd be glad to hear if they do.**

Two later findings extend this. The same sign twist gives bounded, non-trivial companions for several other mock theta functions (fifth-order F₀, F₁, eighth-order T₀); for the third-order ψ it gives a known function, −σ\*/2 of Andrews–Dyson–Hickerson. And a census of 6,327 double sums found one mock-theta-like double sum, **D(q) = Σ_{m,n≥0} q^{m²−mn+n²} / ((−q;q)ₘ (−q;q)ₙ)**, which grows exactly like Ramanujan's third-order f(q) but has golden-ratio asymptotics at q = 1. It is now explained by a proved identity: **D = θ₄·f + 4R**, with R = Σ_{d>k≥0} (−1)^{d+k+1} q^{d²}/(−q;q)ₖ a sum of theta tails.

**Interactive pages:** [quasimondo.github.io/mock-theta-census](https://quasimondo.github.io/mock-theta-census/)
- [Explorer](https://quasimondo.github.io/mock-theta-census/explorer-v4.html): mock theta functions on the unit disk, with a radial probe toward any root of unity, including ψ₁\* and ψ₀\* and the double sum D
- [Fingerprints](https://quasimondo.github.io/mock-theta-census/fingerprints-v3.html): the classifier on 51 known functions, both censuses, the sign twists, and the double sum D
- [Formula sheet](docs/formulas.md): every formula in LaTeX, Unicode, plain text and Mathematica

All results here are **numerical evidence and conjectures** unless marked as proved. Computations used Arb ball arithmetic (python-flint) at 50+ significant digits, and exact integer arithmetic where possible.

## Main results

1. **A validated classifier.** Along q = e^{−t}ζ, fit log|g| = c/t + α log t + β + γt (+ exponentially small terms) at 19 roots of unity ζ, and test whether κ = c·m²/4π² is rational and whether the fit is exact.
   - On a labeled catalog of 50 known functions, it classifies every one correctly. The one sum of two modular forms with different q-offsets is the expected exception.
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
   - Borel singularity at π²/60, with Stokes data π·K² = (5 ± 3√5)/40 and subleading coefficients such as (13 − 6√5)/60.
   - These match the exponential asymptotics of R₀, R₁ at q → −1 exactly.
   - Radial limits at every root of unity are finite. At even order m the series terminates; at odd order m the terms repeat with a factor −2 every m steps, and the limit is one third of the first m terms. **Proved** in `docs/addendum-v3.md` §1 (from (II): f₁\* converges at the root with ratio −½, θ₄R₁ → 0, and reversing the finite sum gives the stated one; the analytic estimates are a sketch).
5. **Transformation behavior.** Under x ↦ x/(2x+1) (Γ₀(2)), ψ₁\* and ψ₀\* transform *together*, to ~10⁻⁹ on held-out roots. Unrelated pairs reach only 10⁻⁴ to 10⁻¹, and the Kontsevich–Zagier control reaches 10⁻¹⁰. But the best-fit weight is ≈ 0.508, not ½, so this is not a classical weight-½ quantum modular form with constant multipliers. `docs/addendum-v3.md` §5 explains why: under x ↦ 1/x the pair does not go to itself but to a **new finite sum W(ξ) = Σ(−1)ʲ ξ^{j²+j}(ξ;ξ²)ⱼ at the dual root**. As X → ∞, ψ₁\*(e^{2πi/X}) = C̄·t^{−1/2}·e^{−S̄/t}·W(e^{−πiX})·(1 + s̄₁t + …) + O(1), with t = −2πi/X, S̄ = π²/3 + i·Cl₂(π/3), C̄ = √(π/√3)·e^{5πi/12}, s̄₁ = 1/4 + i√3/54 (verified on 56 residue classes to 6·10⁻²¹, independently rechecked). This is quantum modularity of the Kashaev-invariant / Garoufalidis–Zagier type, with 2·Cl₂(π/3), the volume of the figure-eight knot complement, as the growth rate. The full matrix (a third, decaying sector) is still open.

6. **Sign twists: a pattern beyond ψ\*.** Inserting (−1)ⁿ into the other catalog mock theta functions gives more companions that are bounded at all 19 roots and are not a known function at ±q: the 5th-order F₀, F₁ and the 8th-order T₀, and the family F(z) at z = i, −ω. The twist of the 3rd-order ψ is also bounded, but it is a known function: −σ\*(q)/2, the Andrews–Dyson–Hickerson companion of σ(q) (ℚ(√6), Maass forms, quantum modular). Three more twists turn out to be false theta functions. None of the dense bounded twists becomes sparse after multiplying by an eta quotient, so none of them is a single Hecke-type double sum over such a product (`census3/hecke/`). F₀\*, F₁\* do follow the ψ\* pattern, with a different ingredient (proved, `docs/addendum-v3.md` §4): (III) F₀\* + Σₙ≥₁ q^{n²}(q;q²)ₙ = θ₄(q²)(−q;q)∞·X(q²) and (IV) qF₁\* − Σₙ≥₀ q^{(n+1)²}(q;q²)ₙ = θ₄(q²)(−q;q)∞·Y(q²), with X = Σ(−1)ⁿq^{n²}/(q;q)₂ₙ, Y = Σ(−1)ⁿq^{(n+1)²}/(q;q)₂ₙ₊₁. Details, the full table and first checks: [`docs/twists.md`](docs/twists.md).
7. **The double sum d1151, and an identity for it.** The rank-2 census left one candidate,

     D(q) = Σ_{m,n≥0} q^{m²−mn+n²} / ((−q;q)ₘ (−q;q)ₙ).

   - **Observed:** the growth fingerprint of Ramanujan's third-order f(q) (κ = 1/24 at all 12 even-order roots), and at q = 1 an exact expansion whose large-order behaviour carries golden-ratio constants (π²/60, π·K² = 8(5+√5)/25, s₁ = (35+6√5)/300; numerical).
   - **Proved:** D = θ₄(q)·f(q) + 4R(q) with R = Σ_{d>k≥0} (−1)^{d+k+1} q^{d²}/(−q;q)ₖ. Each diagonal m − n = d of D is ±f plus a finite sum of 1/(−q)ₖ (a three-line telescoping lemma), and the signs are those of θ₄. The same argument gives an identity for the whole family Σ q^{m²−mn+n²}/((zq)ₘ(q/z)ₙ) in terms of the rank generating function, a theta series and tail sums. So D is not a new mock theta function in the strict sense; whether the general-z identity is known is open.

   Details: [`docs/d1151.md`](docs/d1151.md) (how it was found), [`docs/d1151-identity.md`](docs/d1151-identity.md) (the identity and its proof). Census results: `census2/results/`.

## Open questions

- Are ψ₀\*, ψ₁\*, the family F(z) = Σ zⁿq^{n(n+1)/2}(−q;q)ₙ (which satisfies F(z) = 1 + zqF(zq) + zq²F(zq²)), or the identities in (3) known?
- Is the right structure a matrix-valued holomorphic quantum modular form built from (1, ψ₀\*, ψ₁\*) and θ₄R₀, θ₄R₁?
- Other values of the twist z (`census3/`, all z of order ≤ 12): at a root of odd order m, F(z) grows when z^m is closer than a₀ ≈ 0.10737 of a turn to 1 and is bounded otherwise; at even order it is always bounded. So F(z) is bounded at every root only for z = −1, ±i, −ω, −ω̄. The threshold a₀ is where the real part of a dilogarithm saddle-point value V changes sign; the predicted κ matches the fitted values to 4–5 digits at q = 1 (heuristic, `census3/README.md`).
- What are the non-modular pieces of F₀\*, F₁\* (the analogue of θ₄·R₀, θ₄·R₁ for ψ\*)?
- Is the general-z identity for Σ q^{m²−mn+n²}/((zq)ₘ(q/z)ₙ) (`docs/d1151-identity.md`, section 3) known?

## Layout

| Folder | Contents |
|---|---|
| `census/` | Evaluator and labeled catalog (`catalog.py`), per-root probe (`probe.py`), generator and dedupe (`search.py`), staged search (`stage.py`, `filt.py`, `dclass.py`), q→1 dilogarithm pre-filter (`dilog.py`), identification tools (`ident.py`, `shadow.py`, `companions.py`, …), plus all result data (`cands.json`, `labels.json`, `stageC_summary.json`, `*.jsonl`) |
| `census2/` | Rank-2 (double-sum) census: generator (`gen2.py`), evaluators (`sum2_fast.py`, `sum2.py`), staged runner (`run2.py`), report (`report2.py`); results of the first run in `census2/results/`, the analysis of d1151 in `census2/d1151/`; sign-twist scripts are in `census/twist.py` and `census/twist_checks/` |
| `psi_star/` | Quantum-modularity tests with Kontsevich–Zagier control (`qmf*.py`); exact identities, expansion at q = −1 and large-order analysis (`s1`–`s11`) |
| `census3/` | One-parameter families F(z), g₂(x,q), g₃(x,q) at all roots of unity of order ≤ 12 (`zfam.py`, `zreport.py`, `zextra.py`, results and `README.md`); `hecke/`: sparsity test of the bounded twists against eta quotients and an f_{a,b,c} matcher |
| `docs/` | `d1151.md` (the double sum d1151, how it was found), `d1151-identity.md` (its identity and proof), `d1151-brief.md` (the brief that led to it); `twists.md` (sign twists of the catalog functions); GitHub Pages site (`index.html`, `explorer-v4.html`, `fingerprints-v3.html`); `formulas.md` (every formula in LaTeX, Unicode, plain text and Mathematica, with equivalent forms and coefficients), `handoff-v1.md` (full account of method and results), `addendum-v2.md` (identities, proofs, resurgence analysis), `addendum-v3.md` (odd-root rule, identities (III), (IV), resurgence at q = i and q → 1, the law with W), `psi-star-brief.md`, `psi-star-brief-v2.md` (the briefs that led to them) |
| `psi_star/formulas_check.py` | Checks every equivalent form in `docs/formulas.md` to 50 digits |
| `psi_star/s12_odd_limits.py` | Radial limits of ψ₁\*, ψ₀\* at odd-order roots vs. the closed form (first period)/3 |
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
```

The full census takes a few hours on 2 cores. `stage.py` is resumable, and `census/labels.json` holds the final label for every candidate, so nothing needs re-running to inspect results. Scripts expect to run from their own folder.

## How this was made

This was an exploratory project by Mario Klingemann, an artist and creative technologist, not a professional mathematician. The computations and much of the code were developed in conversation with AI assistants (Anthropic's Claude Opus 5.5 and Claude Fable 5.2). Every claim in the documents is backed by a script in this repository. Corrections are very welcome.

## How to cite

Mario Klingemann, *mock-theta-census: a numerical census of Ramanujan-style q-series and the sign-twisted fifth-order functions ψ₀\*, ψ₁\**, Zenodo, 2026. https://doi.org/10.5281/zenodo.24443194 (all versions; each release also has its own version DOI on Zenodo)

GitHub's "Cite this repository" button (from `CITATION.cff`) gives the same reference in BibTeX and APA.

## License

Code: MIT (see `LICENSE`). Data and documents: CC BY 4.0.
