# mock-theta-census

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24443194.svg)](https://doi.org/10.5281/zenodo.24443194)

A complete numerical census of short Ramanujan-style q-series, classified as **modular**, **mock theta**, **false/quantum type** or **neither** by their radial asymptotics at roots of unity. Along the way it found one pair that does not appear in the literature we checked:

  **ψ₁\*(q) = Σₙ≥0 (−1)ⁿ q^{n(n+1)/2} (−q;q)ₙ**, **ψ₀\*(q) = Σₙ≥0 (−1)ⁿ q^{(n+1)(n+2)/2} (−q;q)ₙ**

These are Ramanujan's fifth-order mock theta functions ψ₁, ψ₀ with a sign twist (−1)ⁿ. They are bounded at every root of unity tested and terminate at even-order roots. At q = −1 they have a divergent asymptotic series with integer coefficients, whose large-order behavior is controlled by Ramanujan's entire Rogers–Ramanujan function at a = 1. As far as we could find, neither the functions nor the identities below appear in the literature. **We'd be glad to hear if they do.**

Two later findings extend this. The same sign twist gives bounded, non-trivial companions for several other mock theta functions (fifth-order F₀, F₁, third-order ψ, eighth-order T₀). And a census of 6,327 double sums found one mock theta candidate that we have not found in the literature so far, **D(q) = Σ_{m,n≥0} q^{m²−mn+n²} / ((−q;q)ₘ (−q;q)ₙ)**. Where it grows it behaves exactly like Ramanujan's third-order f(q); at q = 1 its asymptotic series carries the same golden-ratio constants as ψ\* at q = −1.

**Interactive pages:** [quasimondo.github.io/mock-theta-census](https://quasimondo.github.io/mock-theta-census/)
- [Explorer](https://quasimondo.github.io/mock-theta-census/explorer-v3.html): mock theta functions on the unit disk, with a radial probe toward any root of unity, including ψ₁\* and ψ₀\*
- [Fingerprints](https://quasimondo.github.io/mock-theta-census/fingerprints-v1.html): the classifier's results for 51 known functions, the census, and what ψ₁\*, ψ₀\* turned out to be
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

   with f₀\* = Σ(−1)ⁿq^{n²}/(−q;q)ₙ, f₁\* = Σ(−1)ⁿq^{n²+n}/(−q;q)ₙ, θ₄ = Σ(−1)ⁿq^{n²}, R₀ = Σ(−1)ⁿq^{n²}/(q;q)ₙ, R₁ = Σ(−1)ⁿq^{n²+n}/(q;q)ₙ. These are the z = −1 analogues of the classical bilateral-series identities for f₀ + 2ψ₀ and f₁ + 2ψ₁.
4. **Asymptotics at q = −1** (numerical, 15–18 digits).
   - ψ₁\*(−e^{−t}) = Σ cₖ tᵏ/k! with integer cₖ = 1, −1, 27, −1237, 106587, −14667181, …
   - Borel singularity at π²/60, with Stokes data π·K² = (5 ± 3√5)/40 and subleading coefficients such as (13 − 6√5)/60.
   - These match the exponential asymptotics of R₀, R₁ at q → −1 exactly.
   - Radial limits at every root of unity are finite. At even order m the series terminates; at odd order m the terms repeat with a factor −2 every m steps, and the limit is one third of the first m terms (checked at all odd orders m ≤ 15 to 10⁻¹⁴; not proved).
5. **Transformation behavior.** Under x ↦ x/(2x+1) (Γ₀(2)), ψ₁\* and ψ₀\* transform *together*, to ~10⁻⁹ on held-out roots. Unrelated pairs reach only 10⁻⁴ to 10⁻¹, and the Kontsevich–Zagier control reaches 10⁻¹⁰. But the best-fit weight is ≈ 0.508, not ½, so this is not a classical weight-½ quantum modular form with constant multipliers. Our reading: this fits the resurgent / matrix-valued quantum modularity framework of Garoufalidis–Zagier, as a sign-twisted sibling of the analysis of f₀ near q = −1 in C. Wheeler's thesis (Bonn 2023, Ch. 7).

6. **Sign twists: a pattern beyond ψ\*.** Inserting (−1)ⁿ into the other catalog mock theta functions gives more companions that are bounded at all 19 roots and are not a known function at ±q: the 5th-order F₀, F₁, the 3rd-order ψ and the 8th-order T₀, and the family F(z) at z = i, −ω. Three more twists turn out to be false theta functions. F₀\*, F₁\* do not follow the ψ\* pattern of identities (I), (II) in any simple way. Details, the full table and first checks: [`docs/twists.md`](docs/twists.md).
7. **A double-sum candidate, d1151.** The rank-2 (double-sum) census left one mock theta candidate,

     D(q) = Σ_{m,n≥0} q^{m²−mn+n²} / ((−q;q)ₘ (−q;q)ₙ).

   - **Growth:** where it grows it has exactly the fingerprint of Ramanujan's third-order f(q): κ = 1/24 at all 12 even-order roots, bounded at the odd-order ones.
   - **Not a known function:** it is not a rational linear combination of the third-order functions, the sign twists and standard theta functions (checked to q⁹⁰).
   - **Asymptotics at q = 1:** its exact expansion has integer coefficients (times k!). Their growth carries the golden-ratio constants seen for ψ\* at q = −1: Borel singularity at π²/60, π·K² = 8(5+√5)/25 and s₁ = (35+6√5)/300, all numerical.

   Details and proof status: [`docs/d1151.md`](docs/d1151.md). Census results: `census2/results/`.

## Open questions

- Are ψ₀\*, ψ₁\*, the family F(z) = Σ zⁿq^{n(n+1)/2}(−q;q)ₙ (which satisfies F(z) = 1 + zqF(zq) + zq²F(zq²)), or the identities in (3) known?
- Is the right structure a matrix-valued holomorphic quantum modular form built from (1, ψ₀\*, ψ₁\*) and θ₄R₀, θ₄R₁?
- What happens at other values of the twist z? F(i) and F(−ω) are bounded at all 19 roots (`docs/twists.md`), F(ω) is not.
- What are the non-modular pieces of F₀\*, F₁\* (the analogue of θ₄·R₀, θ₄·R₁ for ψ\*)?
- Is D (d1151) known? Is there an identity expressing it through third-order functions, theta functions and Rogers–Ramanujan-type series?

## Layout

| Folder | Contents |
|---|---|
| `census/` | Evaluator and labeled catalog (`catalog.py`), per-root probe (`probe.py`), generator and dedupe (`search.py`), staged search (`stage.py`, `filt.py`, `dclass.py`), q→1 dilogarithm pre-filter (`dilog.py`), identification tools (`ident.py`, `shadow.py`, `companions.py`, …), plus all result data (`cands.json`, `labels.json`, `stageC_summary.json`, `*.jsonl`) |
| `census2/` | Rank-2 (double-sum) census: generator (`gen2.py`), evaluators (`sum2_fast.py`, `sum2.py`), staged runner (`run2.py`), report (`report2.py`); results of the first run in `census2/results/`, the analysis of d1151 in `census2/d1151/`; sign-twist scripts are in `census/twist.py` and `census/twist_checks/` |
| `psi_star/` | Quantum-modularity tests with Kontsevich–Zagier control (`qmf*.py`); exact identities, expansion at q = −1 and large-order analysis (`s1`–`s11`) |
| `docs/` | `d1151.md` (the double-sum candidate); `twists.md` (sign twists of the catalog functions); GitHub Pages site (`index.html`, `explorer-v3.html`, `fingerprints-v1.html`); `formulas.md` (every formula in LaTeX, Unicode, plain text and Mathematica, with equivalent forms and coefficients), `handoff-v1.md` (full account of method and results), `addendum-v2.md` (identities, proofs, resurgence analysis), `psi-star-brief.md` (focused statement of the open problem) |
| `psi_star/formulas_check.py` | Checks every equivalent form in `docs/formulas.md` to 50 digits |
| `psi_star/s12_odd_limits.py` | Radial limits of ψ₁\*, ψ₀\* at odd-order roots vs. the closed form (first period)/3 |

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
