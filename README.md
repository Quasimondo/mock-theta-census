# mock-theta-census

A complete numerical census of short Ramanujan-style q-series, classified as **modular**, **mock theta**, **false/quantum type** or **neither** by their radial asymptotics at roots of unity. Along the way it found one unexplained pair:

  **ψ₁\*(q) = Σₙ≥0 (−1)ⁿ q^{n(n+1)/2} (−q;q)ₙ**, **ψ₀\*(q) = Σₙ≥0 (−1)ⁿ q^{(n+1)(n+2)/2} (−q;q)ₙ**

These are Ramanujan's fifth-order mock theta functions ψ₁, ψ₀ with a sign twist (−1)ⁿ. They are bounded at every root of unity tested and terminate at even-order roots. At q = −1 they have a divergent asymptotic series with integer coefficients, whose large-order behavior is controlled by Ramanujan's entire Rogers–Ramanujan function at a = 1. As far as we could find, neither the functions nor the identities below appear in the literature. **We'd be glad to hear if they do.**

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
5. **Transformation behavior.** Under x ↦ x/(2x+1) (Γ₀(2)), ψ₁\* and ψ₀\* transform *together*, to ~10⁻⁹ on held-out roots. Unrelated pairs reach only 10⁻⁴ to 10⁻¹, and the Kontsevich–Zagier control reaches 10⁻¹⁰. But the best-fit weight is ≈ 0.508, not ½, so this is not a classical weight-½ quantum modular form with constant multipliers. Our reading: this fits the resurgent / matrix-valued quantum modularity framework of Garoufalidis–Zagier, as a sign-twisted sibling of the analysis of f₀ near q = −1 in C. Wheeler's thesis (Bonn 2023, Ch. 7).

## Open questions

- Are ψ₀\*, ψ₁\*, the family F(z) = Σ zⁿq^{n(n+1)/2}(−q;q)ₙ (which satisfies F(z) = 1 + zqF(zq) + zq²F(zq²)), or the identities in (3) known?
- Is the right structure a matrix-valued holomorphic quantum modular form built from (1, ψ₀\*, ψ₁\*) and θ₄R₀, θ₄R₁?
- What happens at other values of the twist z, e.g. roots of unity i or ω?

## Layout

| Folder | Contents |
|---|---|
| `census/` | Evaluator and labeled catalog (`catalog.py`), per-root probe (`probe.py`), generator and dedupe (`search.py`), staged search (`stage.py`, `filt.py`, `dclass.py`), q→1 dilogarithm pre-filter (`dilog.py`), identification tools (`ident.py`, `shadow.py`, `companions.py`, …), plus all result data (`cands.json`, `labels.json`, `stageC_summary.json`, `*.jsonl`) |
| `psi_star/` | Quantum-modularity tests with Kontsevich–Zagier control (`qmf*.py`); exact identities, expansion at q = −1 and large-order analysis (`s1`–`s11`) |
| `docs/` | `formulas.md` (every formula in LaTeX, Unicode, plain text and Mathematica, with equivalent forms and coefficients), `handoff-v1.md` (full account of method and results), `addendum-v2.md` (identities, proofs, resurgence analysis), `psi-star-brief.md` (focused statement of the open problem) |
| `psi_star/formulas_check.py` | Checks every equivalent form in `docs/formulas.md` to 50 digits |

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

## License

Code: MIT (see `LICENSE`). Data and documents: CC BY 4.0.
