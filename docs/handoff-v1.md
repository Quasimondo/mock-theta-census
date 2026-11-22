# Mock theta search: handoff summary

Handoff from a Claude Opus 5.5 session (20–21 Nov 2026) working with Mario Klingemann.
Everything below was computed in that session. The code and all raw results are in `mock-theta-search-v4.zip`; file names in this document refer to that archive. Numbers are quoted from the result files. Where something is a conjecture, an inference or unverified, it says so.

---

## 0. TL;DR for the next agent

1. **Built and validated a numerical classifier** that sorts q-series into **modular / mock / false-or-quantum / neither** from their radial asymptotics at 19 roots of unity, using 50+ digit ball arithmetic (python-flint/Arb). On 50 labeled known functions it classified every one correctly, except one predictable miss (a sum of two modular forms with different q-offsets).
2. **Ran a complete census** of one fixed "Ramanujan-like" template: 9,936 Eulerian formulas, 8,869 distinct functions. Result: **25 mock, 28 modular, 1,357 bounded at all 19 roots (quantum/false type)**; all others rejected (irrational growth exponent or non-theta growth).
3. **All 25 mock hits are accounted for.** 23 are known functions, rescalings, simple combinations, or specializations of universal mock theta functions. That includes a new-looking but probably elementary identity: **H(x) = (1−x²)·g₂(x,q) + 1 − x**. The other 2 decompose as **½(fifth-order ψ) + ½(sign-alternating companion ψ*)**.
4. **The open object is the pair ψ₁*, ψ₀*:**
   ψ₁*(q) = Σₙ≥0 (−1)ⁿ q^{n(n+1)/2} (−q;q)ₙ,  ψ₀*(q) = Σₙ≥0 (−1)ⁿ q^{(n+1)(n+2)/2} (−q;q)ₙ.
   They are bounded at all tested roots. They are not Hecke-type (multiplying by simple products never gives a sparse series). At even-order roots they terminate to algebraic finite sums, the Kontsevich–Zagier "strange function" pattern. Under x ↦ x/(2x+1) they **transform together, to roughness ~10⁻⁹**, but with **unclean parameters (weight ≈ 0.508)**. Best guess: an incomplete vector-valued quantum modular form, or a higher-depth (logarithmic) one. **This is the most promising thing to resolve.**
5. **Recommended next steps** (§9): resolve the (ψ₁*, ψ₀*) structure (find a third component, or test for log/higher-depth terms); check the literature for ψ₁*, ψ₀* and for the H(x) identity; ask a specialist (Bringmann, Rolen, Mortenson, Zwegers, Lovejoy). Mario is considering an open-data release and possibly an experimental-mathematics paper with a specialist coauthor.

---

## 1. Background and goal

- Mario (artist and creative technologist, long-running interest in "interestingness") asked whether Ramanujan's legacy contains anything worth revisiting computationally. The conversation narrowed to mock theta functions and to the question of which of the infinitely many mock theta functions are "Ramanujan-like" (short Eulerian formula).
- Plan: (1) build a classifier and validate it on known functions; (2) enumerate a natural template of short Eulerian series; (3) classify everything; (4) identify the hits; (5) investigate whatever stays unexplained.
- An interactive explorer (f(q) on the disk plus a radial-limit test) and a results page were published as Claude artifacts for Mario. The agent doesn't need them.

## 2. Conventions

- Approach a root ζ = e^{2πi h/m} (reduced, h/m ∈ [0, ½]) along the ray q = e^{−t}ζ, t → 0⁺. Since τ = h/m + it/2π exactly, modular transformations give exact laws in t.
- **Fit model** (per root): log|g(e^{−t}ζ)| = c/t + α log t + β + γ t (+ exponentially small terms).
  - κ := c·m²/(4π²) is the cusp exponent; for modular and mock forms it is rational.
  - −α is the weight (a half-integer for modular and mock forms).
  - γ is the q-power offset: g = q^{−γ}·(modular object). It recovers known normalizations, e.g. γ = −1/24 for f.
- **Per-root classes:**
  - **E** (exact): fit residual < 10⁻²⁵ and κ within 10⁻²⁰ of a rational with denominator ≤ 5000.
  - **Q** (bounded with a leftover): |κ| < 10⁻⁴ but the residual isn't exponentially small. A genuine power series in t remains: the quantum or "strange" signature.
  - **X**: anything else (irrational κ, or no exact law).
- **Function-level classes:** all E → modular. E and Q mixed, no X → mock. All Q → false/quantum. Any X → neither.
  - Caveat: this implements **Ramanujan's definition** (bounded difference from a theta function at each root, no single theta function working everywhere), **not** the modern one. A mock function plus anything bounded everywhere still classifies as mock. See §7.
- Roots used (19): 0/1, 1/2, 1/3, 1/4, 1/5, 2/5, 1/6, 1/7, 2/7, 3/7, 1/8, 3/8, 1/10, 3/10, 1/12, 5/12, 1/14, 3/14, 5/14. Fingerprint strings follow exactly this order.

## 3. Numerical method (`catalog.py`, `probe.py`)

- **Evaluation:** Arb complex balls (python-flint). Eulerian sums Σ q^{e(n)} Pₙ are formed term by term as exp(e(n)·log q + log Pₙ), with log Pₙ accumulated as a **sum of logs**.
  - Repeated complex-ball multiplication suffers the wrapping effect (radii grow exponentially, ×1.2–1.3 per step). Accumulating logs avoids it.
  - `slog` handles the negative real axis (any consistent branch works, since everything is exponentiated). Without it, ζ = −1 failed with spurious "precision" errors; that's fixed.
- **Adaptive precision:** start at 256 bits and double until the relative radius of |g| < 2⁻¹⁵⁰. The cap is 4,096 bits in the search stages and 16,384–65,536 bits in the full tests.
- **t-grid:** 8 points, t = u/m² with u geometric in [0.004, 0.02].
  - **Second pass A:** if κ > 10⁻⁴ and the fit isn't exact, move closer: t_max = c/150 and t_min = t_max/5, so that bounded parts are suppressed by e^{−150}. Without this, mock functions with small κ were misclassified.
  - **Second pass B** (full mode only): if |κ| ≤ 10⁻⁴ and the fit isn't exact, refit with t/8. An exact fit there means modular-bounded; otherwise Q.
- **"FAST" mode** (env `FAST=1`, used in search filtering): skips pass B. A Q in FAST mode means "bounded at this root", without deciding exact-bounded versus quantum.
- **Least squares** in mpmath at 60 digits.
- **Coefficients** (`search.py`): exact integer truncated power series (numpy convolution). Catalog functions use a numerical Cauchy integral on |q| = ½ (128 points), rounded to integers.

## 4. Round 1: validation on known functions (`results2.json`, `summary.json`)

- **Catalog** (`catalog.py`, 51 functions):
  - 15 modular: θ₃, θ₄, η, b(q) = θ₄/(−q;q)∞, the ternary theta x²+y²+10z², 1/(q;q)∞ via the Durfee sum, Rogers–Ramanujan G and H, Zagier's modular Nahm sums, Euler (−q;q)∞, Göllnitz–Gordon, and (G−H)/q labeled "mixed".
  - 27 mock: third order f, φ, ψ, χ, ω, ν, ρ; fifth f₀, f₁, F₀, F₁, φ₀, φ₁, ψ₀, ψ₁; seventh F₀, F₁, F₂; second μ; sixth φ, ψ; eighth S₀, T₀; tenth φ, ψ, X, χ.
  - 2 false theta functions, σ (Andrews–Dyson–Hickerson), and 6 non-modular rank-1 Nahm sums (Zagier's classification).
  - Third- and fifth-order formulas are cross-checked by identities to 58 digits: 2φ(−q) − f(q) = f(q) + 4ψ(−q) = b(q); χ₀ = 2F₀ − φ₀(−q); χ₁ = 2F₁ + φ₁(−q)/q. Formulas for orders 2, 6, 7, 8 and 10 came from memory; several were later confirmed in the literature (§6), the rest are not independently checked.
- **Result:**

  | Label | Classified correctly |
  |---|---|
  | Modular | 14/14 |
  | Mock | 27/27 |
  | Non-modular Nahm sums | 6/6 |
  | False theta | 2/2 |

  - σ → "false/quantum", consistent with Zagier's description of it as a quantum modular form.
  - (G−H)/q → "neither". Expected: a sum of two modular forms with different offsets fits no single law.
- **Recovered without being told:** offsets q^{−1/24} for f, q^{−1/60} for f₀ and G, q^{−1/168} for the seventh-order F₀, q^{11/60} for H; every mock function has weight ½. For the Nahm sums outside Zagier's list, measured κ at q→1 matches the Rogers dilogarithm prediction:

  | Sum | κ measured | κ predicted |
  |---|---|---|
  | A=4 | 0.0127913286199 | 0.0127913286192 |
  | A=6 | 0.0107523461996 | 0.0107523461981 |
  | A=3 (base q²) | 0.0071739622451 | 0.0071739622449 |
  | A=3/2 (base q⁴) | 0.0045941000030 | 0.0045941000030 |

- **Observation:** the set of bounded roots distinguishes families. Third-order f is bounded exactly at roots of odd order. Each seventh-order function is bounded at a different subset of the 7th and 14th roots: F₀ at 2/7, 3/7, 3/14, 5/14; F₁ at 1/7, 2/7, 1/14, 5/14; F₂ at 1/7, 3/7, 1/14, 3/14.

## 5. Round 2: census of one template (`search.py`, `stage.py`, `filt.py`, `dclass.py`)

**Template** (complete enumeration, not random sampling): Σₙ≥0 s^n q^{(An²+Bn)/2} ∏_{k=1}^{n} F₁(k)^{p₁} [F₂(k)^{p₂}] with
- A ∈ {1, 2, 3, 4}, B ∈ [−A, A+2] with A+B even (18 pairs), s = ±1;
- 8 factor types: (1±qᵏ), (1±q²ᵏ), (1±q^{2k−1}), (1−qᵏ+q²ᵏ), (1+qᵏ+q²ᵏ);
- powers p ∈ {−2, −1, 1}; one factor or two distinct ones (276 product patterns).

Total 9,936 formulas, deduplicated to **8,869 distinct functions**. Canonical form: coefficients up to an additive constant, a monomial factor qᵛ, overall sign, and q → −q. q → q² rescalings were not deduplicated at generation; they were caught afterwards. The generator rediscovered 31 catalog functions, e.g. f, φ, ψ, χ, ω, ν, f₀, f₁, F₀, F₁, φ₀, ψ₁, seventh-order F₀, S₀, tenth-order φ and ψ, G, H, σ.

**Funnel** (FAST mode at the early roots; full 19-root test for the finalists):
- q→1: 3,587 out. q→−1, e^{2πi/3}, i, e^{iπ/3}: further rejections.
- 53 finalists → **25 mock, 28 modular** (all confirmed at 19 roots).
- 1,387 bounded at the first five roots → tested at the remaining 14: **1,357 bounded at all 19**, 1 irrational at an order-5 root, 29 with non-theta growth at order-12 roots (or order 8/10). These 29 were rerun at high precision: the fit stays poor and κ differs between 1/12 and 5/12, so they're genuinely out. 21 of 26 order-12 rejects have (1−qᵏ+q²ᵏ) in the denominator.
- **Census statement:** within this template, at these 19 roots, there are exactly 25 mock theta functions and 28 modular forms.
  - Caveat: "bounded at 19 roots" doesn't mean "bounded at all roots"; first blow-ups at orders 9, 11 or 16+ aren't excluded.
- The 28 modular finalists look like known identities (Slater-type; κ denominators 168/336 suggest the Rogers–Selberg mod-7 family). **Not checked against Slater's list.**

**Pre-filter discovered afterwards** (`dilog.py`, `prefilter_eval.py`). For sign-free series, the q→1 exponent follows from the formula alone:
Φ(x) = −A x²/2 + Σ_f (p_f/d_f)·[Li₂(−σ_f e^{−d_f x}) − Li₂(−σ_f)], c = max_{x≥0} Φ(x), κ = c/4π².
Cyclotomic factors are decomposed: 1−qᵏ+q²ᵏ = (1+q³ᵏ)/(1+qᵏ), 1+qᵏ+q²ᵏ = (1−q³ᵏ)/(1−qᵏ).
- It agrees with the measured κ to 9–10 digits (55/60 in a spot check; the rest are bounded cases).
- On all 4,605 sign-free candidates it would have removed 1,886 of the 2,237 q→1 rejects with **zero false exclusions**. The only two "kept" functions it rejected were the catalog's own non-modular Nahm controls, correctly.
- Not yet extended to ζ = −1 or to (−1)ⁿ series.

**Survival rates** (priors for designing the next template, not rules):

| Feature | Survival |
|---|---|
| No (−1)ⁿ sign | 78/4,605 (1.7%) |
| (−1)ⁿ sign | 6/4,264 (0.14%) |
| A = 2 | 39/1,866 |
| A = 3 | 5/2,628 |
| No denominator factors | 67% end up bounded everywhere |

## 6. Identification of all 25 mock hits (`stageC_summary.json`)

Formulas below are the representatives in generator notation: q^{(An²+Bn)/2} times ∏_{k≤n} factors.

| id | formula | κ (blow-up roots) | identification |
|---|---|---|---|
| c1882, c2360, c2371, c2768, c2986, c5597, c6850, c7734, c7735, c7954, c8234, c8839 | various | various | rescalings q→q² or q³ of catalog functions (φ₁₀, ψ₁, ψ₁₀, ν, f, f₀, ψ, φ₀, f₁), plus μ (2nd) and X (10th) themselves |
| c6630, c7999 | — | — | q→q² of c1656 and c2230 |
| c2230 | Σ q^{n²}(q;q)ₙ/(q³;q³)ₙ | 1/24 | sixth-order γ (Andrews–Hickerson; Lost Notebook). Hickerson–Mortenson give γ = 2m(q,q³,−1) + m(q,q³,−q) |
| c3467 | Σ q^{n(n+2)}(−q;q²)ₙ/(−q²;q²)ₙ | 1/16 | eighth-order S₁ (Gordon–McIntosh) |
| c2437 | Σ q^{n²}(−q;q²)ₙ/(q;q²)ₙ | 1/64 | (V₀+1)/2, eighth order |
| c1106 | Σ q^{n(n+3)/2}(−q;q)ₙ | 1/240 | ψ₀/q (fifth order); missed by canon only because of constant stripping |
| c3326 | Σ q^{n²+2n}(−q;q²)ₙ | 1/120 & 1/480 | φ₁/q (fifth order) |
| c1656 | Σ q^{n(n−1)}/(−q;q)ₙ² | 1/24 | −2f + modular form (power-series test, misfit 10⁻²⁷) |
| c7882 | Σ q^{2n²}/(−q;q)₂ₙ | 1/48 | at q²: (φ(q)+φ(−q))/2, the even part of third-order φ |
| c76 | Σ q^{n(n−1)/2}(−q;q)ₙ²/(−q³;q³)ₙ | 1/144 (blows up at orders 1, 5, 7) | **= H(−ω)**, see below |
| c85 | Σ q^{n(n−1)/2}(q²;q²)ₙ/(q³;q³)ₙ | 1/16 (blows up **only at e^{±2πi/3}**) | **= H(ω)**, see below. An earlier note saying "only at ±i" was a misreading, since corrected |
| c8438 | Σ q^{2n²+n}(−q;q)₂ₙ | 1/240 | **= ½ψ₁ + ½ψ₁*** (even-index half of ψ₁'s series) |
| c9542 | Σ q^{2n²+3n}(−q;q)₂ₙ | 1/240 | **= (½ψ₀ + ½ψ₀*)/q** |

**The H(x) identity** (`ident.py`). Define H(x,q) = Σₙ q^{n(n−1)/2}(−q;q)ₙ/((xq;q)ₙ(q/x;q)ₙ). Since 1∓qᵏ+q²ᵏ = (1±ωqᵏ)(1±ω²qᵏ), with ω = e^{2πi/3}, c76 = H(−ω) and c85 = H(ω). Found numerically:

  **H(x,q) = (1 − x²)·g₂(x,q) + 1 − x**,  g₂(x,q) = Σₙ (−q;q)ₙ q^{n(n+1)/2}/((x;q)_{n+1}(q/x;q)_{n+1}) (Gordon–McIntosh universal mock theta function).

- Verification: exact rational arithmetic to q¹⁹ at x ∈ {2, 3, 5, ½, −2, 7/3}, and Arb at 300 bits at three q values (including |q| = 0.9) for x = ±ω: |difference| ≲ 10⁻⁸⁶.
- **Not proven.** It looks like a telescoping identity and may be known or elementary.
- Hence c76 = (1−ω²)·g₂(−ω,q) + 1 + ω and c85 = (1−ω²)·g₂(ω,q) + 1 − ω. Even-order mock theta functions are known to be related to g₂ (Gordon–McIntosh), so these are known territory. Whether g₂(±ω,q) has been studied as a specific function is unchecked. An explicit Appell–Lerch form follows from the known one for g₂; the exact published formula for g₂ wasn't retrieved.

**The splits** (exact by construction: even-n terms = ½(Σ + Σ(−1)ⁿ)). Confirmed independently by the power-series test at ζ = −1 and i: coefficients ½ and ½ to 10⁻²⁴, misfit 10⁻²⁵. That explains why the earlier fit against ψ₀ and ψ₁ alone stalled at 3·10⁻⁵.

**The power-series ("shadow") test** (`shadow.py`). At a root where a function is bounded, fit g(e^{−t}ζ) ≈ Σₖ aₖ tᵏ + b₁t^{−1/2} + b₂t^{−1} + b₃t^{−3/2}. The negative powers absorb modular pieces with κ = 0.
- Two functions differ by a modular form if and only if their (a₁…a₅) vectors are proportional across the bounded roots.
- Works well at 0/1, 1/3, 1/5: φ(−q) = ½f + modular to 10⁻³² and c1656 = −2f + modular to 10⁻²⁷.
- Near-singular roots need very small t (T0 = 2·10⁻⁵/m²) and the t^{−k/2} terms. At 1/2 and 1/4 it was less reliable before those fixes.

## 7. The open object: ψ₁*, ψ₀*

- **Definitions:** ψ₁*(q) = Σₙ≥0 (−1)ⁿ q^{n(n+1)/2}(−q;q)ₙ (candidate c830); ψ₀*(q) = Σₙ≥0 (−1)ⁿ q^{(n+1)(n+2)/2}(−q;q)ₙ (c1382, as ψ₀*/q).
- **Coefficients:**
  - ψ₁*: 1, −1, −1, 1, 1, 1, 0, −1, −1, −2, 0, 0, 0, 2, 2, 1, 1, 1, −1, −1, −2, −2, −2, −2, −1, −1, 1, 2, 2, 2, …; dense (93% nonzero), slowly growing, max |a| = 26 up to n = 120.
  - ψ₀*/q: 1, 0, −1, −1, 0, 1, 1, 1, 1, −1, −1, −1, −2, −1, 0, 0, 1, 2, 2, 2, 1, 1, 0, −1, −1, −3, −3, −2, −2, −2, …
- **Facts established:**
  1. Both are bounded at all 19 roots (class Q everywhere).
  2. **No relatives:** across the 8,869 template functions, and across ~90,000 further formulas (12 factor types incl. (1±q³ᵏ), (1±q^{3k−1}); powers ±1, ±2; A ∈ 0..6; `companions.py`), nothing satisfies g(±q)·qᵏ = a·ψ₁* + b·ψ₀* + c·ψ₁ + d·ψ₀ + e (exact to q⁵⁹) beyond trivial rewrites and the two halves. The only other relation is the elementary Σ(−1)ⁿ q^{n(n−1)/2}(−q;q)ₙ = 1 − ψ₁* − ψ₀*.
  3. **Not σ**, not a simple Rogers–Fine-type series (Σ(±1)ⁿqⁿ(q;q)ₙ, (−q;q)ₙ, (q;q²)ₙ, (±q²;q²)ₙ and 1/(−q;q)ₙ were tried), and not a false theta function (coefficients too dense).
  4. **Not Hecke-type with a simple multiplier:** ψ₁*, ψ₀* times (q;q)∞, (q;q)∞², (q²;q²)∞, (q;q)∞/(−q;q)∞, (q;q)∞(q²;q²)∞, (q;q)∞²/(q²;q²)∞, 1/(−q;q)∞ or (q²;q²)∞/(−q;q)∞ never becomes sparse (max |c| 14–16, ~93% nonzero). Control: ψ₁·(q²;q²)∞ has max |c| = 2 and 26/60 nonzero.
  5. **Strange-function behavior:** at even-order roots (−ζ;ζ)ₙ vanishes for n ≥ m/2, so the series terminates. Radial limits converge to these finite sums with error ∝ t: ψ₁*(−1) = 1, ψ₁*(i) = 2 − i, ψ₁*(e^{iπ/3}) = 1 − 2√3 i, ψ₁*(e^{iπ/4}) ≈ −3.5355339 − 4.1213203i, ψ₁*(e^{3iπ/4}) ≈ 3.5355339 − 0.1213203i, ψ₁*(e^{iπ/5}) ≈ −8.2082039 + 1.5388418i, ψ₁*(e^{iπ/6}) ≈ −4.4641016 + 12.196152i.

- **Quantum-modularity test** (`qmf*.py`). φ(x) = e^{2πiαx}·F(e^{2πix}) at rationals, with exact finite sums. Cocycle h(x) = φ(x) − Σⱼ εⱼ (cx+d)^{−k} φⱼ(γx). Multipliers are fitted by closed-form least squares on neighbor differences. **Roughness** = Σ|Δh|²/Σ|Δφ|² over sorted neighboring rationals; smooth h makes it scale like (spacing)². Group element γ: x ↦ x/(2x+1) ∈ Γ₀(2), chosen because it preserves even denominators (x ↦ −1/x doesn't). Points: x ∈ (0.05, 2), even denominators ≤ 40–80.

  | Test | Roughness |
  |---|---|
  | Control, Kontsevich F = Σ(q;q)ₙ under S: weight 3/2, α = 1/24, ε = −i | 2·10⁻¹⁰ (wrong params ≥ 6·10⁻³; median 0.54) |
  | Same control under x/(2x+1): weight 3/2, α = 1/24, ε = e^{iπ/6} | 1.1·10⁻¹⁰; scales 8·10⁻⁹ → 1.3·10⁻¹¹ for maxden 24 → 64 |
  | ψ₁* alone, best over k ∈ {0,…,5/2}, α ∈ j/240 | 0.089 (flat in α, \|ε\| = 0.68) |
  | Pair (ψ₁*, ψ₀*), k = ½, α(ψ₁*) = 1/6, α(ψ₀*) = 3/40, separate fits for m ≡ 0 and 2 mod 4 | ~1.3·10⁻⁶, coefficients ≈ (0.43869+0.98266i, −0.13194+0.37248i), identical across classes; plateaus at ~10⁻⁶ (6·10⁻⁴, 2·10⁻⁵, 2·10⁻⁶, 1.5·10⁻⁶, 1·10⁻⁶ for maxden 24, 36, 48, 64, 80) |
  | Adding ψ₁, ψ₀ (mock) terminating values as components | no improvement; coefficients ≈ 0.01 |
  | Pair, continuous optimization (Nelder–Mead over α₁, α₀, k), target ψ₁* | 1.2·10⁻⁹ at α(ψ₁*) = 0.179493, α(ψ₀*) = −0.218409 (also found as 0.7816), k = 0.507778, coefficients ≈ (0.3883+1.0124i, −0.3417+0.2030i) |
  | Same, target ψ₀* | 0.9·10⁻⁹ at α(ψ₁*) = 0.388115, α(ψ₀*) = −0.012957, k = 0.507942 |
  | Held-out validation of the continuous fit | m ≡ 2 mod 4: 1.2·10⁻⁹; m ∈ 66..90: 1.1·10⁻⁸ (ψ₀* target: 1.9·10⁻⁹ and 8.7·10⁻⁹) |
  | Negative controls, same optimizer: (Σ(−1)ⁿq^{n²}(−q;q)ₙ, Σ(−1)ⁿq^{n²+n}(−q;q)ₙ) | 10^{−4.26} |
  | (Σq^{n(n+1)}(−q;q)ₙ, Σ(−1)ⁿq^{(3n²+n)/2}(−q;q)ₙ) | 10^{−1.03} |

- **Interpretation (inference):**
  - (ψ₁*, ψ₀*) come close to a vector-valued quantum modular law of weight ≈ ½ under Γ₀(2). The law is real: it validates out of sample and controls are far worse.
  - But the parameters aren't clean, and the ψ₁* normalizations from the two target fits are inconsistent (0.179 vs 0.388), so the model is misspecified. Candidates:
    - (a) a missing third component;
    - (b) a cocycle with logarithmic terms, i.e. higher-depth quantum modular forms (Bringmann–Kaszian–Milas). A log term would mimic a shifted weight;
    - (c) the wrong group or element, e.g. Γ(2), or using the T-relation jointly;
    - (d) a mixed object.
  - Only one Γ₀(2) generator was tested. Consistency with T: x ↦ x+1 wasn't imposed jointly.
- **Consequence:** c8438 and c9542 are mock in Ramanujan's sense automatically, since ψ* is bounded. Whether they are mock in the modern sense depends on what ψ₀* and ψ₁* are; probably not, since ψ* is not mock-like.

## 8. Things that went wrong, now fixed (avoid repeating)

- Python-flint string constructors (`acb("0.31")`) produce low-precision balls; construct from exact rationals.
- `theta4` sign was off by one with n0 = 1; fixed with init = −1.
- The wrapping effect in repeated complex-ball products: use log accumulation.
- The branch cut at ζ = −1: use `slog`.
- Mock functions with small κ looked "irrational" until the adaptive second pass moved closer to the rim.
- Pure-q-offset terms (e.g. q^{−1/60} for G) need the γ t term in the fit model.
- Canon dedupe strips constants, so ψ₀/q-type shifts with a constant term were missed (c1106, c3326). q→q² rescalings are not deduplicated.
- The OEIS search endpoint blocks automated access (WebFetch fails); Mario searched manually. A truncated search for c76 matched A218085 by coincidence (base-4 digit sums; not related).

## 9. Suggested next steps (priority order)

1. **Pin down (ψ₁*, ψ₀*):**
   - test a Γ₀(2) generating set jointly, with T-consistency imposed (α's tied to the T-action);
   - add a log-cocycle term (h = φ − Σε(cx+d)^{−1/2}φ(γx) − Σδ(cx+d)^{−1/2}log(cx+d)·φ(γx)) and check whether the weight snaps to ½;
   - search for a third component among the 1,357 bounded functions with a (−q;q)ₙ factor, which terminate at even roots;
   - compute the strange power series at t → 0 and compare with the Taylor expansion of the fitted cocycle.
2. **Literature check:** ψ₁*, ψ₀* (Andrews–Berndt *Lost Notebook* vols. on partial theta functions; Bringmann–Folsom–Rhoades on quantum modular forms from q-series; Hikami on quantum invariants; Kontsevich–Zagier strange identity). Also the identity H(x) = (1−x²)g₂(x,q) + 1 − x, and whether g₂(±ω,q) is named.
3. **Prove H(x)** (likely telescoping), and derive explicit Appell–Lerch forms of c76 and c85 from the known g₂ expansion (Hickerson–Mortenson, arXiv:1208.1421).
4. **Engineering:**
   - set Pool size to os.cpu_count() (fixed at 2 now);
   - wire `dilog.py` into generation;
   - add q→q² to the canonical form;
   - fold the power-series test into the classifier, so "mock + bounded-everywhere" mixtures are flagged;
   - a widened template (adding (1±q³ᵏ) and A ≤ 6: ~23k–40k formulas) needs ~6–12 h on 2 cores, or ~1 h on 16.
5. **Possible publication** (Mario's call): open code and data (GitHub + Zenodo) first. Then, with a specialist coauthor, an experimental-mathematics note (method + census + whatever ψ* turns out to be); arXiv math.NT → *Experimental Mathematics* / *Research in Number Theory*. AI use must be disclosed; AI can't be an author.

## 10. File inventory (`mock-theta-search-v4.zip`)

- **Core:**
  - `catalog.py`: evaluator and labeled catalog; `python3 catalog.py` runs the identity checks.
  - `probe.py`: per-root fit and classification. Env: FAST, PCAP.
  - `analyze.py`: round-1 confusion table.
  - `fmt.py`: pretty-printer for probe lines.
- **Search:**
  - `search.py`: generator, exact coefficients, canonical dedupe → `cands.json` (8,869 representatives with id, formula, number of forms, catalog match, first coefficients).
  - `stage.py`: resumable parallel probing.
  - `filt.py`, `dclass.py`: survivor filters.
  - Results: `stageA/B1/B2/B3/C/D1–D4.jsonl`; `search_filter_results.jsonl` (A+B concatenated); `quantum_followup.jsonl` (D concatenated).
  - Summaries: `stageC_summary.json` (53 finalists with fingerprints), `labels.json` (final label per candidate), `groupsB3.json`.
- **Analysis:**
  - `dilog.py`, `prefilter_eval.py`, `dilog_pred.json`: q→1 pre-filter.
  - `patterns.py`: labels and feature tables.
  - `shadow.py`, `shtest.py`, `span.py`, `span2.py`: power-series test.
  - `ident.py`: H, g₂, g₃ at ±ω.
  - `relatives.py`, `relatives.json`, `companions.py`, `companions.json`: companion searches.
  - `qmf.py` (core), `qmf_scan.py`, `qmf_vec.py`, `qmf_vec4.py`, `qmf_split.py`, `qmf_scale.py`, `qmf_opt.py`, `qmf_val.py`, `qmf_negctl.py`: quantum-modularity tests.
- **Requirements:** python ≥ 3.10, `pip install python-flint mpmath numpy scipy`. All scripts expect to run from the archive directory.
