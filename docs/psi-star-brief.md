# Brief: the (ψ₁*, ψ₀*) problem: find the missing component

**For:** an agent continuing a computational investigation with Mario Klingemann.
**Status:** numerical evidence only. Everything here was computed in a previous session; scripts are in `psi-star-tests-v1.zip`. Background on the wider project is in `mock-theta-handoff-v1.md` (optional).

## 1. Where this comes from (one paragraph)

A census of ~8,900 short Eulerian q-series found 25 mock theta functions (in Ramanujan's sense: near each root of unity the function differs from a theta function by something bounded). All but two are known. The two leftovers split exactly as

- Σ q^{2n²+n}(−q;q)₂ₙ = ½ψ₁(q) + ½ψ₁*(q)
- Σ q^{2n²+3n}(−q;q)₂ₙ = (½ψ₀(q) + ½ψ₀*(q))/q

where ψ₀, ψ₁ are Ramanujan's fifth-order mock theta functions and the **companions** are

  **ψ₁*(q) = Σₙ≥0 (−1)ⁿ q^{n(n+1)/2} (−q;q)ₙ**,  **ψ₀*(q) = Σₙ≥0 (−1)ⁿ q^{(n+1)(n+2)/2} (−q;q)ₙ**

(ψ₀, ψ₁ are the same series without the (−1)ⁿ.) Whether the two leftovers are mock in the modern sense depends entirely on what ψ₀*, ψ₁* are. **The question: what kind of object is the pair (ψ₁*, ψ₀*)? It looks like an incomplete vector-valued quantum modular form, and we are looking for the missing component (or the missing log term).**

## 2. What is known about ψ₁*, ψ₀*

- **Coefficients:**
  - ψ₁*: 1, −1, −1, 1, 1, 1, 0, −1, −1, −2, 0, 0, 0, 2, 2, 1, 1, 1, −1, −1, −2, −2, −2, −2, −1, −1, 1, 2, 2, 2, …
  - ψ₀*/q: 1, 0, −1, −1, 0, 1, 1, 1, 1, −1, −1, −1, −2, −1, 0, 0, 1, 2, 2, 2, 1, 1, 0, −1, −1, −3, −3, −2, −2, −2, …
  - Dense (~93% nonzero), slowly growing (max |a| = 26 for n ≤ 120), wave-like sign pattern.
- **Bounded radially at all 19 tested roots of unity** (orders 1–8, 10, 12, 14). The leftover is a genuine power series in t (q = e^{−t}ζ), not exponentially small: the "strange/quantum" signature.
- **Strange-function behavior:** at a root ζ of **even** order m, the factor (−ζ;ζ)ₙ vanishes for n ≥ m/2, so the series terminates. The radial limit equals the finite sum (error ∝ t, checked). Values:
  - ψ₁*(−1) = 1
  - ψ₁*(i) = 2 − i
  - ψ₁*(e^{iπ/3}) = 1 − 2√3 i
  - ψ₁*(e^{iπ/4}) ≈ −3.5355339 − 4.1213203i
  - ψ₁*(e^{3iπ/4}) ≈ 3.5355339 − 0.1213203i
  - ψ₁*(e^{iπ/5}) ≈ −8.2082039 + 1.5388418i
  - ψ₁*(e^{iπ/6}) ≈ −4.4641016 + 12.196152i

  At odd-order roots they don't terminate (but stay bounded).
- **Negative results (don't repeat):**
  - No relatives: among ~100,000 short Eulerian series, nothing equals a·ψ₁* + b·ψ₀* + c·ψ₁ + d·ψ₀ + e (to q⁵⁹, also with q → −q and shifts by q^k), except trivial rewrites. The one elementary relation is Σ(−1)ⁿ q^{n(n−1)/2}(−q;q)ₙ = 1 − ψ₁* − ψ₀*.
  - Not σ(q) (Andrews–Dyson–Hickerson), not a false theta function, not a simple Rogers–Fine-type series.
  - **Not Hecke-type with a simple multiplier:** ψ*·P for P ∈ {(q;q)∞, (q;q)∞², (q²;q²)∞, (q;q)∞/(−q;q)∞, (q;q)∞(q²;q²)∞, (q;q)∞²/(q²;q²)∞, 1/(−q;q)∞, (q²;q²)∞/(−q;q)∞} is never sparse. Control: ψ₁·(q²;q²)∞ is sparse, max |c| = 2.

## 3. The transformation test (the main evidence)

**Setup** (`qmf.py` and friends):
- φ_F(x) := e^{2πiα x}·F(e^{2πi x}) at rationals x, computed exactly with the terminating sums (even denominators only).
- Group element γ = [[1,0],[2,1]], i.e. x ↦ x/(2x+1) ∈ Γ₀(2). It preserves even denominators; S: x ↦ −1/x does not.
- Cocycle: h(x) = φ_target(x) − Σⱼ εⱼ·(2x+1)^{−k}·φⱼ(γx), with the εⱼ fitted by closed-form complex least squares.
- **Roughness** R = Σ|Δh|²/Σ|Δφ_target|² over neighboring sorted rationals, x ∈ (0.05, 2), denominators ≤ maxden. Smooth h ⇒ R ∝ spacing², roughly maxden⁻⁴.

**Calibration** (Kontsevich F(q) = Σ (q;q)ₙ, known quantum modular, weight 3/2, φ = q^{1/24}F):
- Under S, the scan finds weight 3/2, α = 1/24, ε = −i: R = 2·10⁻¹⁰. All wrong parameters give ≥ 6·10⁻³.
- Under the Γ₀(2) map: same parameters, ε = e^{iπ/6}, R = 1.1·10⁻¹⁰.
- R scales 8·10⁻⁹ → 5·10⁻¹⁰ → 9·10⁻¹¹ → 1.3·10⁻¹¹ for maxden 24 → 36 → 48 → 64 (smooth).

**Results for ψ*:**

| Model | Best R |
|---|---|
| ψ₁* alone (k ∈ {0,…,5/2}, α ∈ j/240) | 0.089, flat in α, \|ε\| = 0.68: fails |
| Pair: ψ₁*(x) vs {ψ₁*(γx), ψ₀*(γx)}, k = ½, α(ψ₁*) = 1/6, α(ψ₀*) = 3/40 (grid best) | 1.3·10⁻⁶. Coefficients ε ≈ (0.43869+0.98266i, −0.13194+0.37248i), identical when fitted separately on m ≡ 0 and m ≡ 2 mod 4 (γ swaps these classes). **Plateaus**: 6·10⁻⁴, 2·10⁻⁵, 2·10⁻⁶, 1.5·10⁻⁶, 1·10⁻⁶ for maxden 24, 36, 48, 64, 80 |
| Pair + mock values ψ₁, ψ₀ (terminating sums) as extra columns | no improvement; their ε ≈ 0.01 |
| Pair, continuous Nelder–Mead over (α₁, α₀, k), target ψ₁* | 1.2·10⁻⁹ at α₁ = 0.179493, α₀ = −0.218409 (≡ 0.7816), **k = 0.507778**, ε ≈ (0.3883+1.0124i, −0.3417+0.2030i) |
| Same, target ψ₀* | 0.9·10⁻⁹ at α₁ = 0.388115, α₀ = −0.012957, **k = 0.507942**, ε ≈ (0.3453−0.1970i, 0.6832+0.8422i) |
| Held-out check of the continuous fits | m ≡ 2 mod 4: 1.2·10⁻⁹ (ψ₀*: 1.9·10⁻⁹); denominators 66–90: 1.1·10⁻⁸ (ψ₀*: 8.7·10⁻⁹) |
| Negative controls (same optimizer, unrelated pairs of (−q;q)ₙ-type series) | 10^{−4.3} and 10^{−1.0} |

**Reading:**
- The pair transforms together, far better than unrelated pairs, and the fit generalizes to held-out roots.
- But it isn't clean:
  - At the natural weight ½ the residual plateaus at ~10⁻⁶ instead of falling like a smooth cocycle.
  - The free fit drifts to weight ≈ 0.508.
  - The fitted normalization of ψ₁* is inconsistent between the two target fits (0.179 vs 0.388).

  That's what a misspecified model looks like.

## 4. Hypotheses for what's missing (in order of how I'd test them)

1. **A logarithmic cocycle term (higher-depth quantum modularity).** A term ∝ (cx+d)^{−k}·log(cx+d)·φ(γx) is approximately what a small weight shift k → k+δ mimics: (cx+d)^{−δ} ≈ 1 − δ log(cx+d). Test: fix k = ½ with rational α's and add columns (2x+1)^{−1/2}·log(2x+1)·φⱼ(γx). If R drops to ~10⁻¹⁰ and keeps scaling with maxden, that's strong evidence of depth-two behavior, as in Bringmann–Kaszian–Milas's higher-depth false theta functions.
2. **A third component.** Candidates that terminate at even roots, i.e. contain (−q;q)ₙ:
   - other members of the family Σ (±1)ⁿ q^{(An²+Bn)/2}(−q;q)ₙ (A ∈ {1,2,3}, small B);
   - series with (−q;q)₂ₙ or (−q;q²)ₙ;
   - the 1,357 "bounded everywhere" functions from the census (`labels.json` "quantum" in the full archive).

   Test each as an extra column with k = ½ and rational α ∈ j/120; look for R ≪ 10⁻⁹ **with clean α, k and multiplier entries of modulus-1 type.**
3. **Wrong group element or missing T-consistency.** Impose the T-action jointly: φⱼ(x+1) = e^{2πiαⱼ}φⱼ(x) ties the α's. Fit two generators of Γ₀(2) simultaneously (e.g. T and x ↦ x/(2x+1)), or try Γ(2) or Γ₀(4).
4. **Odd-order roots.** The test currently uses only even denominators, where sums terminate. Radial limits at odd-order roots exist (bounded) and could be computed numerically (slow, needs small t). Including them would allow S itself.

## 5. Reproducing the core computation (no archive needed)

```python
import mpmath as mp
mp.mp.dps = 40
def psi_star(h, m, shift):   # shift=0 -> psi1*, shift=1 -> psi0*; requires m even
    z = mp.e**(2j*mp.pi*mp.mpf(h)/m); s = 0; P = 1
    for n in range(m):
        e = (n*(n+1))//2 if shift == 0 else ((n+1)*(n+2))//2
        s += (-1)**n * z**e * P
        P *= (1 + z**(n+1))
        if abs(P) < mp.mpf(10)**-30: break
    return s
```

For the transformation fit, see `qmf_vec.py` (grid), `qmf_opt.py` (continuous), `qmf_val.py` (held-out), `qmf_negctl.py` (controls), `qmf_scale.py` (density scaling). `qmf.py` holds the shared helpers and the Kontsevich control. Needs `mpmath numpy scipy`. The full pipeline (radial-limit probe with Arb) is in the main archive if radial values at odd-order roots are needed (`catalog.py` and `probe.py`; `python-flint`).

## 6. Deliverable we'd like

Either of these, stated as a conjecture with the numerical evidence:
- "(ψ₁*, ψ₀*, [X]) is a vector-valued quantum modular form of weight ½ on Γ₀(2) with multiplier matrix M"; or
- "ψ₁*, ψ₀* are depth-two quantum modular forms with cocycle …".

Ideally with clean rational α's and an identifiable M. Plus a literature check of whether ψ₁*, ψ₀* (or the two halves) appear anywhere: Andrews–Berndt, *Ramanujan's Lost Notebook* (partial theta chapters); Bringmann–Folsom–Rhoades; Hikami; Kontsevich–Zagier strange identity; Bringmann–Kaszian–Milas on higher-depth false theta functions.
