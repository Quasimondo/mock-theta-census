# Stress test of the proofs for ψ₁\*, ψ₀\* (v1)

27 Nov 2026. Written by Claude Opus 5.5 at the request of Mario Klingemann, following `docs/psi-star-stress-test-brief-v1.md`. Subject: `docs/psi-star-theorems-v1.md` (section references below are to its headings). Scripts: `psi_star/stress/` (one per block, fixed seeds, `README.md`, complete output in `run_log.txt`).

## Summary

**No displayed inequality, identity or uniform-constant claim was violated on its stated domain.** Every claim of Priority 1 (S3, S4, S5) holds, most with a wide margin; so do those of Priority 2 and 3 that could be tested. Every apparent violation met on the way was traced to my own numerics (branch of a logarithm, truncation, cancellation, an under-resolved oscillatory tail, a contour that is not admissible for the e^{w}-version at large |t|). Each was reproduced and resolved before moving on; they are listed under *Numerical notes* in the README.

What the attack does produce:

- **The unstated constants of Parts 4–5, measured.** |H| ≤ C_U·e^{−0.185|s|} on the whole closed sector with **C_U = e^{−0.51} ≈ 0.60**, the same for |t| = 0.4 down to |t| = 0.001 (so t₀ = 0.4 works numerically); |H| ≤ C_ψ·e^{ψ(ρ)/|t|} with C_ψ = e^{−0.12} … e^{−0.14}. With this C_U the proof's C₂ in Corollary R is 48, and Theorem 4 holds with A = 0.37 for B = 240 (A is attained at K = 0, i.e. A ≥ |P·J|).
- **The margins are much larger than proved.** The true descent rate of Proposition 1 is κ_true = inf(−ψ/ρ) ≥ 0.97 on the sector (proved: 0.185), because the bound |Y| ≤ π/2 is loose (max |Y| = 0.37). The sharp constant in Lemma G / Corollary G is 1.49 (proved: L = 3.72). So the "thin margin 0.175 − 0.149" of §0.5 item 5 is an artefact of the bounds; the true margin is about 0.97 − 0.06.
- **Findings** (none breaks a theorem): explicit values for the constants left implicit (F1); a hypothesis missing in a contour move (v19 §4, e^{w}-version, F2); a gap in the wording of the analyticity argument (Corollary of v16, F3); the unstated ε(a) of Theorem 7, which can be taken as 60° (F4); the ranges X ≥ 4 / X ≥ 5 of Theorems 9, 11, which numerically are not needed (X = 2, 3 work, F5); and some smaller points, among them the numerical remark after Lemma F, which is right only near real x. Details below.

Labels used: **holds** (worst ratio < 1, with how close), **holds, different constant**, **violated** (none), **could not test**.

Precision: double precision throughout unless stated; quadrature of analytic integrands by the trapezoid rule (exponentially convergent) or composite Gauss–Legendre; Φ_b by products (both sides of the real axis) and by Faddeev's integral, cross-checked (agreement 5·10⁻¹² relative; against mpmath at 40 digits: products 2·10⁻¹⁴, Faddeev 5·10⁻¹¹).

## Table

"ratio" = LHS/RHS of the claimed inequality (holds if < 1); for identities the relative residual. Where a sup is reported, it is over random samples plus grids including the edges and corners of the domain, then local optimisation (Nelder–Mead) from the worst samples where indicated.

### Priority 1

| id | where | domain tested | samples | worst ratio / value (where) | verdict |
|---|---|---|---|---|---|
| S3.1 | [v15 §1.2] (i), both halves | δ = 10°, θ ∈ [−210°, 10°] (735 values; 30 each within 2° of the ends, 41 within 1° of −90°), ρ ∈ [10⁻³, 100] | 2.94·10⁶ | max \|ψ/ρ − formula\| = 8.7·10⁻¹³ | holds (identity) |
| S3.2 | (ii) \|Y\| ≤ π/2 | same | 2.94·10⁶ | max \|Y\| = 0.370 (ratio 0.236) | holds, far from sharp |
| S3.3 | (ii) λ_α ≤ X ≤ ln 2 | same | 2.94·10⁶ | max X/ln 2 = 0.9999 (ρ → 0); min (X − λ_α) = 8.2·10⁻³ | holds; upper bound sharp as ρ → 0 |
| S3.4 | (iii) ψ ≤ −κρ − cρ² sin 2δ | same | 2.94·10⁶ | min of the margin per unit ρ: 0.785 (θ = 10°, ρ → 0) | holds, wide margin |
| S3.5 | (iv) κ ≥ 0.185 on [−210°, 10°] | θ grid, then bisection | 38001 | min κ = 0.185115 at θ = 10° (formula κ ≥ 0.185 exactly on [−215.7°, 10.0°]; κ > 0 on [−220.6°, 22.4°]) | holds, no margin (by design) |
| S3.6 | ψ(ρ) ≤ −0.185ρ, ρ ∈ [10⁻³, 100] | as S3.1 | 2.94·10⁶ | true margin κ_true(θ) = inf(−ψ/ρ) ≥ 0.971 (θ = 10°, ρ → 0); beyond: κ_true > 0 for θ ∈ [−260°, 40°], < 0 at θ = 50° | holds, wide margin |
| S3.7 | (needed for the Corollary of v16, see F3) κ for δ ≠ 10° | δ = 5, 8, 9, 11, 12, 15°, θ ∈ [−210°, 10°] | 2201 per δ | min κ = 0.333, 0.245, 0.215, 0.155, 0.125, 0.036 | holds |
| S4.1 | [v14 §2.2] Lemma 1 / [v15 §2.2] Lemma 1′: \|ρ_K\| ≤ C(K,η)\|t₁\|^{K−1}/cos(arg t₁) | random ζ₀ (Re ζ₀ ∈ [−5, η/2], \|Im ζ₀\| ≤ 3π) and t₁ (\|t₁\| ∈ [0.004, 0.5], 35 % with \|arg t₁\| within 10° of 90°), hypotheses checked exactly on the ray; η ∈ {¼, ½, 1, 2}; a ∈ {½, 1, random}; K = 1…4 | 500 | C(K,η) estimates: K=1: ≤ 0.023; K=2: ≤ 0.0014; K=3: ≤ 0.0015; K=4: ≤ 2.3·10⁻⁴. No growth as \|t₁\| → 0 (max over \|t₁\| < 0.02 is 20× smaller) or as arg t₁ → ±90° | holds, explicit constants |
| S4.2 | Lemma A | θ ∈ [−210°, −150°] ∪ [−30°, 10°], ρ ∈ [0, 60] both halves; \|w\| ≤ η/4 clause random | 61 000 + 30 000 | min dist/η = 1.0000 (equality at ρ = 0, θ = −30°); clause: min dist/(3η/4) = 1.006 | holds, sharp |
| S4.3 | Lemma 3′ and [v16 §4]: \|H\| ≤ C·e^{−κ\|s\|} | θ: 41 values on the closed sector incl. −209.9, −151…−149, −91…−89, −90 ± 10⁻³, −31…−29, 9.9, 10; \|t\| = 0.4, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, and 0.003, 0.001 at five worst directions; 1100 contour points per (θ, \|t\|) for ρ ∈ [10⁻⁵, 40] | ≈ 3·10⁵ | sup log\|H\| + 0.185\|s\| = −0.506 (\|t\| = 0.4) … −0.510 (\|t\| ≤ 0.02), always at θ = 10°, ρ ≈ −0.45\|t\|. **C_U = 0.60, bounded as \|t\| → 0** | holds; constant found |
| S4.4 | Lemma 3′: \|H\| ≤ C·e^{ψ/\|t\|} | same | ≈ 3·10⁵ | sup log\|H\| − ψ/\|t\| = −0.117 (\|t\| = 0.4) → −0.141 (\|t\| → 0), at θ = 10°; at θ = −90° the value is −½ ln 2 = −0.347 for every \|t\| | holds; C_ψ = 0.89 |
| S4.5 | same, ψ₀\* integrand (extra e^{w}) | same | ≈ 3·10⁵ | sup log\|H₀\| + 0.185\|s\| = −0.64 … −0.51; sup log\|H₀\| − ψ/\|t\| − ρ = −0.16 … −0.14 | holds |
| S4.6 | beyond the sector | θ = −240 … −212, 12 … 35°, \|t\| = 0.05, 0.01 | 24 | log C_U ≤ −0.06 up to θ = 30°, +0.09 at 35° (still bounded); −0.31 at −240° | — |
| S5.1 | [v17 §1] Proposition 2 (shift identity) | 40 random t (θ ∈ [−210°, 10°], \|t\| ∈ [0.008, 0.35]), both versions, N ∈ {1, 2, random up to the hypothesis bound, ≤ 60} | 234 | worst: 0.16 in units of 10⁻¹⁶·Σ\|terms\|, i.e. at rounding level, at θ = −179.25°, \|t\| = 0.059, N = 37, where \|q\| > 1 and the τₙ reach 10¹⁷; spot checks without cancellation (θ = −210, −90, −89.9, −60, 10°) agree to 10⁻¹²…10⁻¹⁴ | holds |
| S5.1′ | hypothesis (N − ½)\|t\|cos γ < π cos α | θ = −90, −60, 0°, \|t\| = 0.3, N around the bound | 12 | holds up to the last admissible N (10⁻¹³), fails from the first N past it (2·10⁻³ … 19): sharp | holds, sharp |
| S5.2 | [v17 §2] Lemma G | ω with \|Re ω\| ≤ 8, \|Im ω\| ≤ 14 and segment distance ≥ 1 checked exactly; Nelder–Mead from the 8 worst | 146 068 | max log(\|g\|/(½e^{L\|ω\|})) → 0 only as ω → 0 (equality at 0); sup Re h(ω)/\|ω\| = 1.49 (L = 3.72); inf \|1 + e^{ω′}\| = 0.632 = 1 − 1/e (c₃ = 0.582) | holds; sharp L ≈ 1.5 |
| S5.3 | Corollary G | θ: 48 values, \|t\| = 0.01, 0.005, 0.002, all N ≤ 0.04/\|t\|, 801 contour points | ≈ 10⁶ | max of log\|D_N\| − (−νN + Lc₁\|s\|) = −0.141; with c₁ replaced by c = N\|t\| up to 0.3 still negative; true constant in place of L: 1.49 | holds |
| S5.4 | Corollary R | θ = −210 … 10° (9), \|t\| = 0.01, 0.008, all N ≤ 0.04/\|t\|, both versions | 160 | with the measured C_U: C₂ = 48; actual \|P·J − Σ_{n<N}\| / (C₂e^{−νN}) ≤ 6·10⁻³; ∫\|H\|\|D_N\|\|ds\| = 0.18 … 0.49 against the proof's 26.7 (N = 1) | holds; the integral is bounded |
| S5.5 | [v17 §3] Lemma T | n ≤ 300, 720 points on \|t\| = 0.3/n and 100 inside, both versions | 2.2·10⁵ | max ratio 0.727 (n = 1, t = −0.3); radius 0.5/n still 0.96; fails for 0.8/n | holds |
| S5.6 | coefficients aₖ, \|aₖ\| ≤ A₀k!B₀ᵏ | k < 46, ψ₁ and ψ₀, 130 digits | 92 | a₀, a₁, a₂ = 1/3, 1/9, 5/162 ✓; afterwards compared with v13's exact values: 1.7·10⁻⁷¹; max \|aₖ\|/(A₀k!B₀ᵏ) = 0.10; \|aₖ/k!\|^{1/k} → 0.262 (1/\|S\| = 0.2905) | holds |
| S5.7 | [v17 §4] Theorem 4, B = 240 | 11 directions incl. −210, −209.9, 9.9, 10°; \|t\| = 0.01, 0.05, 0.3; K ≤ 40; both versions | 2700 | A needed = 0.37 (at K = 0); effective B ≤ 0.294 (≈ 1/\|S\| = 0.2905) | holds; A = 0.37 |
| S5.8 | [v17 §5] Theorem 5 | θ₀ = −120°, −80°; t = e^{iθ₀}z with \|z\| = 0.97R cos φ, φ = −85 … 85° (near the edge of C_R); R = 0.01 (stated), 0.1, 0.3 (beyond); Borel transform by series + Padé (checked against a second Padé) | 28 | \|P·J − 𝒮_{θ₀}[Φ]\|/\|P·J\| ≤ 1.6·10⁻¹³ (5·10⁻¹¹ at \|t\| = 8·10⁻⁴, limited by the product sum) | holds, also at R = 0.3 |

### Priority 2

| id | where | domain tested | samples | worst ratio / value (where) | verdict |
|---|---|---|---|---|---|
| S2.1 | [v16 §1] Lemma Φ (i) | \|x\| ∈ [0.01, ¼], arg x ∈ (−89°, 89°) incl. ±60, ±89; w in the strip | 318 | Faddeev / products − 1: 4.7·10⁻¹² | holds |
| S2.2 | Lemma Φ (ii) at real b | X = 4 … 100 | 40·5 | functional equation (b): 6·10⁻¹²; (1/b), integral itself for π < Im w < π(1+x) vs shifted value: 7·10⁻¹⁴; inversion: 7·10⁻¹⁵; Φ_b(0)²: 6.5·10⁻¹² | holds |
| S2.3 | [v16 §2] Lemma S, defect/\|t\|^{2K−1} | K = 2, 3; ε₁ = π/2, π/6, 0.1; arg x ∈ {0, ±10, ±30, ±45, ±60°}; \|t\| = 1.5 … 0.05; Re w ∈ [−40, 4], Im w on 17 lines incl. ±(π − ε₁); then 50-digit products down to \|t\| = 0.015 at the worst point | ≈ 4·10⁴ | K = 2: 6.1·10⁻⁴, 0.018, 2.6; K = 3: 8·10⁻⁵, 0.021, 78 (for ε₁ = π/2, π/6, 0.1); worst always at arg x = ±60°, Im w = ±(π − ε₁); saturates as \|t\| → 0 (2.2 and 65 at \|t\| = 0.015) | holds |
| S2.4 | Lemma F | \|t\| = π/2 … 0.01, arg x = ±60, ±59.9, ±45, ±30, ±20, ±10, ±5, 0°; Re w ∈ [−40, 0], 43–81 lines incl. Im w = ±π; Nelder–Mead | ≈ 10⁵ | sup = **0.511** (\|x\| = ¼, arg x = ±60°, w = −0.441 ± iπ); 0.432 as \|x\| → 0 at arg x = ±60°; 0.3466 at real x | holds with C₁ = 0.52; the remark "0.347" holds only near real x (smaller points) |
| S2.5 | [v16 §3] Lemma D | margin on \|β\| ≤ 30°; E on the contour for Im x ≠ 0 by products | 60001; 4·10·1200 | margin 0.383022 at β = −30° (no slack, as stated); ρ ≥ 4.8875; \|E\| at the double floor for \|t\| ≤ 0.2, 1.5·10⁻¹³ = e^{−11.8/\|t\|} at \|t\| = 0.4 (doc: e^{−11/\|t\|}) | holds |
| S2.6 | [v19 §2] Lemma B (i) | 0 ≤ arg x ≤ 60°, \|t\| = π/2 … 0.02, closed strip, Re w ∈ [−40, 40] | ≈ 3·10⁴ | C = e^{0.40} (\|x\| = ¼), e^{0.37} (\|t\| = 0.02) | holds |
| S2.7 | Lemma B (ii) | ε₁ = π/6, π/2 | ≈ 2·10⁴ | C = 0.090, 0.043 (sup of \|ratio − 1\|/\|t\|) | holds |
| S2.8 | Lemma B (iii) | δ₁ = 30, 60, 80°, δ₁ + arg x ≤ 90°, ρ ∈ [0, 60] | ≈ 1.2·10⁴ | sup of log(1/\|Φ_b\|) − Re(Li₂/t) − N ln 2 = −0.036 < 0; without 2^N the sup is ln 2 exactly: the bound holds with 2 instead of 2^N | holds |
| S6.1 | [v19 §5] closed formulas on Γ₁, Γ₂, Γ₃ | random points, ρ ≤ 40, η ∈ [0, π], r ≤ 30 | 6500 | 2·10⁻¹³ … 0; Im V − 0.433ρ² + 0.01 ≥ 9.8·10⁻³; Re V ≤ −π²/3 (equality at ρ = 0); Im Li₂ ≥ 0 for ρ ≤ 2π; \|Li₂\| ≤ 0.0043 < 0.01 for ρ ≥ 2π | holds |
| S6.2 | Γ₂: g(η) − g(2π/3) ≤ −0.015(η − 2π/3)² | τ = 0, 0.1, 0.2 | 3·1000 | max of g − g(2π/3) + 0.015(…)² = −2.9·10⁻⁷ (attained next to the saddle, where both sides vanish) | holds |
| S6.3 | (L1) on Γ, 0 ≤ τ ≤ 0.2 | Γ sampled with 5600 points | 15 τ + bisection | c_eff(τ) = 0.313 (τ = 0), 0.163 (τ = 0.2); holds for **τ < 0.3085 (arg x < 17.15°)**, fails beyond | holds; range larger than proved |
| S6.4 | [v19 §4] monotonicity of Im V on the horizontal line | a random in (0, 1) and a = 10⁻³, 10⁻², 0.5, 0.99, 0.999 | 305 | sign pattern never violated; valley margins → 0 only as a → 0 (1.7·10⁻⁴ at a = 10⁻³) | holds |
| S6.5 | [v19 §3] Theorem 6, (exact/asymptotic − 1)/\|t\| | growing: arg x = 0, 2, 5, 8, 11, 11.3°; decaying: 0, 10, 30, 45, 60°; \|t\| = 0.5 … 0.01; exact integrals by quadrature on the descent contours (checked against L: 10⁻¹¹ … 2·10⁻⁵ at \|t\| = 0.5 … 0.8) | 140 | growing 0.122, its b-version 0.252, decaying 0.089, its b-version 0.268; flat in \|t\| | holds (uniform, incl. arg x = 0 and the edges) |
| S6.6 | §0.2 Theorem III forms (P·I, and q^{1/2}P·I(b ± 1/(2b))/P·I(±1/(2b))) | arg x = 0, 5, 11°, \|t\| = 0.2, 0.05, 0.02 | 36 | ≤ 0.35 | holds |
| S6.7 | [v19 §3] Theorem 7, ε(a) | 17 values of a incl. 0.01, 0.05, 0.95, 0.99; ratio at \|t\| = 0.2 … 0.02, arg x = 0 and 30° | 17 + 15·8 | M = 0 for every a tested: Re V has its maximum at w₁ along the line, so (L1) holds for all τ ≥ 0 and ε(a) = 60° (limit of Lemma B); (exact/asymptotic − 1)/\|t\| ≤ 0.09. For a = 0.01 and 0.99 the valley margins (0.017, 0.028) are too small for my quadrature: integral not computed | holds; ε(a) found (F4) |
| S6.8 | Lemma L | its conclusion on the actual integrals = S6.5–S6.7 | — | O(t) uniform | holds (statement itself not tested, see below) |
| S7.1 | [v20 §1] Lemma R1, block bound | X even ≤ 40, m ≤ 300, δ ∈ [10⁻⁷, 1/X] | 4000 | \|B_m\|/min(1 + Cδ, (m+1)Cδ) ≤ 0.74; step \|B_m − (1 − c^X)\| ≤ X²3^{X−1}δ: 0.019; 1 − c^X ≤ min(1, (m+1)X²δ): 1.0000 (sharp); in all samples \|B_m\| ≤ 1, so C = X² would do | holds |
| S7.2 | Lemma R1, the limit | X = 4, 6, 10, 20; δ = 10⁻², 10⁻³, 10⁻⁴ | 12 | error ∝ δ (X = 4: 0.10, 0.013, 0.0014) | holds |
| S7.3 | Lemma R2, each inequality | δ ∈ [10⁻⁴, 0.1] log-uniform, p on and inside \|1 − p\| = (3/2)(1 − \|p\|) | 400 | \|1 − pᵏ\|/(kKδ) ≤ 1.0000 (k = 1 on the edge); (2j−1)!!(Kδ)ʲ/min(Kδ, √2e^{−j}) ≤ 1.0000; tail bound: 0.032 | holds |
| S7.4 | Lemma R2, the bracket | δ ∈ [0, 0.5] | bisection | 0.0766 at δ → 0; positive for **δ < 0.104** (with the e^{−1/K+δ} that the argument gives: δ < 0.172); positive at δ → 0 for any K < 1.763 | holds ("tends to") |
| S7.5 | Lemma R2, W, W₃ → 1 | edge of the region and the path of v20 (X = 6) | 10 | \|W − 1\|, \|W₃ − 1\| ∝ δ (1.5·10⁻³ at δ = 10⁻³ on the edge); on the path \|1 − p\|/(1 − \|p\|) = 1.0000 | holds |
| S7.6 | [v21 §1] Lemma R2′ | X = 5, 7, 11, 21; δ = 0.1 … 10⁻³; on the path and off it (φ = 4δ², δ, 3δ) | 40 | φ/δ² ≤ 0.071 (→ 1/(πX)); min \|(p;p²)ₙ\|/(9/8)ⁿ = 1.69 on the path, 1.13 at φ = 3δ; \|(p;p²)ₙ\|e^{0.3/δ} ≥ 23; log G − bound ≤ −0.153; δ^{1/2}e^{π²/(8δ)}\|ψ(p)\| = 3.28 (constant); \|P_±\|/bound ≤ 3.7·10⁻⁶; F_a + ⅓, F_b − ⅓, W − ⅓, pW₃ + ⅓ ∝ δ | holds |
| S7.7 | [v20 §1] Lemma R3 | X = 4, 5, 8; ε = 10⁻², 10⁻³, 10⁻⁴ | 27 | \|value(ε) − value(0)\|/\|value(0)\| ∝ ε for J, I(1/(2b)), I(−1/(2b)) | holds |
| S7.8 | [v20 §2] Theorem 9, [v21 §2] Theorem 11 | X = 2 … 41 (21 values), ψ₁ and ψ₀; X ≤ 8 with the integrals over L, X > 8 over the descent contours | 42 | relative error ≤ 1.2·10⁻¹² for X ≤ 5 and ≤ 6·10⁻¹² for X ≥ 9; 7·10⁻¹¹ … 8·10⁻⁹ at X = 6 … 8 (cancellation over L); **X = 2 (x = ½) and X = 3: 4·10⁻¹⁴ (ψ₁\*), 1·10⁻¹², 3·10⁻¹³ (ψ₀\*)**; X = 4: 4·10⁻¹³ | holds; also outside the stated range (F5) |
| S7.9 | Theorems 10, 12 | same | 21 | X·\|ψ₁/lead − 1\| ≤ 3.7 (X = 3), ≤ 2.5 for X ≥ 5: O(1/X) | holds |

### Priority 3

| id | where | domain tested | samples | worst ratio / value (where) | verdict |
|---|---|---|---|---|---|
| S1.1 | [v12 §1] Theorem C | random \|u\| ≤ 3, \|q\| ∈ [0.1, 0.9], Re x ∈ [−½, ½], arg x ≤ 120°, δ ∈ [¼, ¾](90° − β); quadrature on L | 150 (21 ill-conditioned skipped) | max \|I_quad/I_C − 1\| = 3·10⁻⁷ (conditioning of D(u) up to 10⁸ allowed); a skipped case (D cancels by 4·10¹³) redone at 50 digits: 1·10⁻¹³ | holds |
| S1.2 | (B1), both recursions | same, by quadrature | ≈ 130 | residual / size of the terms: 4·10⁻⁹ (p⁻¹ can be 10³⁰) | holds |
| S1.3 | [v12 §4.1] G(w), F(w) | random w ∈ H | 80 | 3·10⁻⁹, 2.4·10⁻⁷ | holds |
| S1.4 | [v12 §4.2] steps (X), (Y) | K(λ) by quadrature | 80 | 4.7·10⁻⁸, 4.5·10⁻⁸ | holds |
| S1.5 | Theorem A (A↑), (A↓) | J by quadrature, residue sums Ψ, Ψ̌, dual sums with I from Theorem C | 120 (61 skipped: D(u ∓ m/b) cancels by > 10⁸) | 6.5·10⁻¹⁰ (both routes) where the residue sum and the dual sum do not cancel; where they are ≈ 10²⁰ and cancel to \|J\| ≈ 10⁻³, the residual is 5·10⁻¹⁰ of their size | holds |
| S1.6 | (B2), (B3) | same | 120 | 5·10⁻⁶, 6.5·10⁻¹⁰, (B3): 1.8·10⁻¹⁴ | holds |
| S1.7 | §0.2 Theorem I, four formulas, all integrals by quadrature | random x as S1.1 | 60 | ψ₁\* 3.1·10⁻¹¹, f₁\* 3.4·10⁻¹³, ψ₀\* 4.6·10⁻¹³, f₀\* 3.8·10⁻¹³ | holds |
| S1.8 | link between Part 1 and Parts 3–5 | P·J(1/(2b)) over L (v12) vs e^{t/8}∫H ds over the s-contour (v14–v17); P = 2^{−1/2}e^{π²/(12t)+t/6}/b = (−q;q)∞/(b q^{1/8}(p;p²)∞) | 60 | 2.7·10⁻¹³ (both versions); 2·10⁻¹⁵ | holds |
| SA.1 | [v2 §2] (I), (II); double-sum forms; [v4 §1] (V), (VI); G(−q⁵), G(−1/q) | exact integer power series | to q⁴⁰⁰ | all coefficients zero | holds |
| SA.2 | [v3 §1] step (c): \|(−q;q)_{mJ}\| ≥ 1.4ᴶ if r^{Jm²} ≥ ½ | m = 3 … 11, 1 − r = 10⁻², 10⁻³, 10⁻⁴ | 15 | min ratio 1.25 | holds |
| SA.3 | step (d): θ₄R₁ → 0 | m = 3, 5, 7; m²t = ½ … 1/32 | 15 | m²t·log\|θ₄R₁\| = −3.16 … −3.13 at m²t = 1/32, tending to −π²/3 = −3.29 (bound used: −(π²/4 − 1)) | holds |
| SA.4 | the odd-root rule | all primitive roots of order m = 3, 5, 7, 9; 1 − r = 10⁻² … 10⁻⁴ | 3·12 | error ∝ 1 − r (m = 3: 3.3·10⁻⁵ at 10⁻⁴); A₀ = 0 for m = 3 | holds |

## Findings

Ordered by how much of the chain (§0.3) depends on the statement. None of them is a violation; none breaks a theorem.

**F1 (uniform constants of Parts 4–5; §0.5 items 3 and 5).** Not an error: the constants the proofs leave implicit, measured.

- C_U = e^{−0.51} = 0.60 in |H| ≤ C_U·e^{−0.185|s|} on the closed sector −210° ≤ θ ≤ 10°, for every |t| from 0.4 down to 0.001; the maximum sits at θ = 10° (the end with κ = 0.185115), at ρ ≈ −0.45|t|. The ψ₀\* integrand gives the same up to 0.13 in the logarithm. So t₀ = 0.4 works, and nothing degenerates as |t| → 0 or at the edges. Reproduction: `s4_uniform.py`, part (c).
- With this C_U, C₂ of Corollary R is 48, and the factorial bound of Theorem 4 holds with A = 0.37 for B = 240 (A is |P·J| at K = 0; for K ≥ 1 the ratio is tiny).
- The proved margins are far below the true ones: κ_true ≥ 0.97 (Proposition 1 proves 0.185 because |Y| ≤ π/2 while max |Y| = 0.37), and the constant of Lemma G is 1.49, not L = 3.72. The proof is not affected (it needs only κ₀ > L·c₁), but a referee worried about "0.175 vs 0.149" can be told that the real gap is wide. Sharpening |Y| (for instance |Y| ≤ arctan of the maximal |Im z|/(1 + Re z) on the relevant set) would make Corollary R robust.

**F2 (Part 6, [v19 §4], "Moving the contour", for I(b + 1/(2b)) at real x).** For the kernel-free integral with the extra factor e^{w} the horizontal line through w₁ is an admissible contour only if the integrand still decays on the right, i.e. (2η₁ + 2πa) + |t| < 0; for a = ½ this is |t| < π/3 (X > 6 at x = 1/X). The text proves the move for I(a/b) ("on the right e^{r(2η + 2πa)/|t|}") and for the e^{w}-version only states (L2) ("at most e^{|w|} elsewhere"). Since Theorem 6 is asymptotic (|t| → 0) this is a missing hypothesis, not an error. Minimal reproduction: x = ¼ (X = 4), a = ½: on the line w = r − 2πi/3 the integrand of I(b + 1/(2b)) grows like e^{r(1 − (π/3)/|t|)} = e^{0.33 r}; before I added a guard, the truncated line integral differed from the value over L by a factor 10¹³ (and by 6·10⁵⁵ at X = 2, 10⁻² at X = 6, 4·10⁻⁷ at X = 7, where the decay is slow). `s6_state.I_decay` now returns NaN when the condition fails. Anyone checking Theorems 9, 11 numerically at small X must take the b-versions over L.

**F3 (Corollary of v16, analyticity of P·J).** The proof says "the integrand is analytic in b by Lemma Φ, the contour angle is fixed". It is not: γ = β + 10° moves with t. The argument needs, near each t₀, one contour that works for all t near t₀, i.e. descent for δ = γ(t₀) − β(t) slightly different from 10°, and — for Im x ≤ 0, where P·J is *defined* by the δ = 10° contour — that the integral does not depend on the admissible contour (stated only for Im x > 0, v14 §2.1; it follows from the bounds of Lemma E type, but is not written). Both are true: Proposition 1 (iii) with δ ∈ [5°, 15°] still gives κ ≥ 0.036 on the whole sector (S3.7), and the bounds of S4 hold on the rotated contours as well. A sentence would close it.

**F4 (Theorem 7: ε(a) not given; §0.5 item 6).** On the horizontal line through w₁ the quantity M of the proof (sup of (Re V − Re V(w₁))/min(…)) is 0 for every a tested in (0, 0.99]: Re V has its maximum at w₁ along the line, so (L1) holds for every τ ≥ 0 and one can take **ε(a) = 60°**, the range of Lemma B. An analytic reason for part of the range: Re V″ = −2 + Re(1/(1 + e^{w})) < 0 on the line when η₁ ≥ −5π/6, which is a ≤ 0.711; for a > 0.711 Re V is not concave on the whole line but still maximal at w₁ (numerically). The constants of (L1) degenerate only as a → 0 (c₀ = 0.0048 at a = 0.01, valley margins → 0) and as a → 1 (η₁ → −π, the line approaches the edge of the strip, where Lemma B (ii) loses its constant). Reproduction: `s6_state.py`, part (d).

**F5 (Part 7, ranges X ≥ 4 and X ≥ 5).** Theorem 9's three-integral formula holds numerically at X = 2 (x = ½) and Theorem 11's at X = 3, for ψ₁\* and ψ₀\*, to 4·10⁻¹⁴ … 10⁻¹² (`s7_roots.py`, part (e); at X = 3, A₀ = 0). The odd-root rule itself holds at m = 3 (SA.4). The restriction comes from the range |x| ≤ ¼ of Part 2, but Lemma R3's own remark says that what is used there is Im(i/x) ≥ 2, which holds at x = ⅓ (3) and x = ½ (2). Lemma R2′ needs nothing about X. So X = 3 (and X = 2) can probably be included at the cost of rechecking Lemma F for Im(i/x) ≥ 2 instead of |x| ≤ ¼ — note that Lemma F's constant grows towards the corner |x| = ¼, arg x = ±60° (S2.4), so this should be checked, not assumed.

**Smaller points (no downstream effect).**

- **Lemma F, numerical remark (§0.5 item 4).** "Numerically the supremum … is 0.347 (attained at w = ±iπ), the same for real x and for the neighbouring directions, and independent of |t|". That is right for real x and |arg x| ≤ 10° (0.3466 at |t| = 0.01, 0.348–0.354 at |x| = ¼; 0.3466 ≈ ½ ln 2), but on the lemma's domain |arg x| ≤ 60° the supremum is 0.511 (|x| = ¼, arg x = ±60°, w = −0.441 ± iπ) and 0.432 as |x| → 0 at arg x = ±60°; it depends on both |t| and arg x there. The lemma itself holds with C₁ = 0.52. Reproduction: `s2_faddeev.py`, part (c).
- **Lemma B (iii).** The factor 2^N is not needed numerically: the bound holds with the factor 2 (sup of log(1/|Φ_b|) − Re(Li₂/t) on the rays is ln 2).
- **(L1) for the growing contour.** It holds for τ < 0.3085, i.e. arg x < 17.15°, not only 11.3°; the 0.2 comes from the chord bound for f. Theorem 6 (growing) could be stated up to 17°.
- **Lemma R2.** The bracket with K = 3/2 is positive only for δ < 0.104; fine for "tends to 0", but with the exponent e^{−1/K+δ} that the argument actually produces the range is δ < 0.172, and any K < 1.763 works.
- **Notation in Theorem C.** In D(u), p^{n²/2} must be read as e^{−πi·n²/(2x)}. With the principal power of p the double series changes sign in its odd-n terms when arg p is reduced, and Theorem C fails numerically. §0.1 could say so.
- **Lemma A, Lemma D, Proposition 1 (iv).** Their constants are sharp at an edge (equality at ρ = 0, θ = −30°; 0.383022 at β = −30°; 0.185115 at θ = 10°); there is no slack to absorb a rounding error.

## Not tested

- **Closing lemma [v12 §4.1]**: the two exponent bounds on the vertical line with the constant c, and the remark that the supremum is "this bound times a factor between 1 and 1.5". Not done for lack of time; its consequences (the Fourier transforms F, G and Theorem C) were tested (S1.1, S1.3).
- **Lemma L** as an abstract statement (its proof is sketched): only its conclusion on the actual integrals (S6.5–S6.7).
- **Lemma E** as a separate statement: only through Proposition 2, which uses it (S5.1).
- **Lemma 4′, the moments, the identification of the coefficients (Theorems 1–3)**: superseded by Theorem 4 (§0.3); the coefficients were checked directly (S5.6) and Theorem 4 itself (S5.7).
- **The Nevanlinna–Sokal theorem** (quoted) and the **figure-eight corollary** (quoted Kashaev asymptotics).
- **Theorem 8 / Theorem IV as a whole**: only its three ingredients (Theorem I, S1.7; Theorem 5, S5.8; Theorem 6 and the §0.2 forms, S6.5–S6.6). A direct test needs ψ₁\* at |q| close to 1 in a sector together with Borel sums, which I did not set up.
- **The bounds (2.1), (2.2) of v12 §2**: used, not tested separately.
- **More than double precision** for Theorems 9, 11 at large X: double precision with the descent contours gives 10⁻⁹ … 10⁻¹²; nothing suggested a need for more.
- **Theorem C, A, B outside the sampled set** (|q| < 0.1 or > 0.9, arg x > 120°, δ near 0 or near 90° − β): there the quadrature on L is too ill-conditioned in double precision; 21 of 150 samples in S1.1 were skipped for that reason (one of them redone at 50 digits).
