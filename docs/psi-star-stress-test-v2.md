# Stress test of the proofs, second round (v2): addenda v22–v26

27 Nov 2026. Written by Claude Opus 5.5 at the request of Mario Klingemann, following `docs/psi-star-stress-test-brief-v2.md`. Subject: `docs/addendum-v22.md` … `addendum-v26.md`, and the parts of `docs/psi-star-theorems-v1.md` that round 1 left untested. Scripts: `psi_star/stress2/` (one per block, fixed seeds, `README.md`; complete output in `run_log.txt`).

## Summary

**No claim was violated on its stated domain.** All of Priority 1 holds:
- Theorem 13 (a), (b);
- Theorem 14;
- Lemmas U and J;
- Theorem 15.

So do the Priority 2 and 3 items I could test. As in round 1, every apparent violation along the way came from my own numerics, and each was resolved before moving on (see *Numerical notes* in the README):
- rounding next to the saddle (spurious positive values of the (L1) ratio on C₄);
- the growth of the kernel on one side of a tilted line;
- a product of several hundred factors in double precision.

What this round adds:

- **(L1) on contour C keeps its sign right to the ends.** The margin is proportional to the volume: sup (H − H(w₀))/min(d², d) ≈ −0.315·Vol(z) as a → 0 and ≈ −0.118·Vol(z) as a → 1. It is negative on all five pieces for every a from 10⁻⁵ to 1 − 10⁻⁵. The bound for 1/Φ_b on C holds with the constant 2 exactly, for every a and every X tested. The worst point is the corner 2πi, where the factor of the functional equation equals 2.
- **ε(a) of Theorem 13 (b), measured.** The largest arg x for which (L1) holds on the tilted contour, optimised over the tilt φ₁: 20° near a = 0.3, 17° at a = ½, down to 0.08° at a = 10⁻³ and 0.02° at a = 1 − 10⁻³. Near the ends it is roughly proportional to the volume. Theorem 13 (b) itself holds up to that angle, with a flat O(t) constant.
- **Lemma J is sharp to 3.5 %.** The worst ratio is |P·J(a/b)|·2√2 sin πa = 0.9655, at a = 0.284, x = 1/16. It tends to 0 linearly at both ends of (0, 1). The bound also holds for x up to 1, beyond the stated x ≤ ¼.
- **Theorem 15: the O(1/X) constant has a limit for every a, but it is not uniform in a.** K(a) = lim X·(F/prediction − 1) stays between 1.15 and 1.6 in modulus over the whole of (0, 1), and tends to about −1.152i at both ends. The X needed to reach that limit, however, grows like 1/Vol(a): K settles once X·Vol(a)/2π is about 10 to 30 (within 1 %). Before that, |K| reaches 650 at a = 10⁻³, X = 800. The statement ("the constants depend on a and degenerate") is correct; it is not uniform in a, and nowhere claims to be.
- **v24 and v25 hold everywhere they could be tested.**
  - (R3) and (R4) hold exactly to q¹⁰⁰⁰.
  - Theorem 17 holds on the unit circle, X = 1 … 100, and J(1/b), J(b + 1/b) are continuous up to the real axis (error ∝ ε).
  - Lemma F holds on [¼, 1] with sup 0.383.
  - Lemma D gives |E| ≤ 0.17, against the claimed 1.2.
  - Theorems 9, 11 and 14 hold at X = 1, 2, 3. In 25-digit arithmetic the agreement is 2·10⁻¹⁷ … 10⁻¹⁵ at X = 1, 2, but only 2·10⁻¹⁴ … 10⁻¹² at X = 3, so there my quadrature did not get beyond double precision.
- **Round 1's untested items:**
  - the closing lemma's constant c is sharp: it is the leading term, approached from below;
  - Theorem 8 holds as a whole, uniformly up to |p| = 0.77;
  - Theorem C holds at 30 digits outside the round-1 sample;
  - the corrections of v26 say what round 1 asked for, and Lemma V's internal claims check out.

Labels as in round 1: **holds** (with how close), **holds, different constant**, **violated** (none), **could not test**.

## Table

"ratio" = LHS/RHS of the inequality (holds if < 1), or the sup of a quantity that must be negative; for identities the relative residual.

### Priority 1

| id | where | domain tested | samples | worst ratio / value (where) | verdict |
|---|---|---|---|---|---|
| T13a.1 | [v22 §1] saddle, V′(w₀) = 0, Im V(w₀) = −Vol; [v22 §2] closed forms on C₁, C₂, C₃ (linear), C₅ | 81 values a ∈ [10⁻³, 1 − 10⁻³] | 81·(1 + 150) | ≤ 8·10⁻¹⁴ | holds |
| T13a.2 | [v22 §2] (L1) at τ = 0 on C₁ … C₅ | the same 81 a (clustered at both ends), plus a = 10⁻⁴, 10⁻⁵, 1 − 10⁻⁴, 1 − 10⁻⁵; 3000–9000 points per piece, clustered at the corners 2πi, iπ, r₀ + iπ and at the saddle, then zoomed at the worst point; Taylor series at the saddle for d < 10⁻³ | ≈ 2·10⁶ | sup < 0 on every piece for every a. Least negative −3.6·10⁻⁴ (a = 0.999). The ratio sup/Vol tends to −0.315 (a → 0, attained on C₁ and C₂ at the corner 2πi) and −0.118 (a → 1, on C₄ and C₅). v22's −0.0036 (a = 0.99) and −0.0095 (a = 0.01) reproduced | holds; margin ∝ Vol |
| T13a.3 | [v22 §2] (L2): \|e^{Q/t}/Φ_b\| ≤ C·\|e^{V/t}\| on each piece | a = 0.01 … 0.99 (7), X = 4 … 256, Faddeev's integral (functional equation once for C₁, C₂) | ≈ 10⁵ | log C = ln 2 exactly on C₁, C₂ (the corner 2πi, where 1 + p·w′ = 2); ½ ln 2 on C₃; ≤ 0.35 on C₄; 0 on C₅; independent of a and X | holds, C = 2 |
| T13a.4 | [v22 §5] −Im V rises along the horizontal line through the growing saddle, both directions | the 81 a plus 10⁻⁴, 1 − 10⁻⁴; s ∈ [10⁻⁵, 60] | 5·10⁵ | min of (Im V(w₀) − Im V(w₀ ± s))/min(s², s) = 4.8·10⁻⁵ > 0 (a = 0.9999); 1.20, 7.36 at s = 2, 8 for a = ½ and 3.07, 16.4 for a = 0.2, as v22 says | holds (margin → 0 as a → 1) |
| T13a.5 | Theorem 13 (a), (exact/asymptotic − 1)/\|t\| | a = 0.02 … 0.98 (11), X = 8 … 128; exact integral over C, checked against the 45° line (10⁻¹³ … 10⁻¹⁵ where that is well conditioned) | 47 | 0.017 (a = 0.98) … 0.13 (a = ½) once X·Vol/2π ≳ 2; larger before that (4.8 at a = 0.02, X = 64, falling to 2.8 at X = 128) | holds; constant locally uniform, as stated |
| T13b.1 | [v22 §3] (L1) at τ = 0 on C′ | φ₁ = 1°, 5°, 20°; 20 values of a | 60 | max sup/Vol = −0.118 for φ₁ = 1°, 5°; **positive for φ₁ = 20° at a ≤ 0.1**, as the proof's choice (3π²/2)·tan φ₁ < Vol/2 anticipates | holds where the proof chooses φ₁ |
| T13b.2 | Im V ≥ −(3π²/2)·tan φ₁ on C₂′ | φ₁ = 1, 5, 20°, 20 a | 2.4·10⁵ | min of Im V + (3π²/2) tan φ₁ = 0.26 > 0 | holds |
| T13b.3 | \|p·w′\| ≤ 1 on C₁′ for arg x ≤ φ₁ | φ₁ = 0.5 … 45°, \|x\| = ¼, 1/20, 1/100 | 6.6·10⁴ | max \|x\|·log\|p w′\| = −7·10⁻¹⁰ (equality at the end point P) | holds |
| T13b.4 | ε(a): largest arg x with (L1) on C′, optimised over φ₁ ∈ [0.5°, 45°] | 20 a | 200 bisections | ε = 0.08° (a = 0.001), 0.39 (0.005), 0.78 (0.01), 1.55 (0.02), 3.9 (0.05), 7.8 (0.1), 16.2 (0.2), 20.1 (0.3), 17.2 (0.5), 9.6 (0.7), 6.2 (0.8), 3.0 (0.9), 1.4 (0.95), 0.25 (0.99), 0.02 (0.999). v22 §6 (φ₁ = 5°): 3° for 0.05 ≤ a ≤ 0.75 and 5° for 0.1 ≤ a ≤ 0.75, reproduced (both also hold at a = 0.8; at a = 0.9 the limit is 2.99°) | found |
| T13b.5 | Theorem 13 (b), (exact/asymptotic − 1)/\|t\| at arg x = ε(a)/2 and ε(a) | a = 0.1, 0.3, 0.5, 0.7, 0.9; \|t\| = 0.2, 0.1, 0.05; exact integral over C′ | 30 | ≤ 0.10, flat in \|t\| and in arg x | holds up to ε(a) |
| T14.1 | [v22 §4] Theorem 14 (and v25 Theorem 18 (c) at X = 2) | X = 2, 4, …, 20, 24, 30, 36, 40, 50 (X = 60 was started but did not finish in an hour); four random a and a = ½ each; J over the 45° line, I((a−1)/b) over C, I(a/b) over the horizontal line through the decaying saddle | 75 | 10⁻¹⁴ (X = 2) … 8·10⁻¹² (X = 50), growing with the conditioning; one outlier 3·10⁻⁵ (X = 36, a = 0.976): the horizontal line was too long and the fallback line has conditioning 10¹³ | holds |
| T14.2 | the same with all three integrals over the 45° line | X = 2 … 12, five a | 25 | 3·10⁻¹⁴ … 4·10⁻⁹ (X = 12, a = 0.1); error ≈ 10⁻¹⁶ × conditioning (up to 8·10⁷) | holds |
| T14.3 | dual series on x = 1/X + iε | X = 4, 8, 20; a = 0.1, ½, 0.9; ε = 10⁻², 10⁻³, 10⁻⁴ | 27 | all three defects → 0 linearly in ε once ε ≪ 1/X²; at ε = 10⁻⁴: ≈ 5·10⁻³ (X = 4), 2·10⁻² (X = 8), 8·10⁻² (X = 20) | holds (limit) |
| T14.4 | odd X (no claim) | X = 3 … 21 | 24 | relative error 4·10⁻³ (X = 21, a = ½) … 2.9 (X = 3, a = 0.1): the even-X formula is not exact at odd X, the gap shrinks with X | — (as expected, no claim) |
| J.1 | [v23 §1] Lemma U | 20 values of x ∈ (0, 1] incl. X = 1, 1.25, 1.5; w ∈ [−60, 60] (Faddeev's integral; inversion relation for Re w > 0) | 1.4·10⁴ | max \|Re log Φ_b\| = 1.0·10⁻¹² | holds |
| J.2 | the contour turn: \|Φ_b\| ≤ C in the left wedge (used for the right wedge) | δ = 10° (as in v23), 45°; \|w\| ≤ 200 | 2·10⁵ | sup log\|Φ_b\| = 0 for δ = 10° (C = 1), 0.446 for δ = 45°; sup log\|1/Φ_b\| = 1.77 (δ = 10°), 8.97 (δ = 45°), both at the smallest x = 0.005: these constants depend on x (they grow as x → 0), which the proof (fixed x) allows | holds |
| J.3 | the turned integral equals J over L | a = 0.05 … 0.7, x = 1/40 … ½ | 4 | J over the 45° line vs over the 10° line vs over ℝ: ≤ 10⁻¹³ | holds |
| J.4 | [v23] \|P·J(a/b)\| ≤ 1/(2√2 sin πa) | 24 values a ∈ [10⁻³, 1 − 10⁻³] × 20 values x ∈ (0, 1] (x ≤ ¼ is the stated range) | 480 | max ratio **0.9655** (a = 0.284, x = 1/16); ratio → 0 linearly as a → 0, 1 (0.018 at a = 0.0054); the limit \|1/(1−2z)\| is approached slowly for small a | holds, sharp to 3.5 %; also for x ∈ (¼, 1] |
| T15.1 | [v23 §2] Theorem 15: K(a, X) = X·(F/prediction − 1) | 18 values a ∈ [10⁻³, 1 − 10⁻³], X = 50 … 6400 (exact sums, 30 + 0.08X digits) | 144 | K converges for every a once X·Vol/2π ≈ 7 … 32 (1 % criterion; not yet at X = 6400 for a = 0.001, 0.003, 0.997, 0.999); limit \|K\| ∈ [1.15, 1.6] (a = ½: 0.2015 − 1.5708i); v23's Richardson values reproduced (a = ¼: −0.2166 − 1.3453i; ⅓; 5/6; 1/20) | holds; O(1/X) constant not uniform in a (F1) |

### Priority 2

| id | where | domain tested | samples | worst ratio / value (where) | verdict |
|---|---|---|---|---|---|
| M.1 | [v24 §1] Proposition M (I(0), I(1/b), I(b), I(b + 1/b)) | random x, \|q\| ∈ [0.05, 0.95] (10 with \|q\| ≥ 0.9), Re x ∈ [−1.5, 1.5], arg b ≤ 70°; own quadrature on a line equivalent to L | 50 (1 skipped) | 2·10⁻¹⁵, 7·10⁻¹², 1·10⁻¹⁴, 4·10⁻¹² | holds |
| M.2 | Theorem 16, the four laws (own q-series) | same | 50 | 4·10⁻¹³ (f₁), 1.5·10⁻¹³ (ψ₁), 5·10⁻¹⁰ (f₀), 8·10⁻¹⁰ (ψ₀); with M₁, M₀ replaced by the Mordell integrals (no state integral): ≤ 1.4·10⁻¹⁴ | holds |
| M.3 | Theorem 17 (own Mordell integral on a turned ray) | same | 50 | M₁: 4·10⁻¹³, M₀: 5·10⁻¹⁰ | holds |
| M.4 | (R3), (R4) (verified, not proved, in v24) | exact power series to q¹⁰⁰⁰ modulo 2⁶¹ − 1 and 10⁹ + 7 | 2 | all coefficients zero | holds (to q¹⁰⁰⁰) |
| M.5 | Theorem 17 on the unit circle | x = 1/X, X = 1 … 20, 25, 30, …, 100 | 29 | ≤ 1.2·10⁻¹¹ (growing slowly with X) | holds |
| M.6 | the continuity Theorem 17 needs on the circle ("not written out") | x = 1/X + iε, X = 3, 5, 8, 20; ε = 10⁻², 10⁻³, 10⁻⁴ | 24 | \|J(x + iε) − J(x)\|/\|J\| ∝ ε (1.4·10⁻⁴ at X = 3, 5.5·10⁻³ at X = 20, ε = 10⁻⁴) | holds |
| R3′.1 | [v25 §1] Lemma F for real x ∈ [¼, 1] and \|x/x₀ − 1\| ≤ 0.1 | 9 values x₀, 8 complex points around each; 2440 w each | 2·10⁵ | sup 0.3705 at x = 1 (real), **0.383** at x = 1 − 0.1i | holds, C₁ = 0.39 |
| R3′.2 | Lemma D for \|t\| ≤ 2π·1.01 | \|t\| = 0.5 … 2π·1.01, \|β\| ≤ 30°, left half of the v16 contour | 4.7·10⁵ | \|η\|/e^{−0.37ρ} ≤ 0.956; Σ\|log(1 + η)\| ≤ 0.170 (claim 1.2); number of factors ≤ 0.21ρ | holds |
| R3′.3 | Theorems 9, 11, 14 at X = 1, 2, 3, 25 digits (own mp Faddeev) | ψ₁\*, ψ₀\*; Theorem 14 at X = 2 for a = ¼, 0.7, 0.05 | 12 | X = 1: 2·10⁻¹⁷, 1·10⁻¹⁵; X = 2: 9·10⁻¹⁷, 7·10⁻¹⁵, Theorem 14: 2·10⁻¹⁶ … 5·10⁻¹⁵; X = 3: 2·10⁻¹⁴, 1·10⁻¹² (quadrature limit: 6·10⁻¹⁴ between two resolutions); v25's values at X = 1 reproduced | holds; more than double only at X = 1, 2 |

### Priority 3

| id | where | domain tested | samples | worst ratio / value (where) | verdict |
|---|---|---|---|---|---|
| P3.1 | [v12 §4.1] closing lemma: sup of the exponent on Re ξ = σ_N against −c·σ_N² | 40 random (x, δ), A = 0, ½; N = 2 … 8 (and 2 … 32 for three cases) | 560 | f_N = sup/(−cσ_N²) = 0.982 … 0.999 at N = 8, → 1 from below (0.9998 at N = 32); the remainder (sup + cσ_N²)/(σ_N + 1) ≤ 0.141 | holds: c is the exact leading constant; the remark "times a factor between 1 and 1.5" is slightly off (F5) |
| P3.2 | Lemma L, abstract: a model with entire V, analytic g, Γ = ℝ | arg t = 0, 30, 60, 80°; \|t\| = 0.2 … 0.01 | 20 | (F/lead − 1)/\|t\| ≤ 0.95 where (L1) holds; at arg t = 80°, where (L1) fails, the formula fails too (10¹⁶) | holds |
| P3.3 | Theorem 8 / Theorem IV (ψ₁\*) as a whole | arg x = 10, 5, 2, 1, 0.3°; \|t\| = 0.4 … 0.03 (\|p\| up to 0.77); ψ₁\* summed at 40 digits; Borel sum by series + Padé | 25 | R/\|t\| = 0.25 … 0.28, flat, also at \|p\| = 0.77 | holds, uniform in the sector |
| P3.4 | Theorem C at 30 digits outside the round-1 sample | \|q\| = 0.06, 0.08, 0.92; arg x = 158°; δ = 1°, 2°, 3°, and δ next to 90° − β | 6 | ≤ 8·10⁻³⁰ | holds |
| P3.5 | [v26 §3] Lemma V | 1999 values of a | 1999 | g = 1 has one solution; λ₋ ∈ [0.500, 0.693] ⊂ (½, 2^{−1/2}); g(λ₋) ≤ 0.083 < 0.16; concavity threshold a = 0.808 (v26: 0.807674) | holds |
| P3.6 | [v26 §2] growth rate of the e^{w}-version on Γ₃ at real x | X = 2, 3, 4, 6 | 4 | slope of log\|integrand\| = 1 − π/\|t\| to 4 decimals | holds |
| P3.7 | v26 against the round-1 findings | reading | — | F2: v26 §2 states |t| < π/3 (decaying) and |t| < π (growing), as asked. F3: v26 §1 writes the argument in y on the fixed line L; the zeros of Φ_b and the kernel poles avoid L for \|β\| ≤ 60°, as stated. F4: Lemma V proves M = 0 for all a, so Theorem 7′ has ε(a) = 60°. The smaller points (Lemma F remark, p^{n²/2}, Lemma B (iii), (L1) on Γ up to 17.15°, Lemma R2's bracket, the constants without slack) are all recorded in §4 | done |

## Findings

None is a violation. Ordered by how much depends on them.

**F1 (Theorem 15, v23 §2: the O(1/X) term is not uniform in a).** For every a the limit K(a) = lim X·(F/prediction − 1) exists. It is bounded on all of (0, 1): |K| between 1.15 and 1.6, tending to about −1.152i at both ends. The approach to that limit, however, starts only once X·Vol(a)/2π is about 10:
- a = 0.01 needs X ≈ 1600;
- a = 10⁻³ has not settled at X = 6400, and there |K| reaches 650 at X = 800.

So the error in "F = prediction·(1 + O(1/X))" is bounded by a constant times 1/X only for X ≳ 60/Vol(a). Before that, the factor e^{Vol·X/2π} is not yet large enough to separate the growing term from the rest. The theorem states that the constants depend on a and degenerate as z → 1, which is correct. A uniform statement would need X·Vol(a) as the parameter, not X. Reproduction: `t15.py`.

**F2 (Theorem 13 (b), v22 §3: ε(a) measured).** The proof's choice (3π²/2)·tan φ₁ < Vol/2 is necessary in spirit: with a tilt of 20°, (L1) at τ = 0 fails for a ≤ 0.1. With the best tilt, ε(a) ranges from 20° (a ≈ 0.3) down to 0.08° (a = 10⁻³) and 0.02° (a = 1 − 10⁻³); the table is under T13b.4. v26 §3 notes that Theorem 13 keeps its unspecified ε(a); these numbers could replace it by "numerically ε(a) ≈ …", or by ε(a) ≥ c·Vol(a) with c ≈ 20°. Near a = 1 the limit is set by C₄/C₅, not by the tilt: for a ≥ 0.9 every φ₁ gives the same ε. Rotating C₄, as v22 §7.3 suggests, is the place to gain there.

**F3 (Lemma J, v23 §1: the constant in the contour turn depends on x).** The proof uses 1/|Φ_b| ≤ C in the left wedge, and |Φ_b| ≤ C there for the right wedge, quoting Lemmas F and D "as in Lemma R3". For fixed x that is right: C = 1 for |Φ_b| and C = e^{1.77} for 1/|Φ_b| with δ = 10°, at x = 0.005. But the constant for 1/|Φ_b| grows as x → 0, because 1/|Φ_b| ≈ e^{−Im Li₂(−e^w)/|t|} is large where Im Li₂ < 0. This does no harm, since Lemma J is stated at fixed x and the bound it proves is uniform. It means the turn is pointwise in x, not uniform; a sentence would make that explicit.

**F4 (Theorem 9, 11, 14 at X = 3: precision).** My 25-digit check reached 2·10⁻¹⁷ (X = 1) and 9·10⁻¹⁷ (X = 2) for ψ₁\*, but only 2·10⁻¹⁴ for ψ₁\* and 10⁻¹² for ψ₀\* at X = 3. That is limited by my outer quadrature (two resolutions differ by 6·10⁻¹⁴), not by the formula. So "more than double precision" was achieved at X = 1, 2 only. Nothing suggested a discrepancy.

**F5 (closing lemma, v12 §4.1: the remark on the factor).** The constant c = π(1 − A)·min(|x|²/Im x, |x| sin 2δ/cos²(δ − β)) is the exact leading constant: sup E/(−c·σ_N²) → 1 as N → ∞. In all samples the minimum is attained by the second value, the point on C. At finite N the ratio is slightly below 1 (0.94 at N = 2, 0.982–0.999 at N = 8), because of the lower-order terms e^{O(σ_N)}, which the proof allows; they are at most 0.141·(σ_N + 1). The remark "the supremum on the line is this bound times a factor between 1 and 1.5" is therefore not right as stated: the factor is below 1 at finite N and tends to 1. The proof is not affected.

**Smaller points.**
- **T13a.3:** the bound for 1/Φ_b on C needs exactly the factor 2 of the functional equation, attained at the corner 2πi (where p·w′ = 1 at even X). So the "at most 2" of v22 §2 (L2) is sharp.
- **T13a.4:** the uphill margin along the horizontal line goes to 0 as a → 1 (4.8·10⁻⁵ at a = 0.9999), as does everything else there.
- **T13a.5:** Theorem 13 (a) has a flat O(t) constant ≤ 0.13 for a ∈ [0.25, 0.98] once X·Vol/2π ≳ 2. For a = 0.02 the constant is still 2.8 at X = 128.
- **R3′.1:** Lemma F's constant grows slowly with x: 0.348 at x = ¼, 0.3705 at x = 1, 0.383 in the neighbourhood of 1.
- **T14.1:** the one outlier (3·10⁻⁵ at X = 36, a = 0.976) is a quadrature fallback with conditioning 10¹³, and neighbouring a at the same X give 3–5·10⁻¹²; it is a quadrature artefact, not a failure of the theorem (a dedicated re-run at higher precision was started but not finished).

## Not tested

- **Theorem 13 at a → 0, 1 with |t| small enough to be in the asymptotic regime** (X·Vol/2π ≳ 8 means X ≳ 2·10⁴ at a = 10⁻³). The (L1) inequalities were tested there; the integrals were not.
- **Theorem 14 for a within about 0.03 of 0 or 1 at large X**: the descent contours become too long for my quadrature, and the 45° line too ill-conditioned.
- **The proofs quoted from the literature in v24** ((Q1), (S1)–(S6), (R1), (R2), (R5), (Q2)). Only their consequences were tested, and (R3), (R4) as series.
- **Theorem 17 inside the framework**: not claimed in v24 either.
- **Theorem 8 for f₁\*, ψ₀\*, f₀\***: only ψ₁\*.
- **More than double precision at X = 3** (F4).
