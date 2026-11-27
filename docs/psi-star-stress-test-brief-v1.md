# Brief: stress test of the proofs for ψ₁\*, ψ₀\* (v1)

27 Nov 2026. For a fresh session with access to this repository. Written by Claude Fable 5.2 at the request of Mario Klingemann.

## 1. Goal

`docs/psi-star-theorems-v1.md` collects the theorems on ψ₁\*, ψ₀\* and their proofs (Parts 1–7 and an appendix; §0 states the results and §0.5 lists the weak points). Nobody has checked the analytic estimates. Most lemmas are **inequalities or exact identities with explicit domains**, so they can be attacked numerically.

Your task: try to **break** every displayed inequality, identity and "uniform constant" claim in that document, with code of your own. A claim that survives a serious attack is not proved by it, but a claim that fails is a real finding. Finding one error is worth more than confirming fifty claims.

This is not a proofreading task and not a request for new mathematics. If you notice a logical gap while working (a step that does not follow, a hypothesis used but not stated), write it down as a finding too.

## 2. Rules

- **Independent code.** Write your own implementations of Φ_b (products, and Faddeev's integral for real b), the integrals, the q-series. Use the scripts in `psi_star/v12` … `v21` only to compare afterwards. The existing `check_*.py` scripts tested the final statements at a few points; you test the intermediate steps on their whole stated domains.
- **Adversarial, not confirmatory.** For each inequality LHS ≤ RHS, maximise LHS/RHS over the stated domain: random sampling with fixed seeds, then local optimisation from the worst samples, and always the edges and corners of the domain (the proofs are most likely to fail there). Report the worst ratio found and where.
- **Stated domain first, then beyond.** Test exactly the domain the lemma states. Then push past it and report where the claim starts to fail; a claim that fails right at the edge of its domain has no margin, and that is worth knowing.
- **Uniformity claims.** Where a proof says "uniformly", "for |t| ≤ t₀" or "for small δ" without numbers, estimate the constant on a grid and report (a) whether it stays bounded as the small parameter goes to 0 and as the direction approaches the edges of the sector, (b) an explicit threshold that works.
- **Before reporting a violation:** reproduce it at higher precision and with a second method (for example a different contour or the product formula instead of the integral). State the minimal counterexample.
- **Do not edit** the addenda or the consolidated document. Report.
- Labels: **holds** (worst ratio < 1, say how close), **holds with a different constant** (give it), **violated** (counterexample), **could not test** (say why).

## 3. Numerical pitfalls already met

- State integrals at real x or near it must be integrated over the descent contours of v19 (Part 6); over the line L they cancel by many orders of magnitude (v20 §4, notes).
- Finite sums at q = e^{2πi/X}: terms reach e^{0.1615·X}; use about 0.07·X extra digits (v18 §2.2).
- W, W₃ near p = −1 have alternating terms of size e^{0.3/δ}, δ = 1 − |p| (v21 §4).
- q-series near |q| = 1: 1/(q;q)ₖ reaches e^{1.645/|t|}; set the precision accordingly.
- mpmath `qp` needs `maxterms` raised for |q| close to 1. For |q| > 1 use Φ_{1/b} (Part 2, Lemma Φ).
- Faddeev's integral in the variable w converges for |Im w| < π(1 + Re x); outside use the functional equation and the inversion relation. A double-precision implementation is in `psi_star/v16/fad.py` (for comparison only).

## 4. What to test, in order of priority

Section references are to the headings of `docs/psi-star-theorems-v1.md` ("[v15 §1.2]" etc.).

**Priority 1: everything else rests on these.**

- **S3 (Part 3, Proposition 1).** The two exact formulas for ψ/ρ; the bounds |Y| ≤ π/2 and λ_α ≤ X ≤ ln 2; the formula for κ; and κ ≥ 0.185 for δ = 10°, −210° ≤ θ ≤ 10°. Test ψ(ρ) ≤ −κ·ρ for ρ from 10⁻³ to 100 over the whole range of θ, densely near the ends and near θ = −90°.
- **S4 (Part 4, uniform bound).** (a) The Euler–Maclaurin lemma [v14 §2.2] and Lemma 1′: the remainder bound C(K, η)·|t₁|^{K−1}/cos(arg t₁), for random ζ₀, t₁ meeting the hypotheses, as |t₁| → 0 and as arg t₁ → ±90°. (b) Lemma A (distances to the points iπ(2m+1)). (c) Lemma 3′ and [v16 §4]: sup over the whole contour of |H|·e^{κ|s|} and of |H|/e^{ψ/|t|}, for θ over the closed sector [−210°, 10°] and |t| from 0.4 down to 0.01. Is it bounded as |t| → 0? What t₀ and C does one actually get? Do this for the ψ₁\* and the ψ₀\* integrand.
- **S5 (Part 5).** The shift identity for random N and t (exact). Lemma G for random ω satisfying the distance hypothesis, and the constant c₃. Corollary G on the contour for θ in the closed sector and N up to 0.04/|t|. The margin 0.175 − 0.149 in Corollary R: with your constants from S4, is the integral really bounded? Lemma T on the full disc |t| ≤ 0.3/n. Theorem 4: the ratio |P·J − Σ_{k<K} aₖtᵏ|/(K!·(240|t|)ᴷ) for K up to 40 at many directions including both ends of the sector. Theorem 5 at θ₀ = −120° and −80° with t near the edge of the disc.

**Priority 2.**

- **S2 (Part 2).** Lemma Φ (integral against both products; functional equations and inversion at real b). Lemma S: defect/|t|^{2K−1} for K = 2, 3 over |Im w| ≤ π − ε₁, any Re w, |arg x| ≤ 60°, |x| ≤ ¼. Lemma F: the supremum of |Λ + Li₂(−e^{w})/t|·e^{−Re w/2} on Re w ≤ 0, |Im w| ≤ π, including the edges and x at the corners of its range. Lemma D: the margin 0.383 and the actual size of E. Lemma B (i)–(iii) [v19 §2]: the constants, for 0 ≤ arg x ≤ 60°.
- **S6 (Part 6).** The closed formulas for Re V, Im V on the three pieces of the growing contour, and (L1) there for 0 ≤ τ ≤ 0.2 (then beyond 0.2: where does it fail?). The monotonicity claims on the horizontal line for random a ∈ (0, 1), including a near 0 and 1. Theorem 6: exact/asymptotic − 1 = O(t) uniformly for arg x in [0°, 11°] including arg x = 0 (real x) and the edge. Theorem 7 for random a: find ε(a). Lemma L is only sketched: test its conclusion on the actual integrals.
- **S7 (Part 7).** Lemma R1: the block bound |B_m| ≤ min(1 + Cδ, (m + 1)Cδ) for random X, m, δ. Lemma R2: each inequality of the proof, and the range of δ for which the bracket with K = 3/2 is positive. Lemma R2′: the lower bounds 9/8 and e^{−0.3/δ}, the bound for G(Q), the bound for ψ(p), all along the path and slightly off it. Lemma R3. Theorems 9 and 11 at more X (X = 4 is on the edge of the range of Part 2; X = 3 is excluded: does the formula hold there anyway?), with more than double precision if you can.

**Priority 3: already checked at several points, extend to random parameters.**

- **S1 (Part 1).** Theorems A, B, C at random complex u (|u| ≤ 3), random x with |q| between 0.1 and 0.9, random admissible δ. The Fourier transforms F(w), G(w) on their half-plane and steps (X), (Y). The two bounds of the closing lemma [v12 §4.1] with the stated constant c.
- **SA (Appendix).** Identities (I), (II), (V), (VI) as power series. The steps (c), (d) of the odd-root rule: the lower bound for |(−q;q)_{mJ}|, and θ₄·R₁ → 0.

## 5. Deliverables

- `psi_star/stress/`: one script per block S1 … S7, SA (fixed seeds, runnable from that folder, stated run times), a `README.md`, and the complete output as `run_log.txt`.
- `docs/psi-star-stress-test-v1.md`: one table row per claim tested: id, where it is in the consolidated document, domain tested, number of samples, worst ratio and where it occurred, verdict. Then a section **Findings** (violations, missing hypotheses, constants that are wrong or unstated, logical gaps), each with a minimal reproduction, ordered by how much of the chain depends on it. Then a section **Not tested** with reasons.
- If a claim is violated, say what fails downstream according to the map in §0.3 of the consolidated document, and stop to report before going on to lower priorities.

## 6. Context

- State of the whole project: `docs/psi-star-handoff-v21.md`. README items 5–7 summarise the results.
- The consolidated document is generated: `psi_star/theorems/assemble.py`.
- Known soft spots are listed in §0.5 of the consolidated document; start there within each block.
