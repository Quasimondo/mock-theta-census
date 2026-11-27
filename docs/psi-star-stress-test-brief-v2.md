# Brief: stress test of the proofs, second round (v2)

27 Nov 2026. For a fresh session with access to this repository. Written by Claude Opus 5.5 at the request of Mario Klingemann, after the first round.

## 1. Where things stand

- **The first round is done.** Brief: `docs/psi-star-stress-test-brief-v1.md`. Report: `docs/psi-star-stress-test-v1.md`. Scripts: `psi_star/stress/`. It tested `docs/psi-star-theorems-v1.md` as it was at commit d2ba3c2, i.e. addenda v2–v21. **No claim was violated.** It measured the implicit constants (F1) and reported four other findings, F2–F5, plus some smaller points. Read the report first: its tables, numerical notes and "Not tested" list are the starting point for this round.
- **Status of the findings of round 1** (as of this brief):
  - F1: constants measured. Nothing to fix.
  - F2: a missing hypothesis in the contour move for the e^{w} version of the decaying integral (v19 §4). It needs (2η₁ + 2πa) + |t| < 0, i.e. |t| < π/3 at a = ½. Not yet fixed in the addenda.
  - F3: the analyticity argument in the Corollary of v16. The contour angle moves with t, and independence of the admissible contour is stated only for Im x > 0. Not yet fixed.
  - F4: ε(a) in Theorem 7 can be taken as 60°. Not yet written into the addenda.
  - F5: X = 2, 3 work. **Resolved by v25** (Lemma R3′ and Theorem 18), which extends the continuity lemma to real x ≤ 1. Round 1 warned that Lemma F's constant must then be checked, not assumed (S2.4).
  - Smaller points of round 1: not yet fixed.
- **New since round 1: addenda v22–v25.** They are untested by round 1.
  - v22: the growing state integral of the z-family for every z, on a five-piece contour (Theorem 13), and F(z) at even roots of unity (Theorem 14).
  - v23: |Φ_b| = 1 on the real line (Lemma U), the bound |P·J(a/b)| ≤ 1/(2√2·sin πa) at real x (Lemma J), and the volume rule at even roots of unity as a theorem (Theorem 15).
  - v24: z = 1 and z = q. The state integrals become S-matrix entries (Proposition M), the formula of the family becomes the fifth-order transformation laws (Theorem 16), and the kernel integral is a Mordell integral (Theorem 17).
  - v25: the root-of-unity theorems for X = 1, 2, 3, via the wider continuity lemma R3′.
- **Only v25 §1 is in the consolidated document** (Part 7). v22–v24 are about the z-family and live only in their addenda. Cite them as "[v22 §3]" and so on, as in round 1.
- **The independent checks so far** are the `check_*.py` scripts in `psi_star/v22` … `v25` and the notes at the top of each addendum. They test final statements at a few points. This round tests the intermediate steps on their full stated domains.

## 2. Rules

The same as round 1 (`docs/psi-star-stress-test-brief-v1.md` §2): own code; adversarial (maximise LHS/RHS, edges and corners first); the stated domain first, then beyond; explicit constants for every "uniformly" or "for small"; reproduce a violation with a second method and higher precision before reporting it; do not edit the addenda or the consolidated document; the same labels. Work on a branch `stress-test-v2`.

## 3. Numerical pitfalls (round 1 and the checks of v22–v25)

- All the pitfalls of round 1 §3 and the *Numerical notes* in `psi_star/stress/README.md`.
- **Small X (large x).** At real x ≥ ⅙ the factor e^{2πuy} can outgrow the Gaussian along L = e^{−10°i}·ℝ inside any practical truncation. Integrate on a steeper line, e.g. e^{−45°i}·ℝ. This is equivalent to L because no pole of the kernel or zero of Φ_b lies between them. Reach Φ_b off the strip by the functional equation Φ_b(Y) = Φ_b(Y − ib)/(1 + e^{2πb(Y − ib/2)}) (used in `psi_star/v25/check_small_X.py`).
- **F2:** horizontal-line contours for the e^{w} versions are admissible only for |t| < π/3 at a = ½. In general check the sign of the decay rate on the right before truncating.
- **v22 contour C:** Li₂(−e^{w}) must be continued exactly as v22 §1 states. On Im w = 2π and Im w = π take the boundary values from inside, and for Re w ≤ 0 use the 2πi-periodic power series. The wrong branch gives spurious violations of (L1).
- **Thin margins:** (L1) on contour C has margins down to −0.0036 (a = 0.99) and −0.0095 (a = 0.01) at τ = 0 (`psi_star/v22/check_thm13.py`). Resolve the grid finely near the saddle and near the corners 2πi, iπ and r₀ + iπ.
- **The notation p^{n²/2}** in Theorem C means e^{−πi·n²/(2x)}, not a principal power (round 1, smaller points).

## 4. What to test, in order of priority

**Priority 1: the new proofs that the volume rule rests on.**

- **T13a [v22, Theorem 13 (a)].**
  - (L1) at real x (τ = 0) on all five pieces C₁ … C₅ for a ∈ (0, 1), densely as a → 0 and a → 1, where the margins vanish together with Vol(z). Report the worst ratio as a function of a, and whether it stays negative down to a = 10⁻³ and up to 1 − 10⁻³.
  - The bound for 1/Φ_b that each piece needs (which lemma, and on which domain). In particular C₁ and C₂ lie outside the strip |Im w| ≤ π, and C₁ runs to Re w → −∞ on Im w = 2π.
  - The claim of v22 §5 that the horizontal ray goes uphill. Two checks have confirmed it at a = ½ and 0.2; test it for all a.
- **T13b [v22, Theorem 13 (b)].** The small tilted sector 0 ≤ arg x ≤ ε(a). Find ε(a) numerically, as round 1 did for Theorem 7 (F4).
- **T14 [v22, Theorem 14].** F(z) at even X equals the three-integral formula. Use random a ∈ (0, 1) and X = 2, 4, 6, …, 60, with own Faddeev code. Also test the limits of the dual series on the path x = 1/X + iε, and odd X if the formula suggests a form there.
- **J [v23, Lemmas U, J].**
  - |Φ_b| = 1 on ℝ for real b ∈ (0, 1].
  - The contour turn of Lemma J. In the right wedge the proof uses |Φ_b(−y)| ≤ C, an *upper* bound for |Φ_b| in the left wedge. Test that bound and its constant.
  - |P·J(a/b)| ≤ 1/(2√2·sin πa) for all a ∈ (0, 1) and all real x ∈ (0, 1]. Report the sharpness of the bound as a → 0, 1.
- **T15 [v23, Theorem 15].** The volume rule at even roots of unity, with its error term. How does the O(1/X) constant behave as a → 0 and a → 1? Is the theorem's uniformity in a, if stated, numerically true?

**Priority 2: z = 1 and the small-X extension.**

- **M [v24].**
  - Proposition M and Theorems 16 and 17 at random x with |q| from 0.05 to 0.95, including near the real x-axis and |x| large. Use own quadrature, own Mordell integral and own q-series.
  - The series identities (R3) and (R4), verified but not proved in v24, to higher order (q¹⁰⁰⁰).
  - Theorem 17 on the unit circle (x = 1/X, X = 1 … 100). The continuity it needs is "not written out" in v24.
- **R3′ [v25 §1].**
  - Lemma F at real x ∈ (¼, 1] and in a small complex neighbourhood: sup of |Λ + Li₂(−e^{w})/t|·e^{−Re w/2}. Round 1 found 0.511 at the corner |x| = ¼, arg x = ±60°.
  - Lemma D's claim |E| ≤ 1.2 for |t| ≤ 2π·1.01.
  - Theorems 9, 11 and 14 at X = 1, 2, 3, with more than double precision.

**Priority 3: what round 1 did not test.**

- The closing lemma of [v12 §4.1], with the stated constant c.
- Lemma L as an abstract statement.
- Theorem 8 / Theorem IV as a whole: ψ₁\* at |q| close to 1 in the sector 0 < arg x ≤ 10°, against Borel sum + dual terms.
- Theorems A, B and C outside the sampled set (|q| < 0.1 or > 0.9, arg x > 120°, δ near 0 or near 90° − β), at multiple precision.
- Re-test F2–F4 and the smaller points of round 1 if the addenda have been corrected by the time you start: check that each fix says what the report asked for.

## 5. Deliverables

- `psi_star/stress2/`: one script per block (T13a, T13b, T14, J, T15, M, R3′, P3), fixed seeds, `README.md`, the complete output as `run_log.txt`.
- `docs/psi-star-stress-test-v2.md`: the same format as round 1 (summary, tables with id / where / domain / samples / worst ratio / verdict, findings ordered by how much depends on them, not tested).
- Stop and report before lower priorities if a claim of Priority 1 is violated.

## 6. Context

- README items 5–7 and the "Status" paragraph summarise the results. `docs/psi-star-handoff-v21.md` is the project state for a fresh session; the file name is older than its content.
- Addenda v22–v25 each start with an independent verification note. Their scripts are in `psi_star/v22` … `v25`.
