# Mock theta search: addendum v25 (the formulas at roots of unity hold for X = 1, 2, 3 as well)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). A small item left open in `addendum-v21.md` §5: the theorems at e^{2πi/X} excluded X ≤ 3 only because of the range of one lemma. Scripts and output are in `psi_star/v25/`.

Labels: **proved** (a complete argument is written here or in the addenda cited; nobody else has checked it), **verified** (numerics).

---

## 0. Summary

1. **Lemma R3′ (proved; §1).** The continuity lemma of v20 (Lemma R3) holds at every real 0 < x₀ ≤ 1, not only for x₀ ≤ ¼, and for every u.
2. **Theorem 18 (proved; §1).** Theorems 9 (v20), 11 (v21) and 14 (v22) hold for every X ≥ 1 of the right parity. The new cases:
   - X = 1 (b = 1, q₀ = 1): lim_{r→1⁻} ψ₁\*(r) = ⅓ = P·[J(½) − ⅓·I(½) − ⅓·I(−½)] with P = 2^{−1/2}·e^{−7πi/24}, and lim_{r→1⁻} ψ₀\*(r) = ⅓ = −P·[J(3/2) − ⅓·I(3/2) − ⅓·I(½)].
   - X = 2 (b = 2^{−1/2}, q₀ = −1): ψ₁\*(−1) = 1 = P·[J(1/(2b)) + I(1/(2b)) − I(−1/(2b))], ψ₀\*(−1) = −1, and for every z = e^{2πia}, 0 < a < 1: 1 = i·P·e^{−πia}·[J(a/b) − I((a−1)/b) + z²·I(a/b)].
   - X = 3 (b = 3^{−1/2}, q₀ = ω): lim_{r→1⁻} ψ₁\*(rω) = 1 = P·[J − ⅓·I(1/(2b)) − ⅓·I(−1/(2b))] and lim_{r→1⁻} ψ₀\*(rω) = 0, that is J(b + 1/(2b)) = ⅓·[I(b + 1/(2b)) + I(b − 1/(2b))].
3. **Numerically (verified):** the formulas at X = 1, 2, 3 hold to 10⁻¹¹ or better in double precision; the radial limits are approached like ε.

---

## 1. Continuity at real x up to x = 1

**Lemma R3′.** Let 0 < x₀ ≤ 1 and let u(x) be continuous and bounded near x₀. The integrals I(u), J(u) of v12 over L = e^{−10°·i}·ℝ, with Faddeev's Φ_b, are continuous in x on {Im x ≥ 0} near x₀.

Proof. The proof of Lemma R3 (v20 §1) needs, for x in a neighbourhood N of x₀ in {Im x ≥ 0}, a constant C with |1/Φ_b| ≤ C on the left half of L and |1/Φ_b| ≤ C·|e^{−w²/(2t)}| on the right half. It took this from Lemmas F and D of v16, which were stated for |x| ≤ ¼. Their proofs use the size of x as follows.

- *Lemma F* (|Λ(w) + Li₂(−e^{w})/t| ≤ C₁·e^{Re w/2} for Re w ≤ 0, |Im w| ≤ π). The contour of Faddeev's integral is moved from ℝ + i0 to Im σ = ½. (i) No pole is crossed if the poles σ = i·n/x of r₁(πxσ) lie above the line: Im(i/x) = Re(1/x) > ½. (ii) On the line, u = πxσ must stay away from the poles iπn (n ≠ 0) of 1/sinh u: for real x the distance is at least π·(1 − x/2). For x ∈ N both hold, with Re(1/x) ≥ 0.9 and distance ≥ 1.5 if N is small. The other steps (|r₁(u)| ≤ C·min(|u|, 1/|u|) on the line, and the integral of min(|x|, 1/(|x|·|σ|²))) do not use the size of x. So Lemma F holds on N with a constant C₁(N).
- *Lemma D* (outside the strip). The factors of the products are 1 + η with |η| ≤ exp(−(2πρ/|t|)·0.383), there are at most ρ of them, and ρ ≥ 4.88. For |t| ≤ 2π·1.01 this gives |η| ≤ e^{−0.37ρ} ≤ 0.17 and |E| ≤ Σ |log(1 + η)| ≤ 1.2·ρ·e^{−0.37ρ} ≤ 1.2. So E is bounded (no longer exponentially small, which is not needed here).

With these two bounds and the inversion relation the proof of Lemma R3 applies word for word; the factor e^{2πuy} is e^{O(ρ)} for bounded u. ∎

**Theorem 18.** (a) Theorem 9 of v20 holds for X = 2. (b) Theorem 11 of v21 holds for X = 1 and X = 3. (c) Theorem 14 of v22 holds for X = 2.

Proof. In the three proofs the restriction on X entered only through Lemma R3 (x₀ = 1/X ≤ ¼); Lemma R3′ removes it. The other ingredients have no restriction: Lemma R1 of v20 is stated for X = 2M with any M ≥ 1; Lemmas R2 (v20) and R2′ (v21) are statements about p → ±1 along the path x = 1/X + iε, for any X; the odd-root rule of v3 is proved for every odd m, and its proof covers m = 1 (ζ = 1: (−1;1)ₙ = 2ⁿ, the terms of f₁\* satisfy term_{n+1} = −½·term_n, and θ₄(e^{−t}) → 0). The values: for X = 1, A = A₀ = 1; for X = 2 the sums stop at n = 0 because 1 + q₀ = 0; for X = 3, A = 1 − ω(1 + ω) + (1 + ω)(1 + ω²) = 3 and A₀ = ω − (1 + ω) + 1 = 0. P = 2^{−1/2}·e^{π²/(12t) + t/6}/b at t = −2πi is 2^{−1/2}·e^{−7πi/24}. ∎

Remarks.

- b = 1 is the self-dual point of Φ_b. There the three integrals are taken with Φ₁, and the identity 1 = P·[3·J(½) − I(½) − I(−½)] says that the divergent series ψ₁\*(1) = Σ (−2)ⁿ, summed as 1/(1 + 2) = ⅓ by the radial limit, is given by the formula of v12 at its most symmetric point.
- At X = 3 the radial limit of ψ₀\* vanishes, so there the kernel integral is exactly one third of the sum of the two state integrals.
- Lemma R3′ also removes the remark on the edge case X = 4 in v20.

---

## 2. Numerical checks (verified; `r1_small.py`, `run_log.txt`)

| check | result |
|---|---|
| Theorems 9, 11 at X = 1, 2, 3 (and 4, 5 for comparison), ψ₁\*: formula minus the exact value | 3·10⁻¹⁵, 1·10⁻¹⁴, 1·10⁻¹³ (8·10⁻¹³, 1·10⁻⁹) |
| the same for ψ₀\* | 4·10⁻¹², 1·10⁻¹², 6·10⁻¹² (1·10⁻¹¹, 2·10⁻¹¹) |
| Theorem 14 at X = 2, a = ¼, ⅓, ½, 0.7, 0.9: formula minus 1 | 1·10⁻¹¹, 8·10⁻¹³, 4·10⁻¹⁵, 8·10⁻¹⁵, 4·10⁻¹⁴ |
| radial limits at q = e^{−2πε}, ε = 10⁻², 10⁻³, 10⁻⁴: \|ψ₁\* − ⅓\|; \|ψ₀\* − ⅓\| | 7.1·10⁻³, 7.0·10⁻⁴, 7.0·10⁻⁵; 5·10⁻⁵, 5·10⁻⁷, 5·10⁻⁹ |
| radial limits at q = ω·e^{−2πε}, the same ε: \|ψ₁\* − 1\|; \|ψ₀\*\| | 8.3·10⁻², 6.4·10⁻³, 6.3·10⁻⁴; 9·10⁻³, 4·10⁻⁵, 4·10⁻⁷ |

Values at X = 1: P·J(½) = −0.005723 − 0.239909i, P·I(½) = −0.300817 − 0.306142i, P·I(−½) = −0.716352 − 0.413586i.

---

## 3. Open

1. Roots of unity e^{2πih/k} with h ≠ 1 (unchanged).
2. Whether the three integrals at b = 1 have closed forms (Φ₁ can be written with the ordinary dilogarithm; from memory, not used here).

---

## 4. Scripts (`psi_star/v25/`)

- `r1_small.py`: the checks of §2 (25 seconds).
- `h1_roots.py` (from v20), `k1_odd.py` (from v21), `n1_zroots.py` (from v22), `fad.py` (from v16), `run_log.txt`, `README.md`.
