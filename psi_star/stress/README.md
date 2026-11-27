# Stress test of the proofs for ψ₁\*, ψ₀\* (v1)

Scripts for the brief `docs/psi-star-stress-test-brief-v1.md`. The report is
`docs/psi-star-stress-test-v1.md`; the complete output of the last run is `run_log.txt`.

Everything here is independent of `psi_star/v*`: Φ_b by the products (both sides of the real axis,
Φ_{1/b} for Im x < 0) and by Faddeev's integral (trapezoid rule on Im σ = ½, reduced to
Re w ≤ 0 by the inversion relation and to |Im w| ≤ π by the functional equation), the dilogarithm
from scipy, the q-series directly, the coefficients aₖ by exact power sums. Existing material was
used only afterwards, for comparison: `psi_star/v17/phi_exact.json` (coefficients) and the values
quoted in the consolidated document.

| file | block | what | run time (2 cores, in parallel with S6) |
|---|---|---|---|
| `core.py` | – | Λ(w) = log Φ_b by products and by Faddeev's integral; contour, G, ψ(ρ), log H | – |
| `s1_kernel.py` | S1 | Theorems A, B, C, Fourier transforms F, G, steps (X), (Y), Theorem I, at random u, x, δ | ~5 min |
| `s2_faddeev.py` | S2 | Lemma Φ, S, F, D, B (i)–(iii); accuracy of `core` against mpmath | ~3 min |
| `s3_descent.py` | S3 | Proposition 1 (i)–(iv), and beyond the sector | ~2 min |
| `s4_uniform.py` | S4 | Lemma 1′, Lemma A, sup of \|H\|e^{κ\|s\|} and \|H\|/e^{ψ/\|t\|} | ~4 min |
| `s5_factorial.py` | S5 | Proposition 2, Lemma G, Corollaries G and R, Lemma T, coefficients, Theorems 4, 5 | ~12 min |
| `s6_state.py` | S6 | heights on Γ₁–Γ₃, (L1), monotonicity, Theorems 6, 7, the forms of §0.2 | ~40 min |
| `s7_roots.py` | S7 | Lemmas R1, R2, R2′, R3; Theorems 9, 11 (and 10, 12) for X = 2 … 41 | ~10 min |
| `sa_appendix.py` | SA | identities (I), (II), (V), (VI) exactly to q⁴⁰⁰; steps (c), (d) of the odd-root rule; the rule itself | ~10 s |
| `run_all.sh`, `merge_log.py` | – | runs all of them (output of each in `out_<block>.txt`) and assembles `run_log.txt` | ~40 min |

Run from this folder: `./run_all.sh`, or `python3 sX_….py` for one block. Fixed seeds
(`numpy.random.default_rng(k)` with k the block number; S1 uses 11). Requirements: numpy, scipy,
mpmath (see `requirements.txt` at the top of the repository).

Numerical notes (all met while writing these):

- Double precision throughout, except where stated (mpmath for the coefficients aₖ at 130 digits,
  for the radial limits at the precision the terms need, for a few reference values).
- Φ_b is a logarithm only up to 2πi. Lemmas S and F are statements about Faddeev's integral itself,
  so S2 evaluates them with the integral (not with a sum of principal logarithms of the factors).
- In D(u) (Theorem C), p^{n²/2} means e^{−πi·n²/(2x)}, not a principal power of p; with principal
  powers the double series is wrong by signs for odd n.
- θ must be kept on the sheet [−210°, 10°]: np.angle returns 150° for θ = −210°.
- Identities that cancel large terms (Proposition 2 at |q| > 1, Theorem A when the residue sum and the
  dual sum are 10²⁰ and cancel, (B1) with p⁻¹ large) are measured against the size of the terms.
- At real x the state integrals over L cancel like e^{0.16·X}; S7 takes them over L only for X ≤ 8
  and over the descent contours of v19 otherwise. The horizontal line of v19 §4 is an admissible
  contour for the e^{w}-version I(b + 1/(2b)) only when |t| < π/3 (X > 6), see Findings.
- Oscillatory tails: along the descent contours the phase of e^{V/t} turns at a rate ∝ |w|/|t|;
  the panels shrink like 1/(1 + |w|/2) (`panels_osc` in `s6_state.py`).
