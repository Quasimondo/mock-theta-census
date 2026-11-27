# Stress test of the proofs, second round (v22–v26)

Scripts for `docs/psi-star-stress-test-brief-v2.md`. Report: `docs/psi-star-stress-test-v2.md`;
complete output of the last run: `run_log.txt` (assembled from `out_<block>.txt` by `merge_log.py`).

Own code throughout. Φ_b comes from the round-1 library `psi_star/stress/core.py` (products on both sides
of the real axis, Faddeev's integral on Im σ = ½ with the inversion relation and the functional equation),
plus a new 25-digit Faddeev integral in `r3p.py`. Nothing is imported from `psi_star/v*`. Existing values
were compared only afterwards (the numbers quoted in v22–v26).

| file | block | what | run time (2 cores) |
|---|---|---|---|
| `common2.py` | – | saddles W, Vol(z), V with the continuation of v22 §1, the contours C and C′, the integrals over C and over tilted lines, F(z; q₀) and the prediction of Theorem 15 | – |
| `t13a.py` | T13a | Theorem 13 (a): closed forms, (L1) on C₁…C₅ for a ∈ [10⁻³, 1 − 10⁻³] (and 10⁻⁵, 1 − 10⁻⁵), the bound for 1/Φ_b on each piece, the uphill claim of v22 §5, the asymptotics | ~7 min |
| `t13b.py` | T13b | Theorem 13 (b): C′, the bound on C₂′, ε(a) by bisection over the tilt φ₁, the asymptotics up to arg x = ε(a) | ~6 min |
| `t14.py` | T14 | Theorem 14 for X = 2 … 50 and random a; the dual series on x = 1/X + iε; odd X | ~40 min |
| `lemma_j.py` | J | Lemma U; the wedge bounds used to turn the contour; \|P·J(a/b)\| against 1/(2√2 sin πa) for a ∈ [10⁻³, 1 − 10⁻³], x ∈ (0, 1] | ~25 min |
| `t15.py` | T15 | K(a, X) = X(F/prediction − 1) for even X up to 6400 and a from 10⁻³ to 1 − 10⁻³ | ~2 min |
| `m_z1.py` | M | Proposition M, Theorems 16, 17 at random x; (R3), (R4) to q¹⁰⁰⁰ exactly (python-flint, two primes); Theorem 17 on the unit circle, X = 1 … 100, and the continuity it needs | ~12 min |
| `r3p.py` | R3′ | Lemma F for x ∈ [¼, 1] and complex neighbourhoods; Lemma D for \|t\| ≤ 2π·1.01; Theorems 9, 11, 14 at X = 1, 2, 3 with 25-digit arithmetic | ~25 min |
| `p3.py` | P3 | closing lemma of v12 §4.1; Lemma L on a model; Theorem 8 as a whole; Theorem C at 30 digits outside the round-1 sample; Lemma V and the e^{w} growth rate of v26 | ~15 min |
| `run_all.sh`, `merge_log.py` | – | everything, two streams | ~1.5 h |

Requirements: numpy, scipy, mpmath, python-flint. Fixed seeds.

Numerical notes (met while writing these, in addition to round 1's):

- **Near the saddle** the difference H(w) − H(w₀) is below rounding for |w − w₀| < 10⁻⁴; a first run of
  (L1) "found" ratios of +10² on C₄. `common2.dV` uses the Taylor series at the saddle for
  |w − w₀| < 10⁻³.
- **The kernel 1/(1 + e^{w/x}) grows** like e^{|w|/|x|} on one side of a tilted line when x is complex;
  the truncation of the line has to include that rate (a first run of Theorem 17 at x = −0.13 + 0.12i
  was off by 3·10⁻⁷).
- **Kernel integrals at complex x** are sensitive to the tilt δ of the line when δ approaches 90° − β, where the
  line runs close to the kernel poles i·b·(n + ½): use δ ≈ (90° − β)/2.
- **Faddeev's integral for Re w > 0** directly (without the inversion relation) loses all digits for
  w ~ 60 at small x (terms e^{Re w/2}); Lemma U must be tested with the inversion relation or for |w| ≤ 10.
- **ψ₁\* near |q| = 1** in double precision: the product (−q;q)ₙ over several hundred factors loses about
  1 % at arg x = 0.3°, |t| = 0.05; `p3.py` sums at 40 digits.
- **Slowly decaying lines** (a near 0 or 1 for the horizontal line of Theorem 7, C₁ for a near 0, C₅ for a near 1):
  the decay rate is proportional to min(a, 1 − a) while the phase turns like r/|t|; `t13a.py` and `t14.py`
  skip or switch to the 45° line there.
