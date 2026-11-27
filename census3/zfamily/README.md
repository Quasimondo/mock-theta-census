# F(z) at roots of unity: growth rates are tetrahedron volumes (2026-11-26)

F(z) = Σ zⁿ q^{n(n+1)/2} (−q;q)ₙ, z = e^{2πia}; ψ₁ = F(1) (mock theta), ψ₁* = F(−1).

- `growth.py`: exact values of F(z) at q = e^{2πi/X} (terminating sum at even order, A_z/(1 − 2z^m) at odd order;
  addendum-v3 §2) for X = 400, 800, 1600, and the slope of log|F| in X.
- `branches.py`: the measured slope equals |Im V|/2π at one stationary point of the saddle function V of
  addendum-v3 §2 (after the finite-X correction ½·log 2/800 from the t^{−1/2} prefactor), for z = −1, i, ω, −ω.
- `volumes.py`: that |Im V| is exactly |D(shape)|, the Bloch–Wigner dilogarithm, i.e. the volume of an ideal
  tetrahedron whose shape solves w² + w = 1/z (up to the tetrahedral symmetries of the shape parameter).

| z | growth rate | tetrahedron volume |
|---|---|---|
| −1 | 0.16153 | 1.014942 (regular) |
| i | 0.11088 | 0.696701 |
| ω | 0.13769 | 0.865129 |
| −ω | 0.07745 | 0.486618 |
| 1 | 0 | 0 (real shape, golden ratio) |

Reading (inference): z deforms the tetrahedron; the mock theta function ψ₁ is the degenerate, zero-volume member.
Only growth rates are tested; whether F(i), F(ω), F(−ω) have the full structure of ψ₁* (dual sum, state integrals,
three sectors) is open. The branch of the saddle that dominates is identified numerically, not derived.
