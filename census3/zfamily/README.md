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
**Correction (27 Nov 2026, from `docs/addendum-v18.md` §2.4, rechecked):** `branches.py` selects the saddle only by matching
|Im V| to the measured rate, and that cannot distinguish the two roots of w² + w = 1/z (D(−W₂) = −D(−W₁)). For z = −1, i, −ω
it picks the lower root, for z = ω the upper one, and the candidates it lists are not all critical points. The rates,
volumes and tetrahedron shapes in the table are correct (the shapes of the two roots are equivalent: −W of one root is
1 + W of the other), but the saddle identification is not. The correct rule, fixed by the full complex asymptotics
(phase and constant), is the root W with Im W > 0, w = Log W, on the branch k = −1 of addendum-v18.

Only growth rates are tested here. `docs/addendum-v18.md` answers the follow-up: F(z) has the full structure of ψ₁*
(exact formula with state integrals and dual series 𝒲, three sectors) for every z, and the selected saddle is the root
W with Im W > 0, w = Log W (verified to 10⁻²⁶ … 10⁻³⁶ including the constant; not yet proved).
