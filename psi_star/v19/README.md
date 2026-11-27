# Scripts of addendum v19 (asymptotics of the state integrals)

Python 3 with numpy and mpmath. Run from inside this folder: `python3 e1_checks.py > run_log.txt` (about 40 seconds).

- `e1_checks.py [a b c d e]`
  - a: the height H = −Im V + τ·Re V on the three pieces of the contour of the growing integral (a′ = −½), τ = 0, 0.1, 0.2, and the closed formulas for Re V, Im V on the pieces;
  - b: the same on the horizontal line through the saddle for the decaying integrals, a = ½ and five other a;
  - c: the inequality ln(2cos(η/2)) ≥ (ln 2/π)(π − 3η/2) on [0, 2π/3];
  - d: 𝓑₄, 𝓑₃, 𝓑₂, 𝓑₁ from Theorem C of v12 against (4/3)^{1/4} e^{∓iπ/12} e^{∓i·Cl₂/t}, at arg x = 0.5°, 5°, 11°;
  - e: P·I(±1/(2b)) and q^{1/2}·P·I(b ± 1/(2b)) against the coefficients stated in Theorem 8.
- `common.py` (from v12).
- `run_log.txt`.
