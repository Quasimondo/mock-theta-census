# Scripts of addendum v26 (corrections after the stress test)

Python 3 with numpy and mpmath. Run from inside this folder: `python3 u1_fixes.py > run_log.txt` (about 2.5 minutes).

- `u1_fixes.py [a b c d]`
  - a: Lemma V (Re V has its maximum at the saddle on the horizontal line, for every a).
  - b: Theorem 7 at arg x = 30°, 45°, 60° against the series of Theorem C.
  - c: the e^{w}-version of the decaying integral on the line Im w = −2π/3 at real x: growth rate on the right, and the line against L.
  - d: analyticity of P·J across the real x-axis (Fourier modes on circles in the t-plane).
- `h1_roots.py` (from v20), `fad.py` (from v16), `common.py`, `a1_formula.py` (from v18), `run_log.txt`.
