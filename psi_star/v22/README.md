# Scripts of addendum v22 (the growing state integral for every z; F(z) at even roots of unity)

Python 3 with numpy and mpmath. Run from inside this folder; `run_all.sh` runs everything and writes `run_log.txt` (about 15 minutes).

- `m1_contour.py`: first attempt at the left part of the contour, one straight ray from iπ at 30°, 60°, 75°, 85°; it fails for small a. Also checks the closed forms of Im V on the top edge of the strip and on the real axis, and the monotonicity on the vertical line through the saddle.
- `m2_contour.py`: the contour of Theorem 13 (ray at height 2π, segment down to iπ, top edge, vertical line, real axis), on the imaginary axis (90°) and tilted (85°), at arg x = 0°, 3°, 5°; reports the supremum of (H − H(w₀))/min(d², d) on each piece, which must be negative.
- `n1_zroots.py`: Theorem 14. The finite sum F(z; e^{2πi/X}) against i·P·e^{−πia}·[J(a/b) − I((a−1)/b) + z^X·I(a/b)] at x = 1/X, in double precision, with Φ_b from `fad.py`; the growing integral is integrated over the tilted contour, the decaying one over the horizontal line through its saddle, J over the line L.
- `fad.py` (from v16): Faddeev's integral in the variable w.
- `run_all.sh`, `run_log.txt`.
