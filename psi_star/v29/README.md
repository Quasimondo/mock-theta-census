# Scripts of addendum v29 (Casoratians of the trivial series with the homogeneous solutions)

Python 3 with mpmath. Run from inside this folder: `python3 y1_casoratian.py > run_log.txt` (about 3 minutes).

- `y1_casoratian.py [a b c] [N]`
  - a: the q-side Casoratian D: its first-order equation, its Laurent series, the sum over shifts, the two identities at z = ∓1, and the test that the bilateral sum is not a single product.
  - b: Rogers' false theta identity for σ₋, exact to order N = 300.
  - c: the dual side: the equation of 𝒲, the recursion of I, K(a) + K(a+1) = −I(a/b), K as a sum over shifts, Theorem Z, and K at z = 1.
- `common.py`, `a1_formula.py`, `c2_z1.py` (from v18), `run_log.txt`.
