# Scripts of addendum v24 (z = 1 and z = q: Rogers–Ramanujan S-matrix, the fifth-order laws, the Mordell integral)

Python 3 with numpy and mpmath. Run from inside this folder: `sh run_all.sh > run_log.txt` (about 5 minutes).

- `s1_series.py [N]`: the series identities used, exact to order N = 300 in integer arithmetic: (R1)–(R4) between the fifth-order functions, (R5) Ramanujan's identity for G and H, (S1)–(S6) of Rogers–Slater type.
- `s2_law.py [a b c d e] [digits]`: at six points x, 30 digits:
  - a: the four state integrals against their closed forms (Proposition M);
  - b: the four laws of Theorem 16 against the q-series;
  - c: the kernel integrals against the Mordell integrals (Theorem 17);
  - d: the quoted law for χ₀, χ₁ on its own;
  - e: the steps of the proof of Proposition M (bilinear forms, transformation formulas, eta identity).
- `s3_real.py`: Theorem 17 on the unit circle, x = 1/X, double precision.
- `common.py`, `a1_formula.py`, `c2_z1.py` (from v18), `n1_zroots.py` (from v22), `fad.py` (from v16), `run_all.sh`, `run_log.txt`.
