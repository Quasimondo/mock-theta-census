# Scripts for addendum v13 (the perturbative series from the kernel integral; the kernel integral against the Borel sum)

Python 3 with mpmath and python-flint (s1 only). Run from inside this folder. `run_log.txt` is the output of the commands listed in it.

- `s1_window.py [N]` — exact coefficients of Phi for psi1*, psi0* from the 1/cosh window (Bernoulli polynomials, exact rationals);
  compares with `q2_series_v3.json` (v3, 151 coefficients, 120 digits) and writes `phi_exact.json` (here N = 300; 40 s).
- `borel.py` — Pade-Laplace Borel sums S_theta of Phi from `phi_exact.json`; the kernel integral P J(x) for x on both sides of the real axis.
- `s2_poles.py [n ...]` — poles of the Pade approximants of the Borel transform (50 s).
- `s3_sectors.py [digits] [|t|] [case] [phi ...]` — P J against the Borel sum for directions phi = arg t (about 1 minute per direction at 30 digits).
- `s3b_lateral.py` — the edge directions against lateral Borel sums, several Pade orders (80 s).
