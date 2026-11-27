# Scripts for addendum v12 (proof of the kernel identity)

Python 3 with mpmath. Run from inside this folder. `run_log.txt` is the complete output of `./run_all.sh 30`.

- `common.py` — Phi_b, the integrals I(u) (kernel-free) and J(u) (with kernel) on L = e^(-i 0.4) R, the q-series, B1..B4, D(u), c, P, kappa.
- `t1_kernel_free.py [digits]` — Theorem C: I(u) = c e^(-pi i u^2/2) D(u); special values against B2, B4; recursion (B1) by quadrature.
- `t2_kernel.py [digits]` — Theorem A, upper and lower route, at four values of u; the four formulas for psi1*, psi0*, f1*, f0*.
- `t3_steps.py [digits]` — the steps of the proof of Theorem C (inversion, Fourier transforms F and G, the two linearisations) and the dual sums of Theorem B.
- `t4_series.py [N]` — the three series identities of (B3), exact to p^N.
- `t5_bounds.py [digits]` — the two bounds of the closing lemma (section 4.1).
- `run_all.sh [digits]` — runs everything (about six minutes at 30 digits on two cores).
