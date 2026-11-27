# Scripts of addendum v18 (the z-family F(z) = Σ zⁿ q^{n(n+1)/2} (−q;q)ₙ)

Python 3 with mpmath. Run from inside this folder. `run_all.sh` runs everything and writes `run_log.txt` (about 14 minutes; almost all of it is `a1_formula.py`).

- `common.py` (from v12): Φ_b, the integrals I(u), J(u) by quadrature, the q-series, D(u), c, P.
- `a1_formula.py [digits]`: Theorem Z, both routes, at x = 0.3 + 0.2i, 1/(3 − 0.4i), 0.11 + 0.03i and nine choices of a (z = i, −ω each with two branches of a; z = −1, 1, ω; one complex a). Also D(a/b) through the function H, and I(a/b) by quadrature against Theorem C.
- `b1_roots.py [digits]`: F(z) at q = e^{2πi/X}, X = 802 … 9602, as a finite sum with 0.07·X extra digits, against the saddle-point formula with the rule "root W in the upper half-plane, w = Log W"; the identity Im V = D(−W) on all branches; Vol(z) near z = 1.
- `b2_dominant.py [digits]`: I((a−1)/b) and I(a/b) from Theorem C against their saddle-point values, and F against the m = 1 dual term, at |x| = 1/99, 1/199, 1/399, arg x = 1°. Needs about 90 digits (the q-series cancel).
- `c1_series.py [N]`: the series identities at z = 1 to order N (integer arithmetic).
- `c2_z1.py [digits]`: the strip form of Theorem Z at z = i, −ω; ψ₁ and f₁ with the fifth-order mock theta functions of the dual nome; the Rogers–Ramanujan forms of D(0), D(1/b).
- `d1_compare.py`: the saddles selected in `census3/zfamily/volumes.py` of the repository, translated to the notation of v18; the ratio |ψ₁*|/J_{4₁}(X/2).
- `run_all.sh`, `run_log.txt`.

Note: the a1 section of `run_log.txt` is from the first run of the session; the other sections were rerun after the last edit of the scripts.
