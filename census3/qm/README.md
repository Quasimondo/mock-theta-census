# q -> 1/q duality test for psi0*, psi1*

- `limits.py`: radial limits at all roots h/m, m <= 15 (72 roots), by Neville extrapolation of g(zeta e^{-t})
  with t = 0.003 j / m^2, j = 1..14. Accurate to ~4e-13 for psi0*, psi1*, sigma; it fails for the
  (q;q^2)_n-denominator series (psi*, F0*, F1*, T0*), where exact finite forms are used instead where known.
- `dual.py`: least-squares search for target(zeta) = zeta^e * sum a_i X_i(conj zeta) + sum b_j Y_j(zeta), e in j/120.

Check: it recovers sigma(zeta) = -sigma*(conj zeta) (Andrews-Dyson-Hickerson / Cohen), residual 1e-13.
Result: no such relation for psi1* or psi0* with X, Y in {psi1*, psi0*, sigma, sigma*, 1} (residual 0.67).

## Full library search (1,357 bounded census functions)

- `limits2.py`: radial limits of the 1,357 functions bounded at all 19 roots, at the 46 roots of order <= 12
  (`limits2.jsonl`). Extrapolation is reliable (relative error < 1e-6) at all 46 roots for 531 functions;
  the others lose precision near even-order roots through cancellation in the Eulerian sum.
- `dual2.py`: the fit T(zeta) = a zeta^e X(zeta') + c zeta^e U(zeta') + b for the 531 (zeta' = conj zeta or zeta,
  U = the other psi*), on all roots and separately on odd and even orders.
- `dual3.py`: the same for every function on the roots where it is reliable, if there are at least 24:
  1,110 functions.

Result: no q -> 1/q partner for psi1* or psi0*. The only fits are same-root ones that are series identities,
which check the search: the functions themselves (ids 830 = psi1*, 1382 = q-shift of psi0*), and
sum (-1)^n q^{n(n-1)/2} (-q;q)_n = 1 - psi1* - psi0* (id 278; true as q-series, checked to q^400).
