# Census 3a: one-parameter families at roots of unity

`zfam.py` probes three families at every z, x = e^{2πik/N} with N ≤ 12 (136 functions,
2,584 probes, 19 roots h/m with h ≤ m/2; for complex coefficients the cells of x and of x̄
together cover all 36 roots on the circle):

- F(z) = Σ zⁿ q^{n(n+1)/2} (−q;q)ₙ (ψ₁ = F(1), ψ₁* = F(−1))
- g₂(x,q) = Σ (−q;q)ₙ q^{n(n+1)/2} / ((x;q)ₙ₊₁ (q/x;q)ₙ₊₁)
- g₃(x,q) = Σ q^{n(n+1)} / ((x;q)ₙ₊₁ (q/x;q)ₙ₊₁)

`zreport.py` prints the fingerprints (`zreport.txt`); `zextra.py` tests roots of order 9, 11, 13 and 15.

## Results

**g₂, g₃: all 90 are E/Q at every root**, as expected for the universal mock theta functions.
At q = 1 the exponent follows a clean law in a = k/N (taking a ≤ 1/2; a and 1 − a agree):

- g₃: κ = (1 − 6a)²/24 for a < 1/6, bounded for a ≥ 1/6
- g₂: κ = (1 − 4a)²/16 for a < 1/4, bounded for a ≥ 1/4

**F(z): the growth depends only on the order m of the root and on z^m.**

| Root of order m | z^m | Cells |
|---|---|---|
| even | any | Q: 540 of 540 |
| odd | exactly 1 | R, κ = 1/240 (as for ψ₁ = F(1) at q = 1, but with an inexact law): 28 of 28 |
| odd | within 1/10 of a turn of 1 | X (irrational κ, dilogarithm type): 26 of 26 |
| odd | 1/9 to 1/8 of a turn from 1 | boundary zone: X, Q and timeouts mixed |
| odd | at least 1/7 of a turn from 1 | Q: 233 of 233 |

So F(z) is bounded at every root of unity only if z^m stays far from 1 for every odd m.
This holds for z of order 2, 4 and 6, i.e. z = −1, ±i, −ω, −ω̄. It fails for every other
order N: the roots of order N itself (N odd) or of order m with z^m near 1 give growth.
F(e^{2πi·5/11}) and F(e^{2πi·6/11}) look bounded at all 19 standard roots, but `zextra.py`
shows them growing at roots of order 9 (X) and 11 (R, κ = 1/240). So the bounded companions
in this family are exactly ψ₁* = F(−1) and the already known F(i), F(−ω) with their conjugates.
