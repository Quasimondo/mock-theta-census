# Mock theta search: addendum v10 (the one-sided law as an exact formula)

Session of 26 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). This tries the missing piece named in `addendum-v9.md` §3: an exact statement for the one-sided sums. Scripts and output are in `psi_star/v10/`.

**Independent check** (Claude Opus, same day, `psi_star/v10/phi_psi1.py`, `check_exact.py`): Φ was recomputed exactly by a different method (as half the series of f₁\*, via a log-cosh expansion) and agrees with v3. At 𝒳 = 50 − 0.5i, 60.3 − 0.4i, 80.4 − 0.5i, ψ₁\* and all terms of the formula were computed directly; ψ₁\* minus the two bracket terms agrees with Φ at optimal truncation to 1.2·10⁻¹³, 1.0·10⁻¹⁵, 2.1·10⁻²⁰, i.e. to the truncation error of Φ itself, while the decaying term is 1.9·10⁻⁴, 5.7·10⁻⁵, 1.8·10⁻⁶. So a multiple λ·W in the decaying sector is excluded to |λ| < 10⁻¹⁰ … 10⁻¹⁵. This does not use a Borel sum; the Borel-sum statement itself and §3 were not rechecked.

Everything here is **verified** (numerics) unless marked **inference**. No proof.

---

## 1. Result

Notation as in v5–v7: 𝒳 in the lower half-plane, q = e^{2πi/𝒳} = e^{−t}, t = −2πi/𝒳, p = e^{−πi𝒳}; W, W₃ from v4; X, N₁ and E±, E±′ from v6; and the two state integrals

  𝓑₂ = X(q)·E₊(p) + i·q^{1/4}·N₁(q)·E₊′(p),  𝓑₄ = X(q)·E₋(p) − i·q^{1/4}·N₁(q)·E₋′(p).

Let κ = θ₄(q)·(p;p²)∞ / (2i·q^{1/4}·(−q;q)∞·ψ(p)), and let Φ(t) be the perturbative series of ψ₁\* at q → 1 (v3). Then

  **ψ₁\*(q) = 𝒮[Φ](t) + (κ/2)·[ 𝓑₂·p·W₃(p) − 𝓑₄·W(p) ]**

where 𝒮[Φ](t) is the **Borel sum of Φ along the direction of t**. This is an identity at finite 𝒳, not an asymptotic statement.

- The two bracket terms are the decaying and the growing sector of v5, now with exact coefficients (state integrals, v7).
- The trivial sector is not just "Φ to all orders": it is the Borel sum of Φ, an actual function.
- The freedom "q̃W₃ up to a multiple of W" of v5 is fixed: with this resummation the multiple is **zero**.

Subtracting (II) gives the companion statement f₁\*(q) = 2𝒮[Φ](t) − κ·[𝓑₂·F_b(p) − 𝓑₄·F_a(p)], and the sum of the two is the exact bilateral formula of v6. So this is the one-sided rule in exact form: the bilateral identity splits term by term, n ≥ 0 sums to n ≥ 0 sums, and the remainder is one Laplace integral.

## 2. Evidence (`exact3.py`)

The Borel sum was computed from the 151 known coefficients of Φ: Padé approximant of the Borel transform (orders 60 and 50), then the Laplace integral along the ray of t.

Five points on horocycles, 𝒳 = 50 − 0.5i, 50.7 − 0.5i, 60.3 − 0.4i, 71.1 − 0.6i, 80.4 − 0.5i:

| | size |
|---|---|
| ψ₁\* and the growing term | 9·10³ … 1.7·10⁶ |
| decaying term | 2·10⁻⁴ … 2·10⁻⁶ |
| ψ₁\* − both terms, and 𝒮[Φ] | 0.333… (they agree in all 22 digits printed) |
| difference | 6·10⁻⁵⁶ … 1.3·10⁻⁵³ |

The two Padé orders agree to 10⁻⁵⁶ or better, and the working precision is 60 digits. So the identity holds to about 50 digits beyond the size of the decaying term, at five different values of p with |p| ≈ 0.15 … 0.28. A multiple λ of W in the decaying sector would show up as λ times 10⁻⁴; it is zero to 10⁻⁴⁹.

## 3. Other directions of t (`exact3_sectors.py`; indicative)

- For 𝒳 = 40 − 10i, 30 − 15i, 20 − 20i, 10 − 20i (arg t from −76° to −27°) the same formula holds, to 5·10⁻⁵⁵, 3·10⁻⁴⁸, 3·10⁻⁴¹, 2·10⁻³³. The loss of digits is the Borel–Padé accuracy at larger |t|.
- Between 𝒳 = 10 − 20i and 5 − 20i the difference becomes non-zero. That is the Stokes ray arg t = −17.1°, the direction of the action S = π²/3 − i·Cl₂(π/3). Beyond it, (ψ₁\* − terms − 𝒮[Φ]) / ((κ/2)·𝓑₂·W) = 1.004, 1.030, 1.001, 0.999998 at four points, within the accuracy of the resummation there.
- So across that ray the Borel sum of Φ jumps by (κ/2)·𝓑₂·W(p): Stokes constant 1, carried by W (inference from four points at limited accuracy). In the central sector the decaying coefficient is p·W₃ + W instead of p·W₃. A second jump is expected at arg t = +17.1°; the last point is consistent with it but I did not resolve it.

## 4. What this means for a proof (inference)

The one-sided law for ψ₁\* is now equivalent to one statement: the function

  𝒯(𝒳) = ψ₁\*(q) − (κ/2)·[𝓑₂·p·W₃(p) − 𝓑₄·W(p)]

is analytic in t in the sector −90° < arg t < −17.1° up to its boundary ray on the real 𝒳-axis, with Φ as its asymptotic series there with factorial bounds. Then it is the Borel sum of Φ by Watson's theorem, and everything in v5 follows. The state integrals of v7 already take care of 𝓑₂, 𝓑₄. What is missing is a direct handle on 𝒯: for instance a closed form of the Borel transform of Φ, which would make 𝒮[Φ] an explicit Mordell-type integral. I did not find one.

Not done: ψ₀\* (same test with its own Φ and the factors e^{±2πi/3}; expected to work, not run); the analogous formula for W and 𝒥 in the bounded directions; the closed form of the Borel transform.

## 5. Scripts

`exact3.py` (the test at five points), `exact3_sectors.py` (nine more directions). mpmath only; under a minute in total.
