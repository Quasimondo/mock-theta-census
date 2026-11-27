# Mock theta search: addendum v11 (the trivial sector is a state integral with a kernel)

Session of 26 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). This answers the question relayed from the other session: is 𝒯 = ψ₁\* − (κ/2)[𝓑₂·p·W₃ − 𝓑₄·W] of `addendum-v10.md` itself a state integral, perhaps with one more elementary factor? Script and output are in `psi_star/v11/`.

**Independent check** (Claude Opus, same day, `psi_star/v11/check_kernel.py`; 𝒯 from the code of `psi_star/v10/check_exact.py`, own Φ_b and contour, δ = 0.15): at x = 0.3 + 0.2i, 1/(3 − 0.4i), 0.11 + 0.03i the kernel integral with the stated prefactor equals 𝒯 to 3·10⁻³⁶, 3·10⁻³⁸, 4·10⁻⁴¹ (40-digit quadrature), while ψ₁\* itself is 1.2 − 0.23i, 1.75 − 0.10i, −1.70 − 0.37i there.

**Yes** (verified numerically at three points; no proof).

## 1. Result

With b² = x, q = e^{2πix}, p = e^{−πi/x}, and Φ_b as in v7:

  **𝒯(x) = (−q;q)∞ / (b·q^{1/8}·(p;p²)∞) · ∫_L e^{−πi·y² + π·y/b} / ( Φ_b(y)·(1 + e^{2πy/b}) ) dy,**

where L is the real line rotated by a small negative angle (L = e^{−iδ}·ℝ, through 0). Compared with the integrals of v7 this has 1/Φ_b instead of Φ_b, Gaussian exponent 1 instead of ½, u = 1/(2b), and the extra factor 1/(1 + e^{2πy/b}), which is the one suggested.

Together with v7 and v10, every term of the three-sector formula is now an explicit integral times a q-series in the dual nome:

  ψ₁\*(q) = [kernel integral above] + (κ/2)·[ 𝓑₂·p·W₃(p) − 𝓑₄·W(p) ],

with 𝓑₂, 𝓑₄ the state integrals of v7. By v10 the kernel integral is also the Borel sum of the perturbative series Φ in the sector −90° < arg t < −17.1°.

## 2. How it was found, and the check (`kernel_integral.py`)

- The kernel has simple poles at y = i·b·(n + ½). At those points 1/Φ_b(y) = ((p;p²)∞/(−q;q)∞)·(−q;q)ₙ, the Gaussian gives q^{(n+½)²/2}, and e^{πy/b} gives i·(−1)ⁿ. So the residues at n ≥ 0 are exactly the terms of ψ₁\*, and the prefactor above is the one that would make "integral = sum of residues at n ≥ 0" read ψ₁\*.
- That naive closing is **not** valid, and the integral is not ψ₁\*: it differs from it by order 1 (ratios off by 0.7 … 1.2). It equals 𝒯 instead. At x = 0.3 + 0.2i, 1/(3 − 0.4i), 0.11 + 0.03i the integral and 𝒯 agree in all 20 digits printed, both at 22-digit and at 30-digit working precision. The growing and decaying terms at these points are of size 0.8 … 2 and 10⁻⁴ … 0.1, so the agreement is not a triviality.
- A second candidate, the half-line integral without the kernel, does not match (difference 0.1).

So the difference between "integral" and "residues at n ≥ 0" is exactly the pair of bilinear terms with W and p·W₃. That is where the one-sided dual series come from (inference): the correct closing of the contour must pick up a second family of contributions organised by the dual nome.

## 3. Status of the one-sided law for ψ₁\*

What is now needed for a proof is one contour computation: show that the kernel integral equals ψ₁\* minus the two bilinear terms. Everything else is in place or standard:

- the kernel integral is analytic in b across the positive real axis, like the integrals of v7, so nothing in it can depend on p asymptotically;
- its expansion at b → 0 is the perturbative series Φ (inference, consistent with v10 to 55 digits);
- 𝓑₂, 𝓑₄ are state integrals (v7).

Not done: the same for ψ₀\* (expected: the same kernel integral with another u), for f\*, and for W, 𝒥 in the bounded directions; the contour computation itself.

## 4. On the remarks from the other session

- The figure-eight comparison fits what is here: the saddle values ±i·Cl₂(π/3) and the Hessian √3/2 of v7 §3 are the volume and torsion data you describe, and the kernel integral above has the form you proposed.
- I have **not** read the Garoufalidis–Gu–Mariño papers, only search snippets of related work. So "this is the pattern of the trivial flat connection of 4₁" is your suggestion plus my numerics, not something I have checked against their statements. In particular I do not know whether their integrals for the trivial connection carry the same kind of kernel.

## 5. Script

`kernel_integral.py [digits]`: ψ₁\*, the two bilinear terms, 𝒯, and the two candidate integrals at three values of x. mpmath only; about a minute at 30 digits.
