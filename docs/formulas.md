# Formula sheet: every object in several notations

The same q-series can be written in many equivalent ways, which makes it hard to tell whether a formula has been published before. This sheet lists each object studied here in several forms (rendered LaTeX, raw LaTeX source, plain text and Mathematica syntax), together with equivalent rewritings and the opening coefficients. **Every q-series identity and equivalent form below is checked numerically to 50 digits by [`psi_star/formulas_check.py`](../psi_star/formulas_check.py).** The integral formulas of §10–§13 involve Faddeev's quantum dilogarithm; they are checked by the scripts named with them (the `check_*.py` scripts of the addenda), and the consolidated documents [`psi-star-theorems-v1.md`](psi-star-theorems-v1.md) and [`zfamily-theorems-v1.md`](zfamily-theorems-v1.md) state them with their hypotheses.

Notation: $(a;q)_n=\prod_{k=0}^{n-1}(1-aq^k)$, so $(-q;q)_n=\prod_{k=1}^{n}(1+q^k)=\dfrac{(q^2;q^2)_n}{(q;q)_n}$, also written $(-q)_n$. Note that $q^{n(n+1)/2}=q^{\binom{n+1}{2}}$.

Keywords: sign-twisted, alternating, (−1)ⁿ inserted, z = −1 specialization, fifth-order mock theta functions ψ₀, ψ₁, Rogers–Ramanujan, entire Rogers–Ramanujan function, bilateral series, strange identity, quantum modular form, resurgence, Stokes constants, q = −1 asymptotics, state integral, Faddeev quantum dilogarithm, Borel summation, ideal tetrahedron, Bloch–Wigner dilogarithm, Kashaev invariant, figure-eight knot, Mordell integral, Rogers–Ramanujan S-matrix, transformation law of the fifth-order mock theta functions.

---

## 1. ψ₁\*: sign-twisted ψ₁

$$\psi_1^{\ast}(q)=\sum_{n\ge 0}(-1)^n q^{n(n+1)/2}(-q;q)_n$$

```latex
\psi_1^{*}(q)=\sum_{n\ge 0}(-1)^n q^{n(n+1)/2}(-q;q)_n
```

- Unicode: ψ₁*(q) = Σₙ≥₀ (−1)ⁿ q^{n(n+1)/2} (−q;q)ₙ
- Plain: `psi1*(q) = sum_{n>=0} (-1)^n q^(n(n+1)/2) (-q;q)_n`
- Mathematica: `Sum[(-1)^n q^(n (n + 1)/2) QPochhammer[-q, q, n], {n, 0, Infinity}]`

Equivalent forms:

$$\psi_1^{\ast}(q)=\sum_{n\ge0}(-1)^n q^{\binom{n+1}{2}}\prod_{k=1}^{n}(1+q^k)=\sum_{n\ge0}(-1)^n q^{\binom{n+1}{2}}\frac{(q^2;q^2)_n}{(q;q)_n}={}_2\phi_2\!\left(\begin{matrix}q,\,-q\\ 0,\,0\end{matrix};q,q\right)=F(-1)$$

```latex
\sum_{n\ge0}(-1)^n q^{\binom{n+1}{2}}\prod_{k=1}^{n}(1+q^k)
\sum_{n\ge0}(-1)^n q^{\binom{n+1}{2}}\frac{(q^2;q^2)_n}{(q;q)_n}
{}_2\phi_2\left(\begin{matrix}q,\,-q\\ 0,\,0\end{matrix};q,q\right)
```

(Basic hypergeometric series in the Gasper–Rahman convention; F(z) is the family in §5.)

Coefficients (q⁰ … q⁴⁰):
`1, -1, -1, 1, 1, 1, 0, -1, -1, -2, 0, 0, 0, 2, 2, 1, 1, 1, -1, -1, -2, -2, -2, -2, -1, -1, 1, 2, 2, 2, 3, 3, 3, 2, 0, -1, 0, -2, -4, -3, -5`

## 2. ψ₀\*: sign-twisted ψ₀

$$\psi_0^{\ast}(q)=\sum_{n\ge 0}(-1)^n q^{(n+1)(n+2)/2}(-q;q)_n$$

```latex
\psi_0^{*}(q)=\sum_{n\ge 0}(-1)^n q^{(n+1)(n+2)/2}(-q;q)_n
```

- Unicode: ψ₀*(q) = Σₙ≥₀ (−1)ⁿ q^{(n+1)(n+2)/2} (−q;q)ₙ
- Plain: `psi0*(q) = sum_{n>=0} (-1)^n q^((n+1)(n+2)/2) (-q;q)_n`
- Mathematica: `Sum[(-1)^n q^((n + 1) (n + 2)/2) QPochhammer[-q, q, n], {n, 0, Infinity}]`

Equivalent forms:

$$\psi_0^{\ast}(q)=q\sum_{n\ge0}(-1)^n q^{n(n+3)/2}(-q;q)_n=-\sum_{n\ge1}(-1)^n q^{\binom{n+1}{2}}(-q;q)_{n-1}=q\,F(-q)$$

```latex
q\sum_{n\ge0}(-1)^n q^{n(n+3)/2}(-q;q)_n
-\sum_{n\ge1}(-1)^n q^{\binom{n+1}{2}}(-q;q)_{n-1}
```

Coefficients (q⁰ … q⁴⁰):
`0, 1, 0, -1, -1, 0, 1, 1, 1, 1, -1, -1, -1, -2, -1, 0, 0, 1, 2, 2, 2, 1, 1, 0, -1, -1, -3, -3, -2, -2, -2, -1, 0, 1, 3, 3, 2, 4, 4, 3, 3`

Elementary relation:

$$1-\psi_1^{\ast}(q)-\psi_0^{\ast}(q)=\sum_{n\ge0}(-1)^n q^{n(n-1)/2}(-q;q)_n$$

## 3. The two census hits containing ψ\*

$$\sum_{n\ge0}q^{2n^2+n}(-q;q)_{2n}=\tfrac12\big(\psi_1(q)+\psi_1^{\ast}(q)\big),\qquad \sum_{n\ge0}q^{2n^2+3n}(-q;q)_{2n}=\frac{\psi_0(q)+\psi_0^{\ast}(q)}{2q}$$

```latex
\sum_{n\ge0}q^{2n^2+n}(-q;q)_{2n}=\tfrac12\big(\psi_1(q)+\psi_1^{*}(q)\big)
\sum_{n\ge0}q^{2n^2+3n}(-q;q)_{2n}=\frac{\psi_0(q)+\psi_0^{*}(q)}{2q}
```

Here ψ₁(q) = Σ q^{n(n+1)/2}(−q;q)ₙ and ψ₀(q) = Σ q^{(n+1)(n+2)/2}(−q;q)ₙ are Ramanujan's fifth-order mock theta functions.

## 4. Identities (I), (II)

$$2\,\psi_0^{\ast}(q)=f_0^{\ast}(q)-\vartheta_4(q)\,R_0(q),\qquad 2\,\psi_1^{\ast}(q)=f_1^{\ast}(q)+\vartheta_4(q)\,R_1(q)$$

```latex
2\,\psi_0^{*}(q)=f_0^{*}(q)-\vartheta_4(q)\,R_0(q),\qquad
2\,\psi_1^{*}(q)=f_1^{*}(q)+\vartheta_4(q)\,R_1(q)
```

with

$$f_0^{\ast}(q)=\sum_{n\ge0}\frac{(-1)^n q^{n^2}}{(-q;q)_n},\quad f_1^{\ast}(q)=\sum_{n\ge0}\frac{(-1)^n q^{n^2+n}}{(-q;q)_n},\quad R_0(q)=\sum_{n\ge0}\frac{(-1)^n q^{n^2}}{(q;q)_n},\quad R_1(q)=\sum_{n\ge0}\frac{(-1)^n q^{n^2+n}}{(q;q)_n},$$

$$\vartheta_4(q)=\sum_{n\in\mathbb Z}(-1)^n q^{n^2}=\frac{(q;q)_\infty}{(-q;q)_\infty}.$$

```latex
f_0^{*}(q)=\sum_{n\ge0}\frac{(-1)^n q^{n^2}}{(-q;q)_n},\quad
f_1^{*}(q)=\sum_{n\ge0}\frac{(-1)^n q^{n^2+n}}{(-q;q)_n},\quad
R_0(q)=\sum_{n\ge0}\frac{(-1)^n q^{n^2}}{(q;q)_n},\quad
R_1(q)=\sum_{n\ge0}\frac{(-1)^n q^{n^2+n}}{(q;q)_n}
```

- Plain: `2 psi0* = f0* - theta4 * R0`, `2 psi1* = f1* + theta4 * R1`
- f₀\*, f₁\* are the fifth-order f₀(q) = Σ q^{n²}/(−q;q)ₙ and f₁(q) = Σ q^{n²+n}/(−q;q)ₙ with (−1)ⁿ inserted.
- R₀, R₁ are Ramanujan's entire Rogers–Ramanujan function A_q(a) = Σ (−a)ⁿ q^{n²}/(q;q)ₙ at a = 1 and a = q. For comparison, A_q(−1) = G(q) and A_q(−q) = H(q) are the Rogers–Ramanujan functions.
- These are the z = −1 analogues of the classical bilateral identities B(f₀) = f₀ + 2ψ₀ and B(f₁) = f₁ + 2ψ₁.

Coefficients (q⁰ … q⁴⁰):
- f₀\*: `1, -1, 1, -1, 2, -2, 1, -1, 2, -3, 2, -1, 3, -4, 2, -2, 5, -5, 3, -4, 5, -5, 5, -4, 6, -9, 6, -5, 9, -8, 7, -9, 9, -11, 11, -9, 13, -14, 11, -13, 16`
- f₁\*: `1, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -2, 2, 0, 1, -3, 2, -1, 2, -2, 2, -2, 1, -3, 4, -1, 2, -4, 2, -3, 4, -3, 5, -4, 3, -6, 5, -4, 6, -7`
- R₀: `1, -1, -1, -1, 0, 0, 1, 1, 2, 1, 2, 1, 1, 0, 0, -2, -1, -3, -3, -4, -3, -5, -3, -4, -2, -3, 0, -1, 3, 2, 5, 5, 9, 7, 11, 9, 13, 10, 13, 9, 12`
- R₁: `1, 0, -1, -1, -1, -1, 0, 0, 1, 1, 2, 2, 2, 2, 2, 1, 1, 0, -1, -2, -2, -4, -4, -5, -5, -6, -5, -6, -4, -4, -3, -2, 1, 1, 4, 5, 8, 9, 12, 12, 15`

Double-sum forms:

$$(-q;q)_\infty\,\psi_1^{\ast}(q)=\sum_{\substack{k\ge0\\ k\ \text{even}}}\ \sum_{n\ge0}\frac{(-1)^n q^{k(k+1)/2+n^2+n+kn}}{(q;q)_k},\qquad (-q;q)_\infty\,\psi_0^{\ast}(q)=\sum_{\substack{k\ge1\\ k\ \text{odd}}}\ \sum_{n\ge0}\frac{(-1)^n q^{k(k+1)/2+n^2+kn}}{(q;q)_k}$$

```latex
(-q;q)_\infty\,\psi_1^{*}(q)=\sum_{k\ge0,\ k\ \mathrm{even}}\ \sum_{n\ge0}\frac{(-1)^n q^{k(k+1)/2+n^2+n+kn}}{(q;q)_k}
(-q;q)_\infty\,\psi_0^{*}(q)=\sum_{k\ge1,\ k\ \mathrm{odd}}\ \sum_{n\ge0}\frac{(-1)^n q^{k(k+1)/2+n^2+kn}}{(q;q)_k}
```

## 5. The family F(z) and its q-difference equation

$$F(z;q)=\sum_{n\ge0}z^n q^{n(n+1)/2}(-q;q)_n,\qquad F(z)=1+zq\,F(zq)+zq^2\,F(zq^2)$$

```latex
F(z;q)=\sum_{n\ge0}z^n q^{n(n+1)/2}(-q;q)_n,\qquad F(z)=1+zq\,F(zq)+zq^2\,F(zq^2)
```

- Mathematica: `F[z_] := Sum[z^n q^(n (n + 1)/2) QPochhammer[-q, q, n], {n, 0, Infinity}]`
- Specializations: ψ₁ = F(1), ψ₁\* = F(−1), ψ₀\* = q·F(−q), and 1 − ψ₁\* − ψ₀\* = F(−1/q).

## 6. Behavior at roots of unity and at q = −1 (numerical)

- **Terminating at even-order roots:** (−ζ;ζ)ₙ = 0 for n ≥ m/2, so the series is a finite sum there. Values: ψ₁\*(−1) = 1, ψ₁\*(i) = 2 − i, ψ₁\*(e^{iπ/3}) = 1 − 2√3 i, ψ₁\*(e^{iπ/4}) ≈ −3.5355339 − 4.1213203i.
- **Radial limits at odd-order roots (observation, verified numerically):** for ζ a primitive m-th root of unity with m odd, write the series as Σ aₙ(q). Then the terms repeat with a factor −2 every m steps, a_{n+m}(ζ) = −2·aₙ(ζ), because (−1)^m = −1, ∏_{k=1}^{m}(1+ζᵏ) = 2 and the ζ-power gains a multiple of m. The radial limit is the value of that divergent geometric series, 1/(1+2) times the first period:

$$\lim_{t\to0^+}\psi_1^{\ast}(\zeta e^{-t})=\frac13\sum_{n=0}^{m-1}(-1)^n\zeta^{n(n+1)/2}(-\zeta;\zeta)_n,\qquad \lim_{t\to0^+}\psi_0^{\ast}(\zeta e^{-t})=\frac13\sum_{n=0}^{m-1}(-1)^n\zeta^{(n+1)(n+2)/2}(-\zeta;\zeta)_n\qquad(m\ \text{odd})$$

```latex
\lim_{t\to0^+}\psi_1^{*}(\zeta e^{-t})=\frac13\sum_{n=0}^{m-1}(-1)^n\zeta^{n(n+1)/2}(-\zeta;\zeta)_n,\qquad
\lim_{t\to0^+}\psi_0^{*}(\zeta e^{-t})=\frac13\sum_{n=0}^{m-1}(-1)^n\zeta^{(n+1)(n+2)/2}(-\zeta;\zeta)_n\qquad(m\ \text{odd})
```

  Checked at all 50 cases (both functions at every primitive root of odd order m ≤ 15, up to complex conjugation) by extrapolating high-precision values along the ray to t = 0. Agreement to 10⁻¹⁴ in every case, the accuracy of the extrapolation; script `psi_star/s12_odd_limits.py`. Examples: ψ₁\*, ψ₀\* → 1/3 as q → 1 (formally Σ(−2)ⁿ = 1/3), and ψ₁\* → 1, ψ₀\* → 0 at the primitive cube roots. Together with the even-order case (where the factor is 0 and the series stops), the limit at every root is (sum of the first m terms)/(1 − ρ), with ρ = −2 for odd m and ρ = 0 for even m. This is the same mechanism as f(1) = Σ 1/4ⁿ = 4/3 for Ramanujan's third-order f. **Proved** in `docs/addendum-v3.md` §1 (the odd-root rule); §11 below gives the same limits as exact integral formulas and their asymptotics.
- **Expansion at q = −e^{−t}**, ψ(−e^{−t}) = Σ cₖ tᵏ/k!, with all cₖ integers:
  - ψ₁\*: `1, -1, 27, -1237, 106587, -14667181, 2954135907, -819380500597`
  - ψ₀\*: `-1, 4, -74, 3370, -286886, 39242914, -7873061654, 2177801727490`
- **Large-order behavior:** aₖ = cₖ/k! ≈ K(−1)ᵏ Γ(k+½) A^{−k−½}(1 + s₁A/(k−½) + …) with A = π²/60.
  - ψ₀\*: π·K² = (5 + 3√5)/40, s₁ = (13 − 6√5)/60, s₂ = (781 − 300√5)/7200
  - ψ₁\*: π·K² = (3√5 − 5)/40, s₁ = −(23 + 6√5)/60, s₂ = (1141 + 420√5)/7200

  The same constants appear in R₀, R₁(−e^{−t}) ~ ρ e^{π²/(60t)}(1 − s₁t + s₂t² − …), with ρ = 2√π·|K|.

## 7. The universal mock theta identity used for two other census hits

$$H(x;q)=\sum_{n\ge0}\frac{q^{n(n-1)/2}(-q;q)_n}{(xq;q)_n\,(q/x;q)_n}=(1-x^2)\,g_2(x;q)+1-x,\qquad g_2(x;q)=\sum_{n\ge0}\frac{(-q;q)_n\,q^{n(n+1)/2}}{(x;q)_{n+1}(q/x;q)_{n+1}}$$

```latex
H(x;q)=\sum_{n\ge0}\frac{q^{n(n-1)/2}(-q;q)_n}{(xq;q)_n\,(q/x;q)_n}=(1-x^2)\,g_2(x;q)+1-x
```

Proved by telescoping (see `docs/addendum-v2.md`). With ω = e^{2πi/3}:

$$\sum_{n\ge0}\frac{q^{\binom n2}(-q;q)_n^2}{(-q^3;q^3)_n}=H(-\omega),\qquad \sum_{n\ge0}\frac{q^{\binom n2}(q^2;q^2)_n}{(q^3;q^3)_n}=H(\omega)$$

```latex
\sum_{n\ge0}\frac{q^{\binom n2}(-q;q)_n^2}{(-q^3;q^3)_n}=H(-\omega),\qquad
\sum_{n\ge0}\frac{q^{\binom n2}(q^2;q^2)_n}{(q^3;q^3)_n}=H(\omega)
```

## 8. The double sum d1151 and its identity

$$D(q)=\sum_{m,n\ge0}\frac{q^{m^2-mn+n^2}}{(-q;q)_m\,(-q;q)_n}=\vartheta_4(q)\,f(q)+4R(q),\qquad R(q)=\sum_{d>k\ge0}\frac{(-1)^{d+k+1}q^{d^2}}{(-q;q)_k}$$

```latex
D(q)=\sum_{m,n\ge0}\frac{q^{m^2-mn+n^2}}{(-q;q)_m\,(-q;q)_n}=\vartheta_4(q)\,f(q)+4R(q),\qquad
R(q)=\sum_{d>k\ge0}\frac{(-1)^{d+k+1}q^{d^2}}{(-q;q)_k}
```

- Unicode: D(q) = Σ_{m,n≥0} q^{m²−mn+n²}/((−q;q)ₘ(−q;q)ₙ) = θ₄(q)·f(q) + 4R(q), R(q) = Σ_{d>k≥0} (−1)^{d+k+1} q^{d²}/(−q;q)ₖ
- Plain: `D(q) = theta4(q)*f(q) + 4*R(q)`, `R(q) = sum_{d>k>=0} (-1)^(d+k+1) q^(d^2) / (-q;q)_k`
- Here f(q) = Σ q^{n²}/(−q;q)ₙ² is Ramanujan's third-order mock theta function and θ₄(q) = Σ_{n∈ℤ} (−1)ⁿq^{n²}. **Proved** (`d1151-identity.md`).

Equivalent forms of R:

$$R(q)=\sum_{d\ge1}(-1)^{d-1}q^{d^2}\sum_{k=0}^{d-1}\frac{(-1)^k}{(-q;q)_k}=\sum_{i\ge1}(-1)^{i-1}q^{i^2}\,\Phi(q^{2i}),\qquad \Phi(z)=\sum_{n\ge0}\frac{z^n q^{n^2}}{(-q;q)_n}$$

```latex
R(q)=\sum_{d\ge1}(-1)^{d-1}q^{d^2}\sum_{k=0}^{d-1}\frac{(-1)^k}{(-q;q)_k}=\sum_{i\ge1}(-1)^{i-1}q^{i^2}\,\Phi(q^{2i}),\qquad
\Phi(z)=\sum_{n\ge0}\frac{z^n q^{n^2}}{(-q;q)_n}
```

The key lemma, for the diagonal sums S_d = Σ_{n≥0} q^{n²+dn}/((−q;q)ₙ(−q;q)_{n+d}) (S₀ = f):

$$S_d(q)+S_{d+1}(q)=\frac{2}{(-q;q)_d}$$

```latex
S_d(q)+S_{d+1}(q)=\frac{2}{(-q;q)_d},\qquad S_d(q)=\sum_{n\ge0}\frac{q^{n^2+dn}}{(-q;q)_n\,(-q;q)_{n+d}}
```

General z (proved), with the rank generating function ρ(z) = Σ q^{n²}/((zq;q)ₙ(q/z;q)ₙ):

$$\sum_{m,n\ge0}\frac{q^{m^2-mn+n^2}}{(zq;q)_m\,(q/z;q)_n}=\rho(z)\sum_{d\in\mathbb Z}z^d q^{d^2}-\Big(1-\frac1z\Big)\sum_{d>k\ge0}\frac{z^{d-k}q^{d^2}}{(zq;q)_k}-(1-z)\sum_{d>k\ge0}\frac{z^{k-d}q^{d^2}}{(q/z;q)_k}$$

```latex
\sum_{m,n\ge0}\frac{q^{m^2-mn+n^2}}{(zq;q)_m\,(q/z;q)_n}=\rho(z)\sum_{d\in\mathbb Z}z^d q^{d^2}
-\Big(1-\frac1z\Big)\sum_{d>k\ge0}\frac{z^{d-k}q^{d^2}}{(zq;q)_k}-(1-z)\sum_{d>k\ge0}\frac{z^{k-d}q^{d^2}}{(q/z;q)_k}
```

- z = −1 gives the identity for D; z = 1 gives the Durfee-rectangle identity Σ q^{n²+dn}/((q)ₙ(q)_{n+d}) = 1/(q;q)∞.
- Keywords: A2 quadratic form, sum of tails, theta tails, rank generating function, Durfee rectangle, telescoping.

Coefficients (q⁰ … q⁴⁰):
- D: `1, 3, -4, 7, -7, 7, -11, 19, -22, 22, -30, 42, -45, 49, -65, 80, -87, 95, -111, 138, -157, 164, -196, 239, -257, 280, -324, 374, -416, 447, -510, 596, -649, 695, -786, 903, -979, 1049, -1187, 1335, -1454`
- R: `0, 1, 0, 0, 0, -1, 1, -1, 1, 0, 1, -2, 2, -1, 1, -2, 2, -2, 2, -4, 4, -2, 3, -5, 4, -2, 4, -6, 6, -5, 5, -7, 8, -5, 5, -9, 8, -7, 9, -11, 11`

Verified by `census2/d1151/verify_R.py`.

## 9. Identities (III)–(VII) and the dual series

**(III), (IV)** (proved, `docs/addendum-v3.md` §4), with F₀\* = Σ(−1)ⁿq^{2n²}/(q;q²)ₙ, F₁\* = Σ(−1)ⁿq^{2n²+2n}/(q;q²)ₙ₊₁, X(q) = Σ(−1)ⁿq^{n²}/(q;q)₂ₙ, Y(q) = Σ(−1)ⁿq^{(n+1)²}/(q;q)₂ₙ₊₁:

$$F_0^{\ast}(q)+\sum_{n\ge1}q^{n^2}(q;q^2)_n=\vartheta_4(q^2)\,(-q;q)_\infty\,X(q^2),\qquad qF_1^{\ast}(q)-\sum_{n\ge0}q^{(n+1)^2}(q;q^2)_n=\vartheta_4(q^2)\,(-q;q)_\infty\,Y(q^2)$$

```latex
F_0^{*}(q)+\sum_{n\ge1}q^{n^2}(q;q^2)_n=\vartheta_4(q^2)\,(-q;q)_\infty\,X(q^2),\qquad
qF_1^{*}(q)-\sum_{n\ge0}q^{(n+1)^2}(q;q^2)_n=\vartheta_4(q^2)\,(-q;q)_\infty\,Y(q^2)
```

**The dual series of ψ\*** and **(V), (VI)** (proved, `docs/addendum-v4.md` §1). With ψ(q) = Σ_{n≥0} q^{n(n+1)/2}, P₊(q) = Σ q^{2k²+k}/(q²;q²)ₖ, P₋(q) = Σ q^{2k²−k}/(q²;q²)ₖ:

$$W(q)=\sum_{j\ge0}(-1)^j q^{j^2+j}(q;q^2)_j,\qquad W_3(q)=\sum_{j\ge0}(-1)^j q^{j^2+3j}(q;q^2)_j$$

$$W(q)+\sum_{n\ge1}\frac{q^{2n^2-n}}{(q;q^2)_n}=\psi(q)\,P_+(q),\qquad qW_3(q)+\sum_{n\ge1}\frac{q^{2n^2-3n+1}}{(q;q^2)_n}=\psi(q)\,P_-(q)$$

```latex
W(q)=\sum_{j\ge0}(-1)^j q^{j^2+j}(q;q^2)_j,\qquad W_3(q)=\sum_{j\ge0}(-1)^j q^{j^2+3j}(q;q^2)_j
W(q)+\sum_{n\ge1}\frac{q^{2n^2-n}}{(q;q^2)_n}=\psi(q)\,P_+(q),\qquad
qW_3(q)+\sum_{n\ge1}\frac{q^{2n^2-3n+1}}{(q;q^2)_n}=\psi(q)\,P_-(q)
```

- Plain: `W(q) = sum_{j>=0} (-1)^j q^(j^2+j) (q;q^2)_j`, `W3(q) = sum_{j>=0} (-1)^j q^(j^2+3j) (q;q^2)_j`
- Mathematica: `W[q_] := Sum[(-1)^j q^(j^2 + j) QPochhammer[q, q^2, j], {j, 0, Infinity}]`
- The sums on the left are F_a(q) = Σ_{n≥1} q^{2n²−n}/(q;q²)ₙ and F_b(q) = Σ_{n≥1} q^{2n²−3n+1}/(q;q²)ₙ, the dual series of f₁\*, f₀\*. W = G(−q) and W₃ = G(−q³) in the family G(z) = Σ zʲ q^{j²} (q;q²)ⱼ.

**The duals of F₀\*, F₁\* and (VII)** (`docs/addendum-v8.md`; (VII) is exact to r³²⁹ but not proved):

$$\mathcal J(r)=\sum_{m\ge0}\frac{r^{2m^2+m}}{(-r^2;r^2)_m},\qquad \mathcal K(r)=2\sum_{n\ge0}r^{(n+1)^2}(-r^2;r^2)_n,\qquad \mathcal J(r)+\mathcal K(r)=\psi(r)\,(r^2;r^4)_\infty\sum_{m\ge0}\frac{r^{m(m+1)/2}}{(r^2;r^2)_m}$$

```latex
\mathcal J(r)=\sum_{m\ge0}\frac{r^{2m^2+m}}{(-r^2;r^2)_m},\qquad
\mathcal K(r)=2\sum_{n\ge0}r^{(n+1)^2}(-r^2;r^2)_n,\qquad
\mathcal J(r)+\mathcal K(r)=\psi(r)\,(r^2;r^4)_\infty\sum_{m\ge0}\frac{r^{m(m+1)/2}}{(r^2;r^2)_m}
```

## 10. Faddeev's quantum dilogarithm, the state integrals and the exact formulas

Variables: Im x > 0, b = √x, q = e^{2πix}, t = −2πix, dual nome p = e^{−πi/x}; powers q^α = e^{2πixα}, p^α = e^{−πiα/x}.

$$\Phi_b(y)=\frac{(-q^{1/2}e^{2\pi by};q)_\infty}{(-p\,e^{2\pi y/b};p^2)_\infty},\qquad I(u)=\int_L\frac{e^{-\pi iy^2+2\pi uy}}{\Phi_b(y)}\,dy,\qquad J(u)=\int_L\frac{e^{-\pi iy^2+2\pi uy}}{\Phi_b(y)\,(1+e^{2\pi y/b})}\,dy$$

$$P=\frac{e^{\pi^2/(12t)+t/6}}{\sqrt2\,b}=\frac{(-q;q)_\infty}{b\,q^{1/8}\,(p;p^2)_\infty},\qquad L=e^{-i\delta}\mathbb R,\ \ 0<\delta<90^\circ-|\arg b|$$

```latex
\Phi_b(y)=\frac{(-q^{1/2}e^{2\pi by};q)_\infty}{(-p\,e^{2\pi y/b};p^2)_\infty},\qquad
I(u)=\int_L\frac{e^{-\pi iy^2+2\pi uy}}{\Phi_b(y)}\,dy,\qquad
J(u)=\int_L\frac{e^{-\pi iy^2+2\pi uy}}{\Phi_b(y)\,(1+e^{2\pi y/b})}\,dy
P=\frac{e^{\pi^2/(12t)+t/6}}{\sqrt2\,b}=\frac{(-q;q)_\infty}{b\,q^{1/8}\,(p;p^2)_\infty}
```

- For real b, Φ_b is given by Faddeev's integral, log Φ_b(y) = ∫_{ℝ+i0} e^{−2iyξ}/(4ξ·sinh(bξ)·sinh(ξ/b)) dξ, and |Φ_b(y)| = 1 for real y (Lemma U, `docs/addendum-v23.md`).
- Mathematica (product form, Im x > 0):
  `phib[y_, x_] := QPochhammer[-Exp[I Pi x] Exp[2 Pi Sqrt[x] y], Exp[2 Pi I x]]/QPochhammer[-Exp[-I Pi/x] Exp[2 Pi y/Sqrt[x]], Exp[-2 Pi I/x]]`
  `stateI[u_, x_, d_: 10 Degree] := Exp[-I d] NIntegrate[Exp[-I Pi (r Exp[-I d])^2 + 2 Pi u r Exp[-I d]]/phib[r Exp[-I d], x], {r, -Infinity, Infinity}]`

**The exact formulas** (proved, `docs/addendum-v12.md`; Theorem I of `psi-star-theorems-v1.md`):

$$\psi_1^{\ast}(q)=P\Big[J\big(\tfrac1{2b}\big)+I\big(\tfrac1{2b}\big)\,p\,W_3(p)-I\big(-\tfrac1{2b}\big)\,W(p)\Big],\qquad f_1^{\ast}(q)=2P\Big[J\big(\tfrac1{2b}\big)-I\big(\tfrac1{2b}\big)F_b(p)+I\big(-\tfrac1{2b}\big)F_a(p)\Big]$$

$$\psi_0^{\ast}(q)=q^{1/2}P\Big[J\big(b+\tfrac1{2b}\big)+I\big(b+\tfrac1{2b}\big)\,p\,W_3(p)-I\big(b-\tfrac1{2b}\big)\,W(p)\Big],\qquad f_0^{\ast}(q)=2q^{1/2}P\Big[J\big(b+\tfrac1{2b}\big)-I\big(b+\tfrac1{2b}\big)F_b(p)+I\big(b-\tfrac1{2b}\big)F_a(p)\Big]$$

```latex
\psi_1^{*}(q)=P\Big[J\big(\tfrac1{2b}\big)+I\big(\tfrac1{2b}\big)\,p\,W_3(p)-I\big(-\tfrac1{2b}\big)\,W(p)\Big]
f_1^{*}(q)=2P\Big[J\big(\tfrac1{2b}\big)-I\big(\tfrac1{2b}\big)F_b(p)+I\big(-\tfrac1{2b}\big)F_a(p)\Big]
\psi_0^{*}(q)=q^{1/2}P\Big[J\big(b+\tfrac1{2b}\big)+I\big(b+\tfrac1{2b}\big)\,p\,W_3(p)-I\big(b-\tfrac1{2b}\big)\,W(p)\Big]
f_0^{*}(q)=2q^{1/2}P\Big[J\big(b+\tfrac1{2b}\big)-I\big(b+\tfrac1{2b}\big)F_b(p)+I\big(b-\tfrac1{2b}\big)F_a(p)\Big]
```

- Unicode: ψ₁\*(q) = P·[J(1/(2b)) + I(1/(2b))·p·W₃(p) − I(−1/(2b))·W(p)]
- The state integrals factorise into q- and p-series (Theorem C of v12): I(u) = c·e^{−πiu²/2}·D(u) with c = 2^{−1/2}·e^{−πi/4}·q^{−1/24}·p^{1/12} and the double series D(u) of `docs/addendum-v12.md` §1; `psi_star/v12/check_thmAC.py` checks it.
- **P·J is the Borel sum** of the perturbative series Φ₁(t) = Σ aₖtᵏ, a₀ = 1/3, a₁ = 1/9, a₂ = 5/162, in the directions −120° … −80° (proved, `docs/addendum-v17.md`; 301 exact coefficients in `psi_star/v13/phi_exact.json`). The Borel transform has its nearest singularities at S = π²/3 ∓ i·Cl₂(π/3).

## 11. The laws near q = 1 and at roots of unity

Constants: Cl₂ = Cl₂(π/3) = 1.0149416064…, S̄ = π²/3 + i·Cl₂, A = (π/√3)^{1/2}, s̄₁ = 1/4 + i√3/54.

**Forward law inside the disk** (proved for 0 < arg x ≤ 10°, `docs/addendum-v19.md` Theorem 8):

$$\psi_1^{\ast}(q)=\mathcal S[\Phi_1](t)+A\,t^{-1/2}\Big[e^{5\pi i/12}e^{-\bar S/t}\,W(p)\,(1+O(t))+e^{-5\pi i/12}e^{-S/t}\,p\,W_3(p)\,(1+O(t))\Big]$$

```latex
\psi_1^{*}(q)=\mathcal S[\Phi_1](t)+A\,t^{-1/2}\Big[e^{5\pi i/12}e^{-\bar S/t}\,W(p)\,(1+O(t))+e^{-5\pi i/12}e^{-S/t}\,p\,W_3(p)\,(1+O(t))\Big]
```

**At the roots of unity e^{2πi/X}** (proved, `docs/addendum-v20.md`, `v21.md`, `v25.md`; x = 1/X, integrals with Faddeev's Φ_b at real b):

$$\psi_1^{\ast}(e^{2\pi i/X})=P\Big[J\big(\tfrac1{2b}\big)+I\big(\tfrac1{2b}\big)-I\big(-\tfrac1{2b}\big)\Big]\ (X\ \text{even}),\qquad \lim_{r\to1^-}\psi_1^{\ast}(re^{2\pi i/X})=P\Big[J\big(\tfrac1{2b}\big)-\tfrac13 I\big(\tfrac1{2b}\big)-\tfrac13 I\big(-\tfrac1{2b}\big)\Big]\ (X\ \text{odd})$$

$$\psi_1^{\ast}(e^{2\pi i/X})=w_X\,3^{-1/4}\Big(\frac X2\Big)^{1/2}e^{\mathrm{Cl}_2(\pi/3)X/2\pi}\,e^{2\pi i/3-i\pi X/6}\Big(1-\frac{2\pi i\,\bar s_1}{X}+O(X^{-2})\Big),\qquad w_X=\begin{cases}1&X\ \text{even}\\ 1/3&X\ \text{odd}\end{cases}$$

```latex
\psi_1^{*}(e^{2\pi i/X})=w_X\,3^{-1/4}\Big(\frac X2\Big)^{1/2}e^{\mathrm{Cl}_2(\pi/3)X/2\pi}\,e^{2\pi i/3-i\pi X/6}\Big(1-\frac{2\pi i\,\bar s_1}{X}+O(X^{-2})\Big)
```

- The leading term is proved (Theorems 10, 12); the coefficient −2πi·s̄₁ of 1/X is numerical (`psi_star/v20/check_thm10.py`, `v21/check_thm12.py`, Richardson extrapolation to 10 digits). For ψ₀\* multiply by e^{2πi/3} and replace s̄₁ by −1/12 + i√3/54.
- Plain: `psi1*(e^(2 pi i/X)) ~ w_X 3^(-1/4) sqrt(X/2) exp(Cl2(pi/3) X/(2 pi)) exp(2 pi i/3 - i pi X/6)`
- Figure-eight knot: |ψ₁\*(e^{2πi/X})| / ⟨4₁⟩_{X/2} → 2/X for even X, with ⟨4₁⟩_N = Σ_{k<N} Π_{j≤k} |1 − e^{2πij/N}|² the Kashaev invariant.

## 12. The z-family: exact formula and the volume rule

Notation: z = e^{2πia}, 0 < a < 1; z̃ = e^{2πia/x}, T = p·z̃ = e^{πi(2a−1)/x}; the dual series 𝒲(T; p) = Σ_{k≥0} (−T)ᵏ p^{k²+k} (p;p²)ₖ (so 𝒲(1; p) = W(p), 𝒲(p²; p) = W₃(p)).

**Theorem Z** (proved, `docs/addendum-v18.md` §1; strip form):

$$F(z;q)=i\,P\,e^{-\pi ia}\Big[J\big(\tfrac ab\big)-I\big(\tfrac{a-1}b\big)\,\mathcal W(T)+p\,T\,I\big(\tfrac ab\big)\,\mathcal W(p^2T)\Big]$$

```latex
F(z;q)=i\,P\,e^{-\pi ia}\Big[J\big(\tfrac ab\big)-I\big(\tfrac{a-1}b\big)\,\mathcal W(T)+p\,T\,I\big(\tfrac ab\big)\,\mathcal W(p^2T)\Big]
```

- Mathematica: `calW[T_, p_] := Sum[(-T)^k p^(k^2 + k) QPochhammer[p, p^2, k], {k, 0, Infinity}]`

**Shape, volume and the predicted growth** (`docs/addendum-v18.md` §2, `v22.md`, `v23.md`):

$$W^2+W=\frac1z,\ \ \operatorname{Im}W>0,\qquad \mathrm{Vol}(z)=|D(-W)|,\qquad D(u)=\operatorname{Im}\mathrm{Li}_2(u)+\arg(1-u)\,\log|u|$$

$$F(z;e^{2\pi i/X})=(1+W)\,W^{1/2}\Big(\frac{\pi}{t(1+2W)}\Big)^{1/2}e^{(V(w_0)+\pi^2/12)/t}\,(1+O(1/X)),\qquad t=-\frac{2\pi i}{X},\ \ w_0=\operatorname{Log}W,\ \ V(w)=-\frac{w^2}2-2\pi i(a-1)w+\mathrm{Li}_2(-e^w)$$

```latex
W^2+W=\frac1z,\quad \operatorname{Im}W>0,\qquad \mathrm{Vol}(z)=|D(-W)|,\qquad D(u)=\operatorname{Im}\mathrm{Li}_2(u)+\arg(1-u)\log|u|
F(z;e^{2\pi i/X})=(1+W)\,W^{1/2}\Big(\frac{\pi}{t(1+2W)}\Big)^{1/2}e^{(V(w_0)+\pi^2/12)/t}\,(1+O(1/X)),\qquad
V(w)=-\frac{w^2}2-2\pi i(a-1)w+\mathrm{Li}_2(-e^w)
```

- Hence lim (1/X)·log|F(z; e^{2πi/X})| = Vol(z)/2π for even X and every z ≠ 1 on the unit circle (proved, Theorem 15 of v23; the error terms written out in Theorem 15′ of v27).
- Values: Vol(−1) = 1.0149416 (= Cl₂(π/3)), Vol(i) = 0.6967011, Vol(ω) = 0.8651292, Vol(−ω) = 0.4866176, ω = e^{2πi/3}. As a → 0: Vol ≈ 2π·ln φ·|a|, φ the golden ratio.
- Mathematica: `vol[a_] := Module[{W = Select[w /. NSolve[w^2 + w == Exp[-2 Pi I a], w], Im[#] > 0 &][[1]]}, Abs[Im[PolyLog[2, -W]] + Arg[1 + W] Log[Abs[W]]]]`
- At even X the finite sum equals the three integrals exactly: Σ_{n<X/2} zⁿ q₀^{n(n+1)/2} (−q₀;q₀)ₙ = i·P·e^{−πia}·[J(a/b) − I((a−1)/b) + z^X·I(a/b)], q₀ = e^{2πi/X} (Theorem 14 of v22). The kernel term is bounded: |P·J(a/b)| ≤ 1/(2√2·sin πa) at real x (Lemma J of v23).

## 13. z = 1: the fifth-order mock theta functions and the Mordell integral

At a = 0 the dual series of §12 become Ramanujan's fifth-order functions (proved, `docs/addendum-v18.md` §3.1):

$$\mathcal W(1/p;p)=\varphi_0(-p),\qquad p\,\mathcal W(p;p)=-\varphi_1(-p),\qquad \varphi_0(q)=\sum_{n\ge0}q^{n^2}(-q;q^2)_n,\quad \varphi_1(q)=\sum_{n\ge0}q^{(n+1)^2}(-q;q^2)_n$$

```latex
\mathcal W(1/p;p)=\varphi_0(-p),\qquad p\,\mathcal W(p;p)=-\varphi_1(-p),\qquad
\varphi_0(q)=\sum_{n\ge0}q^{n^2}(-q;q^2)_n,\quad \varphi_1(q)=\sum_{n\ge0}q^{(n+1)^2}(-q;q^2)_n
```

**The state integrals at z = 1, q** are S-matrix entries times powers (proved given classical identities, `docs/addendum-v24.md` Proposition M), with s₁ = sin(π/5), s₂ = sin(2π/5):

$$I(0)=e^{-\pi i/4}\tfrac{2}{\sqrt5}\,s_2\,q^{-1/60}p^{1/30},\quad I\big(\tfrac1b\big)=e^{-\pi i/4}\tfrac{2}{\sqrt5}\,s_1\,q^{-1/60}p^{19/30},\quad I(b)=e^{-\pi i/4}\tfrac{2}{\sqrt5}\,s_1\,q^{-19/60}p^{1/30},\quad I\big(b+\tfrac1b\big)=-e^{-\pi i/4}\tfrac{2}{\sqrt5}\,s_2\,q^{-19/60}p^{19/30}$$

```latex
I(0)=e^{-\pi i/4}\tfrac{2}{\sqrt5}\,s_2\,q^{-1/60}p^{1/30},\quad I\big(\tfrac1b\big)=e^{-\pi i/4}\tfrac{2}{\sqrt5}\,s_1\,q^{-1/60}p^{19/30},\quad
I(b)=e^{-\pi i/4}\tfrac{2}{\sqrt5}\,s_1\,q^{-19/60}p^{1/30},\quad I\big(b+\tfrac1b\big)=-e^{-\pi i/4}\tfrac{2}{\sqrt5}\,s_2\,q^{-19/60}p^{19/30}
```

**The transformation laws** (Theorems 16, 17 of v24; c± = (π(5 ± √5)/(10t))^{1/2}; F₀, F₁ as in §9 without the sign, F₀(q) = Σ q^{2n²}/(q;q²)ₙ, F₁(q) = Σ q^{2n²+2n}/(q;q²)ₙ₊₁):

$$q^{11/60}f_1(q)=M_1+2c_+\,p^{-1/120}\big(F_0(p)-1\big)-2c_-\,p^{71/120}F_1(p),\qquad 2q^{11/60}\psi_1(q)=-M_1+2c_+\,p^{-1/120}\varphi_0(-p)+2c_-\,p^{-49/120}\varphi_1(-p)$$

$$q^{-1/60}f_0(q)=M_0+2c_-\,p^{-1/120}\big(F_0(p)-1\big)+2c_+\,p^{71/120}F_1(p),\qquad 2q^{-1/60}\psi_0(q)=-M_0+2c_-\,p^{-1/120}\varphi_0(-p)-2c_+\,p^{-49/120}\varphi_1(-p)$$

$$M_1=2i\,q^{11/60}P\,J\big(\tfrac1b\big)=\Big(\frac{60t}{\pi}\Big)^{1/2}L\big(\tfrac25,10t\big),\quad M_0=2i\,q^{29/60}P\,J\big(b+\tfrac1b\big)=\Big(\frac{60t}{\pi}\Big)^{1/2}L\big(\tfrac15,10t\big),\quad L(r,\alpha)=\int_0^\infty e^{-3\alpha x^2/2}\,\frac{\cosh((3r-2)\alpha x)+\cosh((3r-1)\alpha x)}{\cosh(3\alpha x/2)}\,dx$$

```latex
q^{11/60}f_1(q)=M_1+2c_+\,p^{-1/120}\big(F_0(p)-1\big)-2c_-\,p^{71/120}F_1(p)
2q^{11/60}\psi_1(q)=-M_1+2c_+\,p^{-1/120}\varphi_0(-p)+2c_-\,p^{-49/120}\varphi_1(-p)
q^{-1/60}f_0(q)=M_0+2c_-\,p^{-1/120}\big(F_0(p)-1\big)+2c_+\,p^{71/120}F_1(p)
2q^{-1/60}\psi_0(q)=-M_0+2c_-\,p^{-1/120}\varphi_0(-p)-2c_+\,p^{-49/120}\varphi_1(-p)
M_1=\Big(\frac{60t}{\pi}\Big)^{1/2}L\big(\tfrac25,10t\big),\quad M_0=\Big(\frac{60t}{\pi}\Big)^{1/2}L\big(\tfrac15,10t\big),\quad
L(r,\alpha)=\int_0^\infty e^{-3\alpha x^2/2}\,\frac{\cosh((3r-2)\alpha x)+\cosh((3r-1)\alpha x)}{\cosh(3\alpha x/2)}\,dx
```

- With M₀, M₁ written as Mordell integrals these are the classical laws (Zwegers, thesis 2002, Propositions 4.10, 4.13; Gordon–McIntosh 2012). The identification of the kernel integrals with L is Theorem 17 of v24 (proved given the known laws and two series identities (R3), (R4) that hold to q¹⁰⁰⁰). The four laws in the Mordell form are checked by `psi_star/formulas_check.py`, independently of any state integral.
- Mathematica: `mordell[r_, a_] := NIntegrate[Exp[-3 a x^2/2] (Cosh[(3 r - 2) a x] + Cosh[(3 r - 1) a x])/Cosh[3 a x/2], {x, 0, Infinity}]` (for Re a > 0).
