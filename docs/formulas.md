# Formula sheet: every object in several notations

The same q-series can be written in many equivalent ways, which makes it hard to tell whether a formula has been published before. This sheet lists each object studied here in several forms (rendered LaTeX, raw LaTeX source, plain text and Mathematica syntax), together with equivalent rewritings and the opening coefficients. **Every equivalent form below is checked numerically to 50 digits by [`psi_star/formulas_check.py`](../psi_star/formulas_check.py).**

Notation: $(a;q)_n=\prod_{k=0}^{n-1}(1-aq^k)$, so $(-q;q)_n=\prod_{k=1}^{n}(1+q^k)=\dfrac{(q^2;q^2)_n}{(q;q)_n}$, also written $(-q)_n$. Note that $q^{n(n+1)/2}=q^{\binom{n+1}{2}}$.

Keywords: sign-twisted, alternating, (−1)ⁿ inserted, z = −1 specialization, fifth-order mock theta functions ψ₀, ψ₁, Rogers–Ramanujan, entire Rogers–Ramanujan function, bilateral series, strange identity, quantum modular form, resurgence, Stokes constants, q = −1 asymptotics.

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
