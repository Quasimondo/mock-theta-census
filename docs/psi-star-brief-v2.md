# Brief v2: what kind of object are ψ₀\*, ψ₁\*?

**For:** an agent continuing the investigation with Mario Klingemann. It follows on from `addendum-v2.md` (an earlier session that found identities (I), (II) and the resurgence at q = −1). Its §5 "next steps" were never done. This brief adds what has been learned since and asks for those steps, plus a few new ones.
**Repo:** https://github.com/Quasimondo/mock-theta-census. Read `README.md`, `docs/addendum-v2.md`, `docs/twists.md`; everything below is in the repo.
**Labels:** **proved** = proof given and machine-checked; **verified** = exact or high-precision computation; **inference** = a reading of the evidence. Please use the same labels.

## 1. The objects

- ψ₁\*(q) = Σₙ≥0 (−1)ⁿ q^{n(n+1)/2} (−q;q)ₙ, ψ₀\*(q) = Σₙ≥0 (−1)ⁿ q^{(n+1)(n+2)/2} (−q;q)ₙ: the fifth-order ψ₁, ψ₀ with a sign twist.
- F(z) = Σ zⁿ q^{n(n+1)/2} (−q;q)ₙ, with F(z) = 1 + zqF(zq) + zq²F(zq²); ψ₁\* = F(−1), ψ₀\* = qF(−q), F(−1/q) = 1 − ψ₁\* − ψ₀\*.
- f₀\* = Σ(−1)ⁿq^{n²}/(−q)ₙ, f₁\* = Σ(−1)ⁿq^{n²+n}/(−q)ₙ, R₀ = Σ(−1)ⁿq^{n²}/(q)ₙ, R₁ = Σ(−1)ⁿq^{n²+n}/(q)ₙ, θ₄ = Σ(−1)ⁿq^{n²}.

## 2. What is established

- **Proved** (addendum-v2 §2): (I) 2ψ₀\* = f₀\* − θ₄R₀, (II) 2ψ₁\* = f₁\* + θ₄R₁.
- **Verified** (addendum-v2 §3): the perturbative series of ψ\* at q = −1 has integer coefficients; its large-order behaviour is governed by R₀, R₁ (Borel singularity at −π²/60, constants in ℚ(√5)).
- **Verified:** bounded at all 19 standard roots (orders 1–8, 10, 12, 14) and at roots of order 9, 11, 13, 15. Radial limits:
  - at even order m the series terminates;
  - at odd order m the limit is one third of the first m terms, to ~10⁻¹³ at all odd m ≤ 15 (not proved).
- **Verified:** under x ↦ x/(2x+1) (Γ₀(2)), ψ₁\*, ψ₀\* transform together to ~10⁻⁹, at weight ≈ 0.508 (`psi-star-brief.md` §3). Addendum-v2 §4 reads this as an effective fit to a matrix cocycle of Garoufalidis–Zagier/Nahm type, not a classical quantum modular form.

## 3. New since addendum-v2 (negative results; don't repeat)

- **Not Hecke-type over a simple product.** ψ₀\*, ψ₁\* times any (q;q)∞^a (q²;q²)∞^b (q⁴;q⁴)∞^c, with −2 ≤ a, b, c ≤ 2, stays dense: at least 386 of the first 400 coefficients are nonzero. Controls: f₀(q;q)∞ has 160 nonzero, ψ₁(q²;q²)∞ has 143 (`census3/hecke/`).
- **Not in the σ family.** The (−1)ⁿ twist of the third-order ψ turned out to be −σ\*/2 (Andrews–Dyson–Hickerson). ψ\* is different.
- **No q → 1/q partner.** σ(ζ) = −σ\*(ζ̄) holds at every root. No relation ψ\*(ζ) = aζ^e X(ζ̄) + cζ^e (other ψ\*)(ζ̄) + b exists:
  - X ∈ {ψ₀\*, ψ₁\*, σ, σ\*, 1}, e ∈ j/120;
  - X among the 1,357 census functions bounded at all roots (1,110 testable, 46 roots of order ≤ 12).
  - The search does recover σ ↔ σ\* (`census3/qm/`).
- **F(z) at other z** (`census3/`). At a root of odd order m, F(z) grows when z^m is within about 1/9 of a turn of 1. It grows with an irrational κ in general, and with κ = 1/240 and an inexact law (class R) when z^m = 1. At even order F(z) is always bounded. So F(z) is bounded everywhere only for z = −1, ±i, −ω, −ω̄.
- **Other bounded dense twists** with the same status as ψ\* but no identities yet (`docs/twists.md`): F₀\*, F₁\* (fifth order, Σ(−1)ⁿq^{2n²}/(q;q²)ₙ etc.), T₀\* (eighth order), F(i), F(−ω). For F₀\*, F₁\* no rational linear relation with ψ\*, f\*, φ\* was found to q²⁰⁰.

## 4. Questions, in priority order

1. **The matrix law** (addendum-v2 §5.2). Test whether (1, ψ₀\*, ψ₁\*) together with homogeneous solutions of the q-difference equation for F at z = −1 transform as a 3×3 matrix-valued quantum modular object. Candidates for those solutions: θ₄R₀, θ₄R₁, or f₀\*, f₁\*. Use x ↦ x+1 and x ↦ x/(2x+1); the cocycle need only be smooth, not of the form ε(cx+d)^{−k}.
   - **Data:** exact limits at even roots from the terminating sums; numerical limits at odd roots for m ≤ 15 are in `census3/qm/limits.json`.
   - **A clean success looks like:** rational exponents, multiplier entries that are roots of unity or simple algebraic numbers, and a residual that falls with denominator size the way the Kontsevich control does in `psi-star-brief.md` §3, instead of plateauing.
2. **Resurgence at other roots** (addendum-v2 §5.1). Repeat the q = −1 analysis at q = i and q → 1. Predict the Borel singularities from the twisted Nahm/R-equation and compare.
3. **Prove the odd-root rule.** At a root ζ of odd order m, ψ₁\*(ζ) = (1/3)·Σ_{n<m} (−1)ⁿ ζ^{n(n+1)/2} (−ζ;ζ)ₙ, and similarly for ψ₀\*. The terms repeat with factor −2 every m steps, so the sum is formally a geometric series with ratio −2. A proof would say why Abel summation (1/(1−(−2)) = 1/3) gives the radial limit.
4. **Analogues of (I), (II) for F₀\*, F₁\*, T₀\*.** `twists.md` sketches why the bilateral route needs a base-q² version for F₀. A proved identity for any of these would show that ψ\* is a pattern rather than a one-off.
5. **Why the F(z) threshold sits at about 1/9 of a turn.** The growth near q = 1 for z = e^{2πia} should come from a saddle point of a dilogarithm-type function V(z). Find V and the value of a where its real part crosses 0. Check it against the data in `census3/zfam_results.jsonl`: κ changes sign between a = 1/9 and a = 1/7.
6. **Literature**, only if cheap: R(−1) = Σ(−1)ⁿq^{n²}/(q)ₙ in the Garoufalidis–Zagier/Wheeler line, and Andrews on the "entire Rogers–Ramanujan function". Also whether (I), (II) or ψ\* appear in the bilateral-series papers (Bajpai–Kimport–Liang–Ma–Ricci; Mc Laughlin).

## 5. Ground rules

- **Verify every claimed identity** as q-series to at least q³⁰⁰ in exact arithmetic, and say so. Verify every numerical constant to at least 15 digits.
- **Give claims as you computed them,** with the label and the script that produced them. A negative result with a clear test is as useful as a positive one.
- **Deliverables:**
  - `psi-star-addendum-v3.md`, in the style of `addendum-v2.md`: summary first, then one section per question attempted.
  - A zip of the scripts. Python with mpmath, numpy, python-flint; the repo's `census/catalog.py` and `census/probe.py` evaluate any Eulerian series with ball arithmetic.

## 6. Reproduce the basic objects (exact integer series)

```python
N = 300
def qpoch(a, step, sign, n=None, inv=False):        # (sign*q^a; q^step)_n as an integer series
    s = [0]*N; s[0] = 1; k, e = 0, a
    while e < N and (n is None or k < n):
        if inv:
            for i in range(e, N): s[i] -= sign*s[i-e]
        else:
            for i in range(N-1, e-1, -1): s[i] += sign*s[i-e]
        k += 1; e += step
    return s
def tw(expo, fac, nmax=60):                          # sum (-1)^n q^expo(n) fac(n)
    out = [0]*N
    for n in range(nmax):
        e = expo(n)
        if e >= N: break
        for i, c in enumerate(fac(n)[:N-e]): out[i+e] += (-1)**n*c
    return out
psi1s = tw(lambda n: n*(n+1)//2, lambda n: qpoch(1, 1, +1, n))
psi0s = tw(lambda n: (n+1)*(n+2)//2, lambda n: qpoch(1, 1, +1, n))
R0 = tw(lambda n: n*n, lambda n: qpoch(1, 1, -1, n, inv=True))
```
