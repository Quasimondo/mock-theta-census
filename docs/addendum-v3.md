# Mock theta search: addendum v3 (ψ₀\*, ψ₁\*)

Session of 26 Nov 2026 (Claude Fable 5.2 with Mario Klingemann), answering `psi-star-brief-v2.md`. It builds on `addendum-v2.md`; identities (I), (II) and the q = −1 analysis are taken from there. Scripts are in `psi_star/v3/`; `run_log.txt` there is their complete output.

**Independent check** (Claude Opus, same day, own code in `psi_star/v3/check_III_IV.py` and `check_law.py`): (III), (IV) hold exactly to q³²⁹; the law of §5.3 holds for ρ = 0, 1, ½, ⅓, 3/2, 2/5 with X up to ≈ 3,200, the leading ratio tending to 1 and s̄₁ = 1/4 + i√3/54 recovered to 10⁻¹⁰ by Richardson extrapolation; ψ₀\*/ψ₁\* → e^{2πi/3}; the quoted values of W agree. The reversal identity A = B of §1(a) agrees to 10⁻²⁸ at all odd orders ≤ 15. The estimates in §1(c), (d) are a sketch and should be written out in full before the theorem is called proved in public.

Status labels: **proved** = proof given here and machine-checked; **verified** = exact or high-precision computation, no proof; **heuristic** = a formal derivation (saddle point) confirmed numerically; **inference** = my reading of the evidence.

Notation as in the brief. In addition Cl₂(π/3) = 1.0149416064096536… (Clausen function), and for a root of unity ζ of order m, q = ζe^{−t}.

---

## 0. Summary

1. **The odd-root rule is proved** (§1). At a root ζ of odd order m, ψ₁\*(ζ) = ⅓·Σ_{n<m}(−1)ⁿζ^{n(n+1)/2}(−ζ;ζ)ₙ, and likewise for ψ₀\*. The reason: by (II), ψ₁\* = ½f₁\* + ½θ₄R₁; the second term vanishes radially, f₁\* converges at ζ as a geometric series with ratio −½, and reversing the finite sum turns it into yours.
2. **The F(z) threshold is found** (§2, heuristic, verified against your data to 10 digits). Growth at odd roots is e^{V/(m²t)} with an explicit dilogarithm function V of a = arg(z^m)/2π. Re V changes sign at **a\* = 0.107365841812… = 1/9.3139 of a turn**, between 1/10 and 1/9. Your "boundary zone" at 1/9 and 1/8 is slow exponential decay, not growth.
3. **Resurgence at q = i and q → 1** (§3, verified). At q = i the Borel singularity is at −π²/240, with constants in ℚ(√5, i), and the attached series is again that of ∓½θ₄R. At q → 1 the singularities are at **π²/3 ∓ i·Cl₂(π/3)**, exactly the saddle values of §2 at z = −1, with closed-form Stokes constants.
4. **Identities for F₀\*, F₁\*** (§4, proved, checked to q³¹⁹). They have the shape of (I), (II): bilateral sum = θ₄(q²)(−q;q)∞ times a sign-twisted Slater-type sum with (q;q)₂ₙ denominators. That denominator type was outside the earlier search.
5. **A modularity law of Zagier type** (§5, verified to 6·10⁻²¹ on 56 residue classes). For X → ∞, ψ₁\*(e^{2πi/X}) grows like e^{Cl₂(π/3)·X/2π}, and the coefficient of the growing term is a **new finite sum at the dual root of unity**: W(ξ) = Σ(−1)ʲ ξ^{j²+j}(ξ;ξ²)ⱼ at ξ = e^{−πiX}. So under x ↦ 1/x the pair does not go to itself; it goes to W. This is the first row of the matrix law asked for in question 1, not yet the whole matrix.
6. **Literature** (§6): references confirmed; ψ\*, (I)–(IV) and W not found in a short search.
7. **Not done:** T₀\*, F(i), F(−ω); the third (decaying) sector of the matrix law; the converse law for W.

---

## 1. The odd-root rule (question 3): proved

**Theorem.** Let ζ be a primitive m-th root of unity, m odd. Then

  lim_{r→1⁻} ψ₁\*(rζ) = ⅓·A, A = Σ_{n=0}^{m−1} (−1)ⁿ ζ^{n(n+1)/2} (−ζ;ζ)ₙ,
  lim_{r→1⁻} ψ₀\*(rζ) = ⅓·A₀, A₀ = Σ_{n=0}^{m−1} (−1)ⁿ ζ^{(n+1)(n+2)/2} (−ζ;ζ)ₙ.

**Proof.** It uses (II): 2ψ₁\* = f₁\* + θ₄R₁ (and (I) for ψ₀\*), valid for |q| < 1.

*(a) A finite identity.* For odd m, ∏_{k=0}^{m−1}(1+ζᵏ) = 2, so (−ζ;ζ)_{m−1} = 1 and

  (−ζ;ζ)_{m−1−n} = 1/∏_{j=1}^{n}(1+ζ^{−j}) = ζ^{n(n+1)/2}/(−ζ;ζ)ₙ.

Also (m−1−n)(m−n)/2 ≡ n(n+1)/2 mod m and (−1)^{m−1−n} = (−1)ⁿ. Replacing n by m−1−n in A gives

  A = B := Σ_{n<m} (−1)ⁿ ζ^{n²+n}/(−ζ;ζ)ₙ, and in the same way A₀ = B₀ := Σ_{n<m} (−1)ⁿ ζ^{n²}/(−ζ;ζ)ₙ.

*(b) f₁\* at the root.* At q = ζ the terms of f₁\* = Σ(−1)ⁿq^{n²+n}/(−q;q)ₙ satisfy term_{n+m} = −½·term_n, because (−ζ;ζ)_{n+m} = 2(−ζ;ζ)ₙ. So the series converges there, to B/(1+½) = ⅔B.

*(c) Radial continuity of f₁\*.* For q = rζ group the factors of (−q;q)ₙ in blocks of m. Within a block rᵏ varies by a factor 1+O(1−r), and for equal radii s the block product is exactly ∏ₖ(1+sζᵏ) = 1+s^m. Hence |(−q;q)_{mJ}| ≥ ∏_{j≤J}(1+r^{jm²})·(1−c(1−r))^J. If r^{Jm²} ≥ ½ this is ≥ 1.4^J; otherwise r^{n²} ≤ 2^{−J}. Either way the n-th term is bounded by C·θ^{n/m} with θ < 1, uniformly for r near 1. By dominated convergence f₁\*(rζ) → ⅔B.

*(d) θ₄R₁ → 0.* The same block comparison for (q;q)ₙ, where the block product is 1−s^m, gives |(q;q)_{mJ}| ≥ (Q;Q)_J(1−ct)^{mJ} with Q = e^{−m²t}. Then

  |R₁(ζe^{−t})| ≤ C Σ_J Q^{J²−ηJ}/(Q;Q)_J ≤ C′ (−Q^{1/2};Q)∞ ≤ C′ exp((1+ε)/(m²t)),

using J²−ηJ ≥ J²/2 − η²/2 and Euler's Σ Q^{J²/2}/(Q;Q)_J = (−Q^{1/2};Q)∞. On the other hand θ₄ = (q;q)∞²/(q²;q²)∞ vanishes at odd-order roots: |θ₄(ζe^{−t})| ≤ C t^{−1/2} exp(−π²/(4m²t)) (standard, from the transformation of η). Since 1 < π²/4, the product tends to 0.

By (II), ψ₁\*(rζ) → ½·⅔B = ⅓A. The proof for ψ₀\* is the same with (I), f₀\*, R₀. ∎

So Abel summation with ratio −2 gives the right value because the series f₁\*, which has the same value at the root by the reversal (a), genuinely converges there with ratio −½, and ⅔·½ = ⅓.

**Corollary (proved).** ψ₁\*(ζe^{−t}) − ½f₁\*(ζe^{−t}) = O(exp(−(π²/4−1−ε)/(m²t))). The perturbative series of ψ₁\* at an odd root is therefore that of ½f₁\*, which can be computed term by term (used in §3).

**Checks** (`q3_oddroot.py`):
- A = B and A₀ = B₀ at all primitive roots of odd order ≤ 25: largest difference 7.6·10⁻³⁷.
- Identity (II) at complex q near seven odd roots: residual ≤ 4·10⁻⁵⁵.
- f₁\* → ⅔A and ψ₁\* → ⅓A with error proportional to t.
- **Verified, not proved:** the actual size of θ₄R₁ is e^{−π²/(3m²t)} (m²t·log|θ₄R₁| = −3.13 … −3.18 at m²t = 1/32, tending to −π²/3 = −3.29), i.e. R₁ itself decays like e^{−π²/(12m²t)} at odd roots.

---

## 2. Growth of F(z) and the threshold (question 5)

**Saddle point (heuristic).** For z = e^{2πia} and q = e^{−t}, with x = nt, the n-th term of F(z) = Σ zⁿq^{n(n+1)/2}(−q;q)ₙ is

  exp(V(x)/t)·e^{−x/2}·√((1+e^{−x})/2)·(1+O(t)), V(x) = 2πi·a·x − x²/2 + Li₂(−e^{−x}) + π²/12.

The saddle is at w² + w = e^{−2πia}, w = e^{−x} (the root continued from w = (√5−1)/2 at a = 0), and

  F(e^{2πia}; e^{−t}) ≈ (1+w)·√(πw/((1+2w)t))·e^{V/t} + (bounded part).

At a root of odd order m the sum over n splits into residues mod m and the same analysis applies with a = arg(z^m)/2π and t replaced by m²t. In your normalisation κ = Re V/(4π²), depending only on z^m, as you observed.

**Results.**

- **a = 0:** V = π²/60, κ = 1/240. ✓ (your class R value)
- **Threshold:** Re V(a\*) = 0 at **a\* = 0.10736584181229166611**, i.e. 1/9.313949 of a turn.
- **a = ½ (z^m = −1):** V = −π²/3 ± i·Cl₂(π/3). This is the non-perturbative sector of ψ₁\* at odd roots (§3, §5).
- **Bounded case (verified to 5–12 digits at five (z, ζ)):** when Re V < 0 the radial limit is A_z/(1 − 2z^m), with A_z the sum of the first m terms. For z = −1 this is the ⅓ rule.

**Comparison with `census3/zfam_results.jsonl`** (`q5_threshold.py`; all rows of F at odd m, grouped by a):

| a | Predicted κ = Re V/4π² | Fitted κ in the census | Classes in the census |
|---|---|---|---|
| 0 | 0.0041666667 | 0.0041666667 | E 6, Q 1, X 28 |
| 1/12 | 0.0016555805 | 0.0016555805 (12 rows) | X 12 |
| 1/11 | 0.0011785967 | 0.0011785967 (14 rows) | X 14 |
| 1/10 | 0.0005516302 | 0.0005516302 (2 rows) | X 2, timeout 8 |
| 1/9 | bounded, Re V = −0.0117 | −0.00039 … −0.00013 | X 8, Q 3, timeout 1 |
| 1/8 | bounded, Re V = −0.0584 | −0.00042 … −0.00025 | X 2, Q 4, timeout 8 |
| ≥ 1/7 | bounded, Re V ≤ −0.1265 | — | Q in every row (225) |

So the sign change is between 1/10 and 1/9, not between 1/9 and 1/7. At a = 1/9 the saddle term decays only like e^{−0.0117/(m²t)}, which at the classifier's range m²t ∈ [0.004, 0.02] is a factor between e^{−2.9} and e^{−0.6}. It is still visible, and the classifier reads the mixture as X with a small negative κ. That is the "boundary zone".

**Precision.** The leading term, prefactor included, was checked at q = e^{−t} for a = 0 and a = 1/12 with t = 1/800 … 1/12800: the ratio F/prediction extrapolates to 1 within 5·10⁻²² (which fixes V to better than 10⁻²⁰). Which saddle contributes is not derived; it is fixed by this check and by the table.

---

## 3. Resurgence at q = i and q → 1 (question 2)

Perturbative series a₀ + a₁t + a₂t² + … were computed to order 150 with ball arithmetic (`q2_series.py`), and the large-order behaviour fitted on k = 75 … 150 (`q2_largeorder.py`). As a control, the same code reproduces the q = −1 table of addendum-v2 (K²π to 5·10⁻³¹, s₁ to 5·10⁻²⁶).

### 3.1 q = i (verified)

- **Series.** ψ₁\*(ie^{−t}) = (2−i) + (−12+7i)t + (508 − 617i/2)t² + …; k!·aₖ is a Gaussian integer (checked for k ≤ 30): 2−i, −12+7i, 1016−617i, −185424+113359i, … For ψ₀\*: −1+2i, 14−22i, −1036+1666i, …
- **Singularity at −π²/240.** aₖ ≈ (−1)ᵏ Σⱼ cⱼ Γ(k+½−j) A^{−(k+½−j)} with A = π²/240 = π²/(15·4²). Fit residual 9·10⁻⁴⁰. Together with q = −1 this suggests A = π²/(15m²) at every even order m.
- **Constants** (found by PSLQ; differences from the fitted values in brackets):

| | π·c₀² | First coefficient s₁ of the attached series |
|---|---|---|
| ψ₁\* | (5−√5)(1−2i)/80 [6·10⁻³¹] | (175+84√5)/300 − i(5+√5)/25 [1·10⁻²⁵] |
| ψ₀\* | −(5+√5)(1+2i)/80 [5·10⁻³¹] | (84√5−125)/300 − i(5−√5)/25 [4·10⁻²⁶] |

  with c₀(ψ₁\*) = 0.13339439113182168 − 0.08244226762806335i and c₀(ψ₀\*) = −0.13339439113182168 + 0.21583665875988503i. The attached series is 1 + s₁t + …
- **It is the θ₄R term again.** Direct evaluation at t = 1/1500 … 1/48000: ½θ₄R₁(ie^{−t})·√t·e^{−π²/(240t)} → π·c₀(ψ₁\*) (ratio 1 within 10⁻²²), and its first correction equals s₁ to 14 digits. The same holds for −½θ₄R₀ and ψ₀\*.

### 3.2 q → 1 (verified)

- **Series.** By the corollary in §1 it is the series of ½f₁\* (resp. ½f₀\*): ψ₁\*: 1/3, 1/9, 5/162, 5/486, 29/17496, −1937/262440, …; ψ₀\*: 1/3, 0, 1/81, 1/81, 139/8748, … 3^{2k+1}·k!·aₖ is an integer (checked for k ≤ 39).
- **Singularities at S = π²/3 − i·Cl₂(π/3) and its conjugate**, |S| = 3.44287. These are −V at the two saddles of F(−1) (§2, a = ½): the prediction from the saddle analysis is confirmed. The coefficients oscillate with period ≈ 21 in k.

  aₖ ≈ 2·Re Σⱼ cⱼ Γ(k+½−j) S^{−(k+½−j)}. Fit residual 6·10⁻³¹; moving S by 10⁻¹² (relative) raises it to 3·10⁻²⁵, so S is confirmed to about 17 digits.
- **Constants** (closed forms; differences in brackets):

| | c₀ | s₁ |
|---|---|---|
| ψ₁\* | e^{iπ/12}/(2·3^{1/4}·√π) [2·10⁻¹⁷] | 1/4 − i√3/54 [4·10⁻¹⁴] |
| ψ₀\* | e^{−7iπ/12}/(2·3^{1/4}·√π) [3·10⁻¹⁷] | −1/12 − i√3/54 [6·10⁻¹⁴] |

  s₁ is known to 13 digits only, below the 15 asked for.
- **Stokes constant.** c₀ = i·C/(2π), where C = (1+w)√(πw/(1+2w)) = √(π/√3)·e^{−5πi/12} is the saddle prefactor of §2 at w = e^{−2πi/3}. Direct evaluation of ½θ₄R₁ at q = e^{−t}, t = 1/60 … 1/140, gives C to 12 digits and the same s₁ to 9.
- **Remark.** 2·Cl₂(π/3) = 2.0298832… is the hyperbolic volume of the figure-eight knot complement. The pair (π²/3, Cl₂(π/3)) plays the role that (Chern–Simons, volume) plays for a Kashaev invariant; see §5.

---

## 4. Identities for F₀\*, F₁\* (question 4)

Let F₀\* = Σ(−1)ⁿq^{2n²}/(q;q²)ₙ, F₁\* = Σ(−1)ⁿq^{2n²+2n}/(q;q²)ₙ₊₁, and

  X(q) = Σ (−1)ⁿ q^{n²}/(q;q)₂ₙ, Y(q) = Σ (−1)ⁿ q^{(n+1)²}/(q;q)₂ₙ₊₁.

**Identities (proved; verified to q³¹⁹, `q4_twists.py`):**

  (III) F₀\*(q) + Σ_{n≥1} q^{n²}(q;q²)ₙ = θ₄(q²)·(−q;q)∞·X(q²)
  (IV) q·F₁\*(q) − Σ_{n≥0} q^{(n+1)²}(q;q²)ₙ = θ₄(q²)·(−q;q)∞·Y(q²)

The sums on the left are φ₀\*(−q) − 1 and −φ₁\*(−q) in your notation, so (III) is the evaluation of the bilateral sum B₀\* of `twists.md`.

**Proof of (III).** The left side is the bilateral sum Σ_{n∈ℤ}(−1)ⁿq^{2n²}/(q;q²)ₙ, since 1/(q;q²)₋ₙ = (−1)ⁿq^{−n²}(q;q²)ₙ. For all n ∈ ℤ, 1/(q;q²)ₙ = (q^{2n+1};q²)∞/(q;q²)∞, and Euler's expansion gives (q^{2n+1};q²)∞ = Σⱼ(−1)ʲq^{j²+2nj}/(q²;q²)ⱼ. The exponent j²+2nj+2n² is positive definite, so the sums may be exchanged. The inner sum Σₙ(−1)ⁿq^{2n²+2nj} vanishes for odd j (substitute n ↦ −n−j) and equals (−1)ⁱq^{−2i²}θ₄(q²) for j = 2i. This leaves θ₄(q²)/(q;q²)∞ · Σᵢ(−1)ⁱq^{2i²}/(q²;q²)₂ᵢ. ∎ (IV) is the same with 1/(q;q²)ₙ₊₁; only odd j survive.

**Remarks.**
- The prefactor is the eta quotient (q²;q²)∞³/((q;q)∞(q⁴;q⁴)∞).
- X, Y are sign-twisted versions of Slater-type sums such as Σq^{n²}/(q)₂ₙ. They are **not** (−q;q)∞·R₀ or q(−q;q)∞·R₁ (checked). The earlier search missed them because (q;q)₂ₙ was not among its denominator types.
- So ψ\* is a pattern: bounded twist = half of a bilateral sum, and the bilateral sum = theta-type product × a twisted Rogers–Ramanujan/Slater sum. I did not analyse X, Y at roots of unity.

**Not done:** T₀\*, F(i), F(−ω). For T₀\* the summand contains a ratio of two Pochhammer symbols; the method needs a partial-fraction step and does not collapse as above.

---

## 5. A modularity law (question 1)

### 5.1 Exact values at every root

With §1, ψ₁\* and ψ₀\* are exact finite sums at every root of unity: the terminating sum at even order, one third of the first m terms at odd order (`jval.py`). This removes the restriction to even denominators in the earlier tests.

### 5.2 The dual sum

Define, at a root of unity ξ,

  **W(ξ) = Σ_{j≥0} (−1)ʲ ξ^{j²+j} (ξ;ξ²)ⱼ.**

At odd order the sum terminates. At even order its terms repeat with a constant factor and it is summed as a geometric series, exactly as in the odd-root rule (`wval.py`). Examples: W(1) = 1, W(−1) = 1/3, W(∓i) = (2±i)/3, W(e^{−iπ/3}) = e^{iπ/3}, W(e^{−2πi/3}) = (5−√3·i)/2.

### 5.3 The law (verified)

Let X → ∞ through rationals with X ≡ ρ mod 2 for a fixed rational ρ, and t = −2πi/X. Then

  ψ₁\*(e^{2πi/X}) = C̄·t^{−1/2}·e^{−S̄/t}·W(e^{−πiX})·(1 + s̄₁t + O(t²)) + O(1),
  ψ₀\*(e^{2πi/X}) = e^{2πi/3}·(the same, with s̄₁ replaced by its ψ₀\* value),

with S̄ = π²/3 + i·Cl₂(π/3), C̄ = √(π/√3)·e^{5πi/12}, and s̄₁ the complex conjugates of §3.2. Here |e^{−S̄/t}| = e^{Cl₂(π/3)·X/2π} = e^{0.16153·X} and its phase is e^{−iπX/6}.

**Evidence** (`q1_law.py`): for all 56 classes ρ with denominator ≤ 9 (both classes mod 2), X up to 2400, Neville extrapolation in 1/X:
- leading coefficient: largest deviation from C̄·W(e^{−πiρ}) is **6·10⁻²¹**;
- first correction: largest deviation from s̄₁ is 4·10⁻¹⁹, the same for every ρ.

**How W was found** (`q1_search.py`): the 56 limits were compared with 3,744 finite sums Σ εʲ ξ^{(aj²+bj)/2} Pⱼ(ξ) at eight choices of dual root, seven Pochhammer types and their inverses. Exactly one function matched (in two equivalent descriptions), for both ψ₁\* and ψ₀\*. Linear combinations of 1, ψ₀\*, ψ₁\* at related points do not fit (relative residual ≥ 0.79).

### 5.4 Second sector (exploratory, 3 digits)

Subtracting the growing term with ten coefficients of its series (taken from the large-order fit of §3.2) from ψ₁\*(e^{2πi/X}) at X ≈ 60–80 leaves the perturbative series Σaₖtᵏ of §3.2 with coefficient 1, to about 10⁻³, for ρ = 0, 1, 1/2, 1/3, 2/3 (`q1_second.py`). So the dual vector appears to start (1, W, ?).

### 5.5 Reading (inference)

- This is quantum modularity in the sense of Zagier's modularity conjecture for the Kashaev invariant and of Garoufalidis–Zagier: ψ\*(e^{2πi/X}) = Σ_σ Φ^{(σ)}(t)·J^{(σ)}(X), with divergent series Φ^{(σ)} (the three sectors 0, S, S̄ of §3.2) and coefficient functions J^{(σ)} on ℚ/2ℤ. The new information is J^{(S̄)} = W.
- Under x ↦ 1/x the pair (ψ₁\*, ψ₀\*) does **not** transform into itself. It transforms into W, a sum of the same type with (ξ;ξ²)ⱼ in place of (−q;q)ₙ, which terminates at odd instead of even order. This matches the pairing of ψ₀, ψ₁ with φ₀(−q), φ₁(−q) in the classical fifth-order transformation formulas (stated from memory; to be checked against Watson or Gordon–McIntosh).
- It explains the earlier fits. The map x ↦ x/(2x+1) is x ↦ 1/x, a shift by 2, and x ↦ 1/x again; W has period 2, so ψ\* returns to ψ\*, but with a cocycle built from the divergent series Φ. A model with constant multipliers and one weight can only be an effective fit, as addendum-v2 §4 argued.
- ψ₁\* and ψ₀\* have the same dual W in the growing sector (up to e^{2πi/3}); they differ in the series and, presumably, in the third sector.

### 5.6 What is missing for the full matrix

- The coefficient function of the decaying sector e^{−S/t}. It is invisible in this direction; the complex-conjugate direction x ↦ −1/x only swaps S and S̄.
- The converse: asymptotics of W(e^{2πi/Y}). A first look at Y ≈ 300 shows no exponential growth (|W| between 2 and 22, irregular), so the converse is not a mirror image and needs its own analysis.
- A proof. The law should follow from the standard stationary-phase analysis of finite q-hypergeometric sums near q = 1 (as for the figure-eight knot); I did not attempt it.

---

## 6. Literature (question 6)

Confirmed references:
- R₀ = A_q(1), R₁ = A_q(q) for Ramanujan's entire function A_q(z) = Σ(−1)ⁿq^{n²}zⁿ/(q;q)ₙ. G. E. Andrews, "Ramanujan's 'Lost' Notebook VIII: The entire Rogers–Ramanujan function", Adv. Math. 191 (2005) 393–407.
- S. Garoufalidis, D. Zagier: "Asymptotics of Nahm sums at roots of unity", Ramanujan J. 55 (2021), arXiv:1812.07690; "Knots, perturbative series and quantum modularity", SIGMA 20 (2024) 055, arXiv:2111.06645; "Knots and their related q-series", SIGMA 19 (2023) 082, arXiv:2304.09377.
- J. Bajpai, S. Kimport, J. Liang, D. Ma, J. Ricci, "Bilateral series and Ramanujan's radial limits", Proc. AMS 143 (2015) 479–492, arXiv:2202.12141. J. Mc Laughlin, arXiv:1906.11997 (contains Σ_{n∈ℤ}q^{n²}/(−q)ₙ = f₀ + 2ψ₀ as theta products).
- For §1 and §5: J. Lovejoy, "Quantum q-series identities" (Hardy–Ramanujan J.), on identities that hold only at roots of unity; A. Folsom et al., "Antiquantum q-series identities and mock theta functions", arXiv:2507.03824, where truncations at roots of unity of otherwise divergent series appear; K. Hikami, J. Lovejoy, "Torus knots and quantum modular forms", Res. Math. Sci. 2 (2015).

Not found in a short search (four queries): ψ₀\*, ψ₁\*, identities (I)–(IV), the function W, or the law of §5. This is "not found", not "new". The people named in addendum-v2 §5.4 remain the right ones to ask; §5 makes the question concrete enough to put to Garoufalidis, Zagier or Wheeler.

---

## 7. Suggested next steps

1. **Prove the law of §5** by stationary phase, and determine the third coefficient function.
2. **Find the partner of W.** W belongs to the family Σ zʲq^{j²}(−q;q²)ⱼ at a half-way point, as ψ\* does in F(z). Its companions and its own identities (analogues of (I), (II)) are the natural next search.
3. **Repeat §5 for F₀\*, F₁\*** using (III), (IV); test whether their duals are again finite sums of the same kind.
4. **Classifier:** add the class "bounded with a slowly decaying exponential" (§2), and use exact finite sums at odd roots where the rule of §1 or its analogue A_z/(1−2z^m) applies.
5. **T₀\*, F(i), F(−ω):** identities still open.

---

## 8. Scripts (`psi_star/v3/`)

Run from the scripts folder. Requirements: `mpmath`, `python-flint` (only for `q2_series.py`).

- `q3_oddroot.py`: §1 checks.
- `q5_threshold.py [path to zfam_results.jsonl]`: §2, threshold, table, precision check, limit rule.
- `q2_series.py 150` → `q2_series.json` (about 2 min); `q2_largeorder.py`: §3 fits, θ₄R comparison, closed forms; `q2_integrality.py`: arithmetic of the coefficients.
- `q4_twists.py`: (III), (IV) to q³¹⁹ (`qs.py` holds the exact series helpers).
- `jval.py`, `wval.py`: exact values of ψ\* and W at roots of unity.
- `q1_limits.py 9` → `q1_limits.json`; `q1_fit.py`: linear fits against 1, ψ₀\*, ψ₁\* (negative); `q1_search.py`: the search that found W; `q1_law.py 9`: the test of §5.3; `q1_second.py`: §5.4.
- `run_log.txt`: output of all of the above as run for this document.
