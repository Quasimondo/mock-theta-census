# Sign twists of the classical mock theta functions

ψ₁\* and ψ₀\* are the fifth-order ψ₁, ψ₀ with a factor (−1)ⁿ inserted into each term. This note asks whether that is a one-off or a pattern. It applies the same twist to every mock theta function in the validation catalog, and runs the growth classifier at all 19 roots of unity on each twisted function.

Status labels as elsewhere in this repo: **verified** = exact or high-precision computation, no proof; **observation** = what the classifier reports.

## Method

- **Functions:** each of the 27 catalog mock theta functions (orders 2, 3, 5, 6, 7, 8, 10), written as an Eulerian series Σ q^{e(n)} Pₙ(q), gets the twist Σ (−1)ⁿ q^{e(n)} Pₙ(q). Three more series are added from the family F(z) = Σ zⁿ q^{n(n+1)/2} (−q;q)ₙ, at z = i, ω and −ω (ω = e^{2πi/3}). Recall that F(1) = ψ₁ and F(−1) = ψ₁\*.
- **Test:** the same radial-growth test as in the census (`census/probe.py`), at the same 19 roots, 570 probes in all. Per root it gives E (exact theta-type law, rational exponent), Q (bounded, power-series remainder) or X (irrational exponent or no exact law).
- **Trivial cases:** each twisted function's coefficients were compared with the catalog functions up to sign, monomial factor and q → −q. Some twists are trivial: when Pₙ depends only on q² and e(n) ≡ n mod 2, the twist is just the original function at −q.
- **Code and data:** `census/twist.py`, `census/twist_results.jsonl`.

## Results (observation)

Cells in the order of the 19 roots 0/1, 1/2, 1/3, 1/4, 1/5, 2/5, 1/6, 1/7, 2/7, 3/7, 1/8, 3/8, 1/10, 3/10, 1/12, 5/12, 1/14, 3/14, 5/14.

### Bounded at all 19 roots, not a catalog function at ±q

| Twist of | Order | Cells | Coefficients |
|---|---|---|---|
| ψ₁ (= ψ₁\*) | 5 | `QQQQQQQQQQQQQQQQQQQ` | dense |
| ψ₀ (= ψ₀\*) | 5 | `QQQQQQQQQQQQQQQQQQQ` | dense |
| F₀ | 5 | `QQQQQQQQQQQQQQQQQQQ` | dense |
| F₁ | 5 | `QQQQQQQQQQQQQQQQQQQ` | dense |
| ψ | 3 | `QQQQQQQQQQQQQQQQQQQ` | only 0, ±1, but about half nonzero |
| T₀ | 8 | `QQQQQQQQQQQQQQQQQQQ` | dense |
| F(z) at z = i | — | `QQQQQQQQQQQQQQQQQQQ` | dense, Gaussian integers |
| F(z) at z = −ω | — | `QQQQQQQQQQQQQQQQQQQ` | dense, in ℤ[ω] |

### Bounded at all 19 roots, false theta functions

These have coefficients ±1 at quadratic exponents. That is a known and simpler type than the dense cases above (exponents verified to q⁶⁹).

| Twist of | Order | Nonzero coefficients at |
|---|---|---|
| ν | 3 | k(3k−2) and k(3k+2): 0, 1, 5, 8, 16, 21, 33, 40, 56, 65 |
| φ (10th order) | 10 | 5k² − 2k and 5k² + 2k: 0, 3, 7, 16, 24, 39, 51 |
| ψ (10th order) | 10 | 5k² − 6k + 2 and 5k² − 4k + 1: 1, 2, 10, 13, 29, 34, 58, 65 |

### Trivial

| Twist of | Order | Is |
|---|---|---|
| φ | 3 | φ(−q); the twist only changes q to −q |

Its cells are `QQQEQQQQQQEEQQEEQQQ`; it is mock, as φ is.

### Neither (some root with non-theta growth)

f, χ, ω, ρ (3rd order); f₀, f₁, φ₀, φ₁ (5th); F₀, F₁, F₂ (7th); μ (2nd); φ, ψ (6th); S₀ (8th); X, χ (10th); and F(z) at z = ω.

The 5th-order twists f₀\*, f₁\* are in this group. That fits identities (I), (II) in `addendum-v2.md`: there f₀\*, f₁\* are the non-modular pieces that ψ₀\*, ψ₁\* are built from.

## F₀\* and F₁\*: first checks (verified)

F₀\*, F₁\* are the natural next pair after ψ₀\*, ψ₁\*: same order, also a pair, and also bounded at all 19 roots. The route that worked for ψ\* was to extend the sum to negative n and find a modular piece in the resulting bilateral series. For F₀\* that gives:

- **Bilateral sum.** Extending F₀ to negative n gives B₀ = F₀(q) + φ₀(−q) − 1. Its twisted version is B₀\* = F₀\*(q) + φ₀\*(−q) − 1. This follows from 1/(q;q²)₋ₙ = (−1)ⁿ q^{−n²} (q;q²)ₙ.
- **Classical case.** B₀ is E at all six roots tested (orders 1 to 6) with κ = 1/120, the fifth-order value. That is consistent with the classical bilateral series being modular. It is not a single infinite product, which suggests a sum of theta quotients.
- **Twisted case.** B₀\* is X at the same six roots, like θ₄·R₁ in the ψ case.
- **No simple factorization found.** B₀\* is not (up to a factor ±1, ±2 and a monomial) an eta quotient ∏ (q^a;q^a)∞^{e_a} (a ∈ {1, 2, 4}, eₐ ∈ {−2, …, 2}) times a single sum Σ (±1)ⁿ q^{(An²+Bn)/2}/(factor)ₙ (A ≤ 8, six product types). This was checked to q²⁶⁰.
- **No linear relation found.** Neither F₀\* nor F₁\* is a rational linear combination of 1, ψ₀\*, ψ₁\*, f₀\*, f₁\* (at q and q²) and φ₀\*, φ₁\* (at ±q). This was checked to q²⁰⁰.
- **Scripts:** `census/twist_checks/` (`F1.py` builds the series; `F2.py`–`F5.py` are the searches).

So F₀\*, F₁\* are not a rearrangement of the ψ\* story. A guess, not verified: F₀ is built on (q;q²)ₙ, so the bilateral argument needs a base-q² version of the partial-fraction step, and the missing factor would be a series outside the library above. The method that identified R₀, R₁ for ψ\* (the exact expansion at q = −1 and its large-order behaviour) is the obvious next step.

## Caveats

- "Bounded at 19 roots" is a classification by the radial-growth test, not a proof of boundedness. Nor does it say what kind of object a function is (quantum modular form, Habiro-type element, or something else).
- The comparison with catalog functions rules out the simplest equivalences only: sign, monomial factor and q → −q. A twisted function could still equal a known function in a less obvious way, for example after a change of variable, or through an identity with theta-function coefficients like (I) and (II).
- F(i) and F(−ω) have non-real coefficients. The classifier looks only at |F|.
