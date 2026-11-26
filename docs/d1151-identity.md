# The identity for d1151: D = θ₄·f + 4R, with R a sum of theta tails

This identifies the piece R = (D − θ₄·f)/4 left open in [`d1151.md`](d1151.md), section 6, for the double sum D(q) = Σ_{m,n≥0} q^{m²−mn+n²}/((−q;q)ₘ(−q;q)ₙ). It was worked out in a separate session (Claude Fable 5.2) from a brief listing the facts and the failed searches ([`d1151-brief.md`](d1151-brief.md)), and checked independently afterwards (Claude Opus 5.5): the proof by hand, and all parts of `census2/d1151/verify_R.py` (about 5 s).

Every claim carries one of the labels **proved**, **verified exactly to q^N**, **numerical**, **heuristic**.

Notation: (−q)ₙ = ∏_{k=1}^{n}(1 + qᵏ), θ₄ = Σ_{n∈ℤ}(−1)ⁿq^{n²}, f = Σ q^{n²}/(−q)ₙ².

## 1. Result

**Theorem (proved; also verified exactly to q³¹⁹ against an independent computation of R).**

  R(q) = Σ_{d > k ≥ 0} (−1)^{d+k+1} q^{d²} / (−q)ₖ

Equivalent forms of the same double sum:

- **Theta tails.** R = Σ_{k≥0} (−1)ᵏ Tₖ / (−q)ₖ, where Tₖ = Σ_{d>k} (−1)^{d−1} q^{d²} = (θ₄⁽ᵏ⁾ − θ₄)/2 and θ₄⁽ᵏ⁾ = Σ_{|j|≤k} (−1)ʲ q^{j²} is the truncated theta series.
- **Partial sums.** R = Σ_{d≥1} (−1)^{d−1} q^{d²} P_d, where P_d = Σ_{k=0}^{d−1} (−1)ᵏ/(−q)ₖ.
- **Two-variable f₀.** R = Σ_{i≥1} (−1)^{i−1} q^{i²} Φ(q^{2i}), where Φ(z) = Σₙ zⁿ q^{n²}/(−q)ₙ is the two-variable f₀. (Verified exactly to q³¹⁹.)

Consequences:

- **Divisibility by 4 is proved.** R has integer coefficients because each 1/(−q)ₖ does.
- **Sum-of-tails form of D (proved; verified to q⁵⁹).**
  D = θ₄/(−q)∞ + 2 Σ_{k≥0} (−1)ᵏ [ θ₄⁽ᵏ⁾/(−q)ₖ − θ₄/(−q)∞ ]
  The same statement without the theta factors is f = 1/(−q)∞ + 2 Σ_{k≥0} (−1)ᵏ [1/(−q)ₖ − 1/(−q)∞]. So D is f with the truncated theta series inserted term by term.

## 2. Proof

**Step 1: diagonals.** Put m = n + d. Then m² − mn + n² = n² + dn + d², so by the symmetry m ↔ n

  D = S₀ + 2 Σ_{d≥1} q^{d²} S_d,  S_d := Σ_{n≥0} q^{n²+dn} / ((−q)ₙ (−q)_{n+d}),  S₀ = f.

**Step 2: lemma.** S_d + S_{d+1} = 2/(−q)_d for all d ≥ 0. (Proved; verified exactly for d = 0…15 to q³¹⁹.)

Proof: let uₙ = q^{n²+dn} / ((−q)_{n−1} (−q)_{n+d}), with the convention 1/(−q)₋₁ = 1 + q⁰ = 2. Then

  uₙ − uₙ₊₁ = q^{n²+dn} [(1 + qⁿ)(1 + q^{n+d+1}) − q^{2n+d+1}] / ((−q)ₙ (−q)_{n+d+1})
       = q^{n²+dn} (1 + qⁿ + q^{n+d+1}) / ((−q)ₙ (−q)_{n+d+1}).

The n-th term of S_d is q^{n²+dn}(1 + q^{n+d+1}) over the same denominator, and the n-th term of S_{d+1} is q^{n²+dn}·qⁿ over it. Their sum is uₙ − uₙ₊₁. The series telescopes to u₀ = 2/(−q)_d. ∎

**Step 3: solve.** With S₀ = f the lemma gives (−1)ᵈ S_d = f − 2 P_d.

**Step 4: sum.** D = Σ_{d∈ℤ} q^{d²} S_{|d|} = Σ_{d∈ℤ} (−1)ᵈ q^{d²} (f − 2P_{|d|}) = θ₄·f + 4 Σ_{d≥1} (−1)^{d−1} q^{d²} P_d. ∎

All rearrangements are of formal power series in q: the coefficient of any fixed power involves finitely many terms.

So θ₄·f appears because every diagonal of D is ±f plus a finite sum of 1/(−q)ₖ, and the signs ± are exactly those of θ₄.

## 3. The same identity for general z

Let D_z = Σ_{m,n≥0} q^{m²−mn+n²} / ((zq)ₘ (q/z)ₙ) and let ρ(z) = Σ q^{n²}/((zq)ₙ(q/z)ₙ) be the rank generating function. The same telescoping gives, for S_d(z) = Σₙ q^{n²+dn}/((q/z)ₙ (zq)_{n+d}):

  S_{d+1}(z) = z·S_d(z) − (z − 1)/(zq)_d  (proved)

and therefore

  D_z = ρ(z)·Σ_{d∈ℤ} zᵈ q^{d²} − (1 − 1/z) Σ_{d>k≥0} z^{d−k} q^{d²}/(zq)ₖ − (1 − z) Σ_{d>k≥0} z^{k−d} q^{d²}/(q/z)ₖ

(proved; verified exactly for z = 2, −3, 5/7 to q²⁷).

- z = −1 is the theorem above: ρ(−1) = f and Σ(−1)ᵈq^{d²} = θ₄.
- z = 1 is the Durfee-rectangle identity Σ q^{n²+dn}/((q)ₙ(q)_{n+d}) = 1/(q)∞: both tail sums vanish.
- Other roots of unity give siblings of D: Jacobi theta at z, times the rank function at z, plus theta tails. At z = i the rank function is the third-order φ(q) and the theta factor is θ₄(q⁴).
- The telescoping needs the diagonal exponent to be exactly n² + dn. Among forms am² + bmn + cn² that forces a = c = 1, b = −1. This is why the A₂ form is the special one.

## 4. What R is, and what it is not

- **R is not θ₄ times a Rogers–Ramanujan-type series.** It is a sum of tails of θ₄ weighted by (−1)ᵏ/(−q)ₖ. This agrees with all the negative searches in `d1151.md`, section 6.
- **Fifth-order flavour.** Φ satisfies q^{s−1} Φ(qˢ) = Φ(q^{s−2}) + Φ(q^{s−1}) − 2 (proved by the same (1 + qⁿ) shift; verified for s = 2…11). Since Φ(1) = f₀ and Φ(q) = f₁, each Φ(q^{2i}) is a combination of 1, f₀, f₁ with Laurent-polynomial coefficients. The homogeneous equation is solved by (−1)ˢ Σ q^{n²+sn}/(q)ₙ, i.e. the Rogers–Ramanujan functions for s = 0, 1. The coefficients have unboundedly negative powers of q, so R is an infinite, not a finite, combination. That is why no linear relation with f₀, f₁, G, H was found.
- **D is not a mock theta function in the strict sense.** D = θ₄·f + 4R is a mixed mock modular piece plus a non-modular piece. Its fingerprint equals f's because θ₄·f dominates: with κ = 1/24 against 1/60, R is smaller by a factor e^{−π²/(10 m² t)} near an even-order root of order m (taking κ = 1/60 for R at all even roots, as the classifier finds; derived below only for m = 2).

## 5. Asymptotics toward q = −1

Put q = −e^{−t}. Then

  R(−e^{−t}) = √(2/(5−√5)) · θ₃(e^{−t}) · e^{π²/(60t)} · (1 + σ₁ t + O(t²)),  σ₁ = −1/15 − (5+2√5)/4 = −(79+30√5)/60.

- **Exponent π²/60 and the constant:** heuristic derivation by a saddle point, confirmed numerically to 12 digits.
- **σ₁:** numerical (11 digits, then PSLQ with coefficients ≤ 79).

Sketch of the saddle point. At q = −q′ the denominator (−q)ₖ becomes (q′;q′²)_{⌈k/2⌉} (−q′²;q′²)_{⌊k/2⌋} and Tₖ becomes −Σ_{d>k} q′^{d²}. With v = q′ᵏ the terms are e^{π²/(24t)}·exp(−ψ(v)/t) times a smooth factor, where ψ = log²v + Li₂(v²)/4. The saddle is at v⁴ = 1 − v², i.e. v² = (√5−1)/2, where ψ = π²/40 by Li₂(w) = π²/10 − log²w at w = (√5−1)/2. Hence π²/24 − π²/40 = π²/60. Even and odd k form two branches whose ratio at the saddle is v²/(1−v) ≠ 1, so they do not cancel.

Two remarks for the asymptotic route of `d1151.md`:

- The factor t^{−1/2} (weight ½) comes from the Gaussian sum over k. It does not require a θ₄ factor: the denominators carry no power of t.
- The leading term and the −1/15 in σ₁ are exactly those of θ₄(q)·G(q⁴), since G(e^{−4t}) = √(2/(5−√5))·e^{π²/(60t)}·e^{−t/15}·(1 + exponentially small). R departs from θ₄(q)·G(q⁴) at order t, by −(5+2√5)/4. These are the power-series corrections seen by the classifier.

## 6. Literature

A short search (two queries) found neither D nor this decomposition. Treat this as "not found", not as "new". (The search for D itself is in `d1151.md`, section 2.)

- Nearest in spirit: q-hypergeometric double sums that are mock theta functions (Lovejoy–Osburn, arXiv:1204.2431 and arXiv:1502.00387), and Hecke-type double sums (Hickerson–Mortenson, arXiv:1208.1421). Their double sums come from Bailey pairs; D is not obviously of that form.
- The form of R belongs to the sums-of-tails literature: Zagier's identity for η, and Andrews–Jiménez-Urroz–Ono, "q-series identities and values of certain L-functions", Duke Math. J. 108 (2001). Those sum tails of products; R sums tails of a theta series with weights 1/(−q)ₖ.
- The ingredients are classical: the diagonal sums S_d(z) are Durfee-rectangle versions of the rank generating function, and the lemma is a contiguous relation of the kind found in Fine's book. A specialist may well recognize section 3.

## 7. Open, in order of likely value

1. Ask a specialist whether the general-z identity of section 3 is known.
2. Asymptotics of R at the other even-order roots and at q = 1 directly from the theorem (the Borel singularity π²/60 at q = 1 should come from the complex saddles of the same ψ).
3. σ₂, σ₃ at q = −1 for a cross-check with the parallel asymptotic route.
4. The siblings D_z at z = i and at cube roots of unity: are they in the census, and do they show the same fingerprint pattern?
