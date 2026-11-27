# Mock theta search: addendum v20 (ψ₁\* and ψ₀\* at even roots of unity: exact formula and asymptotics)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann). Item 2 of `psi-star-addendum-v19.md` §8: radial limits of the exact formulas, and the growth at q = e^{2πi/X} as a theorem. Scripts and output are in `psi_star/v20/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v20/check_thm10.py`, own code): ψ₁\*(i) = 2 − i exactly, and Fable's three parts at X = 4 add up to it. Theorem 10 against the exact finite sums at X = 200 … 3200 (up to 340 digits): X·(ratio − 1) converges, and Richardson extrapolation in 1/X gives 0.2015332627 − 1.570796327i for ψ₁\* and 0.2015332627 + 0.5235987757i for ψ₀\*. These are −2πi·s̄₁ with s̄₁ = 1/4 + i√3/54 (π√3/27 = 0.2015332627, π/2), and −2πi·(−1/12 + i√3/54) for ψ₀\*. So the second-order term asked for in §5.3 is the first correction s̄₁ of v5 at t = −2πi/X, and for ψ₀\* it is −1/12 + i√3/54, as already seen in the v19 check. Theorem 9 itself (the integrals at real b) was not recomputed; the lemmas were read.

Labels: **proved** (a complete argument is written here or in the addenda cited; nobody else has checked it), **verified** (numerics), **quoted**.

---

## 0. Summary

1. **Theorem 9 (proved; §1–§2).** Let X ≥ 4 be even, q₀ = e^{2πi/X}, x = 1/X, b = X^{−1/2}, t = −2πi/X. Then the finite sums

   ψ₁\*(q₀) = Σ_{n<X/2} (−1)ⁿ·q₀^{n(n+1)/2}·(−q₀;q₀)ₙ,  ψ₀\*(q₀) = Σ_{n<X/2} (−1)ⁿ·q₀^{(n+1)(n+2)/2}·(−q₀;q₀)ₙ

   satisfy

   ψ₁\*(q₀) = P·[ J(1/(2b)) + I(1/(2b)) − I(−1/(2b)) ],  ψ₀\*(q₀) = e^{πi/X}·P·[ J(b + 1/(2b)) + I(b + 1/(2b)) − I(b − 1/(2b)) ],

   with P = 2^{−1/2}·e^{π²/(12t) + t/6}/b and the three integrals of v12 taken with Faddeev's Φ_b at real b. So at an even root of unity the three sectors are three integrals, and the dual series have disappeared (W = p·W₃ = 1 there).
2. **Theorem 10 (proved; §3).** As X → ∞ through even integers,

   ψ₁\*(e^{2πi/X}) = 3^{−1/4}·(X/2)^{1/2}·e^{Cl₂(π/3)·X/2π}·e^{2πi/3 − iπX/6}·(1 + O(1/X)),  ψ₀\*(e^{2πi/X}) = e^{2πi/3}·(the same)·(1 + O(1/X)).

   This is the volume rule for z = −1 with its constant and phase.
3. **Corollary (proved, given the known asymptotics of the Kashaev invariant; quoted).** |ψ₁\*(e^{2πi/X})| / ⟨4₁⟩_{X/2} = (2/X)·(1 + O(1/X)): the comparison with the figure-eight knot found in the other session.
4. **How.** Three elementary lemmas: the radial limit of the q-series (the Pochhammer symbols vanish to increasing order at q₀, and a block of X consecutive factors is nearly 1 − y^X); the limit W(p), p·W₃(p) → 1 as p → 1 in a Stolz angle; continuity of the integrals at real x. Then Theorem 6 of v19 at real x.
5. **Numerically (verified):** the formula of Theorem 9 holds to 10⁻⁸ … 10⁻¹³ (double precision) for X = 4 … 60; in Theorem 10, X·(ratio − 1) tends to 0.2015 − 1.571i for ψ₁\* and 0.2017 + 0.524i for ψ₀\*.
6. **Not covered:** odd X, other roots of unity e^{2πih/k}, and z ≠ −1.

---

## 1. Three lemmas

Throughout X = 2M is even, q₀ = e^{2πi/X}, so q₀^M = −1.

**Lemma R1 (radial limit of the q-series).** As ρ → 1⁻: ψ₁\*(ρq₀) → Σ_{n<M} (−1)ⁿ q₀^{n(n+1)/2} (−q₀;q₀)ₙ, and likewise for ψ₀\*.

Proof. Let q = ρq₀, δ = 1 − ρ, X·δ ≤ 1. It suffices to show Σ_{n≥M} |q|^{n(n+1)/2}·|(−q;q)ₙ| → 0; the terms with n < M are continuous.

- *Blocks.* For m ≥ 0 let B_m = Π_{j=1}^{X} (1 + q^{mX+j}). Put c = ρ^{(m+1)X} and c_j = ρ^{mX+j}, so 0 ≤ c_j − c ≤ X·δ. Since Π_{j=1}^{X} (1 + c·q₀ʲ) = 1 − c^X (X even), expanding the product gives |B_m − (1 − c^X)| ≤ (2 + Xδ)^X − 2^X ≤ X²·3^{X−1}·δ. And 0 ≤ 1 − c^X ≤ min(1, (m + 1)·X²·δ). So with C = X²·(1 + 3^{X−1}): |B_m| ≤ min(1 + Cδ, (m + 1)·Cδ).
- *M ≤ n < X.* (−q;q)ₙ contains 1 + q^M = 1 − ρ^M ≤ Mδ; the other factors are at most 2. These finitely many terms tend to 0.
- *n ≥ X.* Write n = BX + r, 0 ≤ r < X, B ≥ 1. Then |(−q;q)ₙ| ≤ 2^X·Π_{m<B} |B_m|. Let B₁ = ⌊1/(2Cδ)⌋. For B ≤ B₁ every factor (m + 1)Cδ is at most ½, so the product is at most min(Cδ, 2^{−B}); the sum over B ≤ B₁ tends to 0. For B > B₁ the product is at most 2^{−B₁}·e^{Cδ·B}, and |q|^{n²/2} ≤ e^{−δB²X²/2}; the sum over B of e^{CδB − δB²X²/2} is O(δ^{−1/2}), so this part is O(δ^{−1/2}·2^{−1/(2Cδ)}) → 0. Each B stands for X values of n. ∎

**Lemma R2 (the dual series at p = 1).** As p → 1 with |p| < 1 and |1 − p| ≤ (3/2)·(1 − |p|): W(p) → 1 and W₃(p) → 1.

Proof. W(p) − 1 = Σ_{j≥1} (−1)ʲ p^{j²+j} (p;p²)ⱼ. Let δ = 1 − |p|, K = 3/2. |1 − pᵏ| ≤ k·Kδ and |1 − pᵏ| ≤ 1 + |p|ᵏ. Let j₁ = ⌊1/(2Kδ)⌋.

- j ≤ j₁: |(p;p²)ⱼ| ≤ (2j − 1)!!·(Kδ)ʲ ≤ min(Kδ, √2·e^{−j}) (each factor (2i + 1)Kδ ≤ 1, and (2j − 1)!! ≤ √2·(2j/e)ʲ). The sum over these j tends to 0.
- j > j₁: |(p;p²)ⱼ| ≤ √2·e^{−j₁}·Π_{i≥j₁} (1 + |p|^{2i+1}) ≤ √2·e^{−j₁}·exp(|p|^{2j₁+1}/(1 − |p|²)) ≤ C·exp(−(1/δ)·[1/(2K) − e^{−1/K + 2δ}/(2 − δ)]). For K = 3/2 the bracket tends to 1/3 − 0.257 = 0.077 > 0. With Σ_j |p|^{j²} = O(δ^{−1/2}) this part tends to 0.

W₃ has the extra factors p^{2j} of modulus < 1. ∎

**Lemma R3 (continuity at real x).** Let 0 < x₀ ≤ ¼. The integrals I(u), J(u) of v12 over L = e^{−10°·i}·ℝ, with Faddeev's Φ_b, for u ∈ {±1/(2b), b ± 1/(2b)}, are continuous in x on {Im x ≥ 0} near x₀.

Proof. Φ_b(y) is analytic in (b, y) off its poles and zeros, and the zeros stay below L. By Lemmas F and D of v16 and the inversion relation, on L one has |1/Φ_b| ≤ C on the left half and |1/Φ_b| ≤ C·|e^{−w²/(2t)}| on the right half, with C uniform near x₀ (|t| is bounded away from 0 there, so the terms Li₂/t are bounded). The Gaussian factor e^{−w²/(2t)} of the integrand has modulus exp(−ρ²·sin 20°/(2|t|)) on L, and the remaining factors are e^{O(ρ)}. Dominated convergence. (For X = 4 the point x₀ = ¼ is on the edge of the range |x| ≤ ¼ of v16; the condition actually used there is Im(i/x) ≥ 2, which holds in a neighbourhood of ¼.) ∎

**The path.** For x = 1/X + iε, ε → 0⁺: q = q₀·e^{−2πε} is radial, and p = e^{−πi/x} has |p| = exp(−πεX²/(1 + ε²X²)) and arg p = πε²X³/(1 + ε²X²) (mod 2π; X even). So p → 1 and |1 − p|/(1 − |p|) → 1: Lemma R2 applies.

---

## 2. Theorem 9

**Theorem 9.** For even X ≥ 4 the two identities of the summary hold.

Proof. For Im x > 0, v12 gives (with (κ/2)·𝓑₂ = P·I(1/(2b)), (κ/2)·𝓑₄ = P·I(−1/(2b)), −i·q^{1/4}·(κ/2)·𝓑₁ = q^{1/2}·P·I(b + 1/(2b)), i·q^{1/4}·(κ/2)·𝓑₃ = q^{1/2}·P·I(b − 1/(2b)), from Theorem C there)

  ψ₁\*(q) = P·[ J(1/(2b)) + I(1/(2b))·p·W₃(p) − I(−1/(2b))·W(p) ],

  ψ₀\*(q) = q^{1/2}·P·[ J(b + 1/(2b)) + I(b + 1/(2b))·p·W₃(p) − I(b − 1/(2b))·W(p) ].

The integrals with the product formula for Φ_b and any admissible δ are those with Faddeev's integral and δ = 10° (v16 Lemma Φ; v14 §2.1), and P = 2^{−1/2}·e^{π²/(12t) + t/6}/b. Let x = 1/X + iε, ε → 0⁺, and use Lemmas R1, R2, R3. ∎

Remarks.

- f₁\*, f₀\* have no such limit: their terms have (−q;q)ₙ in the denominator, and on the other side F_a(p), F_b(p) have poles at p = 1.
- The formula holds for every even X ≥ 4, not only asymptotically: at X = 4, ψ₁\*(i) = 2 − i is the sum of P·J = 0.276 − 0.146i, P·I(1/(2b)) = −0.539 − 0.092i and −P·I(−1/(2b)) = 2.263 − 0.763i.

---

## 3. Theorem 10

**Theorem 10.** As X → ∞ through even integers, the two asymptotic formulas of the summary hold.

Proof. In Theorem 9, by Theorem 6 of v19 at real x (arg x = 0 is included there):

- −P·I(−1/(2b)) = A·e^{5πi/12}·t^{−1/2}·e^{−S̄/t}·(1 + O(t)), A = (π/√3)^{1/2}, S̄ = π²/3 + i·Cl₂ (v19 §6.1);
- P·I(1/(2b)) = O(|t|^{−1/2}·e^{−Cl₂/|t|});
- |P·J(1/(2b))| ≤ C (v16, Theorem 3 with K = 0; in fact P·J = 1/3 + O(1/X)).

At t = −2πi/X: t^{−1/2} = (X/2π)^{1/2}·e^{iπ/4} and e^{−S̄/t} = e^{−i·S̄·X/2π} = e^{Cl₂·X/2π}·e^{−iπX/6}. So the leading term is (π/√3)^{1/2}·(X/2π)^{1/2}·e^{5πi/12 + iπ/4}·e^{Cl₂X/2π − iπX/6} = 3^{−1/4}·(X/2)^{1/2}·e^{2πi/3}·e^{Cl₂X/2π − iπX/6}, and the other two terms are smaller by e^{−Cl₂X/2π}. For ψ₀\*: −e^{πi/X}·P·I(b − 1/(2b)) = e^{2πi/3}·(−P·I(−1/(2b)))·(1 + O(t)) by Theorem 6. ∎

**Corollary (figure-eight).** The Kashaev invariant of 4₁ is ⟨4₁⟩_N = Σ_{k<N} Π_{j≤k} |1 − e^{2πij/N}|², and ⟨4₁⟩_N = 3^{−1/4}·N^{3/2}·e^{N·Vol(4₁)/2π}·(1 + O(1/N)) with Vol(4₁) = 2·Cl₂(π/3) (**quoted** from memory: Andersen–Hansen 2006; I did not look it up). With Theorem 10 at N = X/2:

  |ψ₁\*(e^{2πi/X})| / ⟨4₁⟩_{X/2} = (2/X)·(1 + O(1/X)).

Both have the growth rate Cl₂(π/3)/2π per unit of X and the torsion factor 3^{−1/4}; they differ by one power of N. (Verified directly in v18 §2.5: the ratio times X/2 is 0.99887, 0.99943, 0.99972 at X = 802, 1602, 3202.)

---

## 4. Numerical checks (verified; `h1_roots.py`, `run_log.txt`)

| check | result |
|---|---|
| Theorem 9, ψ₁\*: finite sum against P·J + P·I(+) − P·I(−) at x = 1/X, X = 4, 6, 8, 12, 16, 24, 40, 60 | relative difference 2·10⁻⁹, 1·10⁻⁹, 5·10⁻¹⁰, 7·10⁻¹¹, 9·10⁻¹², 1·10⁻¹², 4·10⁻¹², 5·10⁻¹² |
| Theorem 9, ψ₀\*, the same X | 5·10⁻¹², 3·10⁻¹⁰, 3·10⁻⁸, 2·10⁻¹³, 3·10⁻¹³, 8·10⁻¹³, 3·10⁻¹², 5·10⁻¹² |
| Lemma R1 at X = 8: \|ψ₁\*(q₀e^{−2πε}) − finite sum\| for ε = 10⁻², 10⁻³, 10⁻⁴ | 1.47, 0.68, 0.11 (slow, as the proof suggests) |
| Lemma R2 on the same path: \|W(p) − 1\|, \|p·W₃(p) − 1\| | 0.016, 0.10, 0.018 and 0.87, 0.24, 0.037; \|1 − p\|/(1 − \|p\|) = 1.0023, 1.00003, 1.0000 |
| Theorem 10: X·(exact/asymptotic − 1) at X = 100 … 3200 | ψ₁\*: 0.197 − 1.588i … 0.2014 − 1.5713i; ψ₀\*: 0.208 + 0.538i … 0.2017 + 0.5240i (converging) |

Notes.

- Everything at real x is double precision, with Φ_b from Faddeev's integral (`fad.py` of v16). The integral for J is taken over L. The integrals for I are taken over the descent contours of v19, because over L they cancel by many orders of magnitude (the first attempt over L gave errors of size 10⁴ at X = 8 for that reason).
- For ψ₀\* with X ≤ 8 the decaying integral is taken over L instead: the extra factor e^{w} beats the decay on the horizontal line when |t| > π/3.
- The limits in the last row look like 0.2015 − iπ/2 and 0.2015 + iπ/6; I did not derive the second-order term.

---

## 5. Open

1. Odd X and general roots of unity e^{2πih/k}: these correspond to x → h/k, where the dual nome tends to another root of unity; the formula of v12 still holds for Im x > 0, but the limits of W, W₃ and the size of the integrals are different. Not looked at.
2. The same for the z-family (needs the growing integral for a ≠ ½, v19 §6.3, and a bound for J(a/b) at real x).
3. The second-order term of Theorem 10.
4. The remaining items of v19 §8.

---

## 6. Scripts (`psi_star/v20/`)

- `h1_roots.py [a b c]`: a: Theorem 9 (40 seconds); b: the radial limits; c: Theorem 10.
- `fad.py` (from v16), `run_log.txt`, `README.md`.
