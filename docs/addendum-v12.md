# Mock theta search: addendum v12 (the contour argument: the kernel identity is proved)

Session of 26 Nov 2026 (Claude Fable 5.2 with Mario Klingemann), started from `HANDOFF.md` (state after v11). This does the task named there: derive, by a correct contour argument, the identity

  kernel integral = ψ₁\*(q) − (κ/2)·[𝓑₂·p·W₃(p) − 𝓑₄·W(p)].

Scripts and output are in `psi_star/v12/`.

**Independent check** (Claude Opus, same day, `psi_star/v12/check_thmAC.py`: own Φ_b, own quadrature on L, own double series): Theorem C, I(u) = c·e^{−πiu²/2}·D(u), holds at x = 0.3 + 0.2i for a generic u = 0.3 − 0.2i and for u = 1/(2b), b + 1/(2b), to 1–3·10⁻³¹; Theorem A (A↑) holds at the generic u and at u = 1/(2b) to 6·10⁻³² … 2·10⁻³¹ (at x = 0.45 + 0.3i both also hold, there to 10⁻¹⁷ with a lower-precision constant). The functional equation Φ_b(y + i/b) = Φ_b(y)/(1 + p⁻¹w) behind (B1) and the logic of steps (A), (B), (C) were read and are consistent; the analytic estimates of §2 and §4.1 were not checked line by line. So: the identities are confirmed independently, the proof has been read but not refereed.

Labels: **proved** (a complete argument is written in this addendum; nobody else has checked it), **verified** (numerics), **inference**. Every intermediate identity of the proof was also checked numerically (§6).

**Later remarks:** `addendum-v26.md` §4 (the powers p^{n²/2} in Theorem C mean e^{−πin²/(2x)}); `addendum-v27.md` §4 (the remark on the factor in the closing lemma of §4.1, reworded).

---

## 0. Summary

1. **The identity of v11 is proved** (§1–§4). The mechanism has three steps, each elementary:
   - **(A)** Turn the left end of the contour across the kernel poles y = i·b·(n + ½), n ≥ 0. The residues give ψ₁\*. On the new contour |e^{2πy/b}| > 1, so the kernel is a convergent geometric series in e^{−2πy/b}, and each term is a kernel-free integral I(u − m/b). This is the "second family of contributions" that v11 asked for.
   - **(B)** The integrals I(u − m/b) satisfy a three-term recursion in m with coefficients that are powers of p. Summing it gives exactly I(u)·p·W₃(p) − I(u − 1/b)·W(p). This is where W and W₃ come from.
   - **(C)** The kernel-free integral I(u) is, after two Fourier transforms of Φ_b, the integral of v7. So I(1/(2b)) and I(−1/(2b)) are 𝓑₂ and 𝓑₄ times an elementary factor, and that factor is κ/(2P).
2. **The same argument gives the other three cases** (§5, proved): ψ₀\* (with 𝓑₁, 𝓑₃), and, closing on the other side, f₁\* and f₀\* (with F_a, F_b in place of W, W₃). The difference of the two closings is identity (I)/(II) of v2 combined with (V), (VI) of v4 and the inversion formulas of v6, so the whole set of exact formulas is consistent.
3. **Everything holds for general u** (Theorems A, B, C): a two-variable version, in which W, W₃ become one function 𝒲(t; p).
4. **v7 is now proved** as well: the closing of the contour there is justified by an explicit estimate (§4.1).
5. **Not settled here:** that P·J is the Borel sum of the perturbative series (v10, verified to 55 digits, no proof), and the asymptotic laws of v5. They now need only the small-b asymptotics of two kinds of explicit integrals (§7).
6. **Relation to the paper you sent** (Garoufalidis–Gu–Mariño–Wheeler, arXiv:2111.04763): its state integral for the trivial connection of 4₁ has the same kernel, and its Theorem 7 has the same shape. Its proof does not carry over, for a concrete reason (§7).

---

## 1. Notation and statements

x in the upper half-plane, b = √x with β = arg b ∈ (0, π/2); q = e^{2πix}, p = e^{−πi/x}, q̃ = p²; c_b = i(b + 1/b)/2. For a point y write ξ = y/b, v = e^{2πby} = e^{2πxξ}, w = e^{2πy/b} = e^{2πξ}. As in v7,

  Φ_b(y) = (−q^{1/2}v; q)∞ / (−p·w; p²)∞.

Fix δ with 0 < δ < π/2 − β and let L = e^{−iδ}·ℝ. For u ∈ ℂ put

  g_u(y) = e^{−πi·y² + 2π·u·y} / Φ_b(y),  **I(u) = ∫_L g_u(y) dy,  J(u) = ∫_L g_u(y) / (1 + e^{2πy/b}) dy.**

The integral of v11 is J(1/(2b)). Further

- Ψ(z) = Σ_{n≥0} q^{n(n+1)/2}·zⁿ·(−q;q)ₙ, Ψ̌(z) = Σ_{n≥0} q^{n²+n}·zⁿ/(−q;q)ₙ, so ψ₁\* = Ψ(−1), ψ₀\* = q·Ψ(−q), f₁\* = Ψ̌(−1), f₀\* = Ψ̌(−1/q);
- 𝒲(t) = Σ_{k≥0} (−t)ᵏ·p^{k²+k}·(p;p²)ₖ, so W(p) = 𝒲(1), W₃(p) = 𝒲(p²);
- 𝒢₀(t) = Σ_{k≥0} t^{−k}·p^{2k²+k}/(p;p²)ₖ, 𝒢₁(t) = Σ_{k≥0} t^{−k}·p^{2k²+3k}/(p;p²)ₖ₊₁, so 𝒢₀(1) = F_b(p) − F_a(p), 𝒢₁(1) = F_a(p)/p;
- D(u) = Σ_{k,n≥0} (−1)^{kn+n}·e^{πiu(bk + n/b)}·q^{k²/4}·p^{n²/2+n} / ((q;q)ₖ·(p²;p²)ₙ);
- c = 2^{−1/2}·e^{−iπ/4}·q^{−1/24}·q̃^{1/24};
- P = (−q;q)∞/(b·q^{1/8}·(p;p²)∞) (prefactor of v11), κ as in v10, 𝓑₁…𝓑₄ as in v6.

**Theorem A (kernel integral = one-sided sum + dual sum of kernel-free integrals).** For every u ∈ ℂ, with both sums over m absolutely convergent,

  (A↑) J(u) = −i·b·q^{1/8}·e^{πiub}·((p;p²)∞/(−q;q)∞)·Ψ(e^{2πiub}) − Σ_{m≥1} (−1)ᵐ·I(u − m/b),

  (A↓) J(u) = +i·b·q^{1/8}·e^{−πiub}·((p;p²)∞/(2(−q;q)∞))·Ψ̌(e^{−2πiub}) + Σ_{m≥0} (−1)ᵐ·I(u + m/b).

**Theorem B (the dual sums).** For every u ∈ ℂ:

  (B1) I(u) = p^{−1}·e^{2πiu/b}·[ I(u + 1/b) + p^{−1}·I(u + 2/b) ], and likewise I(u) = q^{1/2}·e^{2πiub}·[ I(u + b) + q^{1/2}·I(u + 2b) ].

  (B2) With t = p·e^{2πiu/b}:

  - Σ_{m≥1} (−1)ᵐ·I(u − m/b) = I(u)·[𝒲(t/p²) − 1] − (t/p)·I(u + 1/b)·𝒲(t),
  - Σ_{m≥0} (−1)ᵐ·I(u + m/b) = I(u)·𝒢₀(t) − I(u + 1/b)·𝒢₁(t).

  (B3) If t = 1, that is u ∈ 1/(2b) + b·ℤ:

  - Σ_{m≥1} (−1)ᵐ·I(u − m/b) = I(u)·p·W₃(p) − I(u − 1/b)·W(p),
  - Σ_{m≥0} (−1)ᵐ·I(u + m/b) = I(u)·F_b(p) − I(u − 1/b)·F_a(p).

**Theorem C (the kernel-free integral).** For every u ∈ ℂ,

  I(u) = c·e^{−πi·u²/2}·D(u).

In particular I(1/(2b)) = c·p^{1/8}·𝓑₂, I(−1/(2b)) = c·p^{1/8}·𝓑₄, I(b + 1/(2b)) = −i·c·q^{−1/4}·p^{1/8}·𝓑₁, I(b − 1/(2b)) = +i·c·q^{−1/4}·p^{1/8}·𝓑₃, and c·p^{1/8} = κ/(2P).

**Corollary (the four exact formulas).**

- ψ₁\*(q) = P·J(1/(2b)) + (κ/2)·[ 𝓑₂·p·W₃(p) − 𝓑₄·W(p) ]   (the identity of v11)
- f₁\*(q) = 2P·J(1/(2b)) − κ·[ 𝓑₂·F_b(p) − 𝓑₄·F_a(p) ]
- ψ₀\*(q) = P·q^{1/2}·J(b + 1/(2b)) − i·q^{1/4}·(κ/2)·[ 𝓑₁·p·W₃(p) + 𝓑₃·W(p) ]
- f₀\*(q) = 2P·q^{1/2}·J(b + 1/(2b)) + i·q^{1/4}·κ·[ 𝓑₁·F_b(p) + 𝓑₃·F_a(p) ]

All of §1 is **proved** (§2–§5) and **verified** (§6).

---

## 2. Proof of Theorem A

**Bounds for 1/Φ_b.** Two forms are used.

- (2.1) 1/Φ_b(y) = (−p·w; p²)∞/(−q^{1/2}v; q)∞. Where |w| ≤ R and |q^{1/2}v| ≤ r < 1, it is bounded by (−|p|R; |p|²)∞/(r; |q|)∞.
- (2.2) The inversion relation Φ_b(y)·Φ_b(−y) = Φ_b(0)²·e^{πi·y²}, Φ_b(0)² = q^{1/24}·q̃^{−1/24} (classical: Jacobi triple product and the theta transformation; checked in §6), gives 1/Φ_b(y) = e^{−πi·y²}·Φ_b(−y)/Φ_b(0)², and Φ_b(−y) = (−q^{1/2}/v; q)∞/(−p/w; p²)∞ is bounded where |1/v| ≤ R and |p/w| ≤ r < 1.

So g_u = O(|e^{−πiy² + 2πuy}|) in regions of the first kind, which include the left sector {Re(by) ≤ 0, Re(y/b) ≤ 0} (arg y from π/2 + β to 3π/2 − β), and g_u = O(|e^{−2πiy² + 2πuy}|) in regions of the second kind, which include the right sector (arg y from −π/2 + β to π/2 − β). Since |e^{−πiy²}| = exp(π·|y|²·sin(2·arg y)), both ends of L lie in directions of Gaussian decay, and I(u), J(u) converge for all u.

**Poles.** The kernel has simple poles at yₙ = i·b·(n + ½), n ∈ ℤ, that is ξ = i(n + ½), with residue −b/(2π). 1/Φ_b has poles at −c_b − i·b·k − i·n/b (k, n ≥ 0), all below L, with Re ξ ≤ −Im x/(2|x|²). The kernel poles with n ≥ 0 are above L, those with n < 0 below.

**Upper route.** Choose ε > 0 with e^{2π|Re x|ε}·|q|^{1/2} < 1 and let L′ be the path that comes down the line ξ = ε + iτ (τ from +∞ to 0) and continues along b·ε + e^{−iδ}·s (s ≥ 0). On the strip 0 ≤ Re ξ ≤ ε, Im ξ ≥ 0 we have |w| ≤ e^{2πε} and |q^{1/2}v| ≤ e^{2π|Re x|ε}|q|^{1/2} < 1, so (2.1) applies there and in the left sector. Close the region between the left ray of L and L′ by the segment Im ξ = N (N a positive integer). On that segment w = e^{2π·Re ξ} > 0, so |1/(1 + w)| ≤ 1, and arg y stays in [π/2 + β − o(1), π − δ], so |g_u| ≤ C·e^{−c·N²}. Between the two right rays there are no poles and Gaussian decay. Hence

  J(u) = 2πi·Σ_{n≥0} Res_{yₙ} + ∫_{L′} g_u/(1 + w) dy.

At yₙ: w = −1, v = q^{n+½}, so 1/Φ_b = (p;p²)∞·(−q;q)ₙ/(−q;q)∞, and e^{−πiyₙ² + 2πuyₙ} = q^{(n+½)²/2}·e^{2πiub(n+½)}. The residue sum is the first term of (A↑).

On L′ we have Re ξ ≥ ε, so 1/(1 + w) = Σ_{m≥1} (−1)^{m−1}·w^{−m} with Σ|w|^{−m} ≤ 1/(e^{2πε} − 1). Since ∫_{L′}|g_u||dy| < ∞, dominated convergence gives ∫_{L′} g_u/(1 + w) = Σ_{m≥1} (−1)^{m−1}·∫_{L′} g_u·w^{−m}, absolutely convergent. Finally g_u·w^{−m} = g_{u−m/b} has no poles between L′ and L, and the same bounds as above let the contour go back: ∫_{L′} g_{u−m/b} = I(u − m/b). This is (A↑). ∎

**Lower route.** Mirror image. Choose ε < Im x/(2|x|²), so that the strip −ε ≤ Re ξ ≤ 0, Im ξ ≤ 0 contains no pole of 1/Φ_b and |p/w| ≤ |p|e^{2πε} < 1, |1/v| ≤ e^{2π|Re x|ε} there; (2.2) applies in the strip and in the right sector. L″ runs along −b·ε + e^{i(π−δ)}·s (from s = ∞) and then down the line ξ = −ε − iτ. On the closing segments Im ξ = −N again w > 0. The contour is moved clockwise, so J(u) = −2πi·Σ_{n<0} Res + ∫_{L″}. At y = −i·b·(j + ½): v = q^{−j−½}, (−q^{−j}; q)∞ = 2·q^{−j(j+1)/2}·(−q;q)ⱼ·(−q;q)∞, which gives the terms q^{j²+j}/(−q;q)ⱼ of Ψ̌ and the factor ½. On L″, |w| ≤ e^{−2πε} and 1/(1 + w) = Σ_{m≥0} (−1)ᵐwᵐ. This is (A↓). ∎

Remark. Subtracting the two routes shows that neither one-sided sum alone is the integral: (A↑) − (A↓) says residue sum(n ≥ 0) + residue sum(n < 0) = Σ_{m∈ℤ} (−1)ᵐ·I(u + m/b), a bilateral dual sum.

---

## 3. Proof of Theorem B

**(B1).** From Φ_b(y + i/b) = Φ_b(y)/(1 + p^{−1}w) one gets g_u(y + i/b) = p^{−1}·e^{2πiu/b}·w·(1 + p^{−1}w)·g_u(y). There is no pole between L and L + i/b and the decay is Gaussian in the strip, so ∫_L g_u(y + i/b) dy = I(u); and w·g_u = g_{u+1/b}. The second relation follows in the same way from Φ_b(y + ib) = Φ_b(y)/(1 + q^{1/2}v). ∎

**(B2).** Put aₘ = I(u − m/b). (B1) at u − m/b reads aₘ = t·p^{2m−2}·aₘ₋₁ + t·p^{2m−3}·aₘ₋₂ for all m ∈ ℤ. By Theorem A, Σ_{m≥1}|aₘ| < ∞, so Ã(z) = a₀ + Σ_{m≥1} aₘzᵐ is continuous on |z| ≤ 1, and the recursion gives

  Ã(z) = t·z·(1 + pz)·Ã(p²z) + a₀ + t·p^{−1}·z·a₋₁.

Iterating (the remainder t^K·z^K·p^{K(K−1)}·(−pz;p²)_K·Ã(p^{2K}z) tends to 0, because Ã is bounded on the disk):

  Ã(z) = Σ_{k≥0} tᵏzᵏp^{k(k−1)}·(−pz; p²)ₖ·[ a₀ + t·p^{2k−1}·z·a₋₁ ].

At z = −1 this is a₀·𝒲(t/p²) − (t/p)·a₋₁·𝒲(t), which is the first formula. For the second put bⱼ = I(u + j/b) = a₋ⱼ, with Σ_{j≥0}|bⱼ| < ∞ again by Theorem A; the recursion solved for the lowest index is bⱼ = −p·bⱼ₋₁ + t^{−1}p^{2j−1}·bⱼ₋₂, so B(z) = Σ_{j≥0} bⱼzʲ satisfies B(z) = b₀ + b₁z/(1 + pz) + t^{−1}p³z²·B(p²z)/(1 + pz), and iteration at z = −1 gives b₀·𝒢₀(t) − b₁·𝒢₁(t). ∎

**(B3).** For t = 1: 𝒲(p^{−2}) = Σ (−1)ᵏp^{k²−k}(p;p²)ₖ = 1 − W + p·W₃ (shift k → k + 1 and use (p;p²)ₖ₊₁ = (p;p²)ₖ(1 − p^{2k+1})), and a₁ = a₀ + p^{−1}a₋₁. So a₀·(p·W₃ − W) − p^{−1}a₋₁·W = a₀·p·W₃ − a₁·W. In the same way 𝒢₀(1) = F_b − F_a, 𝒢₁(1) = F_a/p, and b₀(F_b − F_a) − p^{−1}b₁F_a = I(u)·F_b − I(u − 1/b)·F_a. ∎

So W, W₃, F_a, F_b are the four generating sums of one three-term recursion, read in its two directions.

---

## 4. Proof of Theorem C

The integral I(u) cannot be closed in either half-plane (see the remark after the closing lemma). The integral of v7 can, and I(u) is a Fourier transform away from it.

### 4.1 Closing lemma (this also completes v7)

Let C be the contour of v7 with both angles equal to δ: from ∞·e^{i(π−δ)} through 0 to ∞·e^{iδ}. Consider ∫_C Φ_b(y)·e^{−Aπi·y² + 2πλy} dy with A = ½ (any λ), or with A = 0 and 2πλ = 2πi·w where Im(w·e^{−iδ}) < 0 and Re(w·b) > 0.

**Claim.** The integral equals 2πi times the sum of the residues at all poles c_b + i·b·k + i·n/b (k, n ≥ 0), summed first over k, then over n.

Proof. In the variable ξ the poles are at Re ξ = (n + ½)·Im x/|x|², Im ξ = k + ½ + (n + ½)·Re x/|x|². Close C by the vertical line Re ξ = σ_N := N·Im x/|x|², which runs between the pole columns n = N − 1 and n = N, and by a horizontal line Im ξ = T₀ + K (K integer, T₀ fixed so that the line meets no pole with n < N).

- On Re ξ = σ_N, |w| = |p|^{−2N}, and |(−p·w; p²)∞| ≥ Π_j |1 − |p|^{2j+1−2N}| = |p|^{−N²}·(|p|;|p|²)_N·(|p|;|p|²)∞. So |1/(−p·w; p²)∞| ≤ C·|p|^{N²} = C·exp(−π·σ_N²·|x|²/Im x), uniformly on the whole line.
- For all y, |(−q^{1/2}v; q)∞| ≤ C·exp(π·(σ·Re x − τ·Im x)₊²/Im x), ξ = σ + iτ (triple product: (−a·r^{1/2}; r)∞ ≤ C_r·exp((log a)²/(2·log(1/r)))).
- |e^{−Aπiy²}| = exp(Aπ·Im(xξ²)) = exp(Aπ·[Im x·(σ² − τ²) + 2στ·Re x]).

Adding the exponents on the vertical line, with τ\* = σ_N·Re x/Im x: for τ ≥ τ\* the total is at most −π(1 − A)·σ_N²·|x|²/Im x − πA·Im x·(τ − τ\*)²; for τ ≤ τ\* (down to the point where the line meets C) it is −π(1 − A)·Im(y²). There arg y ∈ [δ, π/2 − β], and Im(y²) = Im(xξ²) is concave in τ, so its minimum is at an end of the segment: the total is at most −c·σ_N² with c = π(1 − A)·min(|x|²/Im x, |x|·sin 2δ/cos²(δ − β)). (The second value belongs to the point on C; `t5_bounds.py` confirms that the supremum on the line is this bound times a factor between 1 and 1.5.) The linear factor costs e^{O(σ_N + |τ|)} for A = ½, and for A = 0 it supplies the decay e^{−2πτ·Re(wb)} as τ → ∞. So the integral over the vertical line tends to 0 as N → ∞ whenever A < 1. On the horizontal line with N fixed, 1/(−p·w;p²)∞ is bounded by periodicity in ξ → ξ + i, |v| is bounded, and the integrand is O(e^{−cK²}) for A = ½, O(e^{−cK}) for A = 0. Let K → ∞, then N → ∞. ∎

Consequences:

- **v7 is proved**: I_½(v) := ∫_C Φ_b(s)·e^{−πi·s²/2 + 2πvs} ds = −i·b·((q)∞/(q̃)∞)·e^{−πi·c_b²/2 + 2πv·c_b}·Σ_{k,n} (−1)^{kn}·Zᵏq^{k²/4}/(q)ₖ·Z̃ⁿq̃^{n²/4}/(q̃)ₙ, Z = i·e^{2πivb}·q^{1/4}, Z̃ = −i·e^{2πiv/b}·q̃^{1/4}, for all v. (Residues as computed in v7.)
- **Fourier transform of Φ_b**: F(w) := ∫_C Φ_b(y)·e^{2πiwy} dy = −i·b·((q)∞/(q̃)∞)·e^{2πi·c_b·w}/Φ_b(−w − c_b). The residue sum is two q-binomial series. Both sides are analytic in w on the half-plane H = {Im(w·e^{−iδ}) < 0}, where the integral converges, so the formula holds on all of H.
- **Fourier transform of 1/Φ_b**: G(w) := ∫_L e^{2πiwy}/Φ_b(y) dy = (i/b)·((q̃)∞/(q)∞)·e^{−2πi·c_b·w}·Φ_b(w + c_b) on H. Same proof, closing downward with the lines Re(by) = N·Im x, on which |1/(−q^{1/2}v; q)∞| ≤ C·|q|^{N²/2}; the roles of the two products are exchanged.

These are the classical Fourier transform formulas for the quantum dilogarithm (Faddeev–Kashaev–Volkov). I did not compare normalisations with the literature; the proof above is included so that complex b is covered.

Remark. The estimate needs A < 1 for Φ_b·e^{−Aπiy²}. (As far as I remember, without having re-read it, this is the condition B > A of Garoufalidis–Kashaev for B = 1.) I(u) is, by (2.2), an integral of Φ_b·e^{−2πiy²}, so A = 2: it cannot be closed, which is why v11 found that "closing naively" fails and why step (C) needs a different tool.

### 4.2 Two Gaussian linearisations

Let Σ_η = {e^{iδ}(ρ − iη) : ρ ∈ ℝ} with η > 0 large, and K(λ) = ∫ Φ_b(z)·e^{πi·z² + 2πλz} dz over any line of direction e^{iδ} below all poles of Φ_b (left end in the left sector at angle π + δ, right end at angle δ; both are directions of Gaussian decay).

**(X)** For all y, e^{−πiy²} = e^{−iπ/4}·∫_{Σ_η} e^{πiσ² + 2πiσy} dσ. Insert this into I(u). For σ ∈ Σ_η and y = r·e^{−iδ} ∈ L, |e^{2πi(σ − iu)y}| = e^{2πr(η + Re(u·e^{−iδ}))}, independent of ρ; |e^{πiσ²}| ≤ C·e^{−π·sin 2δ·ρ² + 2πη|ρ|}; |1/Φ_b| is bounded on the left ray and O(e^{−π·sin 2δ·r²}) on the right ray. For η > −Re(u·e^{−iδ}) the double integral converges absolutely, so by Fubini

  I(u) = e^{−iπ/4}·∫_{Σ_η} e^{πiσ²}·G(σ − iu) dσ = e^{−iπ/4}·(i/b)·((q̃)∞/(q)∞)·e^{πi(iu − c_b)² + 2πi·c_b²}·K(−u − 2i·c_b)

(substitute z = σ − iu + c_b).

**(Y)** For all s, e^{−πis²/2} = √2·e^{−iπ/4}·∫_{Σ_η} e^{2πiσ² + 2πiσs} dσ. Insert this into I_½(v). On the left ray of C the factor |e^{2πi(σ − iv)s}| decays, independently of ρ; on the right ray s = r·e^{iδ} the three Gaussian exponents add up to −π·sin 2δ·[ρ² + (|ρ| − r)²] plus linear terms. Again Fubini applies, and with F from §4.1 and the inversion relation,

  I_½(v) = √2·e^{−iπ/4}·(−i·b)·((q)∞/(q̃)∞)·Φ_b(0)^{−2}·e^{2πi(iv − c_b)² − 2πi·c_b²}·K(−2v − i·c_b).

**Combining.** The two K's coincide for v = u/2 + i·c_b/2 (the two lines are parallel, below all poles, with decay between them). Dividing, and inserting the formula for I_½:

  I(u) = (i/(√2·b))·((q̃)∞/(q)∞)·Φ_b(0)²·e^{πi·c_b²}·e^{−πi·u²/2}·D(u),

because for v = u/2 + i·c_b/2 one has Z = e^{πiub}, Z̃ = −p·e^{πiu/b}, which turns the double sum of v7 into D(u). With (q)∞/(q̃)∞ = q^{−1/24}q̃^{1/24}/√(−ix), √(−ix) = e^{−iπ/4}b, Φ_b(0)² = q^{1/24}q̃^{−1/24} and e^{πi·c_b²} = −i·q^{−1/8}q̃^{1/8}, the constant is c. ∎

**Special values.** At u = ±1/(2b), b ± 1/(2b) the parity split of D(u) in k gives 𝓑₂, 𝓑₄, 𝓑₁, 𝓑₃ exactly as in v7 (the map v = u/2 + i·c_b/2 sends these u to u₂, u₄, u₁, u₃ of v7). And c·p^{1/8} = κ/(2P): from the definitions, κ/(2P) = (b/(4i·q^{1/8}))·θ₄(q)·(p;p²)∞²/((−q;q)∞²·ψ(p)), and by the eta transformations θ₄(q)/(−q;q)∞² · (p;p²)∞²/ψ(p) = (q)∞/(−q)∞³ · (p;p²)∞³/(p²;p²)∞ = 2^{3/2}·q^{1/12}·p^{5/24}/√(−ix). Hence κ/(2P) = 2^{−1/2}·e^{−iπ/4}·q^{−1/24}·p^{5/24} = c·p^{1/8}. So κ of v10 is exactly the constant the proof produces.

---

## 5. The corollary, and consistency

- **ψ₁\*.** (A↑) at u = 1/(2b): e^{πiub} = i, e^{2πiub} = −1, so the residue term is ψ₁\*/P. By (B3) and Theorem C the dual sum is c·p^{1/8}·[𝓑₂·p·W₃ − 𝓑₄·W] = (κ/(2P))·[…]. This is the identity of v11. ∎
- **f₁\*.** (A↓) at the same u: residue term f₁\*/(2P), dual sum (κ/(2P))·[𝓑₂F_b − 𝓑₄F_a].
- **ψ₀\*, f₀\*.** u = b + 1/(2b): e^{2πiub} = −q, e^{πiub} = i·q^{1/2}, t = 1 again; I(u), I(u − 1/b) are ∓i·c·q^{−1/4}p^{1/8} times 𝓑₁, 𝓑₃.
- **Consistency.** Subtracting the first two formulas: 2ψ₁\* − f₁\* = κ·[𝓑₂·(p·W₃ + F_b) − 𝓑₄·(W + F_a)] = κ·ψ(p)·[𝓑₂·P₋(p) − 𝓑₄·P₊(p)] by (V), (VI). By the inversion formula of v6 (i·q^{1/4}N₁ = ½(𝓑₂E₋ − 𝓑₄E₊)) the right side is θ₄R₁. So the contour argument reproduces identity (II) of v2; the same with (I) for ψ₀\*, f₀\*. Nothing new, but it ties the two closings, (V), (VI), v2 and v6 together.
- **tanh form.** Since 1/(1 + w) = ½(1 − tanh(πy/b)), J = ½I − ½Z with Z(u) = ∫_L tanh(πy/b)·g_u dy, and the first formula becomes ψ₁\* = −(P/2)·Z(1/(2b)) + (κ/2)·[𝓑₂·(p·W₃ + ½) − 𝓑₄·W].

---

## 6. Numerical checks (verified; 30-digit working precision, `run_log.txt`)

Four points x = 0.3 + 0.2i, 1/(3 − 0.4i), 0.11 + 0.03i, 0.45 + 0.3i (|q| from 0.15 to 0.83, |p| from 7·10⁻⁴ to 0.28).

| statement | script | worst deviation |
|---|---|---|
| I(±1/(2b)) = c·p^{1/8}·𝓑₂,₄ (relative) | `t1_kernel_free.py` | 6·10⁻³⁰ |
| c·p^{1/8} = κ/(2P) (relative) | `t1_kernel_free.py` | 8·10⁻³¹ |
| Theorem C at four further u (two generic, b ± 1/(2b)) | `t1_kernel_free.py` | 1.2·10⁻²⁸ |
| recursion (B1) by quadrature, m = 1, 2 | `t1_kernel_free.py` | 6·10⁻²⁴ at three points; 5·10⁻¹⁵ and 2·10⁻¹² at x = 0.11 + 0.03i |
| recursion (B1) with I from its series, m = 0…3 | `t3_steps.py` | 3·10⁻²⁵ |
| (A↑) and (A↓) at u = 1/(2b), b + 1/(2b) and two generic u (absolute, terms of size 0.1…10) | `t2_kernel.py` | 1.6·10⁻²⁹ |
| the four formulas of the corollary (absolute, values of size 1) | `t2_kernel.py` | 2·10⁻³⁰ |
| inversion relation, Φ_b(0)², (q)∞/(q̃)∞ | `t3_steps.py` | 2·10⁻³⁰ |
| Fourier transforms F(w), G(w) | `t3_steps.py` | 3·10⁻²⁸ |
| steps (X), (Y) | `t3_steps.py` | 2·10⁻²⁴ (one case; otherwise 3·10⁻²⁷ or better) |
| (B2), both directions, at four u | `t3_steps.py` | 2·10⁻²⁹ |
| 𝒲(p^{−2}) = 1 − W + p·W₃, 𝒢₀(1) = F_b − F_a, p·𝒢₁(1) = F_a | `t4_series.py` | exact to p³⁰⁰ |
| the two bounds of the closing lemma on the lines Re ξ = σ_N, N = 1…6 (20 digits) | `t5_bounds.py` | hold; ratio to the stated bound ≤ 1.46 |

Notes on the table.

- At x = 0.11 + 0.03i the quadrature test of the recursion loses digits (|p| = 7·10⁻⁴; I take this to be cancellation in the quadrature, since the same recursion with I from its series holds to 3·10⁻²⁵ at that point).
- In (X), (Y) the two ratios printed share the quadrature of K, and they deviate by the same amount, so the deviation is the accuracy of that quadrature.
- In `t2_kernel.py` and in part (d) of `t3_steps.py` the kernel-free integrals are evaluated by the series of Theorem C, not by quadrature; J(u) itself is always a quadrature.

---

## 7. What this settles, what it does not, and the paper

**Settled.** The three-sector formula of v10 with the kernel integral of v11 in the place of 𝒯 is a theorem, for ψ₁\*, ψ₀\*, f₁\*, f₀\*. The status line "nothing about the one-sided sums is proved" of the handoff no longer holds.

**Not settled.**

1. P·J(1/(2b)) = Borel sum of the perturbative series Φ (v10: verified to 55 digits). Needed: J is analytic in b in a sector around the positive axis (inference: the definition of J only uses Φ_b, which continues to real b, and the kernel poles stay off L), and its expansion at b → 0 is Φ with factorial bounds on a wide enough sector. Not written.
2. The asymptotic laws of v5 for the one-sided sums. With the theorem they follow from item 1 and the stationary-phase expansion of 𝓑ⱼ (v7 §3, a derivation). Not written.
3. The Stokes jump of v10 §3, identity (VII), W, 𝒥 in the bounded directions: untouched. For W, W₃ the natural guess is now the same construction with q and p exchanged (inference).

**The paper (arXiv:2111.04763v2).** I read §2.6–2.7 (pp. 19–22: eqs. (69)–(84), Theorem 7 and its proof) and the passage on twist knots around eq. (252); not the rest.

- Their state integral for the trivial connection of 4₁, eq. (73), is ∫ tanh(πy/b)·Φ_b(−y + i/(2b))/Φ_b(y − i/(2b)) dy. The kernel is the one found in v11: 1/(1 + e^{2πy/b}) = ½(1 − tanh(πy/b)).
- The same phenomenon appears there: the residues of the kernel in the upper half-plane give the original one-sided sum (for them the Kashaev invariant), but the paper notes that the integrand has "significant contributions from infinity in the upper half plane", so the integral is a different function. That is what v11 observed.
- Their Theorem 7 writes the integral as q-series times q̃-series, including one-sided sums of Appell–Lerch type in q̃, plus ½ of a bilinear term. The tanh form in §5 above has the same shape (inference): one-sided dual sums and ½·𝓑₂.
- Their eq. (75), Z = s(Φ^{(σ₀)}) − (i/2)τ^{−3/2}·s(Φ^{(σ₂)}), is stated there as an upgrade of a conjecture. The analogue here would be −(P/2)·Z(1/(2b)) = 𝒯 − (κ/4)·𝓑₂, which is v10 in tanh form (verified, not proved).
- Their pair, the Habiro series and the inverted Habiro series, corresponds to the pair ψ₁\* / f₁\* here: upper and lower kernel poles, (−q;q)ₙ against 1/(−q;q)ₙ (inference).
- **Why their proof does not carry over.** Theorem 7 is proved by closing the contour downward and summing residues at poles of order two and three (the kernel poles in the lower half-plane sit on dilogarithm poles). Here the lower kernel poles are separate from the poles of 1/Φ_b (offset i/(2b)), and the integrand cannot be closed in either direction (§4.1, remark). The geometric series of step (A) and the Fourier step (C) replace the residue computation, and the one-sided dual series come from a recursion, not from higher-order poles.

So this case has the structure of their trivial-connection row, in a variant with a single dilogarithm and a Gaussian factor. I did not find it in the pages I read, and I have not checked whether it follows from results elsewhere in that paper or in Wheeler's thesis.

---

## 8. Scripts (`psi-star-addendum-v12-scripts.zip`)

- `common.py`: Φ_b, the integrals I(u), J(u) by mpmath quadrature on L, the q-series, 𝓑₁…𝓑₄, D(u), c, P, κ.
- `t1_kernel_free.py [digits]`: Theorem C and the recursion.
- `t2_kernel.py [digits]`: Theorem A (both routes) and the corollary.
- `t3_steps.py [digits]`: the steps of §4 and the sums of Theorem B.
- `t4_series.py [N]`: the three series identities of (B3) as exact power series.
- `t5_bounds.py [digits]`: the bounds of the closing lemma.
- `run_all.sh [digits]`: runs all five and writes `run_log.txt`. mpmath only; about six minutes at 30 digits on two cores.
