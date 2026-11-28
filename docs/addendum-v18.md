# Mock theta search: addendum v18 (the z-family F(z): exact formula, the volume rule, and the limit z → 1)

Session of 27 Nov 2026 (Claude Fable 5.2 with Mario Klingemann), answering the brief "the z-family F(z) = Σ zⁿ q^{n(n+1)/2} (−q;q)ₙ". Scripts and output are in `psi_star/v18/`.

**Independent check** (Claude Opus, 27 Nov 2026, `psi_star/v18/check_rr.py`, `check_saddle.py`, `check_thmZ.py`): the four Rogers–Ramanujan-type identities of §3.1 hold exactly to q³²⁹; the mock-theta reductions at a = 0 (𝒲(1/p) = φ₀(−p) etc.) follow from the definitions. The saddle rule of §2.2 with the full one-loop prefactor was tested against exact values of F at e^{2πi/X}, X = 802, 1602, 3202, for z = i, e^{2πi/7}, −ω, e^{2πi/√2}: same W and k = −1 as in the table, |ratio − 1| falling like 1/X and extrapolating to 10⁻¹⁰. Theorem Z (strip form) at x = 0.3 + 0.2i and 0.45 + 0.3i, z = i and −ω, with own I (Theorem C series) and own J (quadrature): 7·10⁻³² … 5·10⁻³¹. The statement of §2.4 that `census3/zfamily/branches.py` picks the lower root for z = −1, i, −ω and the upper one for z = ω was rechecked and is correct; that README now says so.

Labels: **proved** (a complete argument exists in these addenda; nobody else has checked it), **verified** (numerics), **derivation** (a computation whose analytic steps are not justified), **inference**.

Updated the same day after `census3/zfamily/` appeared in the repository (commit bbac2c2): §2.4 compares with it, §2.5 with the figure-eight comparison of `psi_star/knot/`.

**Later remarks:** `addendum-v30.md` gives the matrix form of Theorem Z (Theorem 20: state integrals = dual frame⁻¹ · Jacobi factor · q-frame) and the behaviour of the q-frame under Γ₀(2) (Theorem 21), which is the precise form of §4.

---

## 0. Summary

1. **Your observation is right (proved; §1).** Theorems A–C of v12 at u = a/b give an exact formula for F(z), z = e^{2πia}, for every a ∈ ℂ. In the strip 0 < Re a < 1 its natural form is

   F(z) = i·P·e^{−πia}·[ J(a/b) − I((a−1)/b)·𝒲(T) + p·T·I(a/b)·𝒲(p²T) ],  T = p·z̃ = e^{πi(2a−1)/x},  z̃ = e^{2πia/x}.

   It is the formula for ψ₁\* with W → 𝒲(T), W₃ → 𝒲(p²T), 𝓑₄ → I((a−1)/b), 𝓑₂ → I(a/b). Checked at z = i, −ω and five other values of z, at three points x, to 10⁻²⁷ … 10⁻³¹.
2. **The dual functions.** 𝒲(T; p) = Σ (−T)ᵏ p^{k²+k} (p;p²)ₖ is the same kind of series as F: with 𝔉(z, y; q) = Σ zⁿ q^{n(n+1)/2} (y;q)ₙ one has F(z; q) = 𝔉(z, −q; q) and 𝒲(T; p) = 𝔉(−T, p; p²). The Jacobi variable transforms as a → a/x. The coefficients I are bilinear in the function H(·; q) of v6 and H(·; p²).
3. **The volume rule is a saddle-point statement (derivation, verified to 10⁻²⁶ … 10⁻³⁶; §2).** The potential of the state integral I((a+k)/b) is V(w) = −w²/2 − 2πi(a+k)w + Li₂(−e^{w}); its critical points are W = e^{w} with W² + W = 1/z; and at every critical point, on every branch k, Im V = D(−W) (Bloch–Wigner; proved). The selected saddle is: **the root W in the upper half-plane, with w = Log W principal**; it belongs to the branch a + k ∈ (−1, 0), that is to the term I((a−1)/b). With this rule the full leading asymptotics at q = e^{2πi/X}, including the constant, is

   F(z) ≈ (1 + W)·W^{1/2}·(π/(t(1 + 2W)))^{1/2}·exp((V(w) + π²/12)/t),  t = −2πi/X.

   The other root gives the decaying sector and belongs to I(a/b).
4. **z → 1 (proved formulas, the interpretation is inference; §3).** Nothing breaks: the formula is entire in a. At a = 0 the dual functions become Ramanujan's fifth-order mock theta functions, 𝒲(1/p) = φ₀(−p), p·𝒲(p) = −φ₁(−p), and on the other route F₀(p), F₁(p); the state integrals I(0), I(1/b) become explicit Rogers–Ramanujan products in q and in p^{1/2}; the volume goes to 0 like 2π·ln((1+√5)/2)·|a|. So at z = 1 the formula is a transformation law of Watson's type for ψ₁ and f₁, with P·J(0) in the place of the Mordell integral.
5. **Quantum Jacobi form (inference; §4):** not in the scalar sense of Bringmann–Folsom; it has the shape of one row of a matrix-valued holomorphic quantum modular form with a Jacobi variable, in the sense of Garoufalidis–Zagier.
6. **Not done:** a proof that the saddle of item 3 is the one the contour passes (steepest descent is not carried out); Borel summability of the trivial sector for a ≠ ½; the identification of P·J(0) with a classical Mordell integral.

---

## 1. The exact formula for F(z)

Notation of v12: Im x > 0, b = √x, q = e^{2πix}, p = e^{−πi/x}; I(u), J(u), 𝒲, 𝒢₀, 𝒢₁, D(u), c, P as there. Put u = a/b, so that e^{2πiub} = z = e^{2πia}, and

  ζ = e^{πia},  ζ̃ = e^{πia/x},  z̃ = ζ̃² = e^{2πia/x},  T = p·z̃.

(T is the variable called t in v12; from v13 on, t = −2πix.) Ψ of v12 is F, and Ψ̌(v) = F̌(v) := Σ q^{n²+n} vⁿ/(−q;q)ₙ.

**Theorem Z (proved: v12 Theorems A, B, C at u = a/b).** For every a ∈ ℂ:

- (Z↑) F(z) = i·P·e^{−πia}·[ J(a/b) + I(a/b)·(𝒲(T/p²) − 1) − (T/p)·I((a+1)/b)·𝒲(T) ]
- (Z↓) F̌(1/z) = −2i·P·e^{πia}·[ J(a/b) − I(a/b)·𝒢₀(T) + I((a+1)/b)·𝒢₁(T) ]
- (strip form) F(z) = i·P·e^{−πia}·[ J(a/b) − I((a−1)/b)·𝒲(T) + p·T·I(a/b)·𝒲(p²T) ]
- I(a/b) = c·e^{−πia²/(2x)}·D(a/b), with D(a/b) = H_e(ζ; q)·H(−ζ̃p; p²) + H_o(ζ; q)·H(ζ̃p; p²),

where H(y; q) = Σ yᵐ q^{m²/4}/(q;q)ₘ (v6) and H_e, H_o are its even and odd parts in y. Replacing a by a ± 1 replaces (ζ, ζ̃) by (−ζ, ζ̃·p^{∓1}).

Proof. (Z↑), (Z↓) are (A↑), (A↓) with (B2), solved for the residue term. Strip form: Σ_{m≥1} (−1)ᵐ I(u − m/b) = −I(u′) − Σ_{m≥1} (−1)ᵐ I(u′ − m/b) with u′ = u − 1/b, and (B2) at u′, whose dual variable is T·p². The formula for D is its definition, split by the parity of k. ∎

Remarks.

- **Why the strip.** |T| = |p|^{1−2·Re a} (for real a). For 0 < a < 1 all terms of 𝒲(T) and 𝒲(p²T) beyond the first are small, of relative size |p|^{3−2a} and |p|^{5−2a}; and the kernel e^{2πas}/(1 + e^{2πs}) of J(a/b) is a window, decaying on both sides. Every a gives a valid formula; a and a + 1 give two formulas for the same F(z), related by the recursion (B1).
- **a = ½** (z = −1): T = 1, and the strip form is the formula of v11/v12 for ψ₁\*.
- **The dual function is of the same family.** 𝔉(z, y; q) := Σ zⁿ q^{n(n+1)/2} (y;q)ₙ. Then F(z; q) = 𝔉(z, −q; q) and 𝒲(T; p) = 𝔉(−T, p; p²). The second argument is −q = q·e^{πi} on one side and (p²)^{1/2} on the other: the two half-periods, exchanged by x → −1/x. The Jacobi variable goes from a to a/x (z̃ = e^{2πia/x}), as for Jacobi forms. 𝒢₀, 𝒢₁ are the corresponding series with the Pochhammer symbol in the denominator.
- **Bilateral identity.** Subtracting (Z↓) from (Z↑): F(z) + F̌(1/z)/(2z) = i·P·e^{−πia}·Σ_{m∈ℤ} (−1)ᵐ I((a+m)/b). At z = −1 this is identity (II) of v2.

---

## 2. Saddles and the volume

### 2.1 Potential and critical points (derivation; the identity for Im V is proved)

In the variable w = 2π·b·y (v13) the integrand of I((a+k)/b) is e^{V_k(w)/t}·(1 + O(t)) where the dual factor of 1/Φ_b is negligible, with

  V_k(w) = −w²/2 − 2πi·(a + k)·w + Li₂(−e^{w}),  t = −2πix.

- **Critical points.** V_k′ = −w − 2πi(a + k) − log(1 + e^{w}) = 0 gives, for W = e^{w}: W² + W = 1/z, for every k. This is the shape equation of the brief. Conversely, for a root W and any integer k, w_k = −Log(1 + W) − 2πi(a + k) is a logarithm of W and a critical point of V_k.
- **Heights (proved, for real a).** With principal Li₂ and Log: Im V_k(w_k) = Im Li₂(−W) − Re w·(Im w + 2π(a + k)) = Im Li₂(−W) + ln|W|·arg(1 + W) = D(−W), the Bloch–Wigner function, for every k. The two roots satisfy W₁ + W₂ = −1, so −W₂ = 1 − (−W₁) and D(−W₂) = −D(−W₁). On |q| = 1 (t = −i|t|) the size of e^{V/t} is e^{−Im V/|t|}. So **every critical point on every branch has height +Vol(z) or −Vol(z)**, Vol(z) = |D(−W)|; +Vol belongs to the root in the upper half-plane.
- **The same potential from the sum.** At q = e^{2πi/X} the n-th term of F is e^{(V₀(w) + π²/12)/t}·(W(1 + W)/2)^{1/2}·(1 + O(t)) at w = −t·n ∈ [0, iπ]; Poisson summation over n produces the branches k, and these are the terms I((a − m)/b) of Theorem A (k = −m).

### 2.2 The selected saddle (verified)

**Rule.** W = the root of W² + W = 1/z with Im W > 0; w = Log W (so 0 < Im w < π: the saddle lies in the strip over the summation segment); k is then fixed by the critical point equation. For 0 < a < 1 this gives k = −1 in all cases tried.

**Prediction** (one-loop saddle point): F(z; e^{2πi/X}) = (1 + W)·e^{w/2}·(π/(t(1 + 2W)))^{1/2}·e^{(V_k(w) + π²/12)/t}·(1 + O(1/X)), principal square root.

Check (`b1_roots.py`): F computed as the finite sum (n < X/2) for X = 802, 1602, …, 9602, divided by the prediction, extrapolated in 1/X.

| z | W | Vol(z) | Vol/2π | \|ratio − 1\| extrapolated |
|---|---|---|---|---|
| −1 | −0.5 + 0.8660i | 1.01494160640965363 | 0.161532973610 | 1.7·10⁻³¹ |
| i | −1.3002 + 0.6248i | 0.69670113910389213 | 0.110883430146 | 1.4·10⁻³⁴ |
| ω | −1.0707 + 0.7587i | 0.865129197895580713 | 0.137689588258 | 3.7·10⁻³³ |
| −ω | 0.4736 + 0.4448i | 0.486617604148760301 | 0.0774475970958 | 2.7·10⁻³⁶ |
| −i | 0.3002 + 0.6248i | 0.69670113910389213 | 0.110883430146 | 5.9·10⁻³⁵ |
| e^{2πi/7} | −1.5114 + 0.3865i | 0.421067509918694248 | 0.067014975579 | 6.9·10⁻³³ |
| e^{2πi·0.1} | −1.5654 + 0.2758i | 0.2986420034137574 | 0.0475303510582 | 5.1·10⁻²⁶ |
| e^{2πi·0.9} | 0.5654 + 0.2758i | 0.2986420034137574 | 0.0475303510582 | 3.1·10⁻²⁶ |
| e^{2πi·0.03} | −1.6133 + 0.0842i | 0.0906064183087769 | 0.0144204593497 | 2.3·10⁻¹⁴ |
| e^{2πi/√2} | 0.1884 + 0.7001i | 0.789592547465243954 | 0.125667556958 | 3.6·10⁻³⁴ |

The four rates of the brief (0.16153, 0.11088, 0.13769, 0.07745) are reproduced, and so is the constant and the phase. For small volume the other contributions are only e^{−Vol·X/2π} smaller, which limits the extrapolation (last but one row).

A numerical trap: the terms of the finite sum reach e^{0.1615·X} whatever z is, so the sum needs about 0.07·X extra digits. With fixed precision every z seems to grow at the rate 0.150, which is rounding noise.

### 2.3 Which term of the exact formula carries it (verified)

For 0 < a < 1 and x → 0 at 1° from the real axis (`b2_dominant.py`, |x| = 1/99, 1/199, 1/399; z = −1, i, ω, −ω):

- I((a−1)/b), computed exactly from Theorem C, divided by its saddle-point value at the selected saddle: 1 + O(t), the deviation halves with |x| (for example z = i: 3.7·10⁻³, 1.8·10⁻³, 9.2·10⁻⁴ in modulus).
- I(a/b) divided by the saddle-point value at the **other** root (lower half-plane, w = Log W, branch k = 0): also 1 + O(t). Its size relative to I((a−1)/b) is 10⁻⁸ … 10⁻⁵⁷.
- F(z) + i·P·e^{−πia}·I((a−1)/b) is small against F: relative size from 6·10⁻⁴ (z = −ω, |x| = 1/99) down to 10⁻²² (z = ω, |x| = 1/399). For z = i it is 1/(1 − 2z) up to 0.06, 0.009, 0.004: the trivial sector, whose leading term is the sum of the geometric series Σ (2z)ⁿ continued analytically. For the other z the remainder is larger than 1; its size is consistent with the next dual term, −I((a−1)/b)·(𝒲(T) − 1) (for z = −1 it is |p|² times the leading term, as it should be; not checked separately for the others).

So the three sectors of v10 exist for every z in the strip: J(a/b) is the trivial one, I((a−1)/b) the growing one (shape W in the upper half-plane, volume +Vol), I(a/b) the decaying one. **The tetrahedron-volume rule is the statement that the term m = 1 of the dual sum of Theorem A is dominated by its critical point W, w = Log W.** What is missing for a proof is the steepest-descent argument for I((a−1)/b) (or the asymptotics of H(·; q) at q → 1, since I is bilinear in H).

### 2.4 Comparison with `census3/zfamily/` (verified; `d1_compare.py`)

That folder uses V(a, x) = 2πi·a·x − x²/2 + Li₂(−e^{−x}) + π²/12 with x = −log w + 2πik and w² + w = 1/z. This is the potential above with x = −w(here), root w(there) = W, and the branch put into x instead of into a.

- **The rates agree** (0.16153, 0.11088, 0.13769, 0.07745), and so does the identification with the Bloch–Wigner function.
- **The selected branch there is underdetermined, and for three of the four z it is the other root.** `branches.py` picks the candidate whose |Im V| matches the measured rate. For each z there are two candidates that match exactly, one for each root: the true critical points of V(a, ·). In the notation here both are on the branch k = 0, and they are the saddles of I(a/b). The other candidates are not critical points and have other values of |Im V| (for z = −1 all candidates tie, because |W| = 1). The table in `volumes.py` has the lower root for z = −1, i, −ω and the upper root for z = ω. |Im V| cannot tell the two roots apart, since D(−W₂) = −D(−W₁).
- **What decides** is the full complex asymptotics (§2.2): phase, constant and sign of Im V. They select the upper root with w = Log W, which is a critical point on the branch k = −1, not k = 0. In the notation of the folder: the root w with Im w > 0, x = −Log w (principal, k = 0), and a replaced by a − 1 in V (the stationarity condition holds with the shift 2πi·j, j = −1, mentioned in a comment of `branches.py`).
- The open question in that README, whether F(i), F(ω), F(−ω) have the full structure of ψ₁\* (dual sum, state integrals, three sectors): yes, §1 and §2.3.

### 2.5 The figure-eight comparison (verified; uses the definition of the Kashaev invariant only)

With J(N) = Σ_{k<N} Π_{j≤k} |1 − e^{2πij/N}|², the Kashaev invariant of 4₁: |ψ₁\*(e^{2πi/X})|/J(X/2)·(X/2) = 0.99887, 0.99943, 0.99972 at X = 802, 1602, 3202, so the ratio tends to 2/X as the other session found. The formula of §2.2 at z = −1 gives |ψ₁\*| ≈ 3^{−1/4}·(X/2)^{1/2}·e^{Cl₂(π/3)·X/2π} (ratio to the exact value 1.00025, 1.00013, 1.00006). Together: J(N) ≈ 3^{−1/4}·N^{3/2}·e^{N·Vol(4₁)/2π} with Vol(4₁) = 2·Cl₂(π/3), which is the known asymptotics of the Kashaev invariant (quoted from memory). So ψ₁\* at X and the figure-eight knot at N = X/2 share volume per unit and torsion 3^{−1/4}, and differ by one power of N.

---

## 3. The limit z → 1

### 3.1 Formulas at a = 0 (proved; series identities exact to order 300; numerics 10⁻³⁰)

At a = 0: T = p, and (`c1_series.py`)

  𝒲(1/p) = φ₀(−p),  p·𝒲(p) = −φ₁(−p),  𝒢₀(p) = F₀(p),  𝒢₁(p) = F₁(p),

with Ramanujan's fifth-order functions φ₀(q) = Σ q^{n²}(−q;q²)ₙ, φ₁(q) = Σ q^{(n+1)²}(−q;q²)ₙ, F₀(q) = Σ q^{2n²}/(q;q²)ₙ, F₁(q) = Σ q^{2n²+2n}/(q;q²)ₙ₊₁. And F(1) = ψ₁(q), F̌(1) = f₁(q) are fifth-order mock theta functions themselves. So (Z↑), (Z↓) read

  ψ₁(q) = i·P·[ J(0) + I(0)·(φ₀(−p) − 1) + I(1/b)·φ₁(−p)/p ],

  f₁(q) = −2i·P·[ J(0) − I(0)·F₀(p) + I(1/b)·F₁(p) ].

The coefficients are modular: I(0) = c·D(0), I(1/b) = c·p^{1/2}·D(1/b), and with r = p^{1/2} = e^{−πi/(2x)}, G and H the Rogers–Ramanujan functions,

  D(0) = [ G(q⁴)·H(−r) + q^{1/4}·H(−q)·H(r) ] / ((q;q²)∞·(−p;p)∞),

  D(1/b) = [ G(q⁴)·G(−r) − q^{1/4}·H(−q)·G(r) ] / ((q;q²)∞·(−p;p)∞).

These rest on four identities of Rogers–Ramanujan type, each checked to order 300: Σ q^{m²}/(q)₂ₘ = G(q⁴)/(q;q²)∞, Σ q^{m²+m}/(q)₂ₘ₊₁ = H(−q)/(q;q²)∞, Σ r^{n²}/(r⁴;r⁴)ₙ = G(r)/(−r²;r²)∞, Σ r^{n²+2n}/(r⁴;r⁴)ₙ = H(r)/(−r²;r²)∞. (They are in Slater's list; I did not look up the numbers.) Note that r is the nome of −1/(4x) and q⁴ that of 4x.

### 3.2 How the structure changes (inference, with the verified facts it rests on)

- **Continuity.** All terms of Theorem Z are entire in a. Nothing degenerates in the formula.
- **The volume closes.** As a → 0 the two roots become real, W → −φ and 1/φ (φ the golden ratio), and Vol(e^{2πia}) = 2π·ln φ·|a| + O(a²) (derivation: dVol/da = 2π·Re w; verified: Vol/a = 3.02021, 3.02317, 3.023539, 3.023543 at a = 0.03, 0.01, 0.001, 0.0001 against 2π·ln φ = 3.023543). The selected saddle w = Log W moves to the end of the summation segment, Im w → π. The growth e^{Vol/|t|} turns off linearly in |a|, with a corner at a = 0.
- **State integrals become modular.** For z ≠ 1 the coefficients I((a−1)/b), I(a/b) grow and decay like e^{±Vol/|t|}. At z = 1 they are the Rogers–Ramanujan products above: bounded up to powers, with real critical values. This is the degenerate member of the family in the sense of the brief: volume 0, and all three sectors of the same order of magnitude.
- **Where the mock theta structure comes from.** From the dual series, not from the state integrals: 𝒲(T; p) and 𝒢(T; p) at T = p^{±1} are fifth-order mock theta functions of the dual nome, as at T = 1, p² they are W, W₃, F_a, F_b. For general a they are two-variable functions of (p, z̃); I do not know whether they are Appell–Lerch sums for general z̃.
- **The window.** For 0 < a < 1 the kernel of J(a/b) decays on both sides and P·J should be the Borel sum of a perturbative series (as proved for a = ½ in v17; not done for other a). At a = 0 the kernel 1/(1 + e^{2πs}) tends to 1 on the left: J(0) = I(0) − J(1/b), and the two choices a = 0, 1 differ by the modular term I(0). So at z = 1 the trivial sector is defined only up to a modular form, as a Mordell integral is.
- **Mordell integral.** The two formulas of §3.1 have the form of Watson's transformation formulas for the fifth-order functions (mock theta function = modular coefficients × mock theta functions of the dual nome + an integral). I did not find Watson's formula for ψ₁ in a form I could quote, and I did not identify i·P·J(0) with a classical Mordell integral. Open.

---

## 4. Is F(z) a quantum Jacobi form? (inference)

From the definitions as I remember them, without having re-read the papers.

- **Bringmann–Folsom** (quantum Jacobi forms): a function of (a, x) on rational points whose failure to transform like a Jacobi form, under x → −1/x with a → a/x and under the elliptic shifts, is a real-analytic function. F does not fit as a scalar: under x → −1/x it is not F(z̃; ·) that appears but 𝒲, with coefficients I that are transcendental and of size e^{±Vol/|t|}; and under the elliptic shift z → q·z it satisfies a second-order q-difference equation (v4, and (B1) of v12) rather than acquiring a factor. At z = ±1 with volume 0 or after projecting out a sector the scalar notion may apply; I did not try.
- **Garoufalidis–Zagier** (matrix-valued holomorphic quantum modular forms, with descendants and a meridian variable): the formula says that F(z; q) minus [state integrals bilinear in q- and p-series] × [the dual series 𝒲] equals i·P·e^{−πia}·J(a/b), which continues analytically across the real x-axis (v16 for a = ½; expected for 0 < a < 1). That is the defining shape of a holomorphic quantum modular form, with (1, 𝒲(T), 𝒲(p²T)) as the dual vector and a/x as the transformed Jacobi variable. So F(z) should be one row of a 3×3 matrix-valued object of that kind, in the family 𝔉(z, y; q).

To make this precise one would have to write the full matrix (the rows for the two solutions of the homogeneous q-difference equation) and check the cocycle property. Not done.

---

## 5. Numerical checks (verified; `run_log.txt`)

| check | script | result |
|---|---|---|
| (Z↑), (Z↓) at x = 0.3 + 0.2i, 1/(3 − 0.4i), 0.11 + 0.03i; z = i (a = ¼ and −¾), −ω (a = 5/6 and −1/6), −1, 1 (a = 0 and 1), ω, and a = 0.3 + 0.1i | `a1_formula.py` | 10⁻²⁹ or better in 50 of 54 cases; worst 3·10⁻²⁷ (a = −¾ at the third point, where \|T\| = 10⁻⁸) |
| D(a/b) through H; I(a/b) by quadrature against Theorem C | `a1_formula.py` | 10⁻³⁰; quadrature 10⁻²⁴ or better |
| strip form at z = i, −ω, two points | `c2_z1.py` | 3·10⁻³⁰ or better |
| saddle rule at roots of unity, ten values of z | `b1_roots.py` | table in §2.2 |
| Im V_k(w_k) = D(−W) and e^{w_k} = W on the branches k = −3 … 3, both roots, ten z | `b1_roots.py` | 4·10⁻⁴⁵ (working precision) |
| I((a−1)/b) and I(a/b) against their saddles; F against the m = 1 term | `b2_dominant.py` | §2.3 |
| series identities of §3.1 | `c1_series.py` | exact to order 300 |
| ψ₁, f₁ formulas with φ₀, φ₁, F₀, F₁; Rogers–Ramanujan forms of D(0), D(1/b), two points | `c2_z1.py` | 3·10⁻³⁰ or better |
| translation of the selected saddles of `census3/zfamily`; figure-eight ratio | `d1_compare.py` | §2.4, §2.5 |

---

## 6. Open

1. Steepest descent for I((a−1)/b): prove that the critical point of §2.2 is the relevant one (this would prove the volume rule). An alternative is the q → 1 asymptotics of H(y; q), a Nahm sum with A = ½ and a Jacobi variable.
2. The trivial sector for general a: v13–v17 for the kernel e^{2πas}/(1 + e^{2πs}), 0 < a < 1. The window argument of v17 should go through; the constants depend on a and degenerate at a → 0, 1.
3. P·J(0) and the Mordell integral of Watson's formulas; whether 𝒲(T; p) is an Appell–Lerch sum for general T.
4. The matrix of §4.
5. The items of v17 §7.

---

## 7. Scripts (`psi_star/v18/`)

- `common.py` (from v12).
- `a1_formula.py [digits]`: Theorem Z at three points and nine cases (12 minutes at 30 digits).
- `b1_roots.py [digits]`: F at roots of unity against the saddle-point formula with the rule of §2.2; Vol(z) near z = 1 (25 seconds).
- `b2_dominant.py [digits]`: the state integrals against their saddles near the real axis (10 seconds at 90 digits).
- `c1_series.py [N]`: exact series identities.
- `c2_z1.py [digits]`: strip form; the case z = 1 (1 minute).
- `d1_compare.py`: comparison with `census3/zfamily` and with the Kashaev invariant of 4₁ (seconds).
- `run_all.sh`, `run_log.txt`, `README.md`.
