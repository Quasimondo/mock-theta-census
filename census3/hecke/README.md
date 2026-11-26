# Hecke-type check of the bounded twists

- `sparse.py`: integer q-series of the bounded twists (ψ₁*, ψ₀*, ψ*, F₀*, F₁*, T₀*) to q⁴⁰⁰, multiplied by every
  eta quotient (q;q)∞^a (q²;q²)∞^b (q⁴;q⁴)∞^c with −2 ≤ a, b, c ≤ 2; prints the sparsest products.
- `ref.py`: the same for known Hecke-type cases (f₀(q)(q;q)∞: 160/400 nonzero; ψ₁(q)(q²;q²)∞: 143/400).
- `hmatch.py`, `hval.py`: search for a single f_{a,b,c}(±q^α, ±q^β; q) (Hickerson–Mortenson notation, a, c ≤ 6)
  equal to a target up to q^s and a constant; `hval.py` checks it finds ψ₁(q)(q;q)∞ = f_{2,4,3}(1, 1; q).

Result: ψ* (twist of the third-order ψ) is sparse as it stands and equals −σ*(q)/2, the Andrews–Dyson–Hickerson
companion of σ(q) (Q(√6), quantum modular). The classifier gives σ(q) itself `QQQQQQQQQQQQQQQQQQQ` as well.
The dense twists stay dense (≥ 386/400 nonzero) under every eta quotient tried, so none of them is a single
Hecke-type double sum over such a product. No f_{a,b,c} with mock-type signs matches ψ*, as expected, since its
Hecke form has the "false" sign pattern.
