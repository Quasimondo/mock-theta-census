"""Exact expansion of d1151(e^{-t}) = sum_{m,n>=0} q^{m^2-mn+n^2} / ((-q;q)_m (-q;q)_n) at q = 1.
log prod_{k<=m} 1/(1+e^{-kt}) = -m log 2 - sum_k L(kt),  L(x) = log((1+e^{-x})/2) = -x/2 + x^2/8 - x^4/192 + ...
sum_k (kt)^j is a Faulhaber polynomial in m, so every t^K coefficient is a polynomial in m, n;
sums against 2^{-m-n} use sum_m m^a 2^{-m} = 2 * Fubini(a)."""
import sympy as sp, sys, json
K = int(sys.argv[1]) if len(sys.argv) > 1 else 40   # 40 terms take about 3 minutes
t, x, m, n = sp.symbols('t x m n')
Lser = sp.series(sp.log((1 + sp.exp(-x)) / 2), x, 0, K + 1).removeO()
Lc = [Lser.coeff(x, j) for j in range(K + 1)]
def faul(j, v):  # sum_{k=1}^{v} k^j
    k = sp.symbols('k')
    return sp.expand(sp.summation(k ** j, (k, 1, v)))
Fm = [faul(j, m) for j in range(K + 1)]; Fn = [faul(j, n) for j in range(K + 1)]
# exponent as series in t (without the -(m+n) log 2 part)
Q = m * m - m * n + n * n
expo = -t * Q - sum(Lc[j] * t ** j * (Fm[j] + Fn[j]) for j in range(1, K + 1))
# exp(expo) truncated at t^K, coefficients polynomials in m,n
ser = [sp.Integer(1)] + [0] * K
P = sp.Poly(sp.expand(expo), t)
e = [sp.expand(P.coeff_monomial(t ** j)) for j in range(K + 1)]
# exp via recursion: E' = expo' E  ->  k E_k = sum_{j=1}^{k} j e_j E_{k-j}
E = [sp.Integer(1)]
for k in range(1, K + 1):
    E.append(sp.expand(sum(j * e[j] * E[k - j] for j in range(1, k + 1)) / k))
fub = {}
def F(a):
    if a not in fub:
        fub[a] = sum(sp.factorial(j) * sp.functions.combinatorial.numbers.stirling(a, j) for j in range(a + 1))
    return fub[a]
S = lambda a: 2 * F(a)  # sum_{m>=0} m^a / 2^m
c = []
for k in range(K + 1):
    poly = sp.Poly(E[k], m, n)
    val = sum(coef * S(a) * S(b) for (a, b), coef in poly.terms())
    c.append(sp.nsimplify(val))
print("d1151(e^{-t}) = sum_k a_k t^k with")
for k, v in enumerate(c):
    print(k, v, float(v))
json.dump([str(v) for v in c], open('a_k.json', 'w'))  # used by largeorder.py
