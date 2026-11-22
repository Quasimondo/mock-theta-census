"""Predict the q->1 growth exponent of a plain (non-alternating) Eulerian series from its formula alone.
log(term_n) ~ (1/t) * Phi(x),  x = n t,  Phi(x) = -A x^2/2 + sum_f p_f/d_f * [Li2(-s_f e^{-d_f x}) - Li2(-s_f)]
c = max_x Phi(x),  kappa = c / (4 pi^2)."""
import mpmath as mp
mp.mp.dps = 30
# factor -> list of (d, sigma, weight) with (1 + sigma q^{d k}) ** weight
DEC = {
 "(1+q^k)": [(1, 1, 1)], "(1-q^k)": [(1, -1, 1)],
 "(1+q^2k)": [(2, 1, 1)], "(1-q^2k)": [(2, -1, 1)],
 "(1+q^(2k-1))": [(2, 1, 1)], "(1-q^(2k-1))": [(2, -1, 1)],
 "(1-q^k+q^2k)": [(3, 1, 1), (1, 1, -1)],
 "(1+q^k+q^2k)": [(3, -1, 1), (1, -1, -1)],
}
def phi(r):
    A, pat = r['A'], r['pat']
    parts = [(d, s, w * p) for f, p in pat for (d, s, w) in DEC[f]]
    def F(x):
        v = -A * x * x / 2
        for d, s, w in parts:
            v += mp.mpf(w) / d * (mp.polylog(2, -s * mp.e ** (-d * x)) - mp.polylog(2, -s))
        return mp.re(v)
    return F
def predict(r):
    F = phi(r)
    xs = [mp.mpf(i) / 50 for i in range(1, 400)]
    vals = [F(x) for x in xs]
    i = max(range(len(xs)), key=lambda k: vals[k])
    best = max(vals[i], 0)
    if vals[i] > 0:
        try:
            x0 = mp.findroot(lambda x: mp.diff(F, x), xs[i])
            best = max(best, F(x0))
        except Exception:
            pass
    return best / (4 * mp.pi ** 2)
