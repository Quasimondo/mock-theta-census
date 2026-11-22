"""Universal mock theta functions and the two-variable family H(x) containing both candidates."""
from flint import acb, arb, ctx
import catalog

def omega():
    return (acb(2) / 3).exp_pi_i()

def Hx(x):
    # H(x,q) = sum q^{n(n-1)/2} (-q;q)_n / ((xq;q)_n (q/x;q)_n)
    def fn(q):
        X = x()
        return catalog.eulerian(q, lambda n: n * (n - 1) // 2,
                                lambda n, Q: (1 + Q(n)) / ((1 - X * Q(n)) * (1 - Q(n) / X)))
    return fn

def g2(x):
    # g2(x,q) = sum (-q;q)_n q^{n(n+1)/2} / ((x;q)_{n+1} (q/x;q)_{n+1})
    def fn(q):
        X = x()
        return catalog.eulerian(q, lambda n: n * (n + 1) // 2,
                                lambda n, Q: (1 + Q(n)) / ((1 - X * Q(n)) * (1 - Q(n + 1) / X)),
                                init=lambda Q: 1 / ((1 - X) * (1 - Q(1) / X)))
    return fn

def g3(x):
    # g3(x,q) = sum q^{n(n+1)} / ((x;q)_{n+1} (q/x;q)_{n+1})
    def fn(q):
        X = x()
        return catalog.eulerian(q, lambda n: n * (n + 1),
                                lambda n, Q: 1 / ((1 - X * Q(n)) * (1 - Q(n + 1) / X)),
                                init=lambda Q: 1 / ((1 - X) * (1 - Q(1) / X)))
    return fn

W = omega
mW = lambda: -omega()
W2 = lambda: omega() ** 2
mW2 = lambda: -(omega() ** 2)
FUN = {
    "F=H(-w)": Hx(mW), "G=H(w)": Hx(W),
    "g2(-w)": g2(mW), "g2(-w2)": g2(mW2), "g2(w)": g2(W), "g2(w2)": g2(W2),
    "g3(-w)": g3(mW), "g3(-w2)": g3(mW2), "g3(w)": g3(W), "g3(w2)": g3(W2),
}
for k, v in FUN.items():
    catalog.BY[k] = dict(fn=v)
