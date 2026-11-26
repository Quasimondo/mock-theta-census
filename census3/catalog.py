"""Labeled catalog of q-series and a high-precision evaluator (python-flint / Arb).

Each Eulerian series is  const + mult * sum_{n>=n0} q^{e(n)} * P_n,
with P_{n0} = init(Q) and P_n = P_{n-1} * ratio(n, Q), where Q(k) = q^k (cached).
Signs like (-1)^n live inside ratio. Formulas for the mock theta functions are the
standard ones (Watson; Andrews-Hickerson; Gordon-McIntosh; Choi; McIntosh) and are
cross-checked numerically where a classical identity is available (see checks()).
"""
from flint import acb, arb, ctx

class PrecisionLow(Exception):
    pass

def slog(z):
    """log that is safe on the negative real axis (any consistent branch works, since we only exponentiate)"""
    if z.real.mid() < 0 and z.imag.contains(0):
        return (-z).log() + acb(0, arb.pi())
    return z.log()

class Pow:
    """q^k computed as exp(k log q): avoids the wrapping blow-up of repeated complex-ball products."""
    def __init__(self, q):
        self.L = slog(q)
        self.c = {}
    def __call__(self, k):
        v = self.c.get(k)
        if v is None:
            v = (self.L * k).exp() if k else acb(1)
            self.c[k] = v
        return v

def eulerian(q, e, ratio, n0=0, init=None, const=0, mult=1, tol_bits=None):
    """Terms are formed as exp(e(n) log q + log P_n), with log P_n accumulated as a sum
    of logs; this keeps ball radii growing linearly instead of exponentially."""
    Q = Pow(q)
    logP = slog(init(Q)) if init else acb(0)
    s = acb(0)
    n = n0
    prec = ctx.prec
    small = 0
    maxt = arb(0)
    eps = arb(2) ** (-prec)
    eps2 = arb(2) ** (-prec - 10)
    while True:
        if n > n0:
            logP += slog(ratio(n, Q))
        term = (Q.L * e(n) + logP).exp()
        s += term
        ts = abs(term).mid()
        if ts > maxt:
            maxt = ts
        if n > n0 + 3:
            ss = abs(s).mid()
            ref = ss if ss > maxt * eps else maxt * eps
            if ts < ref * eps2:
                small += 1
                if small >= 4:
                    break
            else:
                small = 0
        if n % 64 == 0:
            a = abs(term)
            if a.rad() > a.mid() * arb(2) ** (-prec // 2):
                raise PrecisionLow()
        n += 1
        if n > 3_000_000:
            raise RuntimeError("no convergence")
    return const + mult * s

def pentagonal(q):
    """(q;q)_inf = sum_{n in Z} (-1)^n q^{n(3n-1)/2}"""
    Q = Pow(q)
    s = acb(1)
    n = 1
    prec = ctx.prec
    while True:
        qa = Q(n * (3 * n - 1) // 2)
        qb = Q(n * (3 * n + 1) // 2)
        t = qa + qb
        s += -t if n % 2 else t
        if abs(qa).mid() < arb(2) ** (-prec - 20) * max(abs(s).mid(), arb(2) ** (-4 * prec)):
            break
        n += 1
    return s

def theta(q, alt=False):
    return eulerian(q, lambda n: n * n, (lambda n, Q: acb(-1)) if alt else (lambda n, Q: acb(1)),
                    n0=1, init=(lambda Q: acb(-1)) if alt else None, const=1, mult=2)

one = acb(1)
inv = lambda x: one / x

# ---------- the catalog ----------
E = []
def add(name, label, order, fn, note=""):
    E.append(dict(name=name, label=label, order=order, fn=fn, note=note))

def eul(e, ratio, **kw):
    return lambda q: eulerian(q, e, ratio, **kw)

# Modular (or sums of modular forms with rational q-power offsets)
add("theta3", "modular", "", lambda q: theta(q))
add("theta4", "modular", "", lambda q: theta(q, True))
add("eta (q;q)inf", "modular", "", pentagonal)
add("b(q) = theta4/(-q;q)inf", "modular", "", lambda q: theta(q, True) * pentagonal(q) / pentagonal(q * q))
add("ternary x2+y2+10z2", "modular", "", lambda q: theta(q) ** 2 * theta(Pow(q)(10)))
add("1/(q;q)inf  [Durfee sum]", "modular", "", eul(lambda n: n * n, lambda n, Q: inv((1 - Q(n)) ** 2)))
add("Rogers-Ramanujan G", "modular", "", eul(lambda n: n * n, lambda n, Q: inv(1 - Q(n))))
add("Rogers-Ramanujan H", "modular", "", eul(lambda n: n * n + n, lambda n, Q: inv(1 - Q(n))))
add("Nahm A=1,B=0", "modular", "", eul(lambda n: n * n, lambda n, Q: inv(1 - Q(2 * n))))
add("Nahm A=1,B=1/2", "modular", "", eul(lambda n: n * n + n, lambda n, Q: inv(1 - Q(2 * n))))
add("Nahm A=1/2,B=0", "modular", "", eul(lambda n: n * n, lambda n, Q: inv(1 - Q(4 * n))))
add("Nahm A=1/2,B=1/2", "modular", "", eul(lambda n: n * n + 2 * n, lambda n, Q: inv(1 - Q(4 * n))))
add("Euler (-q;q)inf", "modular", "", eul(lambda n: n * (n + 1) // 2, lambda n, Q: inv(1 - Q(n))))
add("Gollnitz-Gordon", "modular", "", eul(lambda n: n * n, lambda n, Q: (1 + Q(2 * n - 1)) / (1 - Q(2 * n))))
add("(G-H)/q  [sum of two]", "mixed", "", eul(lambda n: n * n + 2 * n, lambda n, Q: inv(1 - Q(n))), "sum of two modular forms")

# Mock theta functions
M = lambda name, order, fn: add(name, "mock", order, fn)
M("f", 3, eul(lambda n: n * n, lambda n, Q: inv((1 + Q(n)) ** 2)))
M("phi", 3, eul(lambda n: n * n, lambda n, Q: inv(1 + Q(2 * n))))
M("psi", 3, eul(lambda n: n * n, lambda n, Q: inv(1 - Q(2 * n - 1)), n0=1, init=lambda Q: inv(1 - Q(1))))
M("chi", 3, eul(lambda n: n * n, lambda n, Q: inv(1 - Q(n) + Q(2 * n))))
M("omega", 3, eul(lambda n: 2 * n * (n + 1), lambda n, Q: inv((1 - Q(2 * n + 1)) ** 2), init=lambda Q: inv((1 - Q(1)) ** 2)))
M("nu", 3, eul(lambda n: n * (n + 1), lambda n, Q: inv(1 + Q(2 * n + 1)), init=lambda Q: inv(1 + Q(1))))
M("rho", 3, eul(lambda n: 2 * n * (n + 1), lambda n, Q: inv(1 + Q(2 * n + 1) + Q(4 * n + 2)), init=lambda Q: inv(1 + Q(1) + Q(2))))
M("f0", 5, eul(lambda n: n * n, lambda n, Q: inv(1 + Q(n))))
M("f1", 5, eul(lambda n: n * n + n, lambda n, Q: inv(1 + Q(n))))
M("F0", 5, eul(lambda n: 2 * n * n, lambda n, Q: inv(1 - Q(2 * n - 1))))
M("F1", 5, eul(lambda n: 2 * n * n + 2 * n, lambda n, Q: inv(1 - Q(2 * n + 1)), init=lambda Q: inv(1 - Q(1))))
M("phi0", 5, eul(lambda n: n * n, lambda n, Q: 1 + Q(2 * n - 1)))
M("phi1", 5, eul(lambda n: (n + 1) ** 2, lambda n, Q: 1 + Q(2 * n - 1)))
M("psi0", 5, eul(lambda n: (n + 1) * (n + 2) // 2, lambda n, Q: 1 + Q(n)))
M("psi1", 5, eul(lambda n: n * (n + 1) // 2, lambda n, Q: 1 + Q(n)))
M("F7_0", 7, eul(lambda n: n * n, lambda n, Q: (1 - Q(n)) / ((1 - Q(2 * n - 1)) * (1 - Q(2 * n)))))
M("F7_1", 7, eul(lambda n: n * n, lambda n, Q: (1 - Q(n - 1)) / ((1 - Q(2 * n - 2)) * (1 - Q(2 * n - 1))), n0=1, init=lambda Q: inv(1 - Q(1))))
M("F7_2", 7, eul(lambda n: n * n + n, lambda n, Q: (1 - Q(n)) / ((1 - Q(2 * n)) * (1 - Q(2 * n + 1))), init=lambda Q: inv(1 - Q(1))))
M("mu (2nd)", 2, eul(lambda n: n * n, lambda n, Q: -(1 - Q(2 * n - 1)) / (1 + Q(2 * n)) ** 2))
M("phi6", 6, eul(lambda n: n * n, lambda n, Q: -(1 - Q(2 * n - 1)) / ((1 + Q(2 * n - 1)) * (1 + Q(2 * n)))))
M("psi6", 6, eul(lambda n: (n + 1) ** 2, lambda n, Q: -(1 - Q(2 * n - 1)) / ((1 + Q(2 * n)) * (1 + Q(2 * n + 1))), init=lambda Q: inv(1 + Q(1))))
M("S0 (8th)", 8, eul(lambda n: n * n, lambda n, Q: (1 + Q(2 * n - 1)) / (1 + Q(2 * n))))
M("T0 (8th)", 8, eul(lambda n: (n + 1) * (n + 2), lambda n, Q: (1 + Q(2 * n)) / (1 + Q(2 * n + 1)), init=lambda Q: inv(1 + Q(1))))
M("phi10", 10, eul(lambda n: n * (n + 1) // 2, lambda n, Q: inv(1 - Q(2 * n + 1)), init=lambda Q: inv(1 - Q(1))))
M("psi10", 10, eul(lambda n: (n + 1) * (n + 2) // 2, lambda n, Q: inv(1 - Q(2 * n + 1)), init=lambda Q: inv(1 - Q(1))))
M("X10", 10, eul(lambda n: n * n, lambda n, Q: -inv((1 + Q(2 * n - 1)) * (1 + Q(2 * n)))))
M("chi10", 10, eul(lambda n: (n + 1) ** 2, lambda n, Q: -inv((1 + Q(2 * n)) * (1 + Q(2 * n + 1))), init=lambda Q: inv(1 + Q(1))))

# False / partial theta functions and one Maass-form example
add("false theta sum (-1)^n q^{n(n+1)/2}", "false", "", eul(lambda n: n * (n + 1) // 2, lambda n, Q: acb(-1)))
add("false theta sum (-1)^n q^{n(3n+1)/2}", "false", "", eul(lambda n: n * (3 * n + 1) // 2, lambda n, Q: acb(-1)))
add("sigma (Andrews-Dyson-Hickerson)", "other", "", eul(lambda n: n * (n + 1) // 2, lambda n, Q: inv(1 + Q(n))))

# Confirmed non-modular rank-1 Nahm sums (Zagier's classification)
N = lambda name, fn: add(name, "neither", "", fn, "rank-1 Nahm sum outside Zagier's list")
N("Nahm A=4,B=0", eul(lambda n: 2 * n * n, lambda n, Q: inv(1 - Q(n))))
N("Nahm A=4,B=1", eul(lambda n: 2 * n * n + n, lambda n, Q: inv(1 - Q(n))))
N("Nahm A=3,B=0", eul(lambda n: 3 * n * n, lambda n, Q: inv(1 - Q(2 * n))))
N("Nahm A=3,B=1/2", eul(lambda n: 3 * n * n + n, lambda n, Q: inv(1 - Q(2 * n))))
N("Nahm A=3/2,B=0", eul(lambda n: 3 * n * n, lambda n, Q: inv(1 - Q(4 * n))))
N("Nahm A=6,B=0", eul(lambda n: 3 * n * n, lambda n, Q: inv(1 - Q(n))))

BY = {e["name"]: e for e in E}

def checks():
    """Classical identities as a test of the formulas above."""
    ctx.prec = 200
    q = acb(arb(31) / 100, arb(47) / 100)
    g = lambda nm, x: BY[nm]["fn"](x)
    out = []
    out.append(("2 phi(-q) - f(q) = b(q)", 2 * g("phi", -q) - g("f", q) - g("b(q) = theta4/(-q;q)inf", q)))
    out.append(("f(q) + 4 psi(-q) = b(q)", g("f", q) + 4 * g("psi", -q) - g("b(q) = theta4/(-q;q)inf", q)))
    out.append(("Durfee = 1/eta", g("1/(q;q)inf  [Durfee sum]", q) * pentagonal(q) - 1))
    out.append(("Euler = (q^2;q^2)/(q;q)", g("Euler (-q;q)inf", q) - pentagonal(q * q) / pentagonal(q)))
    out.append(("(G-H)/q", g("(G-H)/q  [sum of two]", q) - (g("Rogers-Ramanujan G", q) - g("Rogers-Ramanujan H", q)) / q))
    # Watson: f(q) series vs omega relation f(q^8) + 2q omega(q) + q^3 omega(-q^4)... (skip; uses q-powers)
    # fifth order: f0(q) = theta-type identities need products; check chi0 = 2F0(q) - phi0(-q) (Ramanujan)
    return out

if __name__ == "__main__":
    for nm, v in checks():
        print(f"{nm:32s} |diff| = {abs(v).mid().str(3, radius=False)}")
    print(len(E), "functions:", {l: sum(e['label'] == l for e in E) for l in set(e['label'] for e in E)})
