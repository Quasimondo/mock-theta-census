"""High-precision evaluation of a rank-2 candidate with Arb balls (python-flint).

Every K-th term is computed exactly as exp(E(m,n) log q + log P1_m + log P2_n + log P12_{m+n}),
with the product logs accumulated as sums of logs (no wrapping blow-up), as in catalog.eulerian.
The terms in between come from their neighbour by one multiplication, which is much cheaper
than an exponential; restarting every K terms keeps the ball radii tight.
"""
import math
from flint import acb, arb, ctx
from catalog import Pow, slog, PrecisionLow
from gen2 import FACT, coupling, expo

def make_fn(c):
    f1, f2, f12 = FACT[c["f1"]], FACT[c["f2"]], (FACT[coupling(c["cp"])] if coupling(c["cp"]) else None)
    s1, s2 = c["s1"], c["s2"]

    def fn(q):
        Q = Pow(q)
        prec = ctx.prec
        eps = arb(2) ** (-prec - 10)
        logs = {1: [acb(0)], 2: [acb(0)], 3: [acb(0)]}
        def plog(which, f, j):
            L = logs[which]
            if f is None:
                return L[0]
            while len(L) <= j:
                k = len(L)
                g = slog(f[1](k, Q))
                L.append(L[-1] + (g if f[3] > 0 else -g))
            return L[j]
        s = acb(0)
        maxt = arb(0)
        # Near a root of unity of order p the terms are not monotone: within a row they can
        # start tiny and grow by many orders of magnitude before they decay. So the stopping
        # rule works on windows of W = 2p consecutive terms (rows): stop only when the window
        # maxima have been decreasing for 3 windows AND are negligible.
        ang = float(q.arg().mid()) / (2 * math.pi) % 1.0
        p = min(range(1, 15), key=lambda d: abs(ang * d - round(ang * d)))
        W = 2 * p
        K = 16
        rows, m = [], 0
        while True:
            win, wins, n = arb(0), [], 0
            row_max = arb(0)
            while True:
                if n % K == 0:
                    # exact restart from the logs every K terms keeps the ball radii tight
                    t = (Q.L * expo(c, m, n) + plog(1, f1, m) + plog(2, f2, n) + plog(3, f12, m + n)).exp()
                    if (s1 < 0 and m % 2) != (s2 < 0 and n % 2):
                        t = -t
                else:
                    # neighbouring term: one multiplication instead of an exponential
                    d = (2 * c["b"] * m + c["c"] * (2 * n - 1) + c["e"]) // 2
                    t = t * (Q(d) if d >= 0 else 1 / Q(-d))
                    g = f2[1](n, Q)
                    t = t * g if f2[3] > 0 else t / g
                    if f12 is not None:
                        g = f12[1](m + n, Q)
                        t = t * g if f12[3] > 0 else t / g
                    if s2 < 0:
                        t = -t
                s += t
                a = abs(t).mid()
                if a > win:
                    win = a
                if a > row_max:
                    row_max = a
                if a > maxt:
                    maxt = a
                if (m + n) % 32 == 0:
                    ab = abs(t)
                    if ab.rad() > ab.mid() * arb(2) ** (-prec // 2):
                        raise PrecisionLow()
                n += 1
                if n % W == 0:
                    wins.append(win); win = arb(0)
                    ref = max(abs(s).mid(), maxt * arb(2) ** (-prec))
                    if len(wins) >= 4 and wins[-1] < ref * eps and wins[-1] < wins[-2] < wins[-3] < wins[-4]:
                        break
                if n > 400000:
                    raise RuntimeError("no convergence")
            rows.append(row_max)
            m += 1
            if m % W == 0 and m >= 4 * W:
                ref = max(abs(s).mid(), maxt * arb(2) ** (-prec))
                w = [max(rows[i:i + W]) for i in range(len(rows) - 4 * W, len(rows), W)]
                if w[-1] < ref * eps and w[-1] < w[-2] < w[-3] < w[-4]:
                    break
            if m > 400000:
                raise RuntimeError("no convergence")
        return s
    return fn


def selftest():
    """Checks the evaluator near roots of unity against two exact closed forms:
    a double sum that is exactly the polynomial 1 - q - q^2 + q^5 (terms cancel massively near the rim),
    and the Andrews-Gordon sum, which equals prod 1/(1-q^n) over n not = 0, +-3 mod 7."""
    from fractions import Fraction
    poly = dict(a=3, b=2, c=2, d=-1, e=4, s1=-1, s2=1, f1="1/(q)", f2="1/(q)", cp="")
    ag = dict(a=2, b=2, c=4, d=0, e=0, s1=1, s2=1, f1="1/(q)", f2="1/(q)", cp="")
    ok = True
    for (h, m, t) in [(2, 5, Fraction(1, 5000)), (1, 2, Fraction(1, 800)), (0, 1, Fraction(1, 200)), (1, 3, Fraction(1, 2000))]:
        for prec in (512, 1024, 2048):
            ctx.prec = prec
            q = (-(arb(t.numerator) / t.denominator)).exp() * ((acb(2 * h) / m).exp_pi_i() if m > 1 else acb(1))
            try:
                v = make_fn(poly)(q)
            except PrecisionLow:
                continue
            if abs(v).rad() < abs(v).mid() * arb(2) ** -100:
                break
        d1 = abs(v - (1 - q - q ** 2 + q ** 5)).mid()
        ctx.prec = 512
        w = make_fn(ag)(q)
        Q = Pow(q)
        L, n = acb(0), 1
        while True:
            qn = Q(n)
            if n % 7 not in (0, 3, 4):
                L -= slog(1 - qn)
            if abs(qn).mid() < arb(2) ** -600:
                break
            n += 1
        P = L.exp()
        d2 = (abs(w - P) / abs(P)).mid()  # relative: the sum is huge near some roots
        good = d1 < arb(10) ** -25 and d2 < arb(10) ** -25
        ok &= good
        print(f"zeta = e^(2 pi i {h}/{m}), t = {t}: polynomial diff {d1.str(2)}, Andrews-Gordon rel. diff {d2.str(2)}  {'OK' if good else 'FAIL'}")
    print("selftest", "passed" if ok else "FAILED")
    return ok


if __name__ == "__main__":
    selftest()
