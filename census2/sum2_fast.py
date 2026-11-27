"""Fast evaluation of a rank-2 candidate: the double sum as ONE polynomial product.

    sum_{m,n} s1^m s2^n q^{E(m,n)} P1_m P2_n P12_{m+n},   2E = a m^2 + 2b mn + c n^2 + d m + e n

The cross term is the only thing that couples m and n, and  mn = ((m+n)^2 - m^2 - n^2)/2,  so

    term(m,n) = A_m * B_n * C_{m+n}
    A_m = s1^m q^{E(m,0) - b m^2/2} P1_m,   B_n = s2^n q^{E(0,n) - b n^2/2} P2_n,   C_k = q^{b k^2/2} P12_k

and the whole double sum is   sum_k C_k * (A conv B)_k.   The convolution is one acb_poly
multiplication inside FLINT (C code, ball arithmetic, built for coefficients of very different
size), so the Python loop over N^2 terms is gone: what is left is O(N) exponentials for the three
one-dimensional sequences.  Every A_m, B_n, C_k is formed exactly as exp(exponent * log q + sum of
logs), as in catalog.eulerian, so there is no restart logic and no wrapping blow-up.

Drop-in for sum2.make_fn:   import sum2_fast as sum2

    python3 sum2_fast.py               closed-form selftest (same two identities as sum2.selftest)
    python3 sum2_fast.py compare 30    old sum2.py against this file on 30 random candidates of cands2.json
"""
import math
import numpy as np
from flint import acb, arb, acb_poly, ctx
from catalog import Pow, slog, PrecisionLow
from gen2 import FACT, coupling, expo

GUARD = 40          # terms below (largest term) * 2^-(prec+GUARD) are left out
CHUNK = 1 << 22     # floats per block in the size scan (keeps the scan under ~100 MB per worker)


class _Seq:
    """One of the sequences A, B, C: values as acb, log-sizes as floats, grown on demand."""
    def __init__(self, Q, e2, f, sign):
        self.Lh = Q.L / 2          # exponents are kept doubled, so half-integers are fine
        self.Q, self.e2, self.f, self.sign = Q, e2, f, sign
        self.logP = acb(0)
        self.v, self.mag = [], []

    def extend(self, n):
        Q, f, v, mag = self.Q, self.f, self.v, self.mag
        for k in range(len(v), n):
            if f is not None and k > 0:
                g = slog(f[1](k, Q))
                self.logP = self.logP + g if f[3] > 0 else self.logP - g
            z = self.Lh * self.e2(k) + self.logP
            x = z.exp()
            v.append(-x if (self.sign < 0 and k % 2) else x)
            lm = float(z.real.mid())
            if not math.isfinite(lm):
                raise PrecisionLow()   # a factor's ball touches zero: the caller raises the precision
            mag.append(lm)


def _scan(a, b, c):
    """row and column maxima of  a[m] + b[n] + c[m+n]  over the rectangle, without storing it"""
    NA, NB = len(a), len(b)
    win = np.lib.stride_tricks.sliding_window_view(c, NB)[:NA]   # win[m] = c[m : m+NB], a view
    rowmax = np.empty(NA)
    colmax = np.full(NB, -np.inf)
    step = max(1, CHUNK // NB)
    for i in range(0, NA, step):
        T = a[i:i + step, None] + b[None, :] + win[i:i + step]
        rowmax[i:i + step] = T.max(axis=1)
        np.maximum(colmax, T.max(axis=0), out=colmax)
    return rowmax, colmax


def make_fn(c):
    f1, f2 = FACT[c["f1"]], FACT[c["f2"]]
    f12 = FACT[coupling(c["cp"])] if coupling(c["cp"]) else None
    s1, s2, b = c["s1"], c["s2"], c["b"]
    E0 = expo(c, 0, 0)
    for (m, n) in ((1, 1), (2, 3), (5, 4), (7, 11)):   # the split below needs exactly this form
        if expo(c, m, n) != expo(c, m, 0) + expo(c, 0, n) - E0 + b * m * n:
            raise ValueError("exponent is not of the form E(m,0) + E(0,n) + b*m*n")
    eA = lambda m: 2 * (expo(c, m, 0) - E0) - b * m * m
    eB = lambda n: 2 * expo(c, 0, n) - b * n * n
    eC = lambda k: b * k * k

    def fn(q):
        prec = ctx.prec
        Q = Pow(q)
        A, B, C = _Seq(Q, eA, f1, s1), _Seq(Q, eB, f2, s2), _Seq(Q, eC, f12, 1)
        # same period estimate as sum2.py: near a root of unity of order p the sizes oscillate with period p
        ang = float(q.arg().mid()) / (2 * math.pi) % 1.0
        p = min(range(1, 15), key=lambda d: abs(ang * d - round(ang * d)))
        W = 2 * p
        NA = NB = max(64, 8 * W)
        drop = (prec + GUARD) * math.log(2)
        while True:
            A.extend(NA); B.extend(NB); C.extend(NA + NB)
            rowmax, colmax = _scan(np.array(A.mag), np.array(B.mag), np.array(C.mag))
            top = rowmax.max()
            cut = top - drop
            M1 = int(np.nonzero(rowmax > cut)[0][-1]) + 1
            N1 = int(np.nonzero(colmax > cut)[0][-1]) + 1
            # a direction is finished when at least 4 windows beyond the kept block are
            # negligible and still falling; otherwise double that direction only
            def done(mx, kept, size):
                if kept + 4 * W > size:
                    return False
                w = [mx[size - (j + 1) * W:size - j * W].max() for j in range(4)]   # last window first
                return w[0] < w[1] < w[2] < w[3]
            okA, okB = done(rowmax, M1, NA), done(colmax, N1, NB)
            if okA and okB:
                break
            if not okA:
                NA *= 2
            if not okB:
                NB *= 2
            if max(NA, NB) > 400000:
                raise RuntimeError("no convergence")
        P = acb_poly(A.v[:M1]) * acb_poly(B.v[:N1])
        s = acb(0)
        Cv = C.v
        for k in range(M1 + N1 - 1):
            s += Cv[k] * P[k]
        # what was left out, as part of the error ball, so that callers that test the radius
        # (probe.evaluate does) also cover the truncation
        tail = (arb(cut) + math.log(64.0 * (M1 + N1))).exp()
        fn.info = dict(M=M1, N=N1, logmax=top / math.log(2))
        return s + acb(arb(0, tail), arb(0, tail))
    return fn


def _q(h, m, t):
    r = (-(arb(t.numerator) / t.denominator)).exp()
    return r * ((acb(2 * h) / m).exp_pi_i() if m > 1 else acb(1))


def selftest():
    """Same two exact closed forms as sum2.selftest: a double sum that is exactly 1 - q - q^2 + q^5
    (massive cancellation near the rim) and the Andrews-Gordon sum = prod 1/(1-q^n), n != 0, +-3 mod 7."""
    from fractions import Fraction
    poly = dict(a=3, b=2, c=2, d=-1, e=4, s1=-1, s2=1, f1="1/(q)", f2="1/(q)", cp="")
    ag = dict(a=2, b=2, c=4, d=0, e=0, s1=1, s2=1, f1="1/(q)", f2="1/(q)", cp="")
    ok = True
    for (h, m, t) in [(2, 5, Fraction(1, 5000)), (1, 2, Fraction(1, 800)), (0, 1, Fraction(1, 200)), (1, 3, Fraction(1, 2000))]:
        for prec in (512, 1024, 2048):
            ctx.prec = prec
            q = _q(h, m, t)
            v = make_fn(poly)(q)
            if abs(v).rad() < abs(v).mid() * arb(2) ** -100:
                break
        d1 = abs(v - (1 - q - q ** 2 + q ** 5)).mid()
        ctx.prec = 512
        q = _q(h, m, t)
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
        d2 = (abs(w - P) / abs(P)).mid()
        good = d1 < arb(10) ** -25 and d2 < arb(10) ** -25
        ok &= good
        print(f"zeta = e^(2 pi i {h}/{m}), t = {t}: polynomial diff {d1.str(2)}, Andrews-Gordon rel. diff {d2.str(2)}  {'OK' if good else 'FAIL'}")
    print("selftest", "passed" if ok else "FAILED")
    return ok


def compare(count=30, seed=0, limit=120):
    """Evaluates random candidates of cands2.json with sum2.py and with this file at the same points.
    'agree' means the two error balls overlap, i.e. the values are equal within their own error bounds."""
    import json, random, signal, time
    from fractions import Fraction
    import sum2
    class Stop(BaseException):
        pass

    def alarm(*a):
        raise Stop()
    signal.signal(signal.SIGALRM, alarm)
    reps = json.load(open("cands2.json"))
    pick = random.Random(seed).sample(reps, min(count, len(reps)))
    pts = [(0, 1, Fraction(1, 200), 256), (1, 2, Fraction(1, 800), 512), (1, 3, Fraction(1, 2000), 512),
           (1, 4, Fraction(1, 3000), 1024), (2, 5, Fraction(1, 5000), 1024)]
    bits = lambda x: float(x.log().mid()) / math.log(2) if x > 0 else float("-inf")
    bad = n = 0
    told = tnew = 0.0
    for c in pick:
        for (h, m, t, prec) in pts:
            ctx.prec = prec
            q = _q(h, m, t)
            t0 = time.time()
            fnew = make_fn(c)
            w = fnew(q)
            t1 = time.time()
            signal.alarm(limit)
            try:
                v = sum2.make_fn(c)(q)
            except PrecisionLow:
                v = "needs more precision"
            except Stop:
                v = f"over {limit} s"
            finally:
                signal.alarm(0)
            t2 = time.time()
            if isinstance(v, str):
                print(f"id {c['id']} at {h}/{m}: old evaluator {v}; new {t1 - t0:.2f} s, relative radius 2^{bits(abs(w).rad() / abs(w).mid()):.0f}")
                continue
            n += 1
            told += t2 - t1
            tnew += t1 - t0
            same = v.overlaps(w)
            bad += not same
            if same and abs(v).rad() > abs(v).mid() * arb(2) ** -20:
                print(f"id {c['id']} at {h}/{m}: agree, but {prec} bits are too few here for either evaluator (heavy cancellation)")
                continue
            print(f"id {c['id']} at {h}/{m}: {'agree' if same else 'DIFFER'}, |old-new|/|old| = 2^{bits(abs(v - w).mid() / abs(v).mid()):.0f}, "
                  f"old {t2 - t1:.2f} s, new {t1 - t0:.2f} s (x{(t2 - t1) / max(t1 - t0, 1e-9):.0f})", flush=True)
    print(f"{n} comparisons, {bad} differ; old {told:.0f} s, new {tnew:.0f} s in total")
    return bad == 0


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "compare":
        compare(int(sys.argv[2]) if len(sys.argv) > 2 else 30)
    else:
        selftest()
