"""Signature 1: growth law along rays q = e^{-t} * zeta toward roots of unity.

For a modular or mock modular form, near the cusp h/m one expects
    log|g(e^{-t} zeta)| = c/t + alpha*log(t) + beta + (exponentially small)
with kappa = c*m^2/(4 pi^2) rational (a cusp exponent) and alpha a half-integer
(minus the weight).  We fit that 3-parameter law in high precision and record
the fit residual and how rational kappa and alpha look.
"""
import sys, json, time, math
from fractions import Fraction
from multiprocessing import Pool
import mpmath as mp
from flint import acb, arb, ctx
import catalog

mp.mp.dps = 60
import os
FAST = os.environ.get('FAST') == '1'
PCAP = int(os.environ.get('PCAP', '65536'))
TFLOOR = float(os.environ.get('TFLOOR', '0.005'))
from math import gcd
ROOTS = [(h, m) for m in [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14] for h in range(0, m // 2 + 1) if gcd(h, m) == 1]
U = [Fraction(15, 1000) * (Fraction(80, 15) ** Fraction(i, 7)) for i in range(8)]  # placeholder, replaced below
# exact rational grid u_i in [0.015, 0.08], roughly geometric
U = [Fraction(x).limit_denominator(100000) for x in [0.004 * (0.02 / 0.004) ** (i / 7) for i in range(8)]]

def evaluate(fn, t, h, m, start=256):
    """start: first precision to try (run() passes what the previous point needed, scaled by 1/t,
    instead of climbing 256, 512, 1024, ... again at every point).
    Afterwards evaluate.lost = bits lost to cancellation at this point."""
    prec = start
    while True:
        ctx.prec = prec
        r = (-(arb(t.numerator) / t.denominator)).exp()
        z = (acb(2 * h) / m).exp_pi_i() if m > 1 else acb(1)
        try:
            g = fn(r * z)
        except catalog.PrecisionLow:
            prec *= 2
            if prec > PCAP:
                return None, prec
            continue
        a = abs(g)
        if a.mid() > 0 and (a.rad() / a.mid()) < arb(2) ** -150:
            y = a.log()
            good = prec if a.rad() == 0 else -float((a.rad() / a.mid()).log().mid()) / math.log(2)
            evaluate.lost = max(0.0, prec - good)
            return mp.mpf(y.mid().str(55, radius=False)), prec
        prec *= 2
        if prec > PCAP:
            return None, prec
evaluate.lost = 0.0

def ratfit(x, D):
    f = Fraction(str(x)).limit_denominator(D)
    return f, abs(mp.mpf(f.numerator) / f.denominator - x)

def probe(args):
    name, h, m = args
    fn = catalog.BY[name]["fn"]
    t0 = time.time()
    hint = []   # [t, bits lost to cancellation] at the last finished point
    def run(ts):
        ys, precs = [], []
        for t in ts:
            # cancellation grows like 1/t: predict the precision this point needs from the last one
            start = 256
            if hint:
                # (about 20 bits of the loss are rounding, not cancellation, and do not scale)
                need = 1.15 * max(0.0, hint[1] - 20) * hint[0] / float(t) + 210
                start = min(PCAP, max(256, 64 * math.ceil(need / 64)))
            try:
                y, p = evaluate(fn, t, h, m, start)
            except Exception as ex:
                return dict(error=f"{ex} at t={float(t):.2e}")
            if y is None:
                return dict(error="precision")
            hint[:] = [float(t), evaluate.lost]
            ys.append(y); precs.append(p)
        T = [mp.mpf(t.numerator) / t.denominator for t in ts]
        A = mp.matrix([[1 / t, mp.log(t), 1, t] for t in T])
        x, _ = mp.qr_solve(A, mp.matrix(ys))
        c, al, be, ga = x[0], x[1], x[2], x[3]
        res = max(abs(ys[i] - (c / T[i] + al * mp.log(T[i]) + be + ga * T[i])) for i in range(len(T)))
        def solve4(idx):
            A4 = mp.matrix([[1 / T[i], mp.log(T[i]), 1, T[i]] for i in idx]); return mp.lu_solve(A4, mp.matrix([ys[i] for i in idx]))[0]
        dc = abs(solve4([0, 1, 2, 3]) - solve4([4, 5, 6, 7]))
        return dict(c=c, al=al, ga=ga, res=res, dc=dc, maxprec=max(precs), tmax=float(ts[-1]))
    fit = run([u / (m * m) for u in U])
    passes = 1
    if "error" not in fit:
        kap0 = fit["c"] * m * m / (4 * mp.pi ** 2)
        if fit["res"] > 1e-25 and abs(kap0) > 1e-4:
            # (also for negative kappa: closer to the rim, the next exponential is suppressed
            # relative to the leading one in the same way)
            # leading exponential is weak: move closer to the rim so that the bounded
            # (non-exponential) part is suppressed by at least e^{-150}
            tmax = Fraction(float(abs(fit["c"]) / 150)).limit_denominator(10 ** 12)
            # census2: for weak growth the ideal second pass lies so close to the rim that a
            # double sum cannot be evaluated in time (cost ~ 1/t^3). Stop at TFLOOR/m^2 and
            # record how strongly the bounded part is still suppressed there.
            floor = Fraction(TFLOOR).limit_denominator(10 ** 6) / (m * m)
            floored = tmax < floor
            if floored:
                tmax = floor
            fit2 = run([tmax * Fraction(x).limit_denominator(10 ** 6) for x in [0.2 * 5 ** (i / 7) for i in range(8)]])
            if "error" not in fit2:
                fit, passes = fit2, 2
                if floored:
                    fit["supp"] = float(abs(fit["c"]) / mp.mpf(tmax.numerator) * tmax.denominator)
        elif fit["res"] > 1e-25 and abs(kap0) <= 1e-4 and not FAST:
            # bounded here: is the leftover exponentially small (modular) or a power series (quantum)?
            fit2 = run([u / (m * m) / 8 for u in U])
            if "error" not in fit2 and fit2["res"] < 1e-25:
                fit, passes = fit2, 2
    if "error" in fit:
        return dict(name=name, h=h, m=m, error=fit["error"])
    c, al, ga, res, dc = fit["c"], fit["al"], fit["ga"], fit["res"], fit["dc"]
    precs = [fit["maxprec"]]
    kappa = c * m * m / (4 * mp.pi ** 2)
    kf, ke = ratfit(kappa, 5000)
    af, ae = ratfit(-al, 4)
    gf, ge = ratfit(ga, 5000)
    return dict(name=name, h=h, m=m, c=float(c), kappa=float(kappa), kappa_rat=str(kf), kappa_err=float(ke),
                gamma=float(ga), gamma_rat=str(gf), gamma_err=float(ge), alpha=float(al), weight_rat=str(af), weight_err=float(ae), resid=float(res),
                dkappa=float(dc * m * m / (4 * mp.pi ** 2)), maxprec=max(precs), secs=round(time.time() - t0, 2), passes=passes, tmax=fit['tmax'],
                supp=fit.get('supp'))

if __name__ == "__main__":
    names = sys.argv[1:] or [e["name"] for e in catalog.E]
    jobs = [(n, h, m) for n in names for (h, m) in ROOTS]
    out = []
    with Pool(2) as p:
        for r in p.imap_unordered(probe, jobs):
            out.append(r)
            print(json.dumps(r), flush=True)
    json.dump(out, open("results.json" if len(sys.argv) == 1 else "results_partial.json", "w"), indent=1)
