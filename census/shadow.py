"""Signature 2: the 'strange' power series at bounded roots.

At a root where a mock theta function stays bounded, g(e^{-t} zeta) ~ a0 + a1 t + a2 t^2 + ...
A modular form contributes only a constant there (plus exponentially small terms), so the
non-constant coefficients (a1, a2, ...) belong to the mock part alone.  Two functions whose
coefficient vectors are proportional across all bounded roots differ by a modular form,
i.e. they share a shadow.  Comparing these vectors tells a genuinely new mock function
from a known one in disguise.
"""
import mpmath as mp
from flint import acb, arb, ctx
from fractions import Fraction

mp.mp.dps = 50
DEG = 7
T0 = 0.00002

def series_at(fn, h, m, prec=512):
    us = [Fraction(x).limit_denominator(10 ** 7) for x in [T0 * (10) ** (i / 15) for i in range(16)]]
    ts, vs = [], []
    for u in us:
        t = u / (m * m)
        p = prec
        while True:
            ctx.prec = p
            try:
                r = (-(arb(t.numerator) / t.denominator)).exp()
                z = (acb(2 * h) / m).exp_pi_i() if m > 1 else acb(1)
                g = fn(r * z)
                if g.real.rad() < arb(2) ** -120 and g.imag.rad() < arb(2) ** -120:
                    break
            except Exception:
                pass
            p *= 2
            if p > 8192:
                return None
        ts.append(mp.mpf(t.numerator) / t.denominator)
        vs.append(mp.mpc(mp.mpf(g.real.mid().str(45, radius=False)), mp.mpf(g.imag.mid().str(45, radius=False))))
    sc = max(ts)
    # power series in t plus possible modular pieces t^(-1/2), t^(-1), t^(-3/2) with zero cusp exponent
    A = mp.matrix([[(t / sc) ** k for k in range(DEG + 1)] + [(t / sc) ** mp.mpf(-0.5), (t / sc) ** -1, (t / sc) ** mp.mpf(-1.5)] for t in ts])
    # complex least squares via real/imag split
    xr, _ = mp.qr_solve(A, mp.matrix([v.real for v in vs]))
    xi, _ = mp.qr_solve(A, mp.matrix([v.imag for v in vs]))
    return [mp.mpc(xr[k], xi[k]) / sc ** k for k in range(DEG + 1)]

def vector(fn, roots):
    """concatenated (a1..a3) over the given bounded roots"""
    out = []
    for h, m in roots:
        a = series_at(fn, h, m)
        if a is None:
            return None
        out += a[1:4]
    return out

def proportional(u, v, tol=1e-12):
    """is u = c v for some complex c?  returns (c, relative misfit)"""
    num = sum(x * mp.conj(y) for x, y in zip(u, v))
    den = sum(abs(y) ** 2 for y in v)
    if den == 0:
        return None, mp.inf
    c = num / den
    mis = mp.sqrt(sum(abs(x - c * y) ** 2 for x, y in zip(u, v)) / max(sum(abs(x) ** 2 for x in u), mp.mpf(10) ** -80))
    return c, mis
