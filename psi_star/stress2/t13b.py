"""T13b: Theorem 13 (b) of v22, the tilted sector 0 <= arg x <= eps(a).

 (a) on C' (tilt phi1): (L1) at tau = 0, and the bound Im V >= -(3 pi^2/2) tan(phi1) on C2';
 (b) eps(a): for each a and phi1 the largest tau with (L1) on C' (bisection), eps = min(phi1, atan tau);
     best over phi1; comparison with the table of v22 §6 (phi1 = 5 deg at arg x = 3, 5 deg);
 (c) the bound on C1': |p w'| <= 1 for arg x <= phi1, directly;
 (d) Theorem 13 (b) itself at arg x = eps(a)/2 and eps(a): (exact/asymptotic - 1)/|t| over C'.
Run:  python3 t13b.py   (about 8 min)
"""
import time
import numpy as np
from math import pi, atan, tan, log
import common2 as c
import core

DEG = pi / 180
PHIS = [0.5, 1, 2, 3, 5, 8, 12, 20, 30, 45]


def sup_L1(a, tau, phi1, n=1500):
    res, w0, v0 = c.heights_L1(a, tau, n=n, phi1=phi1)
    return max(q for q, wb in res.values()), res


def tau_star(a, phi1):
    if sup_L1(a, 0.0, phi1)[0] >= 0:
        return -1.0
    lo, hi = 0.0, 5.0
    if sup_L1(a, hi, phi1)[0] < 0:
        return hi
    for _ in range(30):
        mid = (lo + hi) / 2
        if sup_L1(a, mid, phi1)[0] < 0:
            lo = mid
        else:
            hi = mid
    return lo


def a_list():
    return [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999]


def part_a():
    print("(a) C' at tau = 0: sup of (H - H(w0))/min(d^2, d), and min Im V on C2' + (3 pi^2/2) tan phi1")
    for phi1 in [1, 5, 20]:
        line = "   phi1=%2d:" % phi1
        worst = -1e9
        wb = 1e9
        for a in a_list():
            s, res = sup_L1(a, 0.0, phi1)
            worst = max(worst, s / c.vol(a))
            P = -pi * tan(phi1 * DEG) + 2j * pi
            w = P + np.linspace(0, 1, 4001) * (1j * pi - P)
            iv = c.V(w, a - 1).imag
            wb = min(wb, iv.min() + 1.5 * pi ** 2 * tan(phi1 * DEG))
        print(line + " max over a of sup/Vol = %.4f (claim < 0);  min (Im V + 3pi^2/2 tan phi1) on C2' = %.4f (claim >= 0)"
              % (worst, wb))


def part_b():
    print("\n(b) eps(a): largest tau with (L1) on C' (phi1 = tilt), eps = min(phi1, atan tau) in degrees")
    print("   %7s %8s " % ("a", "Vol") + " ".join("%7s" % ("p=%g" % p) for p in PHIS) + "   best eps  (doc phi1=5: 3 deg / 5 deg)")
    out = {}
    for a in a_list():
        row = []
        for phi1 in PHIS:
            ts = tau_star(a, phi1)
            e = -1.0 if ts < 0 else min(phi1, atan(ts) / DEG)
            row.append(e)
        out[a] = row
        e5 = row[PHIS.index(5)]
        print("   %7.3f %8.5f " % (a, c.vol(a)) + " ".join("%7.2f" % e for e in row)
              + "   %7.2f   %s / %s" % (max(row), "holds" if e5 >= 3 else "fails", "holds" if e5 >= 5 else "fails"))
    return out


def part_c():
    print("\n(c) C1': |p w'| = exp((-|r| cos phi + pi sin phi)/|x|) <= 1 for 0 <= arg x = phi <= phi1 (direct evaluation)")
    worst = -1e9
    for phi1 in PHIS:
        P = -pi * tan(phi1 * DEG)
        for phi in np.linspace(0, phi1, 11):
            for xa in [0.25, 0.05, 0.01]:
                x = xa * np.exp(1j * phi * DEG)
                r = P - np.geomspace(1e-9, 50, 200)
                w = r + 2j * pi
                lg = (-1j * pi / x + w / x).real       # log |p w'| with p = e^{-pi i/x}, w' = e^{w/x}
                worst = max(worst, lg.max() * xa)
    print("   max of |x| log|p w'| over C1', phi <= phi1, |x| in {1/4, 1/20, 1/100}: %.2e (claim <= 0)" % worst)


def part_d(eps):
    print("\n(d) Theorem 13 (b): (exact/asymptotic - 1)/|t| over C' at arg x = eps(a)/2 and eps(a)")
    for a in [0.1, 0.3, 0.5, 0.7, 0.9]:
        row = eps[a]
        k = int(np.argmax(row))
        phi1, e = PHIS[k], row[k]
        line = "   a=%.1f (phi1=%g, eps=%.2f):" % (a, phi1, e)
        for argx in [e / 2, e]:
            vals = []
            for tabs in [0.2, 0.1, 0.05]:
                x = tabs / (2 * pi) * np.exp(1j * argx * DEG)
                try:
                    I, t, b, v0 = c.I_grow_C(a, x=x, phi1=phi1)
                    vals.append("%.3f" % (abs(I / c.asym13(a, t, b) - 1) / tabs))
                except RuntimeError:
                    vals.append("*")
            line += "  arg x=%.2f: %s" % (argx, " ".join(vals))
        print(line)


if __name__ == "__main__":
    T = time.time()
    part_a()
    eps = part_b()
    part_c()
    part_d(eps)
    print("\nrun time %.0f s" % (time.time() - T))
