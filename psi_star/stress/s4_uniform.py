"""S4: Part 4, the uniform bound for the kernel integrand.

(a) Lemma 1' [v15 §2.2] (Euler-Maclaurin with complex step): the remainder
    rho_K of  sum_{n>=0} log(1 + e^{zeta0 - t1 (n+a)})  against the stated expansion,
    for random zeta0, t1 satisfying the hypotheses (distance >= eta from i pi (2m+1)
    along the ray, Re zeta <= eta/2); ratio  |rho_K| cos(arg t1) / |t1|^{K-1}.
(b) Lemma A: distance of the actual paths w - t1 y (left half), -w - t1 y (right half)
    from the lattice, divided by eta = pi min(cos alpha, |cos theta|); and the
    |w| <= eta/4 clause with 3 eta/4.
(c) Lemma 3' and [v16 §4]: sup over the contour of log|H| - psi(rho)/|t| and of
    log|H| + kappa |s| (kappa = 0.185), psi_1^* and psi_0^* integrands, for theta in the
    closed sector [-210, 10] and |t| from 0.4 down to 0.002; then beyond the sector.
Run:  python3 s4_uniform.py     (about 6 min)
"""
import numpy as np
import mpmath as mp
from math import pi, log, cos
import core

DEG = pi / 180
rng = np.random.default_rng(4)


# ------------------------------------------------------------------ (a)
def ray_lattice_dist(z0, t1, mrange=60):
    """min over m of the distance from i pi (2m+1) to the ray z0 - t1 y, y >= 0."""
    m = np.arange(-mrange, mrange)
    P = 1j * pi * (2 * m + 1)
    d = -t1
    y = ((P - z0) * np.conj(d)).real / abs(d) ** 2
    y = np.maximum(y, 0)
    return np.abs(z0 + y * d - P).min()


def em_expansion(z0, t1, a, K):
    u = mp.exp(mp.mpc(z0))
    s = -mp.polylog(2, -u) / t1
    for k in range(1, K + 1):
        s += mp.bernpoly(k, a) * (-t1) ** (k - 1) * mp.polylog(2 - k, -u) / mp.factorial(k)
    return complex(s)


def em_sum(z0, t1, a):
    n_terms = int((max(z0.real, 0) + 40) / t1.real) + 5
    n = np.arange(n_terms) + a
    v = core.log1pexp(z0 - t1 * n)
    import math
    return complex(math.fsum(v.real), math.fsum(v.imag))


def part_a():
    print("(a) Lemma 1' remainder ratio R = |rho_K| cos(arg t1) / |t1|^{K-1}")
    mp.mp.dps = 30
    Ks = [1, 2, 3, 4]
    rows = []
    tried = 0
    while len(rows) < 500 and tried < 200000:
        tried += 1
        eta = rng.choice([0.25, 0.5, 1.0, 2.0])
        T = float(np.exp(rng.uniform(np.log(0.004), np.log(0.5))))
        u = rng.random()
        if u < 0.35:
            phi = rng.choice([-1, 1]) * (90 - 10 ** rng.uniform(-1.3, 1))   # near the edges
        else:
            phi = rng.uniform(-80, 80)
        t1 = T * np.exp(1j * phi * DEG)
        if t1.real < 2.5e-4:
            continue
        z0 = complex(rng.uniform(-5, eta / 2), rng.uniform(-3 * pi, 3 * pi))
        if ray_lattice_dist(z0, t1) < eta:
            continue
        a = rng.choice([0.5, 1.0, float(rng.uniform(0.05, 1))])
        S = em_sum(z0, t1, a)
        R = []
        for K in Ks:
            rem = S - em_expansion(z0, t1, a, K)
            R.append(abs(rem) * cos(phi * DEG) / T ** (K - 1))
        rows.append((eta, T, phi, z0, a, R))
    print("   %d admissible samples (of %d drawn)" % (len(rows), tried))
    for K_i, K in enumerate(Ks):
        for eta in [0.25, 0.5, 1.0, 2.0]:
            sel = [r for r in rows if r[0] == eta]
            if not sel:
                continue
            Rv = np.array([r[5][K_i] for r in sel])
            j = int(np.argmax(Rv))
            small = [r[5][K_i] for r in sel if r[1] < 0.02]
            edge = [r[5][K_i] for r in sel if abs(r[2]) > 85]
            print("   K=%d eta=%.2f  n=%3d  max R=%.3e (|t1|=%.3g, arg=%.2f, zeta0=%.2f%+.2fi)  "
                  "max R(|t1|<0.02)=%.3e  max R(|arg|>85)=%.3e"
                  % (K, eta, len(sel), Rv[j], sel[j][1], sel[j][2], sel[j][3].real, sel[j][3].imag,
                     max(small) if small else float('nan'), max(edge) if edge else float('nan')))
    return rows


# ------------------------------------------------------------------ (b)
def part_b():
    print("\n(b) Lemma A: min distance of the paths to i pi(2m+1), divided by eta")
    worst = (1e9, None)
    worst_small = (1e9, None)
    for th in np.concatenate([np.linspace(-210, -150, 61), np.linspace(-30, 10, 41)]):
        beta, al, ga = core.angles(th)
        t = np.exp(1j * th * DEG)              # only the direction matters
        t1 = t if t.real > 0 else -t
        eta = pi * min(cos(al * DEG), abs(cos(th * DEG)))
        for rho in np.concatenate([[0.0], np.geomspace(1e-3, 60, 300)]):
            for side in (+1, -1):
                w = side * rho * np.exp(1j * al * DEG)
                z0 = w if side < 0 else -w
                d = ray_lattice_dist(z0, t1)
                if d / eta < worst[0]:
                    worst = (d / eta, (th, side * rho))
        for rr in rng.uniform(0, eta / 4, 300):
            z0 = rr * np.exp(1j * rng.uniform(0, 2 * pi))
            d = ray_lattice_dist(z0, t1)
            if d / (0.75 * eta) < worst_small[0]:
                worst_small = (d / (0.75 * eta), (th, z0))
    print("   min dist/eta on the contour paths = %.4f at theta=%.1f, signed rho=%.4g  (claim >= 1)"
          % (worst[0], *worst[1]))
    print("   min dist/(3eta/4) for |w| <= eta/4 = %.4f at theta=%.1f, w=%s  (claim >= 1)"
          % (worst_small[0], worst_small[1][0], worst_small[1][1]))


# ------------------------------------------------------------------ (c)
def sup_on_contour(tabs, th, version=1, rho_max=40.0, npts=500):
    rho = np.concatenate([-np.geomspace(1e-5, rho_max, npts)[::-1], [0.0], np.geomspace(1e-5, rho_max, npts)])
    rho = np.concatenate([rho, rng.uniform(-rho_max, rho_max, 100)])
    lh, s, w = core.logH(rho, tabs, th, version=version, cap=30000)
    psi = core.psi_of(rho, th)
    d1 = lh.real - psi / tabs
    if version == 0:
        d1 = d1 - np.abs(rho)          # the proof bounds |e^w| by e^{|t||s|} = e^{rho}
    d2 = lh.real + 0.185 * np.abs(s)
    i1, i2 = int(np.argmax(d1)), int(np.argmax(d2))
    return d1[i1], rho[i1], d2[i2], rho[i2]


def part_c():
    print("\n(c) sup over the contour of log|H| - psi/|t| (=: log C_psi; for psi_0^* minus rho = |t||s|)"
          " and of log|H| + 0.185|s| (=: log C_U)")
    thetas = sorted(set(list(np.linspace(-210, 10, 23)) + [-209.9, -209, -205, -151, -150, -149, -91, -90.1,
                                                           -90.001, -90, -89.999, -89.9, -89, -31, -30, -29,
                                                           5, 9, 9.9]))
    ts = [0.4, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
    for version in (1, 0):
        print("  version psi_%d^*" % version)
        print("  %8s %10s %9s %10s %9s" % ("|t|", "logC_psi", "at theta", "logC_U", "at theta"))
        for tabs in ts:
            best1, best2 = (-1e9, None), (-1e9, None)
            for th in thetas:
                d1, r1, d2, r2 = sup_on_contour(tabs, th, version)
                if d1 > best1[0]:
                    best1 = (d1, (th, r1))
                if d2 > best2[0]:
                    best2 = (d2, (th, r2))
            print("  %8.3f %10.4f %9.3f %10.4f %9.3f   (rho at max: %.3g, %.3g)"
                  % (tabs, best1[0], best1[1][0], best2[0], best2[1][0], best1[1][1], best2[1][1]))
    print("\n  very small |t| at the worst directions (psi_1^*):")
    for th in [10.0, -210.0, -90.0, -150.0, -30.0]:
        line = "   theta=%7.1f:" % th
        for tabs in [0.01, 0.003, 0.001]:
            d1, r1, d2, r2 = sup_on_contour(tabs, th, 1, npts=300)
            line += "  |t|=%.3f: %.4f / %.4f" % (tabs, d1, d2)
        print(line)
    print("\n  beyond the sector (psi_1^*, |t| = 0.05 and 0.01): log C_psi / log C_U")
    for th in [-240, -230, -225, -220, -215, -212, 12, 15, 20, 25, 30, 35]:
        beta = (th + 90) / 2
        if abs(beta) >= 80:
            continue
        line = "   theta=%6.1f:" % th
        for tabs in [0.05, 0.01]:
            d1, r1, d2, r2 = sup_on_contour(tabs, th, 1, npts=300)
            line += "   |t|=%.2f: %8.4f / %8.4f" % (tabs, d1, d2)
        print(line)


if __name__ == "__main__":
    import time
    t0 = time.time()
    part_a()
    part_b()
    part_c()
    print("\nrun time %.0f s" % (time.time() - t0))
