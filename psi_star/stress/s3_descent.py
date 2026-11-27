"""S3: Part 3, Proposition 1 [v15 §1.2].

Tests, with independent code (scipy's complex dilogarithm, Gauss-Legendre for ell):
  (i)   the two exact formulas for psi/rho, against psi computed from G directly;
  (ii)  |Y| <= pi/2 and lambda_alpha <= X <= ln 2;
  (iii) psi(rho) <= -kappa rho - c rho^2 sin 2delta with the stated kappa;
  (iv)  kappa >= 0.185 for delta = 10 deg, -210 <= theta <= 10, and psi <= -0.185 rho
        for rho in [1e-3, 100]; then beyond the sector.
Run:  python3 s3_descent.py      (about 1 min)
"""
import numpy as np
from math import pi, log
import core

DEG = pi / 180
rng = np.random.default_rng(3)
gl_x, gl_w = np.polynomial.legendre.leggauss(600)
lam_nodes = 0.5 * (gl_x + 1)
lam_wts = 0.5 * gl_w


def ell(w, side):
    """ell = int_0^1 log(1 + e^{-+lambda w}) d lambda (side=+1 right half, -1 left half)."""
    w = np.atleast_1d(w)
    z = np.exp(-side * np.outer(w, lam_nodes))
    return (np.log1p(z) * lam_wts).sum(axis=1)


def lam_alpha(al_deg):
    a = abs(al_deg) * DEG
    if a == 0:
        return 0.0
    if a >= pi / 2:
        return float('nan')
    return log(1 - np.exp(-(pi / 2) / np.tan(a)))


def kappa_formula(theta, delta=10.0):
    beta, al, ga = core.angles(theta, delta)
    return (pi / 2) * np.cos(ga * DEG) - max(log(2), -lam_alpha(al)) * abs(np.sin(ga * DEG))


def theta_grid(lo=-210.0, hi=10.0):
    g = list(np.linspace(lo, hi, 441))
    g += list(lo + np.geomspace(1e-6, 2, 30)) + list(hi - np.geomspace(1e-6, 2, 30))
    g += list(-90 + np.linspace(-1, 1, 41))
    g += list(rng.uniform(lo, hi, 200))
    return np.array(sorted(set(g)))


def rho_grid():
    r = np.geomspace(1e-3, 100, 1500)
    return np.concatenate([r, rng.uniform(1e-3, 100, 500)])


def main():
    delta = 10.0
    sin2d = np.sin(2 * delta * DEG)
    TH = theta_grid()
    RH = rho_grid()
    print("S3  Proposition 1, delta = 10 deg; %d theta values x %d rho values per half" % (len(TH), len(RH)))
    worst_i = 0.0
    maxY, maxX, min_X_minus_lam = 0.0, -1e9, 1e9
    worst_iii = (1e9, None)
    worst_iv = (1e9, None)
    min_kappa = (1e9, None)
    nsamp = 0
    for th in TH:
        beta, al, ga = core.angles(th, delta)
        kap = kappa_formula(th, delta)
        if kap < min_kappa[0]:
            min_kappa = (kap, th)
        for side in (+1, -1):
            r = side * RH
            w = r * np.exp(1j * al * DEG)
            psi = core.psi_of(r, th, delta)
            # (i) exact formula
            l = ell(w, side)
            X, Y = l.real, l.imag
            g, cg, sg = ga * DEG, np.cos(ga * DEG), np.sin(ga * DEG)
            if side > 0:
                f = -(pi - Y) * cg - X * sg - RH * sin2d
                c = 1.0
            else:
                f = -(pi + Y) * cg + X * sg - (RH / 2) * sin2d
                c = 0.5
            err = np.abs(psi / RH - f)
            worst_i = max(worst_i, err.max())
            # (ii)
            maxY = max(maxY, np.abs(Y).max())
            maxX = max(maxX, X.max())
            min_X_minus_lam = min(min_X_minus_lam, (X - lam_alpha(al)).min())
            # (iii) margin of psi <= -kappa rho - c rho^2 sin 2 delta, per unit rho
            m3 = (-kap * RH - c * RH ** 2 * sin2d - psi) / RH
            i3 = np.argmin(m3)
            if m3[i3] < worst_iii[0]:
                worst_iii = (m3[i3], (th, side * RH[i3]))
            # (iv) psi <= -0.185 rho
            m4 = (-0.185 * RH - psi) / RH
            i4 = np.argmin(m4)
            if m4[i4] < worst_iv[0]:
                worst_iv = (m4[i4], (th, side * RH[i4]))
            nsamp += len(RH)
    print("(i)   exact formula: max |psi/rho - formula| = %.2e  (%d samples)" % (worst_i, nsamp))
    print("(ii)  max |Y| = %.6f  (pi/2 = %.6f; ratio %.4f)" % (maxY, pi / 2, maxY / (pi / 2)))
    print("      max X = %.6f  (ln 2 = %.6f; ratio %.4f);  min (X - lambda_alpha) = %.4e"
          % (maxX, log(2), maxX / log(2), min_X_minus_lam))
    print("(iii) min over samples of [(-kappa rho - c rho^2 sin2d) - psi]/rho = %.4e at theta=%.4f, r=%.4g"
          % (worst_iii[0], *worst_iii[1]))
    print("(iv)  min kappa(theta) on [-210,10] = %.6f at theta = %.4f  (claim >= 0.185)" % min_kappa)
    print("      min over samples of [(-0.185 rho) - psi]/rho = %.4e at theta=%.4f, r=%.4g"
          % (worst_iv[0], *worst_iv[1]))
    # true descent margin kappa_true(theta) = inf_rho (-psi/rho) over rho in [1e-3, 100]
    print("\nTrue margin  kappa_true(theta) = min_rho(-psi/rho), rho in [1e-3,100], and the formula kappa(theta):")
    for th in [-260, -250, -240, -230, -220, -215, -212, -210, -200, -150, -120, -90, -60, -30, 0, 5, 10, 12,
               15, 20, 25, 30, 35, 40, 50]:
        beta, al, ga = core.angles(th, delta)
        adm = abs(beta) < 80.0
        kt = 1e9
        for side in (+1, -1):
            psi = core.psi_of(side * RH, th, delta)
            kt = min(kt, (-psi / RH).min())
        print("  theta=%7.1f  beta=%6.1f  gamma=%6.1f  alpha=%6.1f  kappa=%8.4f  kappa_true=%8.4f %s"
              % (th, beta, ga, al, kappa_formula(th, delta), kt, "" if adm else "(delta not admissible)"))
    # where does the formula kappa cross 0.185 and 0 ?
    th_f = np.linspace(-300, 80, 38001)
    kf = np.array([kappa_formula(th) for th in th_f])
    ok = th_f[kf >= 0.185]
    pos = th_f[kf > 0]
    print("\nformula kappa >= 0.185 exactly on theta in [%.3f, %.3f];  kappa > 0 on [%.3f, %.3f]"
          % (ok.min(), ok.max(), pos.min(), pos.max()))
    # analyticity (Corollary of v16) needs a contour that is fixed while t moves: with gamma fixed,
    # delta = gamma - beta changes with theta.  Descent margin for delta in [5, 15] deg:
    print("\nformula kappa(theta, delta) for delta != 10 deg (a fixed contour gamma = gamma(theta0) "
          "has delta = 10 + (theta0 - theta)/2):")
    for d in [5.0, 8.0, 9.0, 11.0, 12.0, 15.0]:
        th_g = np.linspace(-210, 10, 2201)
        kv = np.array([kappa_formula(th, d) for th in th_g])
        print("   delta=%4.1f: min over [-210,10] = %.4f at theta=%.1f" % (d, np.nanmin(kv), th_g[np.nanargmin(kv)]))


if __name__ == "__main__":
    main()
