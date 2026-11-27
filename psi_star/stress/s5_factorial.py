"""S5: Part 5, factorial bounds and Borel summability [v17].

 1. Proposition 2 (shift identity), random N and t, both versions; and past its hypothesis.
 2. Lemma G and the constant c3 (random omega + local optimisation).
 3. Corollary G on the contour, theta in the closed sector, N up to 0.04/|t| (and beyond).
 4. Corollary R: actual |P.J - sum_{n<N}| against C2 e^{-nu N}, C2 from the measured C_U.
 5. Lemma T on the full disc |t| <= 0.3/n, n <= 300 (boundary circle + interior).
 6. Theorem 4: coefficients a_k (independent series computation), then
    |P.J - sum_{k<K} a_k t^k| / (K! (240|t|)^K), K <= 40, many directions incl. both ends.
 7. Theorem 5 at theta0 = -120, -80, with t near the edge of the disc C_R.
Run:  python3 s5_factorial.py     (about 10 min)
"""
import json
import time
import numpy as np
import mpmath as mp
from math import pi, log, cos, sin, factorial
from scipy.optimize import minimize
import core

DEG = pi / 180
rng = np.random.default_rng(5)
C1 = 0.04
C3 = 2 * np.exp(-0.5) * np.sin(0.5)
L = 2 + 1 / C3
NU = log(2) - L * C1


# ------------------------------------------------------------------ the integral P.J
def pj_direct(t, version=1, h=0.02, Rr=None, delta=10.0, th=None):
    """P.J = e^{t/8} int_Gamma H ds  (version 0: e^{t/8 - t/2} int H e^w ds), trapezoid in r,
    s = r e^{-i gamma}.  Returns (value, |integrand| sum for conditioning)."""
    t = complex(t)
    tabs = abs(t)
    if th is None:
        th = np.angle(t) / DEG
        if th > 10 + 1e-9:
            th -= 360          # theta lives in [-210, 10] (sheet of the sector)
    beta, al, ga = core.angles(th, delta)
    if Rr is None:
        Rr = 48.0
    r = np.arange(-Rr, Rr + h / 2, h)
    lh, s, w = core.logH(tabs * r, tabs, th, version=version, cap=60000)
    val = np.exp(lh).sum() * h * np.exp(-1j * ga * DEG)
    pref = np.exp(t / 8) if version == 1 else np.exp(t / 8 - t / 2)
    return pref * val


def tau_n(n, t, version=1):
    t = complex(t)
    j = np.arange(1, n + 1)
    lg = -(n * n + (n if version == 1 else 0)) * t - log(2) - core.log1pexp(-j * t).sum()
    return np.exp(lg)


def logD_N(w, t, N, version=1):
    j = np.arange(N)
    w = np.atleast_1d(w)
    out = -N * N * t - 2 * N * w - core.log1pexp(-w[:, None] - (j[None, :] + 0.5) * t).sum(axis=1)
    if version == 0:
        out = out + N * t
    return out


def shift_rhs(t, N, version=1, h=0.02, Rr=48.0, delta=10.0):
    t = complex(t)
    tabs, th = abs(t), np.angle(t) / DEG
    if th > 10 + 1e-9:
        th -= 360
    beta, al, ga = core.angles(th, delta)
    r = np.arange(-Rr, Rr + h / 2, h)
    lh, s, w = core.logH(tabs * r, tabs, th, version=version, cap=60000)
    integ = np.exp(lh + logD_N(w, t, N, version)).sum() * h * np.exp(-1j * ga * DEG)
    pref = np.exp(t / 8) if version == 1 else np.exp(t / 8 - t / 2)
    res = sum((-1) ** n * tau_n(n, t, version) for n in range(N))
    return res + (-1) ** N * pref * integ, pref * integ


def part1():
    print("1. Proposition 2 (shift identity): |P.J - (sum_{n<N} (-1)^n tau_n + R_N)| / |P.J|")
    worst = (0, None)
    cnt = 0
    for k in range(40):
        th = rng.uniform(-210, 10) if k > 5 else [-210, 10, -90, -150, -30, -90.0][k]
        tabs = float(np.exp(rng.uniform(log(0.008), log(0.35))))
        t = tabs * np.exp(1j * th * DEG)
        beta, al, ga = core.angles(th)
        Nmax = int((pi * cos(al * DEG) / (tabs * cos(ga * DEG))) + 0.5)
        for version in (1, 0):
            pj = pj_direct(t, version)
            for N in sorted(set([1, 2, int(rng.integers(1, max(2, min(Nmax, 60))))])):
                if (N - 0.5) * tabs * cos(ga * DEG) >= pi * cos(al * DEG):
                    continue
                # distance (in s) of the nearest pole of 1/Phi_b from Gamma - iN, and step size
                lam_p = 0.5 + pi * cos(al * DEG) / (tabs * cos(ga * DEG))
                dpole = (lam_p - N) * cos(ga * DEG)
                hh = min(0.02, dpole / 6)
                rhs, rn = shift_rhs(t, N, version, h=hh)
                # conditioning: size of the terms that cancel
                big = sum(abs(tau_n(n, t, version)) for n in range(N)) + abs(rn)
                e = abs(pj - rhs) / max(abs(pj), 1e-16 * big)
                cnt += 1
                if e > worst[0]:
                    worst = (e, (th, tabs, N, version))
    print("   %d tests; worst of |difference| / max(|P.J|, 1e-16 * sum of |terms|) = %.2e"
          " at theta=%.2f |t|=%.4f N=%d version psi_%d^*" % (cnt, worst[0], *worst[1]))
    print("   past the hypothesis (N - 1/2)|t| cos(gamma) >= pi cos(alpha): poles of 1/Phi_b are crossed")
    for th, tabs in [(-90.0, 0.3), (-60.0, 0.3), (0.0, 0.3)]:
        t = tabs * np.exp(1j * th * DEG)
        beta, al, ga = core.angles(th)
        Nc = pi * cos(al * DEG) / (tabs * cos(ga * DEG)) + 0.5
        pj = pj_direct(t, 1)
        line = "   theta=%6.1f |t|=%.2f N_crit=%.2f:" % (th, tabs, Nc)
        for N in [int(Nc) - 1, int(Nc), int(Nc) + 1, int(Nc) + 2]:
            dpole = abs(Nc - N) * cos(ga * DEG)
            rhs, _ = shift_rhs(t, N, 1, h=min(0.02, dpole / 6))
            line += "  N=%d: %.1e" % (N, abs(pj - rhs) / abs(pj))
        print(line)


# ------------------------------------------------------------------ 2. Lemma G
def seg_lattice_dist(om):
    m = np.arange(-30, 30)
    P = 1j * pi * (2 * m + 1)
    y = (P * np.conj(om)).real / abs(om) ** 2
    y = np.clip(y, 0, 1)
    return np.abs(y * om - P).min()


def lemG_ratio(om):
    # log(|g(e^{-om})| / (1/2 e^{L|om|}))
    return (log(2) - 2 * om - core.log1pexp(np.array([-om]))[0]).real - L * abs(om)


def part2():
    print("\n2. Lemma G: |g(e^{-omega})| <= 1/2 e^{L|omega|} when dist([0,omega], i pi(2m+1)) >= 1")
    best = []
    n = 0
    while n < 200000:
        om = complex(rng.uniform(-8, 8), rng.uniform(-14, 14))
        if rng.random() < 0.3:
            om = complex(rng.uniform(-2, 2), rng.uniform(-4, 4))
        n += 1
        if abs(om) < 1e-9 or seg_lattice_dist(om) < 1:
            continue
        best.append((lemG_ratio(om), om))
    best.sort(key=lambda z: -z[0])
    print("   %d admissible samples; max log-ratio %.4f at omega=%s" % (len(best), best[0][0], best[0][1]))

    def f(v):
        om = complex(v[0], v[1])
        if abs(om) < 1e-9 or seg_lattice_dist(om) < 1:
            return 10.0
        return -lemG_ratio(om)
    opt = []
    for r0, om0 in best[:8]:
        res = minimize(f, [om0.real, om0.imag], method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-12})
        opt.append((-res.fun, complex(*res.x)))
    opt.sort(key=lambda z: -z[0])
    print("   after local optimisation: max log-ratio %.4f at omega=%s (claim <= 0)" % opt[0])
    # sharp constant: sup over admissible omega of (log|2g| ) / |omega|
    lam = max((r + L * abs(om)) / abs(om) for r, om in best[:20000] if abs(om) > 0.05)
    print("   sup Re h(omega)/|omega| over samples = %.4f  (L = %.4f)" % (lam, L))
    # c3: inf |1 + e^{omega'}| at distance >= 1 from the lattice
    a = np.linspace(-4, 4, 1601)
    b = np.linspace(-pi, pi, 1257)
    A, B = np.meshgrid(a, b)
    eta = A + 1j * B
    msk = np.abs(eta) >= 1
    val = np.abs(1 - np.exp(eta))[msk]
    k = np.argmin(val)
    print("   inf |1 + e^{omega'}| = %.5f at eta=%s   (c3 = %.5f, claimed true inf 1 - 1/e = %.5f)"
          % (val[k], eta[msk][k], C3, 1 - np.exp(-1)))
    # |h'| <= L on the segment, numerical sup near the contours
    return opt[0]


# ------------------------------------------------------------------ 3. Corollary G
def part3():
    print("\n3. Corollary G: max over the contour of log|D_N(w)| - (-nu N + L c1 |s|)   (claim <= 0)")
    worst = (-1e9, None)
    worst_beyond = {}
    ths = list(np.linspace(-210, 10, 45)) + [-209.9, 9.9, -90.0]
    for tabs in [0.01, 0.005, 0.002]:
        for th in ths:
            beta, al, ga = core.angles(th)
            t = tabs * np.exp(1j * th * DEG)
            rho = np.concatenate([-np.geomspace(1e-4, 80, 400)[::-1], [0], np.geomspace(1e-4, 80, 400)])
            w = rho * np.exp(1j * al * DEG)
            s = w / (1j * t)
            Nmax = int(C1 / tabs)
            for N in range(1, int(0.3 / tabs) + 1, max(1, int(0.3 / tabs) // 60)):
                ld = logD_N(w, t, N).real
                m = ld - (-NU * N + L * C1 * np.abs(s))
                i = int(np.argmax(m))
                if N <= Nmax:
                    if m[i] > worst[0]:
                        worst = (m[i], (tabs, th, N, rho[i]))
                else:
                    c = N * tabs
                    key = round(c, 2)
                    nu_c = log(2) - L * c
                    mm = (ld - (-nu_c * N + L * c * np.abs(s))).max()
                    worst_beyond[key] = max(worst_beyond.get(key, -1e9), mm)
    print("   N <= 0.04/|t|: max = %.4f at |t|=%.3f theta=%.1f N=%d rho=%.3g" % (worst[0], *worst[1]))
    print("   beyond, with c1 replaced by c = N|t| (nu = ln2 - L c):  max over the grid by c:")
    ks = sorted(worst_beyond)
    print("   " + "  ".join("c=%.2f:%.3f" % (k, worst_beyond[k]) for k in ks[::max(1, len(ks) // 12)]))
    # actual exponents: sup_s (log|D_N| ) / |s| as N|t| = c1
    print("   actual constant: sup over contour and theta of (log|D_N| + N ln2) / sum_j |omega_j|, N|t| = 0.04:")
    sup = -1e9
    for th in np.linspace(-210, 10, 45):
        beta, al, ga = core.angles(th)
        tabs = 0.002
        t = tabs * np.exp(1j * th * DEG)
        N = 20
        rho = np.concatenate([-np.geomspace(1e-3, 80, 300)[::-1], np.geomspace(1e-3, 80, 300)])
        w = rho * np.exp(1j * al * DEG)
        om = np.abs(w[:, None] + t * (np.arange(1, N + 1)[None, :] - 0.5)).sum(axis=1)
        v = (logD_N(w, t, N).real + N * log(2)) / om
        sup = max(sup, v.max())
    print("   = %.4f   (Lemma G gives L = %.4f)" % (sup, L))


# ------------------------------------------------------------------ 4. Corollary R
def part4(logCU=-0.505):
    print("\n4. Corollary R: |P.J - sum_{n<N}(-1)^n tau_n| <= C2 e^{-nu N}")
    CU = np.exp(logCU)
    t0 = 0.01
    k0 = 0.185 - t0
    C2 = np.exp(t0 / 2 + C1) * CU * 2 / (k0 - L * C1)
    print("   with C_U = e^{%.3f} (measured, S4) and kappa0 = %.3f: C2 = %.2f;  kappa0 - L c1 = %.4f"
          % (logCU, k0, C2, k0 - L * C1))
    worst = (0, None)
    for th in [-210, -180, -150, -120, -90, -60, -30, 0, 10]:
        for tabs in [0.01, 0.008]:
            t = tabs * np.exp(1j * th * DEG)
            for version in (1, 0):
                pj = pj_direct(t, version)
                for N in range(1, int(C1 / tabs) + 1):
                    res = sum((-1) ** n * tau_n(n, t, version) for n in range(N))
                    rat = abs(pj - res) / (C2 * np.exp(-NU * N))
                    if rat > worst[0]:
                        worst = (rat, (th, tabs, N, version, abs(pj - res)))
    print("   max ratio actual/bound = %.3e at theta=%d |t|=%.3f N=%d psi_%d (|R_N| = %.3e)" % (worst[0], *worst[1]))
    # the integral in the proof, actually
    print("   int |H||D_N||ds| actually (vs the proof's C_U 2/(kappa0 - L c1) e^{-nu N}):")
    for th in [-210, -90, 10]:
        tabs = 0.01
        t = tabs * np.exp(1j * th * DEG)
        beta, al, ga = core.angles(th)
        h = 0.02
        r = np.arange(-48, 48 + h / 2, h)
        lh, s, w = core.logH(tabs * r, tabs, th, cap=60000)
        line = "   theta=%5d:" % th
        for N in [1, 2, 4]:
            I = np.exp(lh.real + logD_N(w, t, N).real).sum() * h
            line += "  N=%d: %.3e (bound %.3e)" % (N, I, CU * 2 / (k0 - L * C1) * np.exp(-NU * N))
        print(line)
    print("   decay well past N = 0.04/|t| (theta=-90, |t|=0.05): |R_N| for N = 1..30 step 5:")
    t = 0.05 * np.exp(-1j * pi / 2)
    pj = pj_direct(t, 1)
    print("   " + "  ".join("%d:%.2e" % (N, abs(pj - sum((-1) ** n * tau_n(n, t) for n in range(N))))
                          for N in range(1, 31, 5)))


# ------------------------------------------------------------------ 5. Lemma T
def part5():
    print("\n5. Lemma T: |tau_n(t)| <= (e^c/2) lambda^n on |t| <= 0.3/n")
    c = 0.3
    mc = 1 + np.exp(-c) * np.cos(c)
    lam = np.exp(c) / mc
    worst = (0, None)
    for version in (1, 0):
        for n in list(range(1, 60)) + list(range(60, 301, 10)):
            phis = np.linspace(0, 2 * pi, 721)[:-1]
            pts = list((c / n) * np.exp(1j * phis))
            pts += list((c / n) * np.sqrt(rng.random(100)) * np.exp(2j * pi * rng.random(100)))
            for t in pts:
                v = abs(tau_n(n, t, version)) / (np.exp(c) / 2 * lam ** n)
                if v > worst[0]:
                    worst = (v, (n, version, t))
    print("   max ratio = %.4f at n=%d psi_%d t=%s (claim <= 1)" % worst[0:1] + worst[1][:2] + (worst[1][2],)
          if False else "   max ratio = %.4f at n=%d version psi_%d, t=%.4g%+.4gi (claim <= 1)"
          % (worst[0], worst[1][0], worst[1][1], worst[1][2].real, worst[1][2].imag))
    # how much room: largest c' such that the same bound holds with radius c'/n (n = 1..60)
    print("   the bound with the stated lambda for radius c'/n, sup over n <= 60 of the ratio:")
    line = "  "
    for cc in [0.3, 0.5, 0.8, 1.0, 1.5]:
        wv = 0
        for n in range(1, 61):
            for t in (cc / n) * np.exp(1j * np.linspace(0, 2 * pi, 361)):
                wv = max(wv, abs(tau_n(n, t)) / (np.exp(c) / 2 * lam ** n))
        line += "  c'=%.1f: %.3g" % (cc, wv)
    print(line)


# ------------------------------------------------------------------ 6. coefficients, Theorem 4
def coeffs(Kmax, version=1, nmax=None, dps=130):
    """a_k = sum_n (-1)^n [t^k] tau_n, by exact power sums and series exponentiation."""
    mp.mp.dps = dps
    if nmax is None:
        nmax = 400 + 14 * Kmax
    # log cosh(v) = sum_m l_m v^m
    lc = mp.taylor(lambda v: mp.log(mp.cosh(v)), 0, Kmax)
    # power sums S_p(n) = sum_{j<=n} j^p, computed incrementally (integers)
    S = [0] * (Kmax + 1)
    a = [mp.mpf(0)] * Kmax
    for n in range(0, nmax + 1):
        if n >= 1:
            for p in range(Kmax + 1):
                S[p] += n ** p
        # log tau_n = -(n+1) log 2 + c1 t - sum_m l_m (t/2)^m S_m(n)
        lin = -mp.mpf(3) / 4 * n * (n + 1) if version == 1 else (-n * n + mp.mpf(n * (n + 1)) / 4)
        ser = [mp.mpf(0)] * Kmax
        if Kmax > 1:
            ser[1] = lin
        for m in range(2, Kmax):
            if lc[m] != 0:
                ser[m] -= lc[m] * S[m] / mp.mpf(2) ** m
        # exp of series with ser[0] = 0
        e = [mp.mpf(0)] * Kmax
        e[0] = mp.mpf(1)
        for k in range(1, Kmax):
            acc = mp.mpf(0)
            for j in range(1, k + 1):
                if ser[j] != 0:
                    acc += j * ser[j] * e[k - j]
            e[k] = acc / k
        sc = (-1) ** n * mp.mpf(2) ** (-(n + 1))
        for k in range(Kmax):
            a[k] += sc * e[k]
    return a


def part6():
    print("\n6. Coefficients and Theorem 4")
    t0 = time.time()
    a1 = coeffs(46, 1)
    a0 = coeffs(46, 0)
    print("   a_0..a_2 (psi_1) = %s, %s, %s   (claim 1/3, 1/9, 5/162)" %
          (mp.nstr(a1[0], 15), mp.nstr(a1[1], 15), mp.nstr(a1[2], 15)))
    print("   a_0..a_2 (psi_0) = %s, %s, %s" % (mp.nstr(a0[0], 15), mp.nstr(a0[1], 15), mp.nstr(a0[2], 15)))
    # compare (afterwards) with the coefficients of v13 (psi_star/v17/phi_exact.json)
    try:
        ex = json.load(open("../v17/phi_exact.json"))
        d1 = max(abs(a1[k] - mp.mpf(mp.fraction(*map(int, ex["psi1"][k].split("/"))))
                     if "/" in ex["psi1"][k] else abs(a1[k] - int(ex["psi1"][k]))) / max(1, abs(a1[k]))
                 for k in range(46))
        key0 = [k for k in ex.keys() if k.startswith("psi0")]
        msg = ""
        if key0:
            d0 = max(abs(a0[k] - mp.mpf(mp.fraction(*map(int, ex[key0[0]][k].split("/"))))
                         if "/" in ex[key0[0]][k] else abs(a0[k] - int(ex[key0[0]][k]))) / max(1, abs(a0[k]))
                     for k in range(46))
            msg = ", psi_0: %.1e" % d0
        print("   comparison with v13's exact coefficients (k < 46): max rel. diff psi_1: %.1e%s" % (d1, msg))
    except Exception as e:
        print("   (comparison with v17/phi_exact.json skipped: %s)" % e)
    # growth: (|a_k|/k!)^{1/k}
    print("   |a_k/k!|^{1/k} at k = 10, 20, 30, 45: " +
          ", ".join("%.4f" % float((abs(a1[k]) / mp.factorial(k)) ** (mp.mpf(1) / k)) for k in [10, 20, 30, 45]) +
          "   (1/|S| = %.4f)" % (1 / abs(pi ** 2 / 3 - 1j * core.CL2)))
    # Lemma T consequence |a_k| <= A0 k! B0^k
    A0, B0 = 3.22, 15.91
    print("   max_k |a_k| / (A0 k! B0^k) = %.3e" % max(float(abs(a1[k]) / (A0 * mp.factorial(k) * B0 ** k)) for k in range(46)))
    print("   (coefficients in %.0f s)" % (time.time() - t0))
    af1 = [complex(v) for v in a1]
    af0 = [complex(v) for v in a0]
    print("   ratio |P.J - sum_{k<K} a_k t^k| / (K! (240|t|)^K): max over K <= 40 (= A needed), and K at the max;")
    print("   B_eff = max over K with remainder above the double floor of (rem/K!)^{1/K}/|t|")
    worstA = (0, None)
    for version, af in ((1, af1), (0, af0)):
        for th in [-210, -209.9, -180, -150, -120, -90, -60, -30, 0, 9.9, 10]:
            line = "   psi_%d theta=%7.1f:" % (version, th)
            for tabs in [0.01, 0.05, 0.3]:
                t = tabs * np.exp(1j * th * DEG)
                pj = pj_direct(t, version)
                partial = 0j
                best, kbest, beff = 0, 0, 0
                for K in range(0, 41):
                    rem = abs(pj - partial)
                    rat = rem / (factorial(K) * (240 * tabs) ** K)
                    if rat > best:
                        best, kbest = rat, K
                    floor = 1e-8 if abs(cos(th * DEG)) < 0.05 else 1e-12   # Faddeev vs products
                    if K >= 1 and rem > 1e3 * floor:
                        beff = max(beff, (rem / factorial(K)) ** (1 / K) / tabs)
                    partial += af[K] * t ** K
                if best > worstA[0]:
                    worstA = (best, (version, th, tabs, kbest))
                line += "  |t|=%.2f: A=%.3f(K=%d) Beff=%.3f" % (tabs, best, kbest, beff)
            print(line)
    print("   largest A needed: %.4f (psi_%d, theta=%.1f, |t|=%.2f, K=%d)" % (worstA[0], *worstA[1]))
    return a1, a0


# ------------------------------------------------------------------ 7. Theorem 5
def part7(a1):
    print("\n7. Theorem 5: P.J = S_theta0[Phi](t) at theta0 = -120, -80, t near the edge of C_R")
    mp.mp.dps = 60
    K = len(a1)
    bc = [a1[k] / mp.factorial(k) for k in range(K)]
    Lp = (K - 1) // 2
    P, Q = mp.pade(bc, Lp, K - 1 - Lp)
    P2, Q2 = mp.pade(bc[:K - 4], (K - 5) // 2, K - 5 - (K - 5) // 2)

    def B(xi, which=1):
        if abs(xi) < 1.5:
            return mp.polyval(bc[::-1], xi)
        if which == 1:
            return mp.polyval(P[::-1], xi) / mp.polyval(Q[::-1], xi)
        return mp.polyval(P2[::-1], xi) / mp.polyval(Q2[::-1], xi)

    def borel(t, th0, which=1):
        e = mp.expj(th0 * mp.pi / 180)
        rate = (e / t).real
        rmax = 80 / rate
        f = lambda r: mp.exp(-r * e / t) * B(r * e, which) * e
        pts = [0, min(1.0, rmax)] + ([rmax] if rmax > 1 else [])
        return complex(mp.quad(f, pts) / t)

    for R in [0.01, 0.1, 0.3]:
        for th0 in [-120.0, -80.0]:
            line = "   R=%.2f theta0=%6.1f:" % (R, th0)
            for phi in [-85, -60, 0, 60, 85]:
                z = 0.97 * R * cos(phi * DEG) * np.exp(1j * phi * DEG)
                t = complex(np.exp(1j * th0 * DEG) * z)
                if not (-210 <= np.angle(t) / DEG <= 10 or -210 <= np.angle(t) / DEG + 360 <= 10):
                    continue
                pj = pj_direct(t, 1)
                bs = borel(mp.mpc(t), th0)
                bs2 = borel(mp.mpc(t), th0, 2) if R > 0.05 else bs
                line += "  phi=%d(|t|=%.3g): %.1e[pade %.0e]" % (phi, abs(t), abs(pj - bs) / abs(pj),
                                                                abs(bs - bs2) / abs(pj))
            print(line)


if __name__ == "__main__":
    T = time.time()
    part1()
    part2()
    part3()
    part4()
    part5()
    a1, a0 = part6()
    part7(a1)
    print("\nrun time %.0f s" % (time.time() - T))
