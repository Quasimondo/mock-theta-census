"""Addendum v17, checks for Theorem 4 (factorial constants) and Theorem 5 (Borel sum).
 a  the bound |D_N(w)| <= exp(-N (ln 2 - L c1) + L c1 |s|) on the contour, N = floor(c1/|t|), at the ends and in the middle of the sector
 b  |tau_n(t)| <= (e^c/2) lambda^n on the circle |t| = c/n  (c = 0.3)
 c  a_k = sum_n (-1)^n [t^k] tau_n  (exact rational series arithmetic, n <= 400) against phi_exact.json, and |a_k| against the bound A0 k! B0^k
 d  Gevrey ratios ( |P J - sum_{k<K} a_k t^k| / (K! |t|^K) )^(1/K)
 e  P J against the Borel sums in the directions -120, -100, -80 degrees
Usage: python3 g2_checks.py [a b c d e]"""
import sys, json, math
import numpy as np
from fractions import Fraction as Fr

C1, CC = 0.04, 0.3
C3 = 2*math.exp(-0.5)*math.sin(0.5); L = 2 + 1/C3
LAM = math.exp(CC)/(1 + math.exp(-CC)*math.cos(CC))
EX = json.load(open("phi_exact.json"))
def a_float(case, k):
    n, d = EX[case][k].split("/"); return float(Fr(int(n), int(d)))

def part_a():
    print(f"(a) sup over the contour (|w| <= 40) of  log|D_N| + N (ln 2 - L c1) - L c1 |s|   (must be <= 0);  L = {L:.4f}, c1 = {C1}")
    rho = np.linspace(-40, 40, 160001)
    for th in (-210, -150, -90, -30, 10):
        row = []
        for T in (0.02, 0.01, 0.004, 0.001):
            N = int(C1/T); t = T*np.exp(1j*np.radians(th)); beta = np.radians((th + 90)/2); al = beta - np.radians(10)
            w = rho*np.exp(1j*al); tot = np.zeros_like(rho)
            for j in range(1, N + 1):
                om = w + t*(j - 0.5)
                l1p = np.where(om.real < 0, -om.real + np.log(np.abs(1 + np.exp(np.minimum(om.real, 0) + 1j*om.imag))), np.log(np.abs(1 + np.exp(-np.maximum(om.real, 0) - 1j*om.imag))))
                tot += -2*om.real - l1p
            val = tot + N*(math.log(2) - L*C1) - L*C1*np.abs(rho)/T
            i0 = len(rho)//2
            row.append(f"|t|={T}, N={N}: {val.max():+.3f} (log|D_N(0)| + N ln 2 = {tot[i0] + N*math.log(2):+.4f})")
        print(f"    theta = {th:5d}:  " + "   ".join(row))
    print("    sup of |h'| = |-2 + 1/(1 + e^omega)| on the strip around the contour (|eps| <= c1), by direction alpha:")
    r = np.linspace(-30, 30, 60001)
    for al_deg in (-70, -40, -10, 20, 40):
        al = np.radians(al_deg); best = 0
        for off in (-C1, 0, C1):
            om = r*np.exp(1j*al) + 1j*off*np.exp(1j*al)
            best = max(best, np.max(np.abs(-2 + 1/(1 + np.exp(om)))))
        print(f"        alpha = {al_deg:4d}: {best:.4f}   (L = {L:.4f})")

def tau_np(n, t, case):
    q = np.exp(-t); den = 2.0*np.ones_like(t)
    for j in range(1, n + 1): den = den*(1 + q**j)
    return q**(n*n + (n if case == "psi1" else 0))/den

def part_b():
    print(f"(b) max over |t| = c/n of |tau_n(t)| / ((e^c/2) lambda^n),  c = {CC}, lambda = {LAM:.4f}   (must be <= 1)")
    ph = np.exp(2j*np.pi*np.arange(4096)/4096)
    for case in ("psi1", "psi0"):
        row = []
        for n in (1, 2, 3, 5, 10, 20, 40, 80):
            m = np.max(np.abs(tau_np(n, CC/n*ph, case)))
            row.append(f"n={n}: {m/(math.exp(CC)/2*LAM**n):.3g}")
        print(f"    {case}:  " + "  ".join(row))

def mul(a, b, K):
    c = [Fr(0)]*K
    for i, x in enumerate(a):
        if x:
            for j in range(K - i): c[i + j] += x*b[j]
    return c
def inv(a, K):
    b = [Fr(0)]*K; b[0] = 1/a[0]
    for k in range(1, K): b[k] = -b[0]*sum(a[i]*b[k - i] for i in range(1, k + 1))
    return b
def expser(m, K):                      # e^(-m t)
    out = [Fr(1)]
    for k in range(1, K): out.append(out[-1]*Fr(-m, k))
    return out

def part_c(K=13, NMAX=400):
    print(f"(c) a_k from the residues: sum_(n <= {NMAX}) (-1)^n [t^k] tau_n, exact rationals, against phi_exact.json;  and |a_k|/(A0 k! B0^k)")
    A0 = math.exp(CC)/(2*(1 - LAM)); B0 = 1/(CC*(1 - LAM))
    for case in ("psi1", "psi0"):
        prod = [Fr(2)] + [Fr(0)]*(K - 1); tot = [Fr(0)]*K
        for n in range(NMAX + 1):
            if n: prod = mul(prod, [x + (1 if i == 0 else 0) for i, x in enumerate(expser(n, K))], K)
            term = mul(expser(n*n + (n if case == "psi1" else 0), K), inv(prod, K), K)
            for k in range(K): tot[k] += (-1)**n*term[k]
        worst = 0; row = []
        for k in range(K):
            num, den = EX[case][k].split("/"); ex = Fr(int(num), int(den))
            worst = max(worst, abs(float(tot[k] - ex)))
            row.append(f"{abs(float(ex))/(A0*math.factorial(k)*B0**k):.1e}")
        print(f"    {case}: max |difference| for k < {K}: {worst:.1e};   a_0..a_3 = " + ", ".join(str(Fr(int(EX[case][k].split('/')[0]), int(EX[case][k].split('/')[1]))) for k in range(4)))
        print(f"          |a_k|/(A0 k! B0^k), k = 0..{K-1}: " + " ".join(row))
    print(f"    A0 = e^c/(2(1 - lambda)) = {A0:.3f},  B0 = 1/(c (1 - lambda)) = {B0:.2f}")
    for case in ("psi1", "psi0"):
        r = [ (abs(a_float(case, k))/math.factorial(k))**(1/k) for k in (50, 100, 150)]
        print(f"    {case}: (|a_k|/k!)^(1/k) at k = 50, 100, 150: " + ", ".join(f"{v:.4f}" for v in r) + "    (1/|S| = 0.2905)")

def part_d():
    from g1_shift import integral, tx
    print("(d) Gevrey ratios (|P J - sum_{k<K} a_k t^k|/(K! |t|^K))^(1/K);  1/|S| = 0.2905")
    Ks = (1, 2, 3, 4, 6, 8, 10, 12, 14, 16)
    for case in ("psi1", "psi0"):
        for th in (-90, -30, -150):
            for T in (0.2, 0.1):
                t, x = tx(T, th); pj = integral(t, x, 0, case); part = 0j; row = []
                for K in range(0, max(Ks) + 1):
                    if K in Ks:
                        rem = abs(pj - part)
                        row.append(f"{(rem/(math.factorial(K)*T**K))**(1/K):.3f}" if rem > 5e-12 else "  -  ")
                    part += a_float(case, K)*t**K
                print(f"    {case}, theta = {th:5d}, |t| = {T}:  K = {Ks}:  " + " ".join(row))
    from mpmath import mp, mpf, mpc, exp, pi, factorial, radians
    import borel
    mp.dps = 30
    Ks = (1, 2, 4, 8, 12, 16, 20, 24)
    for th in (-200, 5):
        for T in ("0.2", "0.1"):
            t = mpf(T)*exp(mpc(0, 1)*radians(th)); x = mpc(0, 1)*t/(2*pi)
            pj = borel.PJ(x, case="psi1"); part = 0; row = []
            for K in range(0, max(Ks) + 1):
                if K in Ks: row.append(f"{float((abs(pj - part)/(factorial(K)*mpf(T)**K))**(mpf(1)/K)):.3f}")
                part += borel.coeff("psi1", K)*t**K
            print(f"    psi1, theta = {th:5d}, |t| = {T} (product formula, 30 digits):  K = {Ks}:  " + " ".join(row))

def part_e():
    from g1_shift import integral, tx
    from mpmath import mp, mpf, mpc, exp, pi, radians
    import borel
    mp.dps = 30
    print("(e) P J against the Borel sums S_theta0[Phi], theta0 = -120, -100, -80 degrees  (Pade [60/60], Laplace integral)")
    for case in ("psi1", "psi0"):
        t, x = tx(0.3, -90); pj = integral(t, x, 0, case)
        row = [f"{abs(complex(borel.borel_sum(case, mpc(t), radians(th0), 60)) - pj):.1e}" for th0 in (-120, -100, -80)]
        print(f"    {case}, real x (|t| = 0.3, theta = -90; double precision integral):  " + "  ".join(row))
    for th in (-200, -150, -30, 5):
        t = mpf("0.3")*exp(mpc(0, 1)*radians(th)); x = mpc(0, 1)*t/(2*pi)
        pj = borel.PJ(x, case="psi1"); row = []
        for th0 in (-120, -100, -80):
            row.append(f"{float(abs(borel.borel_sum('psi1', t, radians(th0), 60) - pj)):.1e}" if abs(th - th0) < 88 else "(outside)")
        print(f"    psi1, |t| = 0.3, theta = {th:5d} (product formula, 30 digits):  " + "  ".join(row))

if __name__ == "__main__":
    parts = sys.argv[1:] or list("abcde")
    for p in parts: {"a": part_a, "b": part_b, "c": part_c, "d": part_d, "e": part_e}[p]()
