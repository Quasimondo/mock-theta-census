"""Addendum v15, checks for Theorem 2 (P J ~ Phi on both sides of the real x-axis), contour angle delta = 10 degrees.
t = |t| e^(i theta), x = i t/(2 pi), b = sqrt(x), beta = arg b, w = i t s on the line e^(i alpha) R, alpha = beta - delta.
t1 = t (Re t > 0, |q| < 1) or -t (Re t < 0, |q| > 1);  q1 = e^(-t1), p1 = e^(-2 pi^2/t1);  w' = e^(2 pi s).
   1/Phi_b(b s) = [ (-p1 w'; p1^2)_inf / (-q1^(1/2) e^w; q1)_inf ]^(+1 or -1)        (+1 for |q| < 1, -1 for |q| > 1)
 (a) for |q| > 1:  Phi_b(y) Phi_b(-y) = e^(-t/24 + pi^2/(6t)) e^(pi i y^2)   (inversion relation of the continued function)
 (b) Euler-Maclaurin on the contour, both sides:
       log[ C0 e^(-w^2/(2t)) / Phi_b(b s) ] = G(w)/t - (1/2) ln 2 + t/24 - sum_{k=2}^K B_k(1/2) (-t)^(k-1) Li_(2-k)(-e^w)/k! + O(t^K),
     C0 = 2^(-1/2) e^(pi^2/(12t) + t/24),  G(w) = -w^2/2 + Li2(-e^w) + pi^2/12   (left half of the contour, Re w < 0);
     defect/|t|^K printed for |t| halving.
 (c) uniform bound:  sup over the contour of |H| / exp(psi(rho)/|t|),  H = C0 e^(t s^2/2) / (Phi_b(b s) 2 cosh(pi s)).
 (d) the statement:  (P J - sum_{k<K} a_k t^k)/(a_K t^K) -> 1.
Usage: python3 w2_checks.py [a b c d]"""
import sys
from mpmath import mp, mpf, mpc, exp, log, pi, sqrt, qp, polylog, bernpoly, factorial, cos, sin, cosh, nstr, arg
from borel import PJ, coeff, I
DELTA = mpf(10)

def setup(th, T):
    t = mpf(T)*exp(I*th*pi/180); x = I*t/(2*pi); b = sqrt(x); upper = t.real > 0
    t1 = t if upper else -t
    return t, x, b, upper, t1, exp(-t1), exp(-2*pi**2/t1)

def invPhi_s(s, t, x, upper, t1, q1, p1):          # 1/Phi_b(b s)
    w = I*t*s; wd = exp(2*pi*s)
    r = qp(-p1*wd, p1*p1)/qp(-exp(-t1/2 + w), q1, maxterms=10**6)
    return r if upper else 1/r

def Gfun(w):
    if w.real <= 0: return -w*w/2 + polylog(2, -exp(w)) + pi**2/12
    return -w*w - polylog(2, -exp(-w)) - pi**2/12

def part_a():
    mp.dps = 30; worst = 0
    for th, T in ((-150, "0.5"), (-200, "0.3"), (-110, "0.8")):
        t, x, b, upper, t1, q1, p1 = setup(th, T)
        for y in (mpc("0.3", "0.2"), mpc("-0.7", "0.1")):
            Ph = lambda yy: 1/invPhi_s(yy/b, t, x, upper, t1, q1, p1)
            worst = max(worst, abs(Ph(y)*Ph(-y)/(exp(-t/24 + pi**2/(6*t))*exp(pi*I*y*y)) - 1))
    print("(a) |q| > 1, inversion relation at 6 points: max defect", nstr(worst, 3))

def part_b():
    mp.dps = 40
    print("(b) Euler-Maclaurin on the contour (delta = 10 deg, left half, rho = |w|), defect/|t|^K, for |t| = 0.2, 0.1, 0.05, 0.025")
    for th in (-200, -150, -110, -70, 0):
        for rho in (mpf("0.5"), mpf("3")):
            rows = []
            for T in ("0.2", "0.1", "0.05", "0.025"):
                t, x, b, upper, t1, q1, p1 = setup(th, T)
                w = -rho*exp(I*(arg(b) - DELTA*pi/180)); s = -I*w/t
                N = int(60/abs(t1.real)) + 50
                lq = sum(log(1 + exp(w - t1*(n + mpf(1)/2))) for n in range(N))
                ld = log(qp(-p1*exp(2*pi*s), p1*p1))
                lhs = pi**2/(12*t) - log(2)/2 + t/24 - w*w/(2*t) + (ld - lq if upper else lq - ld)
                out = []
                for K in (1, 3, 5):
                    em = Gfun(w)/t - log(2)/2 + t/24 - sum(bernpoly(k, mpf(1)/2)*(-t)**(k - 1)*polylog(2 - k, -exp(w))/factorial(k) for k in range(2, K + 1))
                    out.append(abs(lhs - em)/abs(t)**K)
                rows.append(out)
            for j, K in enumerate((1, 3, 5)):
                print(f"    theta = {th:5d}, rho = {nstr(rho, 2):>3}, K = {K}:  " + "  ".join(nstr(r[j], 3) for r in rows))
        sys.stdout.flush()

def part_c():
    mp.dps = 20
    print("(c) sup over the contour (delta = 10 deg; rho in (0, 6], both halves) of |H| / exp(psi(rho)/|t|)")
    grid = [mpf(k)/50 for k in range(1, 51)] + [mpf(k)/10 for k in range(11, 61)]
    for th in (-200, -150, -110, -70, -20, 0, 8):
        row = []
        for T in ("0.4", "0.2", "0.1"):
            t, x, b, upper, t1, q1, p1 = setup(th, T); al = arg(b) - DELTA*pi/180; ga = arg(b) + DELTA*pi/180
            C0 = exp(pi**2/(12*t) + t/24)/sqrt(2); best = 0; where = None
            for sgn in (1, -1):
                for rho in grid:
                    w = sgn*rho*exp(I*al); s = -I*w/t
                    H = C0*exp(t*s*s/2)*invPhi_s(s, t, x, upper, t1, q1, p1)/(2*cosh(pi*s))
                    psi = (exp(-I*th*pi/180)*Gfun(w)).real - pi*rho*cos(ga)
                    r = abs(H)/exp(psi/mpf(T))
                    if r > best: best, where = r, sgn*rho
            row.append(f"|t|={T}: {nstr(best, 4)} (at {nstr(where, 3)})")
        print(f"    theta = {th:5d}:  " + "   ".join(row)); sys.stdout.flush()

def part_d(cases=(("psi1", 0), ("psi1", -150), ("psi1", -200), ("psi0", -130))):
    mp.dps = 25
    print("(d) ( P J - sum_{k<K} a_k t^k ) / (a_K t^K),  |t| = 0.2, 0.1, 0.05   (contour angle delta = 10 deg)")
    for case, th in cases:
        vals = {}
        for T in ("0.2", "0.1", "0.05"):
            t = mpf(T)*exp(I*th*pi/180); x = I*t/(2*pi)
            vals[T] = (t, PJ(x, delta=DELTA*pi/180, case=case))
        for K in ((1, 3, 6, 10) if case == "psi1" else (2, 3, 6, 10)):
            out = []
            for T in ("0.2", "0.1", "0.05"):
                t, v = vals[T]
                out.append(nstr((v - sum(coeff(case, k)*t**k for k in range(K)))/(coeff(case, K)*t**K), 8))
            print(f"    {case}, theta = {th:5d}, K = {K:2d}:  " + "    ".join(out))
        sys.stdout.flush()

if __name__ == "__main__":
    for pt in (sys.argv[1:] or ["a", "b", "c", "d"]): {"a": part_a, "b": part_b, "c": part_c, "d": part_d}[pt]()
