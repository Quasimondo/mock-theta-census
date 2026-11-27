"""Addendum v14, checks for the proof of  P J ~ Phi  on  -90 deg < arg t < -arctan(ln 2/pi) = -12.44 deg.
Contour: tau = i t s real (s = y/b), i.e. the line L of v12 with delta = arg b.   q = e^(-t), p = e^(-2 pi^2/t), theta = arg t.
 (a) the real potential U(tau) = -tau^2/2 + Li2(-e^tau) + pi^2/12:  U(0) = 0, U'(0) = -ln 2, U concave, reflection identity,
     maximum pi^2/60 at tau = ln((sqrt5 - 1)/2);  and Lemma 2:
       e^(-tau^2/(2t)) (-q;q)_inf (-p w;p^2)_inf / ((p;p^2)_inf (-q^(1/2) e^tau;q)_inf)
          = 2^(-1/2) e^(-pi^2/(12t) + t/12) e^(-tau^2/t) (-q^(1/2) e^(-tau);q)_inf / (-p/w;p^2)_inf .
 (b) Euler-Maclaurin with remainder:
       log (-q^(1/2) e^tau; q)_inf = -Li2(-e^tau)/t + sum_{k=1}^K B_k(1/2) (-t)^(k-1) Li_(2-k)(-e^tau)/k! + O(t^K)
       log (-q; q)_inf             =  pi^2/(12 t)   + sum_{k=1}^K B_k(1)   (-t)^(k-1) Li_(2-k)(-1)/k!     + O(t^K)
     the defect divided by |t|^K is printed for |t| halving; it must stay bounded.
 (c) the uniform bound:  |H(tau)| <= C exp( (cos(theta) U(tau) - pi |sin(theta)| |tau|)/|t| )  for all real tau, where
     H = e^(-tau^2/(2t)) (-q;q)_inf (-p w;p^2)_inf / ( (p;p^2)_inf (-q^(1/2) e^tau;q)_inf 2 cos(pi tau/t) ),  w = e^(-2 pi i tau/t);
     sup over a grid of tau of |H| / bound is printed for several |t|; it must stay bounded as |t| -> 0.
 (d) the statement itself:  ( P J - sum_{k<K} a_k t^k ) / (a_K t^K) -> 1  as |t| -> 0  (a_k exact, from v13), psi1 series;
 (e) the same for the psi0 series (P q^(1/2) J(b + 1/(2b))).
Usage: python3 u1_checks.py [part ...]     parts: a b c d e (default all)"""
import sys
from mpmath import mp, mpf, mpc, exp, log, pi, sqrt, qp, polylog, bernpoly, factorial, cos, sin, nstr, diff
from borel import PJ, coeff, I

def U(tau):
    return -tau**2/2 + polylog(2, -exp(tau)).real + pi**2/12

def part_a():
    mp.dps = 30
    print("(a) U(0) =", nstr(U(mpf(0)), 5), "  U'(0) + ln 2 =", nstr(diff(U, mpf(0)) + log(2), 5))
    ts = log((sqrt(5) - 1)/2)
    print("    tau* = ln((sqrt5-1)/2) =", nstr(ts, 12), "  U'(tau*) =", nstr(diff(U, ts), 5), "  U(tau*) - pi^2/60 =", nstr(U(ts) - pi**2/60, 5))
    worst = max(abs(U(mpf(k)/7) - (-(mpf(k)/7)**2 - polylog(2, -exp(-mpf(k)/7)).real - pi**2/12)) for k in range(-60, 61))
    print("    reflection  U(tau) = -tau^2 - Li2(-e^(-tau)) - pi^2/12 :  max defect on [-8.6, 8.6] =", nstr(worst, 3))
    # Lemma 2: the inversion form of the integrand
    worst2 = 0
    for th, T, tau in ((-70, "0.3", "0.8"), (-30, "0.2", "2.5"), (-50, "0.5", "-0.4")):
        t = mpf(T)*exp(I*th*pi/180); q = exp(-t); p = exp(-2*pi**2/t); tau = mpf(tau); w = exp(-2*pi*I*tau/t)
        left = exp(-tau**2/(2*t))*qp(-q, q)/qp(-exp(-t/2 + tau), q, maxterms=10**6)*qp(-p*w, p*p)/qp(p, p*p)
        right = exp(-pi**2/(12*t) + t/12)/sqrt(2)*exp(-tau**2/t)*qp(-exp(-t/2 - tau), q, maxterms=10**6)/qp(-p/w, p*p)
        worst2 = max(worst2, abs(left/right - 1))
    print("    inversion form (Lemma 2), three points: max |left/right - 1| =", nstr(worst2, 3))
    conc = max(diff(U, mpf(k)/7, 2) for k in range(-60, 61))
    slope = max(U(mpf(k)/7) + log(2)*mpf(k)/7 for k in range(-60, 61))
    print("    max U'' on the grid =", nstr(conc, 5), " (U'' = -1 - e^tau/(1+e^tau) < 0);   max of U(tau) + ln2 tau =", nstr(slope, 3), " (<= 0 by concavity)")

def part_b():
    mp.dps = 40
    print("(b) Euler-Maclaurin, defect/|t|^K  (K = 1, 3, 5), for |t| = 0.2, 0.1, 0.05, 0.025")
    for th in (-75, -45, -20):
        for tau in (mpf("-0.6"), mpf("0.25")):
            rows = []
            for T in ("0.2", "0.1", "0.05", "0.025"):
                t = mpf(T)*exp(I*th*pi/180); q = exp(-t)
                # the logs are of products of factors 1 + u, u on a spiral; sum the individual logs to stay on the right branch
                lhs = sum(log(1 + exp(tau - t*(n + mpf(1)/2))) for n in range(int(60/abs(t.real)) + 50))
                lhs0 = sum(log(1 + exp(-t*(n + 1))) for n in range(int(60/abs(t.real)) + 50))
                out = []
                for K in (1, 3, 5):
                    em = -polylog(2, -exp(tau))/t + sum(bernpoly(k, mpf(1)/2)*(-t)**(k - 1)*polylog(2 - k, -exp(tau))/factorial(k) for k in range(1, K + 1))
                    em0 = pi**2/(12*t) + sum(bernpoly(k, 1)*(-t)**(k - 1)*polylog(2 - k, mpf(-1))/factorial(k) for k in range(1, K + 1))
                    out.append((abs(lhs - em)/abs(t)**K, abs(lhs0 - em0)/abs(t)**K))
                rows.append(out)
            for j, K in enumerate((1, 3, 5)):
                print(f"    theta = {th:4d}, tau = {nstr(tau, 3):>5}, K = {K}:  tau-product: " + "  ".join(nstr(r[j][0], 3) for r in rows) + "     (-q;q): " + "  ".join(nstr(r[j][1], 3) for r in rows))
        sys.stdout.flush()

def part_c():
    mp.dps = 25
    print("(c) sup_tau |H(tau)| / exp((cos(theta) U - pi |sin(theta)| |tau|)/|t|)   (grid tau in [-6, 6], step 0.02 near 0, 0.1 outside)")
    grid = [mpf(k)/50 for k in range(-50, 51)] + [mpf(k)/10 for k in range(-60, -10)] + [mpf(k)/10 for k in range(11, 61)]
    for th in (-80, -60, -40, -20, -14):
        row = []
        for T in ("0.4", "0.2", "0.1", "0.05"):
            t = mpf(T)*exp(I*th*pi/180); q = exp(-t); p = exp(-2*pi**2/t); c, s_ = cos(th*pi/180), abs(sin(th*pi/180))
            pref = exp(pi**2/(12*t) + t/24)/sqrt(2)          # (-q;q)_inf/(p;p^2)_inf
            best = 0; where = None
            for tau in grid:
                w = exp(-2*pi*I*tau/t)
                invPhi = qp(-p*w, p*p)/qp(-exp(-t/2 + tau), q, maxterms=10**6)
                H = exp(-tau**2/(2*t))*pref*invPhi/(2*cos(pi*tau/t))
                bd = exp((c*U(tau) - pi*s_*abs(tau))/mpf(T))
                r = abs(H)/bd
                if r > best: best, where = r, tau
            row.append(f"|t|={T}: {nstr(best, 4)} (at tau={nstr(where, 3)})")
        print(f"    theta = {th:4d}:  " + "   ".join(row)); sys.stdout.flush()

def part_d(cases=(("psi1", -60), ("psi1", -30), ("psi1", -80)), label="(d)"):
    mp.dps = 25
    print(label + " ( P J - sum_{k<K} a_k t^k ) / (a_K t^K)   (contour: tau real)")
    for case, th in cases:
        vals = {}; Ts = ("0.4", "0.2", "0.1") if th == -80 else ("0.2", "0.1", "0.05")
        for T in Ts:
            t = mpf(T)*exp(I*th*pi/180); x = I*t/(2*pi)
            vals[T] = (t, PJ(x, delta=mp.arg(sqrt(x)), case=case))
        for K in ((1, 3, 6, 10) if case == "psi1" else (2, 3, 6, 10)):      # a_1 = 0 for the psi0 series
            out = []
            for T in Ts:
                t, v = vals[T]
                part = sum(coeff(case, k)*t**k for k in range(K))
                out.append(nstr((v - part)/(coeff(case, K)*t**K), 8))
            print(f"    {case}, theta = {th:4d}, |t| = {', '.join(Ts)}, K = {K:2d}:  " + "    ".join(out))
        sys.stdout.flush()

if __name__ == "__main__":
    parts = sys.argv[1:] or ["a", "b", "c", "d", "e"]
    part_e = lambda: part_d(cases=(("psi0", -45),), label="(e)")
    for pt in parts: {"a": part_a, "b": part_b, "c": part_c, "d": part_d, "e": part_e}[pt]()
